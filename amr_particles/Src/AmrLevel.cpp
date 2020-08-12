
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

#include "AmrLevel.h"

AmrLevel::AmrLevel () noexcept
{
   BL_PROFILE("AmrLevel::AmrLevel()");
   parent = 0;
   level = -1;
}

AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
		            const DistributionMapping& dm)
: geom(level_geom), grids(ba), dmap(dm) 
{
    BL_PROFILE("AmrLevel::AmrLevel()");
    parent = &papa;
    level = lev;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
    }
    if (level < parent->maxLevel())
    {
        fine_ratio = parent->refRatio(level);
    }

    m_factory.reset(new CellFabFactory());
    state.reset(new CellFabArray(grids, dmap, parent->nComp(), parent->nGrow(), MFInfo(), *m_factory));
    
    // By default whole level is propagated
    propagate_fab.resize(state->local_size());
    propagate_fab.assign(propagate_fab.size(), true);

    if (parent->useFixedCoarseGrids()) constructAreaNotToTag();
}

void
AmrLevel::init (AmrLevel* old) 
{
    const int ncomp = state->nComp();
    const IntVect ngrow = state->nGrowVect();

    // make empty initial cells (including ghost cells)
    for (AmrIter mfi(*this); mfi.isValid(); ++mfi) 
    {
        const Box& box = mfi.growntilebox(ngrow);
        auto cfab = state->array(mfi);
        
        Loop(box, ncomp, [&] (int i, int j, int k, int n) 
        {
            if (!cfab(i,j,k,n))
                cfab(i,j,k,n).reset(new Cell());
        });
    }

    if (old) {
        // Fill level from old data
        this->FillPatch(*old, 0, 0, state->n_comp); 
    }
    
    // connect pointers to local ghost cells
    const FabArrayBase::FB& TheFB = state->getFB(ngrow, geom.periodicity());
    for(auto& tag : *TheFB.m_LocTags) 
    {
        const auto* sfab = &(state->get(tag.srcIndex));
              auto* dfab = &(state->get(tag.dstIndex));
        dfab->copy<RunOn::Host>(*sfab, tag.sbox, 0, tag.dbox, 0, ncomp);
    }	
}

void 
AmrLevel::errorEst (TagBoxArray& tb,
                    int          clearval,
                    int          tagval,
                    Real         time,
                    int          n_error_buf,
                    int          ngrow) 
{
    //Print() << "Tagging cells at level = " << level <<"\n";
    //int tagged = 0;

    // Tag cells for refinement
    for (MFIter mfi(tb); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        auto tagArray  = tb[mfi].array();
        auto dataArray = (*state)[mfi].array();
/*#ifdef _OPENMP
#pragma omp parallel
#endif*/
        ParallelFor(box, [&] (int i, int j, int k)
        {
        // TODO: tag if overlaps with fine areas meeting criteria
            if (i > 5 && j > 5 /*dataArray(i,j,k)->number_of_particles >= 50*/) {
                tagArray(i,j,k) = tagval;
                //tagged++;
            }
            else {
                tagArray(i,j,k) = clearval;
            }
        });
    }
    
    //Print() << "Tagged cells: " << tagged << "\n";
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
}

// Fill an entire CellFabArray by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrLevel::FillFromCoarsePatch (CellFabArray& dest,
                               int       icomp,
                               int       ncomp,
                               int       ngrow)
{
    BL_PROFILE("AmrLevel::FillCoarsePatch()");

    // Must fill this region on crse level and interpolate.
    BL_ASSERT(level != 0);
    BL_ASSERT(ngrow <= mf.nGrow());

    const Box&                 pdomain = Domain();
    const BoxArray&            dest_BA = dest.boxArray();
    const DistributionMapping& dest_DM = dest.DistributionMap();
    AmrLevel&                  clev    = parent->getLevel(level-1);
    const Geometry&            cgeom   = clev.geom;
    MyInterpolater*            mapper  = &my_interpolater;

    BoxArray crseBA(dest_BA.size());
    Box domain_g(pdomain);
    // DISTRIBUTION SHOULD ALREADY BE CORRECT
    
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            domain_g.grow(i,ngrow);
        }
    }
    
    for (int j = 0, N = crseBA.size(); j < N; ++j)
    {
        const Box bx = amrex::grow(dest_BA[j], ngrow) & domain_g;
        crseBA.set(j, mapper->CoarseBox(bx, crse_ratio));
    }
    
    // should contain only data needed by interpolater
    CellFabArray coarse(crseBA, dest_DM, ncomp, 0);
        
    // ProperlyNested checks that cells only have neighbours one level coarser
    if ( level == 1 
        || ProperlyNested(crse_ratio, parent->blockingFactor(level),
                            ngrow, dest_BA.ixType(), mapper))
    {
        // Filling coarse patch
        coarse.FillPatchSingleLevel(IntVect::Zero, *clev.state, icomp, icomp, ncomp, cgeom);
    }
    else
    {
        Abort("level was not properly nested\n");
    }
    
    for (AmrIter mfi(*this); mfi.isValid(); ++mfi)
    {
        auto& sfab = coarse[mfi];         // source data
        auto& dfab = dest[mfi];           // destination data
        const Box dbox = mfi.validbox();
        
        // interpolating to fine level
        mapper->interp(sfab, icomp, dfab, icomp, ncomp, dbox, 
                        crse_ratio, cgeom, geom, RunOn::Cpu);
        
        // clear source cells
        Box cbox = amrex::coarsen(dbox, crse_ratio);
        ParallelFor(cbox, ncomp, 
        [&] (int i, int j, int k, int n) {
            sfab(IntVect(i,j,k), n)->clear();
        });
    }
}

// Creates fine boundary for this and coarse boundary for fine level
void AmrLevel::constructCrseFineBdry(AmrLevel* fineLevel) 
{
    if (level == 0) coarse_boundary.reset();

    // By default all cells are propagated
    propagate_fab.resize(state->local_size());
    propagate_fab.assign(propagate_fab.size(), true);

    if (fineLevel == nullptr) {
        fine_boundary.reset();
        return;
    }

    const int ncomp = state->nComp();
    const int rank = ParallelDescriptor::MyProc();

    // CellFabs overlapping with fine region will not be propagated
    BoxArray grids_to_not_propagate(fineLevel->boxArray());
    grids_to_not_propagate.coarsen(fine_ratio);
    
    // FINE TO THIS/COARSE BOUNDARY

    Box cdomain(Domain());
    const IntVect ngrow = state->nGrowVect();

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            cdomain.grow(i, ngrow[i]);
        }
    }

    // propagated area
    BoxArray crse_ba = complementIn(cdomain, grids_to_not_propagate);
    crse_ba.grow(ngrow);
    BoxList crse_boxes(crse_ba);

    // We will only add local boxes and set which ones to propagate
    BoxList cb_boxes; 
	for (MFIter mfi(*state); mfi.isValid(); ++mfi) 
    {
        Box b = mfi.validbox();
        if (grids_to_not_propagate.contains(b)) 
        {
            propagate_fab[mfi.LocalIndex()] = false;
            cb_boxes.push_back(b);
        }
    }
    cb_boxes.intersect(crse_boxes);
    cb_boxes.simplify();
    
    if (cb_boxes.isEmpty()) 
    {
        fine_boundary.reset();
    }
    else
    {
        BoxArray cb_ba(cb_boxes);
        DistributionMapping dmap { Vector<int>(cb_boxes.size(), rank) };

        CellFabArray* new_coarse = new CellFabArray(cb_ba, dmap, ncomp, 0);
        cb_ba.refine(fine_ratio);
        CellFabArray* new_fine = new CellFabArray(cb_ba, dmap, ncomp, 0);
        
        // fine boundary of this level
        fine_boundary.reset(new_coarse, this, new_fine, fineLevel, fineLevel);
    }

    // THIS/COARSE TO FINE BOUNDARY

    const BoxList fine_boxes(fineLevel->boxArray());
	const DistributionMapping fine_dmap = fineLevel->DistributionMap();
    
    Box fdomain(fineLevel->Domain());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            fdomain.grow(i, ngrow[i]);
        }
    }

    // local boxes at the coarse-to-fine boundary of fineLevel
	BoxList fb_boxes;   
	for (int i = 0; i < fine_boxes.size(); i++)
	{
        if (rank == fine_dmap[i]) 
        {
            Box b = amrex::grow(fine_boxes.data()[i], ngrow);
            b &= fdomain;
            BoxList complement = amrex::complementIn(b, fine_boxes);
            fb_boxes.join(std::move(complement));
        }
	}
    // combine boxes if possible and remove overlaps
    fb_boxes.simplify();
    
    if (fb_boxes.isEmpty()) 
    {
        coarse_boundary.reset();
    } 
    else 
    {
        BoxArray fb_ba(fb_boxes);
        DistributionMapping dmap { Vector<int>(fb_boxes.size(), rank) };

        CellFabArray* new_fine = new CellFabArray(fb_ba, dmap, ncomp, 0);

        fb_ba.coarsen(fine_ratio);
        CellFabArray* new_coarse = new CellFabArray(fb_ba, dmap, ncomp, 0);
        
        // coarse boundary of fine level
        fineLevel->coarse_boundary.reset(new_coarse, this, new_fine, fineLevel, this);
    }
    
    //FabArrayBase::flushCPCache();
}
 
void AmrLevel::constructAreaNotToTag ()
{
    if (level == 0 || !parent->useFixedCoarseGrids() || parent->useFixedUpToLevel()>level)
        return;
        
    // We are restricting the tagging on the finest fixed level
    if (parent->useFixedUpToLevel()==level)
    {
        // We use the next coarser level shrunk by one blockingfactor
        //    as the region in which we allow tagging. 
        // Why level-1? Because we always use the full domain at level 0 
        //    and therefore level 0 in initialba is level 1 in the AMR hierarchy, etc.
        const Vector<BoxArray>& initialba = parent->getInitialBA();
        Box tagarea(initialba[level-1].minimalBox());
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;

        // We disallow tagging in the remaining part of the domain.
        BoxArray tagba = amrex::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;

        BoxArray bxa(parent->Geom(level).Domain());
        BL_ASSERT(bxa.contains(m_AreaNotToTag));
    }

    if (parent->useFixedUpToLevel()<level)
    {
        Box tagarea = parent->getLevel(level-1).getAreaToTag();
        tagarea.refine(parent->refRatio(level-1));
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;
        BoxArray tagba = amrex::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;
    }
}

void
AmrLevel::FillPatch (AmrLevel&  Amrlevel,
                     int        boxGrow,
                     int        icomp,
                     int        ncomp)
{
    BL_ASSERT(boxGrow <= leveldata.nGrow());
    BL_ASSERT(icomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= idx && idx < AmrLevel::state.size());

    CellFabArray& leveldata = *state;
    CellFabArray& old = Amrlevel.getCells();
    const Geometry& geom = Amrlevel.Geom();
    const IndexType& boxType = leveldata.boxArray().ixType();
    const int level = Amrlevel.Level();

    Print() << "Filling Patch Level " << level << "\n";
    if (level == 0)
    {
        leveldata.FillPatchSingleLevel(IntVect::Zero, old, icomp, icomp, ncomp, geom);

        // TODO: should also take averages from fine levels
        if (level != parent->finestLevel()) 
        {
            //leveldata.AverageFromFineLevel();
        }
    }
    else if (level == 1 || 
            ProperlyNested(Amrlevel.crse_ratio,
                           Amrlevel.parent->blockingFactor(level),
                           boxGrow, boxType, &my_interpolater))
    {
        Print() << "    FillFromTwoLevels\n";
        AmrLevel& fine_level = Amrlevel;
        AmrLevel& crse_level = Amrlevel.parent->getLevel(level-1);

        const Geometry& geom_fine = fine_level.geom;
        const Geometry& geom_crse = crse_level.geom;
        
        CellFabArray& coarse = *crse_level.state;
        CellFabArray& fine = *fine_level.state;

        leveldata.FillPatchTwoLevels(coarse, fine, icomp, icomp, ncomp, 
                                     geom_crse, geom_fine, Amrlevel.crse_ratio);
    }
    else Abort("FillPatch failed!!!\n");
}

void
AmrLevel::BoundaryContainer::reset(CellFabArray* new_coarse, AmrLevel* coarse_parent, 
                                   CellFabArray* new_fine, AmrLevel* fine_parent,
                                   AmrLevel* source_level) noexcept
{
    coarse.reset(new_coarse);
    fine.reset(new_fine);

    src = source_level;

    if (coarse) connect_cells(*coarse, coarse_parent);
    if (fine)   connect_cells(*fine, fine_parent);
}

void
AmrLevel::BoundaryContainer::connect_cells(CellFabArray& bndry_data, AmrLevel* parent)
{
    int ncomp = bndry_data.nComp();
    CellFabArray& all_cells = parent->getCells();

    // create tags for copying
    const FabArrayBase::CPC &cpc = bndry_data.getCPC(
        bndry_data.nGrowVect(), all_cells, all_cells.nGrowVect(), parent->geom.periodicity()
    );
    /*const FabArrayBase::CPC &cpc = all_cells.getCPC(
        all_cells.nGrowVect(), bndry_data, IntVect::Zero, parent->geom.periodicity()
    );*/

    // connecting local ghost cells to parent cells
    for (auto &tag : *cpc.m_LocTags)
    {
        const auto* sfab = &(all_cells.get(tag.srcIndex));
              auto* dfab = &(bndry_data.get(tag.dstIndex));
        //Print() << "copy: " << tag.sbox << " -> " << tag.dbox << "\n";
        dfab->copy<RunOn::Host>(*sfab, tag.sbox, 0, tag.dbox, 0, ncomp);
        
        /*auto src_array = all_cells.array(tag.dstIndex);
        auto data_array = data.array(tag.srcIndex);
        
        //if (parent->level > 0) continue;
        Print() << "  " << tag.dbox << " -> " << tag.sbox << "\n";
        Dim3 offset = (tag.sbox.smallEnd() - tag.dbox.smallEnd()).dim3();
        
        Print() << "  " << offset << "\n";
        //continue;
        
        ParallelFor(tag.sbox, ncomp, [&] (int i, int j, int k, int n) 
        {
            // connecting pointers
            auto& p = data_array(i,j,k,n);
            auto& src = src_array(i+offset.x, j+offset.y, k+offset.z, n);
            if(!p) {
                p = src;
                //p->number_of_particles = 11;
                //Print() << "connect " << IntVect{i,j,k} << ", " << p->number_of_particles << "\n";
            }
            //else src = p;
        });*/
    }
} 

void
AmrLevel::FillCoarseToFine()
{
    if (!coarse_boundary.fine) return;
    //Print() << "Fill " << level-1 << "->" << level << "\n";

    CellFabArray* fine     = coarse_boundary.fine.get();
    CellFabArray* coarse   = coarse_boundary.coarse.get();
    AmrLevel* coarse_level = coarse_boundary.src;

    const int ncomp        = fine->nComp();
    const IntVect& ratio   = coarse_level->fineRatio();
    const Geometry& cgeom  = coarse_level->Geom();
    const Geometry& fgeom  = this->Geom();
    MyInterpolater* mapper = &my_interpolater;

    for (MFIter mfi(*fine); mfi.isValid(); ++mfi)
    {
        auto& sfab = (*coarse)[mfi];
        auto& dfab = (*fine)[mfi];
        
        //Print() << "  " << sfab.box() << "->" << dfab.box() << "\n";
        // interpolate from coarse to fine cells
        mapper->interp(sfab, 0, dfab, 0, ncomp, mfi.validbox(), 
                        ratio, cgeom, fgeom, RunOn::Cpu);
    }
}

void
AmrLevel::FillFineToCoarse()
{
    if (!fine_boundary.coarse) return;
    //Print() << "Fill " << level+1 << "->" << level << "\n";
    CellFabArray* fine   = fine_boundary.fine.get();
    CellFabArray* coarse = fine_boundary.coarse.get();
    AmrLevel* fine_level = fine_boundary.src;

    const int ncomp        = coarse->nComp();
    const IntVect& ratio   = this->fineRatio();
    const Geometry& cgeom  = this->Geom();
    const Geometry& fgeom  = fine_level->Geom();
    MyInterpolater* mapper = &my_interpolater;

    for (MFIter mfi(*fine); mfi.isValid(); ++mfi)
    {
        auto& sfab = (*fine)[mfi];
        auto& dfab = (*coarse)[mfi];
        
        //Print() << "  " << sfab.box() << "->" << dfab.box() << "\n";
        mapper->average(dfab, 0, sfab, 0, ncomp, mfi.validbox(), 
                        ratio, cgeom, fgeom, RunOn::Cpu);
    }
}