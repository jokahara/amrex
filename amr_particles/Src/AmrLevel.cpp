
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
    children.resize(state->local_size(), -1);

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
    const FabArrayBase::FB& TheFB = state->FabArrayBase::getFB(ngrow, geom.periodicity());
    for(auto& tag : *TheFB.m_LocTags) 
    {
        const auto& sfab = state->get(tag.srcIndex);
              auto& dfab = state->get(tag.dstIndex);

        dfab.copy<RunOn::Host>(sfab, tag.sbox, 0, tag.dbox, 0, ncomp);
    }
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
    //state->clearThisBD(true);
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
        coarse.FillPatchSingleLevel(*clev.state, icomp, icomp, ncomp, cgeom);
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

void
AmrLevel::post_regrid(int lbase, int new_finest)
{
    // make coarse/fine boundaries
    if (this->level < new_finest) {
        AmrLevel& fine_level = (*parent)[this->level + 1];
        constructCrseFineBdry(&fine_level);
    }
    else {
        constructCrseFineBdry(nullptr);
    }
};

// Creates fine boundary for this and coarse boundary for fine level
void AmrLevel::constructCrseFineBdry(AmrLevel* fineLevel) 
{
    if (level == 0) coarse_boundary.reset();

    // By default all cells are propagated
    children.resize(state->local_size(), -1);

    if (fineLevel == nullptr) {
        fine_boundary.reset();
        return;
    }

    const int ncomp = state->nComp();
    const int rank = ParallelDescriptor::MyProc();

    // CellFabs overlapping with fine region will not be propagated
    BoxArray grids_with_children(fineLevel->boxArray());
    grids_with_children.coarsen(fine_ratio);
    
    // FINE BOUNDARY FOR COARSE/THIS LEVEL

    Box cdomain(Domain());
    const IntVect ngrow = state->nGrowVect();

    // growing periodic domains so that all ghost cells are included
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            cdomain.grow(i, ngrow[i]);
        }
    }

    // the propagated area
    BoxArray crse_ba = complementIn(cdomain, grids_with_children);
    crse_ba.grow(ngrow);
    BoxList crse_boxes(crse_ba);

    // Add only local boxes to boundary box list
    // and set up parent-child connections
    BoxList cb_boxes; 
	for (MFIter it_fine(*fineLevel->state); it_fine.isValid(); ++it_fine) 
    {
        Box cbox = it_fine.validbox().coarsen(fine_ratio);
        for (MFIter it_crse(*state); it_crse.isValid(); ++it_crse) 
        {
            if (cbox == it_crse.validbox()) 
            {
                children[it_crse.LocalIndex()] = it_fine.LocalIndex();
                cb_boxes.push_back(cbox);
                break;
            }
        }
    }
    cb_boxes.intersect(crse_boxes);
    cb_boxes = removeOverlap(cb_boxes);
    // sorts and combines boxes if possible
    cb_boxes.simplify();
    
    if (cb_boxes.isEmpty()) {
        fine_boundary.reset();
    }
    else
    {
        BoxArray cb_ba(cb_boxes);
        // all boxes belong to this rank
        DistributionMapping dmap { Vector<int>(cb_ba.size(), rank) };
        // cells to be filled from fine level
        CellFabArray* new_coarse = new CellFabArray(cb_ba, dmap, ncomp, 0);
        // refine to get the fine area
        cb_ba.refine(fine_ratio);
        CellFabArray* new_fine = new CellFabArray(cb_ba, dmap, ncomp, 0);
        
        // fine boundary of this level (this also connects pointers to level data)
        fine_boundary.reset(new_coarse, this, new_fine, fineLevel, fineLevel);
    }

    // COARSE BOUNDARY OF FINE LEVEL

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
            fb_boxes.join(amrex::complementIn(b, fine_boxes));

        }
	}
    fb_boxes = removeOverlap(fb_boxes);
    fb_boxes.simplify();
    
    if (fb_boxes.isEmpty()) {
        fineLevel->coarse_boundary.reset();
    } 
    else 
    {
        BoxArray fb_ba(fb_boxes);
        // all boxes belong to this rank
        DistributionMapping dmap { Vector<int>(fb_ba.size(), rank) };
        // cells to interpolated to from coarse data
        CellFabArray* new_fine = new CellFabArray(fb_ba, dmap, ncomp, 0);

        // coarsen to get coarse area
        fb_ba.coarsen(fine_ratio);
        CellFabArray* new_coarse = new CellFabArray(fb_ba, dmap, ncomp, 0);
        
        // coarse boundary of fine level
        fineLevel->coarse_boundary.reset(new_coarse, this, new_fine, fineLevel, this);
    }
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

// Tag cells for refinement
void 
AmrLevel::errorEst (TagBoxArray& tb,
                    int          clearval,
                    int          tagval,
                    Real         time,
                    int          n_error_buf,
                    int          ngrow) 
{
#define MAX_PARTICLES 12
    
/*#ifdef _OPENMP
#pragma omp parallel
#endif*/
    for (AmrIter it(*this); it.isValid(); ++it)
    {
        Box box = it.validbox();
        auto tagArray  = tb[it].array();
        auto dataArray = (*this)[it].array();
        
        int childLID = it.childLID();
        if (childLID < 0)
        {
            ParallelFor(box, [&] (int i, int j, int k)
            {
                // Tag cells meeting chosen criteria.
                // Cells are clearval by default.
                if (dataArray(i,j,k)->number_of_particles >= MAX_PARTICLES) {
                    tagArray(i,j,k) = tagval;
                }
            });
        }
        else
        {
            // Cells with children are tagged by default
            tb[it].setVal(tagval, box);
            
            AmrLevel& fineLevel = (*parent)[level+1];

            // If the child cells also have children this CellFab is kept refined.
            if (fineLevel.children[childLID] >= 0) continue;

            const auto& fineArray = fineLevel.state->atLocalIdx(childLID).array();
            
            ParallelFor(box, [&] (int i, int j, int k)
            {
                IntVect iv{i,j,k};
                Box fine_box(Box::TheUnitBox());
                fine_box.shift(iv).refine(fine_ratio);
                
                int sum = 0;
                // iterate through children
                For(fine_box, [&] (int ii, int jj, int kk) 
                {
                    Cell* child = fineArray(ii,jj,kk).get();
                    sum += child->particles.size();
                });

                if (sum < MAX_PARTICLES) {
                    // unrefine cell if requirement is met
                    tagArray(i,j,k) = clearval;
                }
            });
        }
    }
}

void 
AmrLevel::manual_tags_placement (TagBoxArray&           tags,
                                 const Vector<IntVect>& bf_lev)
{
    if (level > 0) return;

    // making sure that the cells at periodic boundaries are not refined
    Box domain = geom.Domain();
    Box domain_g = geom.Domain();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (geom.isPeriodic(i)) {
            domain.grow(i,-1);
            domain_g.grow(i, 1);
        }
    }
    BoxArray boxes = amrex::boxComplement(domain_g, domain);

    tags.setVal(boxes, TagBox::CLEAR);
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
    const IndexType& boxType = leveldata.boxArray().ixType();
    const int level = Amrlevel.Level();

    if (level == 0)
    {
        CellFabArray& old = Amrlevel.getCells();
        const Geometry& geom = Amrlevel.Geom();
        leveldata.FillPatchSingleLevel(old, icomp, icomp, ncomp, geom);
    }
    else if (level == 1 || 
            ProperlyNested(Amrlevel.crse_ratio,
                           Amrlevel.parent->blockingFactor(level),
                           boxGrow, boxType, &my_interpolater))
    {
        AmrLevel& fine_level = Amrlevel;
        AmrLevel& crse_level = Amrlevel.parent->getLevel(level-1);

        const Geometry& fine_geom = fine_level.geom;
        const Geometry& crse_geom = crse_level.geom;
        
        CellFabArray& coarse = *crse_level.state;
        CellFabArray& fine = *fine_level.state;

        MyInterpolater* mapper = &my_interpolater;

        BoxArray area_to_fill = amrex::complementIn(fine_geom.Domain(), fine.boxArray());
        area_to_fill = amrex::intersect(area_to_fill, leveldata.boxArray());

        // coarse should already have same distribution as leveldata
        for (AmrIter it(crse_level); it.isValid(); ++it)
        {
            auto& sfab = coarse[it];

            Box cbox = it.validbox();
            cbox.refine(crse_ratio);
            
            for (AmrIter it2(*this); it2.isValid(); ++it2)
            {
                Box fine_region = it2.validbox() & cbox;
                if (fine_region.isEmpty()) continue;
                
                auto& dfab = leveldata[it2];
                //interpolate to fine fab
                mapper->interp(sfab, icomp, dfab, icomp, ncomp, fine_region, 
                                crse_ratio, crse_geom, fine_geom, RunOn::Cpu);
            }   
        }

        leveldata.FillPatchSingleLevel(fine, icomp, icomp, ncomp, geom);
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

    src = source_level; // where data is filled from

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

    // connecting local ghost cells to parent cells
    for (auto &tag : *cpc.m_LocTags)
    {
        auto const& src = all_cells.array(tag.srcIndex);
        auto const& dst = bndry_data.array(tag.dstIndex);
              
        const auto dlo = amrex::lbound(tag.dbox);
        const auto slo = amrex::lbound(tag.sbox);
        const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
        
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(RunOn::Host, tag.dbox, ncomp, i, j, k, n,
        {
            auto& d = dst(i,j,k,n);
            auto& s = src(i+offset.x,j+offset.y,k+offset.z,n);
            if (d) {
                // pointer has already been connected to different cell
                s = d;
            }
            else d = s;
        });
    }
} 

void
AmrLevel::FillCoarseToFine()
{
    if (!coarse_boundary.fine) return;

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
        
        // interpolate from coarse to fine cells
        mapper->interp(sfab, 0, dfab, 0, ncomp, mfi.validbox(), 
                        ratio, cgeom, fgeom, RunOn::Cpu);
    }
}

void
AmrLevel::FillFineToCoarse()
{
    if (!fine_boundary.coarse) return;
    
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
        
        // average from fine to coarse level
        mapper->average(dfab, 0, sfab, 0, ncomp, mfi.validbox(), 
                        ratio, cgeom, fgeom, RunOn::Cpu);
    }
}