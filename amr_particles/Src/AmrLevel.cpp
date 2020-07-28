
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

#include "AmrLevel.h"
#include "Interpolater.h"

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
		            const DistributionMapping& dm,
                    Real            time)
    : geom(level_geom), grids(ba), dmap(dm), areaToPropagate(ba)
{
    BL_PROFILE("AmrLevel::AmrLevel(dm)");
    level  = lev;
    parent = &papa;

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
    
    state.define(grids, dmap, parent->nComp(), parent->nGrow(), MFInfo(), *m_factory);

    if (parent->useFixedCoarseGrids()) constructAreaNotToTag();
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
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(tb); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            auto tagArray  = tb[mfi].array();
            auto dataArray = state[mfi].array();

            //tag cells for refinement
            ParallelFor(box, [&] (int i, int j, int k)
            {
                if (dataArray(i,j,k).number_of_particles >= 6) {
                    tagArray(i,j,k) = tagval;
                    //tagged++;
                }
                else {
                    tagArray(i,j,k) = clearval;
                }
            });
        }
    }
    //Print() << "Tagged cells: " << tagged << "\n";
}

Long
AmrLevel::countCells () const noexcept
{
    return grids.numPts();
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
}

const BoxArray&
AmrLevel::getEdgeBoxArray (int dir) const noexcept
{
    BL_ASSERT(dir >=0 && dir < AMREX_SPACEDIM);
    if (edge_grids[dir].empty()) {
        edge_grids[dir] = grids;
        edge_grids[dir].surroundingNodes(dir);
    }
    return edge_grids[dir];
}

const BoxArray&
AmrLevel::getNodalBoxArray () const noexcept
{
    if (nodal_grids.empty()) {
        nodal_grids = grids;
        nodal_grids.surroundingNodes();
    }
    return nodal_grids;
}

// Interpolate from coarse level to the valid area in dest.
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

    Box domain_g = pdomain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            domain_g.grow(i,ngrow);
        }
    }
    
    BoxArray crseBA(dest_BA.size());
    
    for (int j = 0, N = crseBA.size(); j < N; ++j)
    {
        const Box& bx = amrex::grow(dest_BA[j],ngrow) & domain_g;
        crseBA.set(j, amrex::coarsen(bx, crse_ratio).grow(ngrow));
    }

    // contains only data needed by interpolater
    CellFabArray coarse(crseBA, dest_DM, ncomp, 0);
    
    // ProperlyNested checks that cells only have neighbours one level coarser
    if ( level == 1 
        || ProperlyNested(crse_ratio, parent->blockingFactor(level),
                            ngrow, dest_BA.ixType()) )
    {
        coarse.FillPatchSingleLevel({0,0,0}, clev.state, icomp, 0, ncomp, cgeom);
    }
    else
    {
        Abort("level was not properly nested\n");
        // TODO: FillPatch should not be necessary for finer levels.
        // FillPatch(clev, coarse, 0, icomp, ncomp);
    }

#ifdef _OPENMP
//#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(coarse); mfi.isValid(); ++mfi)
    {
        auto& sfab = coarse[mfi];       // source data
        auto& dfab = dest[mfi];         // destination data
        const Box dbox = dfab.box();    
        const Box sbox = mfi.validbox();
        // interpolating to dest
        dest.interpolateTo(dfab, dbox, geom, sfab, sbox, cgeom);

        // TODO: clear interpolated data (maybe swap during FillPatch would be in fact better)
        /*ParallelFor(sbox, ncomp, [&] (int i, int j, int k, int n) 
        {
            const IntVect iv(i,j,k);
            Cell &parent = sfab(iv);
            
        });*/
    }

    if (dest.nGrow() > 0) {
        //dest.FillPatchTwoLevels();
        //or dest.FillBoundary(geom.periodicity(), true);
    }   
}

const BoxArray& AmrLevel::getAreaNotToTag () noexcept
{
    return m_AreaNotToTag;
}

const Box& AmrLevel::getAreaToTag () noexcept
{
    return m_AreaToTag;
}

void AmrLevel::setAreaNotToTag (BoxArray& ba) noexcept
{
    m_AreaNotToTag = ba;
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

FillPatchIteratorHelper::FillPatchIteratorHelper (AmrLevel&  Amrlevel,
                                                  CellFabArray& leveldata)
    :
    m_Amrlevel(Amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_Amrlevel.level+1)
{}

FillPatchIterator::FillPatchIterator (AmrLevel& Amrlevel,
                                      CellFabArray& leveldata)
    :
    MFIter(leveldata),
    m_Amrlevel(Amrlevel),
    m_leveldata(leveldata),
    m_ncomp(0)
{}

FillPatchIteratorHelper::FillPatchIteratorHelper (AmrLevel&     Amrlevel,
                                                  CellFabArray& leveldata,
                                                  int           boxGrow,
                                                  Real          time,
                                                  int           icomp,
                                                  int           ncomp)
    :
    m_Amrlevel(Amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_Amrlevel.level+1),
    m_time(time),
    m_growsize(boxGrow),
    m_icomp(icomp),
    m_ncomp(ncomp)
{
    Initialize(boxGrow,time,icomp,ncomp);
}

FillPatchIterator::FillPatchIterator (AmrLevel&     Amrlevel,
                                      CellFabArray& leveldata,
                                      int           boxGrow,
                                      Real          time,
                                      int           icomp,
                                      int           ncomp)
    :
    MFIter(leveldata),
    m_Amrlevel(Amrlevel),
    m_leveldata(leveldata),
    m_ncomp(ncomp)
{
    Initialize(boxGrow,time,icomp,ncomp);
}

static
bool
NeedToTouchUpPhysCorners (const Geometry& geom)
{
    return geom.isAnyPeriodic() && !geom.isAllPeriodic();
}

void
FillPatchIterator::Initialize (int  boxGrow,
                               Real time,
                               int  icomp,
                               int  ncomp)
{
    BL_PROFILE("FillPatchIterator::Initialize");

    BL_ASSERT(icomp >= 0);
    BL_ASSERT(ncomp >= 1);

    m_ncomp = ncomp;

    m_fabs.define(m_leveldata.boxArray(), m_leveldata.DistributionMap(),
		          m_ncomp, boxGrow, MFInfo(), m_leveldata.Factory());

    const Geometry& geom = m_Amrlevel.Geom();

    const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_Amrlevel.level;

    if (level == 0)
    {
        FillFromLevel0(time, icomp, ncomp);
    }
    else
    {
        if (level == 1 || 
        ProperlyNested(m_Amrlevel.crse_ratio,
                        m_Amrlevel.parent->blockingFactor(m_Amrlevel.level),
                        boxGrow, boxType))
        {
            FillFromTwoLevels(time, icomp, ncomp);
        } else {
            static bool first = true;
            if (first) {
                first = false;
                if (ParallelDescriptor::IOProcessor() && amrex::Verbose()) {
                    IntVect new_blocking_factor = m_Amrlevel.parent->blockingFactor(m_Amrlevel.level);
                    new_blocking_factor *= 2;
                    for (int j = 0; j < 10; ++j) {
                        if (ProperlyNested(m_Amrlevel.crse_ratio,
                                            new_blocking_factor,
                                            boxGrow, boxType)) {
                            break;
                        } else {
                            new_blocking_factor *= 2;
                        }
                    }
                    amrex::Print() << "WARNING: Grids are not properly nested.  We might have to use\n"
                                << "         two coarse levels to do fillpatch.  Consider using\n";
                    if (new_blocking_factor < IntVect{AMREX_D_DECL(128,128,128)}) {
                        amrex::Print() << "         amr.blocking_factor=" << new_blocking_factor << "\n";
                    } else {
                        amrex::Print() << "         larger amr.blocking_factor.\n";
                    }
                }
            }

            Print() << "    FillPatchIteratorHelper\n";
            // This initialization sends all data needed by fill
            FillPatchIteratorHelper* fph = 0;
            fph = new FillPatchIteratorHelper(m_Amrlevel, m_leveldata, boxGrow, time, icomp, ncomp);

#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
            for (MFIter mfi(m_fabs); mfi.isValid(); ++mfi)
            {
                fph->fill(m_fabs[mfi],mfi.index());
            }
            delete fph;
        }
    }
    
    // TODO:
    // Call hack to touch up fillPatched data.
    // m_Amrlevel.set_preferred_boundary_values(m_fabs, icomp, 0, ncomp, time);
}

void
FillPatchIterator::FillFromLevel0 (Real time, int icomp, int ncomp)
{
    Print() << "    FillFromLevel0\n";
    BL_ASSERT(m_Amrlevel.level == 0);

    CellFabArray& src = m_Amrlevel.state;
    //src.FillBoundary();

    const Geometry& geom = m_Amrlevel.geom;

    m_fabs.FillPatchSingleLevel(src, icomp, icomp, ncomp, geom);
    //amrex::FillPatchSingleLevel (m_fabs, time, smf, stime, icomp, ncomp, geom, physbcf, icomp);
}

void
FillPatchIterator::FillFromTwoLevels (Real time, int icomp, int ncomp)
{
    Print() << "    FillFromTwoLevels\n";
    int ilev_fine = m_Amrlevel.level;
    int ilev_crse = ilev_fine-1;

    BL_ASSERT(ilev_crse >= 0);

    AmrLevel& fine_level = m_Amrlevel;
    AmrLevel& crse_level = m_Amrlevel.parent->getLevel(ilev_crse);

    const Geometry& geom_fine = fine_level.geom;
    const Geometry& geom_crse = crse_level.geom;
    
    CellFabArray& smf_crse = crse_level.state;
    //smf_crse.FillBoundary();

    CellFabArray& smf_fine = fine_level.state;
    //smf_fine.FillBoundary();

    //m_fabs.FillPatchTwoLevels
}

void
FillPatchIteratorHelper::Initialize (int  boxGrow,
                                     Real time,
                                     int  icomp,
                                     int  ncomp)
{
    BL_PROFILE("FillPatchIteratorHelper::Initialize()");

    BL_ASSERT(icomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[idx].inRange(icomp,ncomp));
    BL_ASSERT(0 <= idx && idx < AmrLevel::state.size());

    m_time         = time;
    m_growsize     = boxGrow;
    m_icomp        = icomp;
    m_ncomp        = ncomp;

    const int         MyProc     = ParallelDescriptor::MyProc();
    auto&             AmrLevels  = m_Amrlevel.parent->getAmrLevels();
    const AmrLevel&   topLevel   = *AmrLevels[m_Amrlevel.level];
    const Box&        topPDomain = topLevel.Domain();
    const IndexType&  boxType    = m_leveldata.boxArray().ixType(); // should be cell centered

    for (int l = 0; l <= m_Amrlevel.level; ++l)
    {
        m_mfid[l].resize(1); 
        m_mfid[l][0] = m_mfcd.RegisterFabArray(&AmrLevels[l]->state);
    }
    for (int i = 0, N = m_leveldata.boxArray().size(); i < N; ++i)
    {
        //
        // A couple typedefs we'll use in the next code segment.
        //
        typedef std::map<int,Vector<Vector<Box> > >::value_type IntAABoxMapValType;

        typedef std::map<int,Vector<Vector<Vector<FillBoxId> > > >::value_type IntAAAFBIDMapValType;

        if (m_leveldata.DistributionMap()[i] != MyProc) continue;
        //
        // Insert with a hint since the indices are ordered lowest to highest.
        //
        IntAAAFBIDMapValType v1(i,Vector<Vector<Vector<FillBoxId> > >());

        m_fbid.insert(m_fbid.end(),v1)->second.resize(m_Amrlevel.level+1);

        IntAABoxMapValType v2(i,Vector<Vector<Box> >());

        m_fbox.insert(m_fbox.end(),v2)->second.resize(m_Amrlevel.level+1);
        m_cbox.insert(m_cbox.end(),v2)->second.resize(m_Amrlevel.level+1);

        m_ba.insert(m_ba.end(),std::map<int,Box>::value_type(i,amrex::grow(m_leveldata.boxArray()[i],m_growsize)));
    }

    BoxList         tempUnfillable(boxType);
    BoxList         unfillableThisLevel(boxType);
    Vector<Box>     unfilledThisLevel;
    Vector<Box>     crse_boxes;
    Vector<IntVect> pshifts(27);

    for (std::map<int,Box>::const_iterator it = m_ba.begin(), 
         End = m_ba.end(); it != End; ++it)
    {
        const int  bxidx = it->first;
        const Box& box   = it->second;

        unfilledThisLevel.clear();
        unfilledThisLevel.push_back(box);

        if (!topPDomain.contains(box))
        {
            unfilledThisLevel.back() &= topPDomain;

            if (topLevel.geom.isAnyPeriodic())
            {
                //
                // May need to add additional unique pieces of valid region
                // in order to do periodic copies into ghost cells.
                //
                topLevel.geom.periodicShift(topPDomain,box,pshifts);

                for (const auto& iv : pshifts)
                {
                    Box shbox = box + iv;
                    shbox    &= topPDomain;

                    // Boxes are assumed cell centered
                    /*if (boxType.nodeCentered())
                    {
                        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                        {
                            if (iv[dir] > 0)
                            {
                                shbox.growHi(dir,-1);
                            }
                            else if (iv[dir] < 0)
                            {
                                shbox.growLo(dir,-1);
                            }
                        }
                    }*/

                    if (shbox.ok())
                    {
                        BoxList bl = amrex::boxDiff(shbox,box);

                        unfilledThisLevel.insert(unfilledThisLevel.end(), bl.begin(), bl.end());
                    }
                }
            }
        }

	    // cells outside physical boundaries are not included in unfilledThisLevel

        bool Done = false;

        Vector< Vector<Box> >&                 TheCrseBoxes = m_cbox[bxidx];
        Vector< Vector<Box> >&                 TheFineBoxes = m_fbox[bxidx];
        Vector< Vector< Vector<FillBoxId> > >& TheFBIDs     = m_fbid[bxidx];

        for (int l = m_Amrlevel.level; l >= 0 && !Done; --l)
        {
            unfillableThisLevel.clear();

            AmrLevel&       theAmrLevel = *AmrLevels[l];
            CellFabArray&   theState    = theAmrLevel.state;
            const Box&      thePDomain  = theAmrLevel.Domain();
            const Geometry& theGeom     = theAmrLevel.geom;
            const bool      is_periodic = theGeom.isAnyPeriodic();
            const IntVect&  fine_ratio  = theAmrLevel.fine_ratio;
            Vector<Box>&    FineBoxes   = TheFineBoxes[l];
            //
            // These are the boxes on this level contained in thePDomain
            // that need to be filled in order to directly fill at the
            // highest level or to interpolate up to the next higher level.
            //
            FineBoxes = unfilledThisLevel;
            //
            // Now build coarse boxes needed to interpolate to fine.
            //
            // If we're periodic and we're not at the finest level, we may
            // need to get some additional data at this level in order to
            // properly fill the CoarseBox()d versions of the fineboxes.
            //
            crse_boxes.clear();

            for (const auto& fbx : FineBoxes)
            {
                crse_boxes.push_back(fbx);

                if (l != m_Amrlevel.level)
                {
                    const Box& cbox = amrex::coarsen(fbx,fine_ratio).grow(1);
		            crse_boxes.back() = cbox;

                    if (is_periodic && !thePDomain.contains(cbox))
                    {
                        theGeom.periodicShift(thePDomain,cbox,pshifts);

                        for (const auto& iv : pshifts)
                        {
                            Box shbox = cbox + iv;
                            shbox    &= thePDomain;

                            // Boxes are assumed cell centere
                            /*if (boxType.nodeCentered())
                            {
                                for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                                {
                                    if (iv[dir] > 0)
                                    {
                                        shbox.growHi(dir,-1);
                                    }
                                    else if (iv[dir] < 0)
                                    {
                                        shbox.growLo(dir,-1);
                                    }
                                }
                            }*/

                            if (shbox.ok())
                            {
                                crse_boxes.push_back(shbox);
                            }
                        }
                    }
                }
            }

            Vector< Vector<FillBoxId> >& FBIDs = TheFBIDs[l];
            Vector<Box>&             CrseBoxes = TheCrseBoxes[l];

            FBIDs.resize(crse_boxes.size());
            CrseBoxes.resize(crse_boxes.size());
            //
            // Now attempt to get as much coarse data as possible.
            //
            for (int i = 0, M = CrseBoxes.size(); i < M; i++)
            {
                BL_ASSERT(tempUnfillable.isEmpty());

                CrseBoxes[i] = crse_boxes[i];

                BL_ASSERT(CrseBoxes[i].intersects(thePDomain));

                /*theState.InterpAddBox(m_mfcd) -> FillFab*/
                
                unfillableThisLevel.catenate(tempUnfillable);
            }

            unfillableThisLevel.intersect(thePDomain);

            if (unfillableThisLevel.isEmpty())
            {
                Done = true;
            }
            else
            {
                unfilledThisLevel.clear();

                unfilledThisLevel.insert(unfilledThisLevel.end(),
                                         unfillableThisLevel.begin(),
                                         unfillableThisLevel.end());
            }
        }
    }
    Cell::transfer_particles = false;
    
    // handles data sends and receives
    m_mfcd.CollectData();
}

void
FillPatchIteratorHelper::fill (CellFab& fab,
                               int      idx)
{
    BL_PROFILE("FillPatchIteratorHelper::fill()");

    BL_ASSERT(fab.box() == m_ba[idx]);

    Vector< Vector<std::unique_ptr<CellFab> > > cfab(m_Amrlevel.level+1);
    Vector< Vector<Box> >&                 TheCrseBoxes = m_cbox[idx];
    Vector< Vector<Box> >&                 TheFineBoxes = m_fbox[idx];
    Vector< Vector< Vector<FillBoxId> > >& TheFBIDs     = m_fbid[idx];
    auto&                                  AmrLevels    = m_Amrlevel.parent->getAmrLevels();
    //
    // Build all coarse fabs from which we'll interpolate and
    // fill them with coarse data as best we can.
    //
    for (int l = 0; l <= m_Amrlevel.level; l++)
    {
        CellFabArray&                      TheState  = AmrLevels[l]->state;
        const Vector<Box>&                 CrseBoxes = TheCrseBoxes[l];
        auto&                              CrseFabs  = cfab[l];
        const Vector< Vector<FillBoxId> >& FBIDs     = TheFBIDs[l];
        const int                          NC        = CrseBoxes.size();

        CrseFabs.resize(NC);

        Box domain_box = amrex::convert(AmrLevels[l]->Geom().Domain(), fab.box().ixType());
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (AmrLevels[l]->Geom().isPeriodic(idim)) {
                int n = domain_box.length(idim);
                domain_box.grow(idim, n);
            }
        }

        for (int i = 0; i < NC; i++)
        {
            BL_ASSERT(CrseBoxes[i].ok());
            CrseFabs[i].reset(new CellFab(CrseBoxes[i],m_ncomp));
	    }

        for (int i = 0; i < NC; i++)
        {
            //
            // Set to special value we'll later check
            // to ensure we've filled the FABs at the coarse level.
            //
            // TODO: TheState.InterpFillFab
        }
    }
    //
    // Now work from the bottom up interpolating to next higher level.
    //
    for (int l = 0; l < m_Amrlevel.level; l++)
    {
        auto&              CrseFabs   = cfab[l];
        AmrLevel&          TheLevel   = *AmrLevels[l];
        CellFabArray&      TheState   = TheLevel.state;
        const Box&         ThePDomain = TheLevel.Domain();
        const int          NC         = CrseFabs.size();
        // TODO: Fillboundary with periodics only should work I think?
        TheState.FillBoundary();
        
        //
        // Interpolate up to next level.
        //
        AmrLevel&           crseAmrLevel  = *AmrLevels[l];
        AmrLevel&           fineAmrLevel  = *AmrLevels[l+1];
        const IntVect&      fine_ratio    = crseAmrLevel.fine_ratio;
        const Vector<Box>&  FineBoxes     = TheFineBoxes[l];
        CellFabArray&       fState        = fineAmrLevel.state;
        const Box&          fDomain       = fineAmrLevel.Domain();
        auto&               FinerCrseFabs = cfab[l+1];
        //const Vector<BCRec>& theBCs       = fState.getBCs();
        const int           NF            = FineBoxes.size();

        for (int ifine = 0; ifine < NF; ++ifine)
        {
            Vector<BCRec> bcr(m_ncomp);
            CellFab    finefab(FineBoxes[ifine],m_ncomp);
            CellFab    crsefab(amrex::coarsen(finefab.box(),fine_ratio).grow(0),m_ncomp);
            
            // TODO: get rid off copies:

            // Fill crsefab from m_cbox via copy on intersect.
            for (int j = 0; j < NC; j++) {
                crsefab.copy<RunOn::Host>(*CrseFabs[j]);
	        }
            
            // TODO: Interpolate up to fine patch.
            // m_map->interp

            // Copy intersect finefab into next level m_cboxes.
            for (int j = 0, K = FinerCrseFabs.size(); j < K; ++j) {
                FinerCrseFabs[j]->copy<RunOn::Host>(finefab);
            }
        }

        CrseFabs.clear();
    }
    //
    // Now for the finest level stuff.
    //
    CellFabArray&      FineState      = m_Amrlevel.state;
    const Box&         FineDomain     = m_Amrlevel.Domain();
    const Geometry&    FineGeom       = m_Amrlevel.geom;
    auto&              FinestCrseFabs = cfab[m_Amrlevel.level];
    //
    // Copy intersect coarse into destination fab.
    //
    for (int i = 0, N = FinestCrseFabs.size(); i < N; ++i) {
        fab.copy<RunOn::Host>(*FinestCrseFabs[i], 0, 0, m_ncomp);
    }
    // TODO: another periodic FillBoundary
    /*if (FineGeom.isAnyPeriodic() && !FineDomain.contains(fab.box()))
    {
        Vector<IntVect> pshifts(27);

        FineGeom.periodicShift(FineDomain,fab.box(),pshifts);

        for (int i = 0, N = FinestCrseFabs.size(); i < N; i++)
        {
            for (const auto& iv : pshifts)
            {
                fab.shift(iv);

                Box src_dst = FinestCrseFabs[i]->box() & fab.box();
                src_dst    &= FineDomain;

                if (src_dst.ok())
                    fab.copy<RunOn::Host>(*FinestCrseFabs[i],src_dst, 0, 0, m_ncomp);

                fab.shift(-iv);
            }
        }
    }*/
    //
    // No longer need coarse data at finest level.
    //
    FinestCrseFabs.clear();
}

static
bool
HasPhysBndry (const Box&      b,
              const Box&      dmn,
              const Geometry& geom)
{
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i))
        {
            if (b.smallEnd(i) < dmn.smallEnd(i) || b.bigEnd(i) > dmn.bigEnd(i))
            {
                return true;
            }
        }
    }

    return false;
}

void
AmrLevel::FillPatch (AmrLevel&  Amrlevel,
		             CellFabArray& leveldata,
                     int        boxGrow,
                     int        icomp,
                     int        ncomp)
{
    BL_ASSERT(boxGrow <= leveldata.nGrow());
    BL_ASSERT(icomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[idx].inRange(icomp,ncomp));
    BL_ASSERT(0 <= idx && idx < AmrLevel::state.size());

    Print() << "Filling Patch with FillPatchIterator\n";
    // fill leveldata with data from Amrlevel

    CellFabArray& old_data = Amrlevel.getData();
    const Geometry& geom = Amrlevel.Geom();
    const IndexType& boxType = leveldata.boxArray().ixType();
    const int level = Amrlevel.level;

    if (level == 0)
    {
        leveldata.FillPatchSingleLevel({0,0,0}, old_data, icomp, icomp, ncomp, geom);
    }
    else if (level == 1 || ProperlyNested(Amrlevel.crse_ratio,
                                        Amrlevel.parent->blockingFactor(Amrlevel.level),
                                        boxGrow, boxType))
    {
        Print() << "    FillFromTwoLevels\n";
        AmrLevel& fine_level = Amrlevel;
        AmrLevel& crse_level = Amrlevel.parent->getLevel(level-1);

        const Geometry& geom_fine = fine_level.geom;
        const Geometry& geom_crse = crse_level.geom;
        
        CellFabArray& coarse = crse_level.state;
        CellFabArray& fine = fine_level.state;
        //smf_crse.FillBoundary();
        //smf_fine.FillBoundary();

        leveldata.FillPatchTwoLevels(coarse, fine, icomp, icomp, ncomp, 
                                     geom_crse, geom_fine, Amrlevel.crse_ratio);
    }
    else Print() << "FillPatch failed!!!\n";

    //FillPatchIterator fpi(Amrlevel, leveldata, boxGrow, time, icomp, ncomp);
    //const CellFabArray& mf_fillpatched = fpi.get_mf();
    //amrex::Copy(leveldata, mf_fillpatched, icomp, icomp, ncomp, boxGrow);
}

// fill an entire CellFabArray by interpolating from the coarser level
// this comes into play when a new level of refinement appears
/*void
AmrLevel::FillCoarsePatch (CellFabArray& fine, int icomp, int ncomp, int ngrow)
{
    int lev = Level();

    Print() << "Interpolate to: " << lev << "\n";
    CellFabArray& coarse = parent->getLevel(lev-1).getData();
    const BoxArray& ba = fine.boxArray();
    const DistributionMapping& dm = fine.DistributionMap();

    const Geometry &fine_geom = Geom();
    const Geometry &coarse_geom = parent->Geom(lev-1);

    IntVect ngrow{0,0,0}; // not adding ghost data for now
    IntVect ratio = parent->refRatio(lev-1);

    const Box &fine_domain = fine_geom.Domain();

    // make an array of coarsened boxes corresponding to fine boxes
    BoxArray ba_crse_patch(ba.size());
    for (int i = 0, N = ba.size(); i < N; ++i)
    {
        Box bx = amrex::grow(ba[i], ngrow);
        bx &= fine_domain;

        ba_crse_patch.set(i, bx.coarsen(ratio));
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ba_crse_patch, dm); mfi.isValid(); ++mfi)
    {
        auto& sfab = coarse[mfi];       // source data
        auto& dfab = fine[mfi];         // destination data
        const Box dbox = dfab.box();    // this box is grown by 1
        const Box sbox = mfi.validbox();

        // copying particles to fine cells
        For(sbox, [&] (int i, int j, int k) 
        {
            const IntVect iv(i,j,k);
            Cell &parent = sfab(iv);
            
            for (uint i = 0; i < parent.number_of_particles; i++)
            {
                auto& particle = parent.particles[i];
                // get the new location corresponding to fine geometry
                IntVect new_iv = fine_geom.CellIndex(particle.data());

                dfab(new_iv).particles.push_back(particle);
                dfab(new_iv).number_of_particles++;
            }

            // clear for now
            parent.particles.clear();
            parent.number_of_particles = 0;
        });
    }
}*/

/*void
AmrLevel::FillPatchAdd (AmrLevel&  Amrlevel,
                        CellFabArray& leveldata,
                        int        boxGrow,
                        Real       time,
                        int        icomp,
                        int        ncomp)
{
    BL_ASSERT(boxGrow <= leveldata.nGrow());

    FillPatchIterator fpi(Amrlevel, leveldata, boxGrow, time, icomp, ncomp);
    const CellFabArray& mf_fillpatched = fpi.get_mf();
    amrex::Add(leveldata, mf_fillpatched, icomp, icomp, ncomp, boxGrow);
}*/