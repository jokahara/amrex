
#include <iostream>
#include <algorithm>

#include <unistd.h>

#include <AMReX_RealBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include "StateData.h"
#include "StateDescriptor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

StateData::StateData () 
    : desc(nullptr),
      new_time{INVALID_TIME,INVALID_TIME},
      old_time{INVALID_TIME,INVALID_TIME},
      arena(nullptr)
{
}

StateData::StateData (const Box&                 p_domain,
                      const BoxArray&            grds,
		              const DistributionMapping& dm,
                      const StateDescriptor*     d,
                      Real                       cur_time,
                      Real                       dt,
                      const FabFactory<CellFab>& factory)
{
    define(p_domain, grds, dm, *d, cur_time, dt, factory);
}

StateData::StateData (StateData&& rhs) noexcept
    : m_factory(std::move(rhs.m_factory)),
      desc(rhs.desc),
      domain(rhs.domain),
      grids(std::move(rhs.grids)),
      dmap(std::move(rhs.dmap)),
      new_time(rhs.new_time),
      old_time(rhs.old_time),
      new_data(std::move(rhs.new_data)),
      old_data(std::move(rhs.old_data)),
      arena(rhs.arena)
{   
}

void
StateData::operator= (StateData const& rhs)
{
    m_factory.reset(rhs.m_factory->clone());
    desc = rhs.desc;
    arena = rhs.arena;
    domain = rhs.domain;
    grids = rhs.grids;
    dmap = rhs.dmap;
    new_time = rhs.new_time;
    old_time = rhs.old_time;

    new_data.reset(new CellFabArray(grids,dmap,desc->nComp(),desc->nExtra(),
                                MFInfo().SetTag("StateData").SetArena(arena),
                                *m_factory));
    CellFabArray::Copy(*new_data, *rhs.new_data, 0, 0, desc->nComp(),desc->nExtra());

    if (rhs.old_data) {
        old_data.reset(new CellFabArray(grids,dmap,desc->nComp(),desc->nExtra(),
                                    MFInfo().SetTag("StateData").SetArena(arena),
                                    *m_factory));
        CellFabArray::Copy(*old_data, *rhs.old_data, 0, 0, desc->nComp(),desc->nExtra());
    } else {
        old_data.reset();
    }
}

void
StateData::define (const Box&                 p_domain,
                   const BoxArray&            grds,
		           const DistributionMapping& dm,
                   const StateDescriptor&     d,
                   Real                       time,
                   Real                       dt,                 
                   const FabFactory<CellFab>& factory)
{
    BL_PROFILE("StateData::define()");
    domain = p_domain;
    desc = &d;
    arena = nullptr;
    grids = grds;
    dmap = dm;
    m_factory.reset(factory.clone());
    //
    // Convert to proper type.
    //
    IndexType typ(desc->getType());
    StateDescriptor::TimeCenter t_typ(desc->timeType());
    if (!typ.cellCentered())
    {
        domain.convert(typ);
        grids.convert(typ);
    }
    if (t_typ == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time - dt;
    }
    else
    {
        new_time.start = time;
        new_time.stop  = time+dt;
        old_time.start = time-dt;
        old_time.stop  = time;
    }
    int ncomp = desc->nComp();

    new_data.reset(new CellFabArray(grids,dmap,ncomp,desc->nExtra(),
                                MFInfo().SetTag("StateData").SetArena(arena),
                                *m_factory));
    old_data.reset();
}

void
StateData::copyOld (const StateData& state)
{
    const CellFabArray& MF = state.oldData();
    
    int nc = MF.nComp();
    int ng = MF.nGrow();
    
    BL_ASSERT(nc == (*old_data).nComp());
    BL_ASSERT(ng == (*old_data).nGrow());
    
    CellFabArray::Copy(*old_data, state.oldData(), 0, 0, nc, ng);
    
    old_time = state.old_time;
}

void
StateData::copyNew (const StateData& state)
{
    const CellFabArray& MF = state.newData();

    int nc = MF.nComp();
    int ng = MF.nGrow();
    
    BL_ASSERT(nc == (*new_data).nComp());
    BL_ASSERT(ng == (*new_data).nGrow());
    
    CellFabArray::Copy(*new_data, state.newData(), 0, 0, nc, ng);

    new_time = state.new_time;
}

void
StateData::reset ()
{
    new_time = old_time;
    old_time.start = old_time.stop = INVALID_TIME;
    std::swap(old_data, new_data);   
}

void
StateData::restart (std::istream&              is,
                    const Box&                 p_domain,
                    const BoxArray&            grds,
                    const DistributionMapping& dm,
                    const FabFactory<CellFab>& factory,
                    const StateDescriptor&     d,
                    const std::string&         chkfile)
{
    desc = &d;
    arena = nullptr;
    domain = p_domain;
    grids = grds;
    dmap = dm;
    m_factory.reset(factory.clone());

    // Convert to proper type.
    IndexType typ(desc->getType());
    if (!typ.cellCentered()) {
        domain.convert(typ);
        grids.convert(typ);
    }

    {
	Box domain_in;
	BoxArray grids_in;
	is >> domain_in;
	grids_in.readFrom(is);
	BL_ASSERT(domain_in == domain);
	BL_ASSERT(amrex::match(grids_in,grids));
    }

    restartDoit(is, chkfile);
}

void 
StateData::restartDoit (std::istream& is, const std::string& chkfile)
{
}

void 
StateData::restart (const StateDescriptor& d,
		    const StateData& rhs)
{
    desc = &d;
    arena = nullptr;
    domain = rhs.domain;
    grids = rhs.grids;
    old_time.start = rhs.old_time.start;
    old_time.stop  = rhs.old_time.stop;
    new_time.start = rhs.new_time.start;
    new_time.stop  = rhs.new_time.stop;
    old_data.reset();
    new_data.reset(new CellFabArray(grids,dmap,desc->nComp(),desc->nExtra(),
                                MFInfo().SetTag("StateData").SetArena(arena),
                                *m_factory));
    new_data->setVal(0.);
}

StateData::~StateData()
{
    desc = nullptr;
}

void
StateData::allocOldData ()
{
    if (old_data == nullptr)
    {
        old_data.reset(new CellFabArray(grids,dmap,desc->nComp(),desc->nExtra(),
                                    MFInfo().SetTag("StateData").SetArena(arena),
                                    *m_factory));
    }
}

BCRec
StateData::getBC (int comp, int i) const noexcept
{
    BCRec bcr;
    amrex::setBC(grids[i],domain,desc->getBC(comp),bcr);
    return bcr;
}

void
StateData::setOldTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        old_time.start = old_time.stop = time;
    }
    else
    {
        amrex::Error("StateData::setOldTimeLevel called with Interval");
    }
}

void
StateData::setNewTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
    }
    else
    {
        amrex::Error("StateData::setNewTimeLevel called with Interval");
    }
}

void
StateData::syncNewTimeLevel (Real time)
{
    Real teps = (new_time.stop - old_time.stop)*1.e-3;
    if (time > new_time.stop-teps && time < new_time.stop+teps)
    {
	if (desc->timeType() == StateDescriptor::Point)
	{
	    new_time.start = new_time.stop = time;
	}
	else
	{
	    new_time.stop = time;
	}
    }
}

void
StateData::setTimeLevel (Real time,
                         Real dt_old,
                         Real dt_new)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time - dt_old;
    }
    else
    {
        new_time.start = time;
        new_time.stop  = time+dt_new;
        old_time.start = time-dt_old;
        old_time.stop  = time;
    }
}

void
StateData::swapTimeLevels (Real dt)
{
    old_time = new_time;
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start += dt;
        new_time.stop  += dt;
   }
    else
    {
        new_time.start = new_time.stop;
        new_time.stop += dt;
    }
    std::swap(old_data, new_data);
}

void
StateData::replaceOldData (CellFabArray&& mf)
{
    old_data.reset(new CellFabArray(std::move(mf)));
}

// This version does NOT delete the replaced data.

void
StateData::replaceOldData (StateData& s)
{
    CellFabArray::Swap(*old_data, *s.old_data, 0, 0, old_data->nComp(), old_data->nGrow());
}

void
StateData::replaceNewData (CellFabArray&& mf)
{
    new_data.reset(new CellFabArray(std::move(mf)));
}

// This version does NOT delete the replaced dat     a.

void
StateData::replaceNewData (StateData& s)
{
    CellFabArray::Swap(*new_data, *s.new_data, 0, 0, new_data->nComp(), new_data->nGrow());     
}
     
void
StateData::FillBoundary (CellFab&       dest,
                         Real           time,
                         const Real*    dx,
                         const RealBox& prob_domain,
                         int            dest_comp,
                         int            src_comp,
                         int            num_comp)
{
    BL_PROFILE("StateData::FillBoundary(dx)");
    BL_ASSERT(dest.box().ixType() == desc->getType());
   
    if (domain.contains(dest.box())) return;

    const Box& bx  = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();

    Vector<int> bcrs;

    Real xlo[AMREX_SPACEDIM];
    BCRec bcr;
    const Real* problo = prob_domain.lo();

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    for (int i = 0; i < num_comp; )
    {
        const int dc  = dest_comp+i;
        const int sc  = src_comp+i;
        Cell*     dat = dest.dataPtr(dc);

        if (desc->master(sc))
        {
            const int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= num_comp)
            {
                //
                // Can do the whole group at once.
                //
                bcrs.resize(2*AMREX_SPACEDIM*groupsize);

                int* bci  = bcrs.dataPtr();

                for (int j = 0; j < groupsize; j++)
                {
                    amrex::setBC(bx,domain,desc->getBC(sc+j),bcr);

                    const int* bc = bcr.vect();

                    for (int k = 0; k < 2*AMREX_SPACEDIM; k++)
                        bci[k] = bc[k];

                    bci += 2*AMREX_SPACEDIM;
                }
                //
                // Use the "group" boundary fill routine.
                //
		        desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcrs.dataPtr(),groupsize);
                i += groupsize;
            }
            else
            {
                amrex::setBC(bx,domain,desc->getBC(sc),bcr);
                desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcr.vect());
                i++;
            }
        }
        else
        {
            amrex::setBC(bx,domain,desc->getBC(sc),bcr);
            desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcr.vect());
            i++;
        }
    }

#ifdef AMREX_USE_GPU
    // Add a synchronize here in case the user code launched kernels
        // to handle the boundary fills.
        Gpu::synchronize();
#endif
}

void
StateData::FillBoundary (Box const&      bx,
                         CellFab&        dest,
                         Real            time,
                         const Geometry& geom,
                         int             dest_comp,
                         int             src_comp,
                         int             num_comp)
{
    BL_PROFILE("StateData::FillBoundary(geom)");
    BL_ASSERT(bx.ixType() == desc->getType());
   
    if (domain.contains(bx)) return;

    Vector<BCRec> bcr(num_comp);

    for (int i = 0; i < num_comp; )
    {
        const int dc  = dest_comp+i;
        const int sc  = src_comp+i;

        if (desc->master(sc))
        {
            const int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= num_comp)
            {
                for (int j = 0; j < groupsize; j++)
                {
                    amrex::setBC(bx,domain,desc->getBC(sc+j),bcr[j]);
                }
                //
                // Use the "group" boundary fill routine.
                //
		        desc->bndryFill(sc)(bx,dest,dc,groupsize,geom,time,bcr,0,sc);
                i += groupsize;
            }
            else
            {
                amrex::setBC(bx,domain,desc->getBC(sc),bcr[0]);
                desc->bndryFill(sc)(bx,dest,dc,1,geom,time,bcr,0,sc);
                i++;
            }
        }
        else
        {
            amrex::setBC(bx,domain,desc->getBC(sc),bcr[0]);
            desc->bndryFill(sc)(bx,dest,dc,1,geom,time,bcr,0,sc);
            i++;
        }
    }
}

void
StateData::RegisterData (MultiFabCopyDescriptor&  CellFabArrayCopyDesc,
                         Vector<CellFabId>&      mfid)
{
    mfid.resize(2);
    mfid[MFNEWDATA] = CellFabArrayCopyDesc.RegisterFabArray(new_data.get());
    mfid[MFOLDDATA] = CellFabArrayCopyDesc.RegisterFabArray(old_data.get());
}

void
StateData::InterpAddBox (MultiFabCopyDescriptor& CellFabArrayCopyDesc,
                         Vector<CellFabId>&     mfid,
                         BoxList*                unfillableBoxes,
                         Vector<FillBoxId>&      returnedFillBoxIds,
                         const Box&              subbox,
                         Real                    time,
                         int                     src_comp,
                         int                     dest_comp,
                         int                     num_comp,
                         bool                    extrap)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == nullptr)
        {
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = CellFabArrayCopyDesc.AddBox(mfid[MFNEWDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
        }
        else
        {
            amrex::InterpAddBox(CellFabArrayCopyDesc,
                                unfillableBoxes,
                                returnedFillBoxIds,
                                subbox,
                                mfid[MFOLDDATA],
                                mfid[MFNEWDATA],
                                old_time.start,
                                new_time.start,
                                time,
                                src_comp,
                                dest_comp,
                                num_comp,
                                extrap);
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = CellFabArrayCopyDesc.AddBox(mfid[MFNEWDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
        }
        else if (old_data != nullptr        &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = CellFabArrayCopyDesc.AddBox(mfid[MFOLDDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
        }
        else
        {
            amrex::Error("StateData::Interp(): cannot interp");
            }
   }     
}

void
StateData::InterpFillFab (MultiFabCopyDescriptor&  CellFabArrayCopyDesc,
                          const Vector<CellFabId>& mfid,
                          const Vector<FillBoxId>& fillBoxIds,
                          CellFab&                 dest,
                          Real                     time,
                          int                      src_comp,
                          int                      dest_comp,
                          int                      num_comp,
                          bool                     extrap)
            {
    BL_PROFILE("StateData::InterpFillFab()");
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == nullptr)
        {
            CellFabArrayCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
        }
        else
        {
            amrex::InterpFillFab(CellFabArrayCopyDesc,
                                fillBoxIds,
                                mfid[MFOLDDATA],
                                mfid[MFNEWDATA],
                                dest,
                                old_time.start,
                                new_time.start,
                                time,
                                src_comp,
                                dest_comp,
                                num_comp,
                                extrap);
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
            CellFabArrayCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
        }
        else if (old_data != nullptr        &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
            CellFabArrayCopyDesc.FillFab(mfid[MFOLDDATA], fillBoxIds[0], dest);
        }
        else
        {
            amrex::Error("StateData::Interp(): cannot interp");
        }
    }
}

void
StateData::getData (Vector<CellFabArray*>& data,
		    Vector<Real>& datatime,
		    Real time) const
{
    data.clear();
    datatime.clear();

    if (desc->timeType() == StateDescriptor::Point)
    {
	BL_ASSERT(new_data != nullptr);
        if (old_data == nullptr)
        {
	    data.push_back(new_data.get());
	    datatime.push_back(new_time.start);
        }
        else
        {
	    const Real teps = (new_time.start - old_time.start)*1.e-3;
	    if (time > new_time.start-teps && time < new_time.start+teps) {
		data.push_back(new_data.get());
		datatime.push_back(new_time.start);
	    } else if (time > old_time.start-teps && time < old_time.start+teps) {
	    	    data.push_back(old_data.get());
		    datatime.push_back(old_time.start);
	    } else {
		data.push_back(old_data.get());
		data.push_back(new_data.get());
		datatime.push_back(old_time.start);
		datatime.push_back(new_time.start);
	    }
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
	    data.push_back(new_data.get());
	    datatime.push_back(time);
        }
        else if (old_data != nullptr        &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
	    data.push_back(old_data.get());
	    datatime.push_back(time);
        }
        else
        {
            amrex::Error("StateData::getData(): how did we get here?");
        }
    }
}

void
StateData::checkPoint (const std::string& name,
                       const std::string& fullpathname,
                       std::ostream&  os,
                       bool           dump_old)
{
}

void
StateData::printTimeInterval (std::ostream &os) const
{
    os << '['
       << old_time.start
       << ' '
       << old_time.stop
       << "] ["
       << new_time.start
       << ' '
       << new_time.stop
       << ']'
       << '\n';
}

StateDataPhysBCFunct::StateDataPhysBCFunct (StateData&sd, int sc, const Geometry& geom_)
    : statedata(&sd),
      src_comp(sc),
      geom(geom_)
{ }

void
StateDataPhysBCFunct::operator() (CellFabArray& mf, int dest_comp, int num_comp, IntVect const& /* */,
                                  Real time, int /*bccomp*/)
{
    BL_PROFILE("StateDataPhysBCFunct::()");

    const Box&     domain      = statedata->getDomain();
    const int*     domainlo    = domain.loVect();
    const int*     domainhi    = domain.hiVect();
    const Real*    dx          = geom.CellSize();
    const RealBox& prob_domain = geom.ProbDomain();

    bool has_bndryfunc_fab = statedata->desc->hasBndryFuncFab();
    bool run_on_gpu = statedata->desc->RunOnGPU() && Gpu::inLaunchRegion();

#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel if (!run_on_gpu)
#endif
#endif
    {
	CellFab tmp;

	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
        CellFab& dest = mf[mfi];
        Array4<Cell> const& desta = dest.array();
        const Box& bx = dest.box();
    
	    bool has_phys_bc = false;
	    bool is_periodic = false;
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            bool touch = bx.smallEnd(i) < domainlo[i] || bx.bigEnd(i) > domainhi[i];
            if (geom.isPeriodic(i)) {
                is_periodic = is_periodic || touch;
            } else {
                has_phys_bc = has_phys_bc || touch;
            }
	    }

	    if (has_phys_bc)
	    {
            if (has_bndryfunc_fab) {
                statedata->FillBoundary(bx, dest, time, geom, dest_comp, src_comp, num_comp);
            } else {
                statedata->FillBoundary(dest, time, dx, prob_domain, dest_comp, src_comp, num_comp);
            }
		
		if (is_periodic) // fix up corner
		{
		    Box GrownDomain = domain;
		    
		    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
		    {
                if (!geom.isPeriodic(dir))
                {
                    const int lo = domainlo[dir] - bx.smallEnd(dir);
                    const int hi = bx.bigEnd(dir) - domainhi[dir];
                    if (lo > 0) GrownDomain.growLo(dir,lo);
                    if (hi > 0) GrownDomain.growHi(dir,hi);
                }
		    }
		    
		    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
		    {
                if (!geom.isPeriodic(dir)) continue;
                
                Box lo_slab = bx;
                Box hi_slab = bx;
                lo_slab.shift(dir, domain.length(dir));
                hi_slab.shift(dir,-domain.length(dir));
                lo_slab &= GrownDomain;
                hi_slab &= GrownDomain;
                
                if (lo_slab.ok())
                {
                    if (run_on_gpu)
                    {
                        tmp.resize(lo_slab,num_comp);
                        Array4<Cell> const& tmpa = tmp.array();
                        const int ishift = -domain.length(dir);
                        amrex::launch(lo_slab,
                        [=] AMREX_GPU_DEVICE (Box const& tbx) noexcept
                        {
                            const Box db = amrex::shift(tbx, dir, ishift);
                            const auto dlo = amrex::lbound(db);
                            const auto tlo = amrex::lbound(tbx);
                            const auto len = amrex::length(db);
                            for (int n = 0; n < num_comp; ++n) {
                                for         (int k = 0; k < len.z; ++k) {
                                    for     (int j = 0; j < len.y; ++j) {
                                        AMREX_PRAGMA_SIMD
                                        for (int i = 0; i < len.x; ++i) {
                                            tmpa(i+tlo.x,j+tlo.y,k+tlo.z,n)
                                                = desta(i+dlo.x,j+dlo.y,k+dlo.z,n+dest_comp);
                                        }
                                    }
                                }
                            }
                        });
                        if (has_bndryfunc_fab) {
                            statedata->FillBoundary(lo_slab, tmp, time, geom, 0, src_comp, num_comp);
                        } else {
                            statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                        }
                        amrex::launch(lo_slab,
                        [=] AMREX_GPU_DEVICE (Box const& tbx) noexcept
                        {
                            const Box db = amrex::shift(tbx, dir, ishift);
                            const auto dlo = amrex::lbound(db);
                            const auto tlo = amrex::lbound(tbx);
                            const auto len = amrex::length(db);
                            for (int n = 0; n < num_comp; ++n) {
                                for (int k = 0; k < len.z; ++k) {
                                    for (int j = 0; j < len.y; ++j) {
                                        AMREX_PRAGMA_SIMD
                                        for (int i = 0; i < len.x; ++i) {
                                            desta(i+dlo.x,j+dlo.y,k+dlo.z,n+dest_comp)
                                                = tmpa(i+tlo.x,j+tlo.y,k+tlo.z,n);
                                        }
                                    }
                                }
                            }
                        });
                    }
                    else
                    {
                        tmp.resize(lo_slab,num_comp);
                        const Box db = amrex::shift(lo_slab, dir, -domain.length(dir));
                        tmp.copy<RunOn::Host>(dest, db, dest_comp, lo_slab, 0, num_comp);
                        if (has_bndryfunc_fab) {
                            statedata->FillBoundary(lo_slab, tmp, time, geom, 0, src_comp, num_comp);
                        } else {
                            statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                        }
                        dest.copy<RunOn::Host>(tmp, lo_slab, 0, db, dest_comp, num_comp);
                    }
                }
                
                if (hi_slab.ok())
                {
                    if (run_on_gpu)
                    {
                        tmp.resize(hi_slab,num_comp);
                        Array4<Cell> const& tmpa = tmp.array();
                        const int ishift = domain.length(dir);
                        amrex::launch(hi_slab,
                        [=] AMREX_GPU_DEVICE (Box const& tbx) noexcept
                        {
                            const Box db = amrex::shift(tbx, dir, ishift);
                            const auto dlo = amrex::lbound(db);
                            const auto tlo = amrex::lbound(tbx);
                            const auto len = amrex::length(db);
                            for (int n = 0; n < num_comp; ++n) {
                                for         (int k = 0; k < len.z; ++k) {
                                    for     (int j = 0; j < len.y; ++j) {
                                        AMREX_PRAGMA_SIMD
                                        for (int i = 0; i < len.x; ++i) {
                                            tmpa(i+tlo.x,j+tlo.y,k+tlo.z,n)
                                                = desta(i+dlo.x,j+dlo.y,k+dlo.z,n+dest_comp);
                                        }
                                    }
                                }
                            }
                        });
                        if (has_bndryfunc_fab) {
                            statedata->FillBoundary(hi_slab, tmp, time, geom, 0, src_comp, num_comp);
                        } else {
                            statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                        }
                        amrex::launch(hi_slab,
                        [=] AMREX_GPU_DEVICE (Box const& tbx) noexcept
                        {
                            const Box db = amrex::shift(tbx, dir, ishift);
                            const auto dlo = amrex::lbound(db);
                            const auto tlo = amrex::lbound(tbx);
                            const auto len = amrex::length(db);
                            for (int n = 0; n < num_comp; ++n) {
                                for         (int k = 0; k < len.z; ++k) {
                                    for     (int j = 0; j < len.y; ++j) {
                                        AMREX_PRAGMA_SIMD
                                        for (int i = 0; i < len.x; ++i) {
                                            desta(i+dlo.x,j+dlo.y,k+dlo.z,n+dest_comp)
                                                = tmpa(i+tlo.x,j+tlo.y,k+tlo.z,n);
                                        }
                                    }
                                }
                            }
                        });
                    }
                    else
                    {
                        tmp.resize(hi_slab,num_comp);
                        const Box db = amrex::shift(hi_slab, dir, domain.length(dir));
                        tmp.copy<RunOn::Host>(dest, db, dest_comp, hi_slab, 0, num_comp);
                        if (has_bndryfunc_fab) {
                            statedata->FillBoundary(hi_slab, tmp, time, geom, 0, src_comp, num_comp);
                        } else {
                            statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                        }
                        dest.copy<RunOn::Host>(tmp, hi_slab, 0, db, dest_comp, num_comp);
                    }
                }
		    }
		}
	    }
	}
    }
}
