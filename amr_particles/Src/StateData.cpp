
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
    : new_time{INVALID_TIME,INVALID_TIME},
      old_time{INVALID_TIME,INVALID_TIME}
{
}

StateData::StateData (const Box&                 p_domain,
                      const BoxArray&            grds,
		              const DistributionMapping& dm,
                      Real                       cur_time,
                      Real                       dt,
                      const FabFactory<CellFab>& factory)
{
    define(p_domain, grds, dm, cur_time, dt, factory);
}

StateData::StateData (StateData&& rhs) noexcept
    : m_factory(std::move(rhs.m_factory)),
      domain(rhs.domain),
      grids(std::move(rhs.grids)),
      dmap(std::move(rhs.dmap)),
      new_time(rhs.new_time),
      old_time(rhs.old_time),
      new_data(std::move(rhs.new_data)),
      old_data(std::move(rhs.old_data))
{   
}

void
StateData::operator= (StateData const& rhs)
{
    m_factory.reset(rhs.m_factory->clone());
    domain = rhs.domain;
    grids = rhs.grids;
    dmap = rhs.dmap;
    new_time = rhs.new_time;
    old_time = rhs.old_time;

    new_data.reset(new CellFabArray(grids,dmap,rhs.new_data->nComp(),rhs.new_data->nGrow(),
                                    MFInfo(), *m_factory));
    amrex::Copy(*new_data, *rhs.new_data, 0, 0, rhs.new_data->nComp(),rhs.new_data->nGrow());

    if (rhs.old_data) {
        old_data.reset(new CellFabArray(grids,dmap,rhs.old_data->nComp(),rhs.old_data->nGrow(),
                                        MFInfo(), *m_factory));
        amrex::Copy(*old_data, *rhs.old_data, 0, 0, rhs.old_data->nComp(),rhs.old_data->nGrow());
    } else {
        old_data.reset();
    }
}

void
StateData::define (const Box&                 p_domain,
                   const BoxArray&            grds,
		           const DistributionMapping& dm,
                   Real                       time,
                   Real                       dt,                 
                   const FabFactory<CellFab>& factory)
{
    BL_PROFILE("StateData::define()");
    domain = p_domain;
    grids = grds;
    dmap = dm;
    m_factory.reset(factory.clone());
    
    new_time.start = time;
    new_time.stop  = time+dt;
    old_time.start = time-dt;
    old_time.stop  = time;

    new_data.reset(new CellFabArray(grids,dmap,ncomp,desc->nExtra(),
                                    MFInfo(), *m_factory));
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
    
    amrex::Copy(*old_data, state.oldData(), 0, 0, nc, ng);
    
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
    
    amrex::Copy(*new_data, state.newData(), 0, 0, nc, ng);

    new_time = state.new_time;
}

void
StateData::reset ()
{
    new_time = old_time;
    old_time.start = old_time.stop = INVALID_TIME;
    std::swap(old_data, new_data);   
}

StateData::~StateData()
{
}

void
StateData::allocOldData ()
{
    if (old_data == nullptr)
    {
        old_data.reset(new CellFabArray(grids, dmap, new_data->nComp(), new_data->nGrow(),
                                        MFInfo(), *m_factory));
    }
}

BCRec
StateData::getBC (int comp, int i) const noexcept
{
    BCRec bcr;
    amrex::setBC(grids[i],domain,bc[comp],bcr);
    return bcr;
}

void
StateData::setOldTimeLevel (Real time)
{
    amrex::Error("StateData::setOldTimeLevel called with Interval");
}

void
StateData::setNewTimeLevel (Real time)
{
    amrex::Error("StateData::setNewTimeLevel called with Interval");
}

void
StateData::syncNewTimeLevel (Real time)
{
    Real teps = (new_time.stop - old_time.stop)*1.e-3;
    if (time > new_time.stop-teps && time < new_time.stop+teps)
    {
	    new_time.stop = time;
    }
}

void
StateData::setTimeLevel (Real time,
                         Real dt_old,
                         Real dt_new)
{
    new_time.start = time;
    new_time.stop  = time+dt_new;
    old_time.start = time-dt_old;
    old_time.stop  = time;
}

void
StateData::swapTimeLevels (Real dt)
{
    old_time = new_time;
    
    new_time.start = new_time.stop;
    new_time.stop += dt;
    
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

// This version does NOT delete the replaced data.
void
StateData::replaceNewData (StateData& s)
{
    CellFabArray::Swap(*new_data, *s.new_data, 0, 0, new_data->nComp(), new_data->nGrow());     
}
     
void
StateData::FillBoundary ()
{
    new_data->FillBoundary();
    old_data->FillBoundary();
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
    new_data->FillBoundary();
    old_data->FillBoundary();
}

void
StateData::RegisterData (CopyDescriptor&     CellCopyDesc,
                         Vector<FabArrayId>& mfid)
{
    mfid.resize(2);
    mfid[MFNEWDATA] = CellCopyDesc.RegisterFabArray(new_data.get());
    mfid[MFOLDDATA] = CellCopyDesc.RegisterFabArray(old_data.get());
}

void
StateData::InterpAddBox (CopyDescriptor&        CellCopyDesc,
                         Vector<FabArrayId>&    mfid,
                         BoxList*               unfillableBoxes,
                         Vector<FillBoxId>&     returnedFillBoxIds,
                         const Box&             subbox,
                         Real                   time,
                         int                    src_comp,
                         int                    dest_comp,
                         int                    num_comp,
                         bool                   extrap)
{
    const Real teps = (new_time.start - old_time.start)*1.e-3;

    returnedFillBoxIds.resize(1);

    if (time > new_time.start-teps && time < new_time.stop+teps)
    {
        returnedFillBoxIds[0] = CellCopyDesc.AddBox(mfid[MFNEWDATA],
                                    subbox, unfillableBoxes, src_comp, dest_comp, num_comp);
    }
    else if (old_data != nullptr        &&
                time > old_time.start-teps &&
                time < old_time.stop+teps)
    {
        returnedFillBoxIds[0] = CellCopyDesc.AddBox(mfid[MFOLDDATA],
                                    subbox, unfillableBoxes, src_comp, dest_comp, num_comp);
    }
    else
    {
        amrex::Error("StateData::Interp(): cannot interp");
    }
}

void
StateData::InterpFillFab (CopyDescriptor&           CellCopyDesc,
                          const Vector<FabArrayId>& mfid,
                          const Vector<FillBoxId>&  fillBoxIds,
                          CellFab&                  dest,
                          Real                      time,
                          int                       src_comp,
                          int                       dest_comp,
                          int                       num_comp,
                          bool                      extrap)
            {
    BL_PROFILE("StateData::InterpFillFab()");

    const Real teps = (new_time.start - old_time.start)*1.e-3;
    if (time > new_time.start-teps && time < new_time.stop+teps)
    {
        CellCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
    }
    else if (old_data != nullptr        &&
                time > old_time.start-teps &&
                time < old_time.stop+teps)
    {
        CellCopyDesc.FillFab(mfid[MFOLDDATA], fillBoxIds[0], dest);
    }
    else
    {
        amrex::Error("StateData::Interp(): cannot interp");
    }
}

void
StateData::getData (Vector<CellFabArray*>& data,
		    Vector<Real>& datatime,
		    Real time) const
{
    data.clear();
    datatime.clear();

    const Real teps = (new_time.start - old_time.start)*1.e-3;

    if (time > new_time.start-teps && time < new_time.stop+teps)
    {
        data.push_back(new_data.get());
        datatime.push_back(time);
    }
    else if (old_data != nullptr &&
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
