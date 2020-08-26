#ifndef CellFabArray_H
#define CellFabArray_H

#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>

#include "Interpolater.h"
#include "Cell.h"

using namespace amrex;
using CellFab = BaseFab< std::shared_ptr<Cell> >;
using CellFabFactory = DefaultFabFactory<CellFab>;

// Global interpolater 
extern MyInterpolater my_interpolater;

class CellFabArray : public FabArray<CellFab>
{
public:
    
    /**
    * \brief Constructs an empty CellFabArray.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    CellFabArray () noexcept {};

    /**
    * \brief Construct a FabArray<CellFab> with a valid region defined by bxs
    * and a region of definition defined by the grow factor ngrow
    * and the number of components nvar.
    */
    CellFabArray (const BoxArray&      bxs,
            const DistributionMapping& dm,
            int                        ncomp,
            int                        ngrow,
	        const MFInfo&              info = MFInfo(),
            const FabFactory<CellFab>& factory = DefaultFabFactory<CellFab>())
    : CellFabArray(bxs,dm,ncomp,IntVect(ngrow),info,factory) {}

    CellFabArray (const BoxArray&      bxs,
            const DistributionMapping& dm,
            int                        ncomp,
            const IntVect&             ngrow,
	        const MFInfo&              info = MFInfo(),
            const FabFactory<CellFab>& factory = DefaultFabFactory<CellFab>())
    : FabArray<CellFab>(bxs,dm,ncomp,ngrow,info,factory) {}

    CellFabArray (const CellFabArray& rhs, MakeType maketype, int scomp, int ncomp);

    //! The destructor -- deletes all FABs in the array.
    ~CellFabArray () {};

    bool is_local(IntVect loc) {
        for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
            if (mfi.validbox().contains(loc))
                return true;
        }
        return false;
    }

    CellFabArray (CellFabArray&& rhs) noexcept
    : FabArray<CellFab>(std::move(rhs)) {}
    CellFabArray& operator= (CellFabArray&& rhs) noexcept;

    CellFabArray (const CellFabArray& rhs) = delete;
    CellFabArray& operator= (const CellFabArray& rhs) = delete;
    
    void EstimateWork (MultiFab& mf);

    /** 
     * Returns a map of tags pointing to boxes which receive data 
     * from remote neighbours during FillBoundary
    */
    const MapOfCopyComTagContainers& get_receive_tags() {
        return *fill_boundary_tags->m_RcvTags;
    }

    const MapOfCopyComTagContainers& get_send_tags() {
        return *fill_boundary_tags->m_SndTags;
    }

    const CopyComTagsContainer& get_local_tags() {
        return *fill_boundary_tags->m_LocTags;
    }

    /* Fills in the boundary regions of each FAB in the FabArray.  
    * If cross=true, corner cells are not filled.
    * If the length of periodic is provided, periodic boundaries are
    * also filled.  Note that FabArray itself does not contains
    * any periodicity information. FillBoundary expects that its 
    * cell-centered version of its BoxArray is non-overlapping.
    */
    void FillBoundary (bool cross = false);
    void FillBoundary (const Periodicity& period, bool cross = false);
    void FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross = false);

    //! Same as FillBoundary(), but only copies ncomp components starting at scomp.
    void FillBoundary (int scomp, int ncomp, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, 
                       bool cross=false, bool enforce_periodicity_only=false);

    void FillBoundary_nowait (bool cross = false);
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, 
                              bool cross=false, bool enforce_periodicity_only=false);
    void FillBoundary_finish ();

    /**
    * This function copies data from src to this FabArray. Each FAB
    * in src is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the valid regions.
    */
    inline void ParallelCopy (CellFabArray& src,
                        int src_comp, int dest_comp, int num_comp,
                        const Periodicity& period = Periodicity::NonPeriodic(),
                        CpOp op = FabArrayBase::COPY,
                        const FabArrayBase::CPC* a_cpc = nullptr)
    { ParallelCopy(src,src_comp,dest_comp,num_comp,0,0,period,op,a_cpc); }

    inline void ParallelCopy (CellFabArray& src,
                        int src_comp, int dest_comp, int num_comp,
                        int src_nghost, int dst_nghost,
                        const Periodicity& period = Periodicity::NonPeriodic(),
                        CpOp op = FabArrayBase::COPY,
                        const FabArrayBase::CPC* a_cpc = nullptr)
    { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op,a_cpc); }

    void ParallelCopy (CellFabArray& src,
                        int src_comp, int dest_comp, int num_comp,
                        const IntVect& src_nghost, const IntVect& dst_nghost,
                        const Periodicity& period = Periodicity::NonPeriodic(),
                        CpOp op = FabArrayBase::COPY,
                        const FabArrayBase::CPC* a_cpc = nullptr);

    /*void ParallelCopy_local (const CPC& thecpc, CellFabArray const& src,
                            int scomp, int dcomp, int ncomp, CpOp op) 
    { PC_local_cpu(thecpc, src, scomp, dcomp, ncomp, op); }*/

    inline void FillPatchSingleLevel (CellFabArray& src,
                               int scomp, int dcomp, int ncomp,
                               const Geometry& geom);


    void FillPatchTwoLevels (CellFabArray& coarse, CellFabArray& fine, 
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             const IntVect& ratio);

    /*void FillCoarseToFine (CellFabArray& coarse,
                           int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           const IntVect& ratio);

    void FillFineToCoarse (CellFabArray& fine,
                           int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           const IntVect& ratio);*/

private:
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

    const FB* fill_boundary_tags;

    #ifdef BL_USE_MPI
    //! Prepost nonblocking receives
    void PostReceives (const MapOfCopyComTagContainers& m_RcvTags,
                    Vector<MPI_Request>& recv_reqs,
                    int icomp, int ncomp, int SeqNum);
    
    void PostSends (const MapOfCopyComTagContainers& m_SndTags,
                    Vector<MPI_Request>& send_reqs,
                    int icomp, int ncomp, int SeqNum);
    #endif

    void FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp);

    void PC_local_cpu (const CPC& thecpc, CellFabArray const& src,
                        int scomp, int dcomp, int ncomp, CpOp op);
                        
};


inline void
CellFabArray::EstimateWork (MultiFab& mf) 
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto& mfab = mf[mfi];
        auto& sfab = get(mfi);
        auto box = mfi.validbox();
        
        ParallelFor(box, [&] (int i, int j, int k)
        {
            IntVect idx{i,j,k};
            mfab(idx) = sfab(idx)->number_of_particles;
        });
    }
}

inline void
CellFabArray::FillBoundary (bool cross)
{
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait(0, nComp(), n_grow, Periodicity::NonPeriodic(), cross);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (const Periodicity& period, bool cross)
{
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait(0, nComp(), n_grow, period, cross);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross)
{
    if ( nghost.max() > 0 ) {
        FillBoundary_nowait(0, nComp(), nghost, period, cross);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, bool cross)
{
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait(scomp, ncomp, n_grow, Periodicity::NonPeriodic(), cross);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait(scomp, ncomp, n_grow, period, cross);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, const IntVect& nghost,
                            const Periodicity& period, bool cross,
                            bool enforce_periodicity_only)
{
    if ( nghost.max() > 0 ) {
        FillBoundary_nowait(scomp, ncomp, nghost, period, cross, enforce_periodicity_only);
        FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary_nowait (bool cross)
{
    FillBoundary_nowait(0, nComp(), nGrowVect(), Periodicity::NonPeriodic(), cross);
}

inline void
CellFabArray::FillBoundary_nowait (const Periodicity& period, bool cross)
{
    FillBoundary_nowait(0, nComp(), nGrowVect(), period, cross);
}

inline void
CellFabArray::FillBoundary_nowait (int scomp, int ncomp, bool cross)
{
    FillBoundary_nowait(scomp, ncomp, nGrowVect(), Periodicity::NonPeriodic(), cross);
}

inline void
CellFabArray::FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    FillBoundary_nowait(scomp, ncomp, nGrowVect(), period, cross);
}

// Fill this with data from src
inline void
CellFabArray::FillPatchSingleLevel (CellFabArray& src,
                                    int scomp, int dcomp, int ncomp,
                                    const Geometry& geom)
{
    BL_PROFILE("FillPatchSingleLevel");
    
    if (this != &src or scomp != dcomp) 
    {
        Periodicity period = geom.periodicity();
        // create tags for cells which are copied
        // ghost cells are not included
        const CPC& cpc = getCPC(IntVect::Zero, src, IntVect::Zero, period);   

        if (ParallelDescriptor::NProcs() > 1)
        {
            Cell::transfer_particles = false;
            ParallelCopy(src, scomp, dcomp, ncomp, period, FabArrayBase::COPY, &cpc);

            // resizing receiving cells
            for (auto const& kv : *cpc.m_RcvTags)
            {
                for (auto const& tag : kv.second)
                {
                    const auto& bx = tag.dbox;
                    auto dfab = this->array(tag.dstIndex);
                    amrex::Loop( bx, ncomp,
                    [&] (int ii, int jj, int kk, int n) noexcept
                    {
                        const IntVect idx{ii,jj,kk};
                        dfab(idx, n+dcomp)->resize();
                    });
                }
            }
        }
        
        Cell::transfer_particles = true;
        ParallelCopy(src, scomp, dcomp, ncomp, period, FabArrayBase::COPY, &cpc);
    }
}

#endif