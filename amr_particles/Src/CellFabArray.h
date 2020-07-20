#ifndef CellFabArray_H
#define CellFabArray_H

#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Periodicity.H>
#include "Cell.h"

using namespace amrex;
using CellFab = BaseFab<Cell>;
using CellFabFactory = DefaultFabFactory<CellFab>;

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

    /**
    * \brief Swap from src to dst including nghost ghost cells.
    * The two CellFabArrays MUST have the same underlying BoxArray.
    * The swap is local.
    */
    /*static void Swap (CellFabArray& dst, CellFabArray& src,
                      int srccomp, int dstcomp, int numcomp, int nghost)
    { Swap(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost)); }
    
    static void Swap (CellFabArray& dst, CellFabArray& src,
                      int srccomp, int dstcomp, int numcomp, const IntVect& nghost);*/

    // Fills ghost cells of this CellFabArray
    void FillBoundary (bool cross = false);
    void FillBoundary (const Periodicity& period, bool cross = false);
    void FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross = false);

    //! Same as FillBoundary(), but only copies ncomp components starting at scomp.
    void FillBoundary (int scomp, int ncomp, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);

    void FillBoundary_nowait (bool cross = false);
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);
    void FillBoundary_finish ();

    // Copies from data from src to this CellFabArray (needed for AMR)
    void ParallelCopy (const CellFabArray& src,
                       int src_comp, int dest_comp, int num_comp,
                       const Periodicity& period = Periodicity::NonPeriodic(),
                       CpOp op = FabArrayBase::COPY)
    { ParallelCopy(src,src_comp,dest_comp,num_comp,0,0,period,op); }

    void ParallelCopy (const CellFabArray& src,
                        int src_comp, int dest_comp, int num_comp,
                        int src_nghost, int dst_nghost,
                        const Periodicity& period = Periodicity::NonPeriodic(),
                        CpOp op = FabArrayBase::COPY)
    { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op); }

    void ParallelCopy (const CellFabArray& src,
                        int src_comp, int dest_comp, int num_comp,
                        const IntVect& src_nghost, const IntVect& dst_nghost,
                        const Periodicity& period = Periodicity::NonPeriodic(),
                        CpOp op = FabArrayBase::COPY,
                        const FabArrayBase::CPC* a_cpc = nullptr);

    void FillPatchSingleLevel (const CellFabArray& src,
                               int scomp, int dcomp, int ncomp,
                               const Geometry& geom)
    { FillPatchSingleLevel(nGrowVect(), src, scomp, dcomp, ncomp, geom); }

    void FillPatchSingleLevel (amrex::IntVect const& nghost,
                               const CellFabArray& src,
                               int scomp, int dcomp, int ncomp,
                               const Geometry& geom);

    void FillPatchTwoLevels (const CellFabArray& coarse, const CellFabArray& fine, 
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             const IntVect& ratio);

private:
    typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

    #ifdef BL_USE_MPI
    //! Prepost nonblocking receives
    void PostReceives (const MapOfCopyComTagContainers& m_RcvTags,
                    Vector<MPI_Request>& recv_reqs,
                    int icomp, int ncomp, int SeqNum);
    
    void PostSends (const MapOfCopyComTagContainers& m_SndTags,
                    Vector<MPI_Request> send_reqs,
                    int icomp, int ncomp, int SeqNum);
    #endif

    void FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp);

    void PC_local_cpu (const CPC& thecpc, CellFabArray const& src,
                        int scomp, int dcomp, int ncomp, CpOp op);
};

inline void
CellFabArray::FillBoundary (bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), n_grow, Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), n_grow, period, cross);
	FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), nghost, period, cross);
	FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, n_grow, Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, n_grow, period, cross);
	FillBoundary_finish();
    }
}

inline void
CellFabArray::FillBoundary (int scomp, int ncomp, const IntVect& nghost,
                             const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, nghost, period, cross);
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

inline void
CellFabArray::FillPatchSingleLevel (amrex::IntVect const& nghost,
                                    const CellFabArray& src,
                                    int scomp, int dcomp, int ncomp,
                                    const Geometry& geom)
{
    BL_PROFILE("FillPatchSingleLevel");

    if (this == &src and scomp == dcomp) {
        FillBoundary(dcomp, ncomp, nghost, geom.periodicity());
    } else {
        ParallelCopy(src, scomp, dcomp, ncomp, IntVect{0}, nghost, geom.periodicity());
    }
    // Boundary conditions need to applied separately
}

inline void
CellFabArray::FillPatchTwoLevels (const CellFabArray& coarse, const CellFabArray& fine, 
                                  int scomp, int dcomp, int ncomp,
                                  const Geometry& cgeom, const Geometry& fgeom,
                                  const IntVect& ratio)
{
	BL_PROFILE("FillPatchTwoLevels");

    IntVect const& nghost = nGrowVect();
/*
	Interpolater* mapper = &cell_cons_interp;

	if (nghost.max() > 0 || getBDKey() != fine->getBDKey())
	{
	    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

	    Box fdomain = fgeom.Domain();
	    fdomain.convert(boxArray().ixType());
	    Box fdomain_g(fdomain);
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (fgeom.isPeriodic(i)) {
                fdomain_g.grow(i,nghost[i]);
            }
	    }

	    const FabArrayBase::FPinfo& fpc 
        = FabArrayBase::TheFPinfo(*fine, mf, fdomain_g, nghost, coarsener,
                                  amrex::coarsen(fgeom.Domain(),ratio), nullptr);

	    if ( ! fpc.ba_crse_patch.empty())
	    {
            CellFabArray coarse_patch = make_mf_crse_patch<CellFabArray>(fpc, ncomp);
            //amrex::mf_set_domain_bndry(coarse_patch, cgeom);

            FillPatchSingleLevel(coarse_patch, cmf, scomp, 0, ncomp, cgeom);

            int idummy1=0, idummy2=0;
            bool cc = fpc.ba_crse_patch.ixType().cellCentered();
            ignore_unused(cc);

#ifdef _OPENMP
#pragma omp parallel if (cc && Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(coarse_patch); mfi.isValid(); ++mfi)
            {
                auto& sfab = coarse_patch[mfi];
                int li = mfi.LocalIndex();
                int gi = fpc.dst_idxs[li];
                auto& dfab = mf[gi];
                const Box& dbx = fpc.dst_boxes[li] & dfab.box();

                //pre_interp(sfab, sfab.box(), 0, ncomp);

                interp(sfab, 0, dfab, dcomp, ncomp, dbx,
                        ratio, cgeom, fgeom, RunOn::Gpu);

                //post_interp(dfab, dbx, dcomp, ncomp);
            }
	    }
	}
*/
	FillPatchSingleLevel(nghost, fine, scomp, dcomp, ncomp, fgeom);
}

#endif