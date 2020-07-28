#ifndef CellFabArray_H
#define CellFabArray_H

#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Interpolater.H>
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
    
    void EstimateWork (MultiFab& mf);

    /** 
     * Returns a map of tags pointing to boxes which receive data 
     * from remote neighbours during FillBoundary
    */
    const MapOfCopyComTagContainers& get_receive_tags() {
        const FB& theFB = *m_TheFBCache.find(getBDKey())->second;
        return *theFB.m_RcvTags;
    }

    const MapOfCopyComTagContainers& get_send_tags() {
        const FB& theFB = *m_TheFBCache.find(getBDKey())->second;
        return *theFB.m_SndTags;
    }

    const CopyComTagsContainer& get_local_tags() {
        const FB& theFB = *m_TheFBCache.find(getBDKey())->second;
        return *theFB.m_LocTags;
    }

    // Fills ghost cells of this CellFabArray
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

    // Copies from data from src to this CellFabArray
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

    inline void FillPatchSingleLevel (CellFabArray& src,
                               int scomp, int dcomp, int ncomp,
                               const Geometry& geom)
    { FillPatchSingleLevel(nGrowVect(), src, scomp, dcomp, ncomp, geom); }

    void FillPatchSingleLevel (amrex::IntVect const& nghost,
                               CellFabArray& src,
                               int scomp, int dcomp, int ncomp,
                               const Geometry& geom);

    void FillPatchTwoLevels (CellFabArray& coarse, CellFabArray& fine, 
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             const IntVect& ratio);

    // Used for interpolating coarse data to finer level.
    void interpolateTo(CellFab& dfab, const Box& dbox, const Geometry& fgeom,
                       const CellFab& sfab, const Box& sbox, const Geometry& cgeom) 
    {
        // copying particles to fine cells
        ParallelFor(sbox, n_comp, [&] (int i, int j, int k, int n) 
        {
            const IntVect iv(i,j,k);
            const Cell &parent = sfab(iv, n);
            
            for (uint i = 0; i < parent.number_of_particles; i++)
            {
                auto& particle = parent.particles[i];
                // get the new location corresponding to fine geometry
                IntVect new_iv = fgeom.CellIndex(particle.data());

                // add particles to new locations
                dfab(new_iv).particles.push_back(particle);
                dfab(new_iv).number_of_particles++;
            }
        });
    }

private:
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
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
                        
    // Swaps vector pointers, instead of copying. Much faster than PC_local_cpu.
    void PC_local_swap_cpu (const CPC& thecpc, CellFabArray& src,
                            int scomp, int dcomp, int ncomp);

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
            mfab(idx) = sfab(idx).number_of_particles;
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
CellFabArray::FillPatchSingleLevel (amrex::IntVect const& nghost,
                                    CellFabArray& src,
                                    int scomp, int dcomp, int ncomp,
                                    const Geometry& geom)
{
    BL_PROFILE("FillPatchSingleLevel");

    Periodicity period = geom.periodicity();
    
    if (this != &src or scomp != dcomp) {
        // tags for cells which are moved / copied
        // ghost cells are ignored for now
        const CPC& cpc = getCPC(nghost, src, nghost, period);   

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
                    dfab(idx, n+dcomp).resize();
                });
            }
        }

        Cell::transfer_particles = true;
        ParallelCopy(src, scomp, dcomp, ncomp, period, FabArrayBase::COPY, &cpc);
    }
    
    // Finishing by filling ghost cells
    /*Cell::transfer_particles = false;
    FillBoundary(dcomp, ncomp, n_grow, period);

    for (auto const& kv : get_receive_tags())
    {
        for (auto const& tag : kv.second)
        {
            const auto& bx = tag.dbox;
            auto dfab = this->array(tag.dstIndex);
            amrex::Loop( bx, ncomp,
            [&] (int ii, int jj, int kk, int n) noexcept
            {
                dfab(ii, jj, kk, n+dcomp).resize();
            });
        }
    }

    Cell::transfer_particles = true;
    FillBoundary(dcomp, ncomp, n_grow, period);*/
}

// Fill this with data from coarse and fine
inline void
CellFabArray::FillPatchTwoLevels (CellFabArray& coarse, CellFabArray& fine, 
                                  int scomp, int dcomp, int ncomp,
                                  const Geometry& cgeom, const Geometry& fgeom,
                                  const IntVect& ratio)
{
	BL_PROFILE("FillPatchTwoLevels");

    IntVect const& ngrow = nGrowVect(); // only grow boundaries?
	Interpolater* mapper = &cell_cons_interp;

	if (ngrow.max() > 0 || getBDKey() != fine.getBDKey())
	{
	    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

	    Box fdomain = fgeom.Domain();
	    fdomain.convert(boxArray().ixType());
	    Box fdomain_g(fdomain);
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (fgeom.isPeriodic(i)) {
                fdomain_g.grow(i, ngrow[i]);
            }
	    }

        // Fill patch info (for ghost cells)
	    const FabArrayBase::FPinfo& fpc 
        = FabArrayBase::TheFPinfo(fine, *this, fdomain_g, ngrow, coarsener,
                                  amrex::coarsen(fgeom.Domain(),ratio), nullptr);
                                  
	    if ( ! fpc.ba_crse_patch.empty())
	    {
            //CellFabArray coarse_patch(fpc.ba_crse_patch, fpc.dm_crse_patch, ncomp, 0);
            //coarse_patch.FillPatchSingleLevel(coarse, scomp, 0, ncomp, cgeom);
            
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(fpc.ba_crse_patch, fpc.dm_crse_patch); mfi.isValid(); ++mfi)
            {
                auto& sfab = coarse[mfi];
                int li = mfi.LocalIndex();
                int gi = fpc.dst_idxs[li];
                auto& dfab = (*this)[gi];
                const Box& dbx = fpc.dst_boxes[li] & dfab.box();
                Print() << "interp: " << sfab.box() << " -> " << dbx << "\n";

                //interpolate to fine box

                //interp(sfab, 0, dfab, dcomp, ncomp, dbx, ratio, cgeom, fgeom, RunOn::Gpu);
            }
	    }
	}

    // fill from fine data
	FillPatchSingleLevel(ngrow, fine, scomp, dcomp, ncomp, fgeom);
}

#endif