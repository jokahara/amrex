
#include <AMReX_FabArray.H>
#include <AMReX_Periodicity.H>
#include "Cell.h"

using namespace amrex;
/**
 * \brief
 * A collection (stored as an array) of FArrayBox objects.
 * This class is useful for storing floating point data on a domain defined by
 * a union of rectangular regions embedded in a uniform index space.y
 * CellFabArray class extends the function of the underlying FabArray class just
 * as the FArrayBox class extends the funtion of BaseFab<Real>.
 * Additional member functions are defined for I/O and simple arithmetic operations on these aggregate objects.
 * This class does NOT provide a copy constructor or assignment operator.
 */

class CellFabArray : public FabArray< BaseFab<Cell> >
{
public:
    typedef BaseFab<Cell> CellFab;
    /**
    * \brief Constructs an empty CellFabArray.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    CellFabArray () noexcept {};

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

    #ifdef BL_USE_MPI
    //! Prepost nonblocking receives
    void PostReceives (const MapOfCopyComTagContainers& m_RcvTags,
                    Vector<MPI_Request>& recv_reqs,
                    int icomp, int ncomp, int SeqNum);
    
    void PostSends (const MapOfCopyComTagContainers& m_SndTags,
                    Vector<MPI_Request> send_reqs,
                    int icomp, int ncomp, int SeqNum);
    #endif
    

private:
    typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

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

void
CellFabArray::FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost,
                                    const Periodicity& period, bool cross 
                                    /*,bool enforce_periodicity_only*/)
{
    bool enforce_periodicity_only = false;

    fb_cross = cross;
    fb_epo   = enforce_periodicity_only;
    fb_scomp = scomp;
    fb_ncomp = ncomp;
    fb_nghost = nghost;
    fb_period = period;

    fb_recv_reqs.clear();

    bool work_to_do;
    if (enforce_periodicity_only) {
	    work_to_do = period.isAnyPeriodic();
    } else {
	    work_to_do = nghost.max() > 0;
    }
    if (!work_to_do) return;

    // contains information on what cells to send
    const FB& TheFB = getFB(nghost, period, cross, enforce_periodicity_only);

    if (ParallelContext::NProcsSub() == 1)
    {
        //
        // There can only be local work to do.
        //
	    int N_locs = (*TheFB.m_LocTags).size();
        if (N_locs == 0) return;
/*#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))
            if (Gpu::inGraphRegion())
            {
                FB_local_copy_cuda_graph_1(TheFB, scomp, ncomp);
            }
            else
#endif
            {
                FB_local_copy_gpu(TheFB, scomp, ncomp);
            }
        }
        else
#endif*/
        {
            FB_local_copy_cpu(TheFB, scomp, ncomp); 
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum = ParallelDescriptor::SeqNum();
    fb_tag = SeqNum;

    const int N_locs = TheFB.m_LocTags->size();
    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        // No work to do.
        return;

    // Post receives
    if (N_rcvs > 0) {
        PostReceives(*TheFB.m_RcvTags, fb_send_reqs,
                    scomp, ncomp, SeqNum);
        //fb_recv_stat.resize(N_rcvs);
    }

    // Post sends
    if (N_snds > 0)
    {
        PostSends(*TheFB.m_SndTags, fb_send_reqs,
                    scomp, ncomp, SeqNum); 
    }

    //FillBoundary_test();  

    // Do the local work.  Hope for a bit of communication/computation overlap.
    if (N_locs > 0)
    {
/*#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
            if (Gpu::inGraphRegion()) {
                FB_local_copy_cuda_graph_n(TheFB, scomp, ncomp);
            }
            else
#endif
            {
                FB_local_copy_gpu(TheFB, scomp, ncomp);
            }
        }
        else
#endif*/
            FB_local_copy_cpu(TheFB, scomp, ncomp);
	}

    //FillBoundary_test();  
#endif
    return;
}


#ifdef BL_USE_MPI
void
CellFabArray::PostSends (const MapOfCopyComTagContainers& m_SndTags,
                        Vector<MPI_Request> send_reqs,
                        int icomp, int ncomp, int SeqNum)
{

    Vector<int>                         send_rank;
    Vector<const CopyComTagsContainer*> send_cctc;
    send_reqs.clear();

    for (auto const& kv : m_SndTags)
    {
        send_rank.push_back(kv.first);
        send_reqs.push_back(MPI_REQUEST_NULL);
        send_cctc.push_back(&kv.second);
    }

    const int nsend = send_rank.size();

    MPI_Comm comm = ParallelContext::CommunicatorSub();
    // Get mpi datatypes:
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
    for (int j = 0; j < nsend; ++j)
    {
        auto const& cctc = *send_cctc[j];

        int number_of_sends = 0;
        for (auto const& tag : cctc)
            number_of_sends += tag.sbox.numPts();
        number_of_sends *= ncomp;

        std::vector<void*> addresses(number_of_sends, NULL);
        std::vector<int> counts(number_of_sends, -1);
        std::vector<MPI_Datatype> datatypes(number_of_sends, MPI_DATATYPE_NULL);

        int offset = 0;
        for (auto const& tag : cctc)
        {
            const Box& bx = tag.sbox;
            auto sfab = this->array(tag.srcIndex);
            
            // TEST: LoopConcurrentOnCpu or for(cell: numPts)
            amrex::Loop( bx, ncomp,
            [&] (int ii, int jj, int kk, int n) noexcept
            {
                const IntVect idx = IntVect(ii,jj,kk);
                const Long cell = bx.index(idx) + bx.numPts()*n + offset;
                std::tie(
                    addresses[cell], 
                    counts[cell], 
                    datatypes[cell]
                ) = sfab(idx, n+icomp).get_mpi_datatype();
            });
            offset += bx.numPts()*ncomp;
        }

        // get displacements in bytes for incoming user data
        std::vector<MPI_Aint> displacements(addresses.size(), 0);
        for (size_t i = 0; i < addresses.size(); i++) {
            displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
        }

        MPI_Datatype send_datatype;
        MPI_Type_create_struct(
            number_of_sends,
            &counts[0],
            &displacements[0],
            &datatypes[0],
            &send_datatype
        );
        
        MPI_Type_commit(&send_datatype);

        const int rank = ParallelContext::global_to_local_rank(send_rank[j]);

        MPI_Isend(
            addresses[0],
            1,
            send_datatype,
            rank,
            SeqNum,
            comm,
            &send_reqs[j]);

        MPI_Type_free(&send_datatype);
    }
    // TODO: Do MPI_Isends outside omp loop
}

void
CellFabArray::PostReceives (const MapOfCopyComTagContainers& m_RcvTags, 
                        Vector<MPI_Request>& recv_reqs,
                        int icomp, int ncomp, int SeqNum)
{

    Vector<int>&                        recv_from = fb_recv_from;
    Vector<const CopyComTagsContainer*> send_cctc;

    recv_from.clear();
    recv_reqs.clear();

    Vector<const CopyComTagsContainer*> recv_cctc;

    for (const auto& kv : m_RcvTags)
    {
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
        recv_cctc.push_back(&kv.second);
    }

    const int nrecv = recv_from.size();

    MPI_Comm comm = ParallelContext::CommunicatorSub();
    for (int j = 0; j < nrecv; ++j)
    {
        const int rank = ParallelContext::global_to_local_rank(recv_from[j]);

        auto const& cctc = *recv_cctc[j];

        int number_of_receives = 0;
        for (auto const& tag : cctc) {
            number_of_receives += tag.dbox.numPts();
        }
        number_of_receives *= ncomp;

        std::vector<void*> addresses(number_of_receives, NULL);
        std::vector<int> counts(number_of_receives, -1);
        std::vector<MPI_Datatype> datatypes(number_of_receives, MPI_DATATYPE_NULL);

        int offset = 0;
        for (auto const& tag : cctc)
        {
            const Box& bx = tag.dbox;
            auto dfab = this->array(tag.dstIndex);

            // TEST: LoopConcurrentOnCpu or for(cell: numPts)
            amrex::Loop( bx, ncomp,
            [&] (int ii, int jj, int kk, int n) noexcept
            {
                const IntVect idx = IntVect(ii,jj,kk);
                const Long cell = bx.index(idx) + bx.numPts()*n  + offset;

                std::tie(
                    addresses[cell], 
                    counts[cell], 
                    datatypes[cell]
                ) = dfab(idx, n+icomp).get_mpi_datatype();
            });
            offset += bx.numPts()*ncomp;
        }
        
        // get displacements in bytes for incoming user data
        std::vector<MPI_Aint> displacements(number_of_receives, 0);
        for (size_t i = 0; i < addresses.size(); i++) {
            displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
        }

        int ret_val;
        MPI_Datatype recv_datatype;
        ret_val = MPI_Type_create_struct(
            number_of_receives,
            &counts[0],
            &displacements[0],
            &datatypes[0],
            &recv_datatype
        );
        
        ret_val = MPI_Type_commit(&recv_datatype);

        ret_val = MPI_Irecv(
            addresses[0],
            1,
            recv_datatype,
            rank,
            SeqNum,
            comm,
            &recv_reqs[j]);

        MPI_Type_free(&recv_datatype);
    }
}
#endif

void
CellFabArray::FillBoundary_finish () 
{
    BL_PROFILE("FillBoundary_finish()");

    if ( n_grow.allLE(IntVect::TheZeroVector()) && !fb_epo ) return; // For epo (Enforce Periodicity Only), there may be no ghost cells.

    n_filled = fb_nghost;

    if (ParallelContext::NProcsSub() == 1) return;
    
#ifdef AMREX_USE_MPI

    const FB& TheFB = getFB(fb_nghost,fb_period,fb_cross,fb_epo);
    const int N_rcvs = TheFB.m_RcvTags->size();
    if (N_rcvs > 0)
    {
        fb_recv_stat.resize(N_rcvs);
        MPI_Waitall(fb_recv_reqs.size(), fb_recv_reqs.dataPtr(), fb_recv_stat.dataPtr());
        // = ParallelDescriptor::Waitall(fb_recv_reqs, fb_recv_stat);
    }
    
    const int N_snds = TheFB.m_SndTags->size();
    if (N_snds > 0) 
    {
        Vector<MPI_Status> stats(N_snds);
        MPI_Waitall(fb_send_reqs.size(), fb_send_reqs.dataPtr(), stats.dataPtr());
        // = ParallelDescriptor::Waitall(fb_send_reqs, stats);
    }
#endif
}

void
CellFabArray::FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;

    bool is_thread_safe = TheFB.m_threadsafe_loc;
    if (is_thread_safe)
    {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            const auto* sfab = &(get(tag.srcIndex));
                  auto* dfab = &(get(tag.dstIndex));
            dfab->template copy<RunOn::Host>(*sfab, tag.sbox, scomp, tag.dbox, scomp, ncomp);
        }
    }
    else
    {
        LayoutData<Vector<FabCopyTag<CellFab> > > loc_copy_tags(boxArray(),DistributionMap());
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            loc_copy_tags[tag.dstIndex].push_back
                ({this->fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
        }
/*#ifdef _OPENMP
#pragma omp parallel
#endif*/
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const auto& tags = loc_copy_tags[mfi];
            auto dfab = this->array(mfi);
            for (auto const & tag : tags)
            {
                auto const sfab = tag.sfab->array();
                const auto offset = tag.offset.dim3();
                amrex::LoopOnCpu(tag.dbox, ncomp,
                [=] (int i, int j, int k, int n) noexcept
                {
                    dfab(i,j,k,n+scomp) = sfab(i+offset.x,j+offset.y,k+offset.z,n+scomp);
                });
            }
        }
    }
}

void 
CellFabArray::ParallelCopy (const CellFabArray& src,
                            int scomp, int dcomp, int ncomp,
                            const IntVect& snghost, const IntVect& dnghost,
                            const Periodicity& period,
                            CpOp op,
                            const FabArrayBase::CPC* a_cpc)
{
    BL_PROFILE("FabArray::ParallelCopy()");
    amrex::Print() << "Warning: Calling ParallelCopy()\n";

    if (size() == 0 || src.size() == 0) return;

    n_filled = dnghost;

    if ((src.boxArray().ixType().cellCentered() || op == FabArrayBase::COPY) 
        && (boxarray == src.boxarray && distributionMap == src.distributionMap)
	    && snghost == IntVect::TheZeroVector() && dnghost == IntVect::TheZeroVector()
        && !period.isAnyPeriodic())
    {
/*#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif*/
        for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

            // avoid self copy or plus
            if (this != &src) {
                auto const sfab = src.array(fai);
                auto       dfab = this->array(fai);
                if (op == FabArrayBase::COPY) {
                    /*AMREX_HOST_DEVICE_PARALLEL_FOR_4D*/
                    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,dcomp+n) = sfab(i,j,k,scomp+n);
                    });
                } else {
                    /*AMREX_HOST_DEVICE_PARALLEL_FOR_4D*/
                    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,dcomp+n) += sfab(i,j,k,scomp+n);
                    });
                }
            }
        }

        return;
    }

    const CPC& thecpc = (a_cpc) ? *a_cpc : getCPC(dnghost, src, snghost, period);

    if (ParallelContext::NProcsSub() == 1)
    {
        // There can only be local work to do.
	    int N_locs = (*thecpc.m_LocTags).size();

        if (N_locs == 0) return;
/*#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            PC_local_gpu(thecpc, src, scomp, dcomp, ncomp, op);
        }
        else
#endif*/
        {
            PC_local_cpu(thecpc, src, scomp, dcomp, ncomp, op);
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum  = ParallelDescriptor::SeqNum();

    const int N_snds = thecpc.m_SndTags->size();
    const int N_rcvs = thecpc.m_RcvTags->size();
    const int N_locs = thecpc.m_LocTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        // No work to do.
        return;
    }

    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int ipass = 0, SC = scomp, DC = dcomp; ipass < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);
        Vector<MPI_Request> recv_reqs;
        Vector<MPI_Request> send_reqs;

        // Post rcvs.
        if (N_rcvs > 0) {
            PostReceives(*thecpc.m_RcvTags, recv_reqs, SC, NC, SeqNum);
        }

        // Post send's
        if (N_snds > 0) {
            PostSends(*thecpc.m_SndTags, send_reqs, SC, NC, SeqNum);
        }

        // Do the local work.  Hope for a bit of communication/computation overlap.
        if (N_locs > 0)
	    {
/*#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion())
            {
                PC_local_gpu(thecpc, src, SC, DC, NC, op);
            }
            else
#endif*/
            {
                PC_local_cpu(thecpc, src, SC, DC, NC, op);
            }
        }

        if (N_rcvs > 0)
        {
            Vector<MPI_Status> stats(N_rcvs);
            MPI_Waitall(recv_reqs.size(), recv_reqs.dataPtr(), stats.dataPtr());
        }
	
        if (N_snds > 0 && !thecpc.m_SndTags->empty()) {
            Vector<MPI_Status> stats(N_snds);
            MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr());
        }

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

    return;

#endif /*BL_USE_MPI*/
}

void
CellFabArray::PC_local_cpu (const CPC& thecpc, CellFabArray const& src,
                             int scomp, int dcomp, int ncomp, CpOp op)
{
    int N_locs = thecpc.m_LocTags->size();
    if (N_locs == 0) return;
    bool is_thread_safe = thecpc.m_threadsafe_loc;

    if (is_thread_safe)
    {
/*#ifdef _OPENMP
#pragma omp parallel for
#endif*/
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];
            if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                // avoid self copy or plus
                const auto* sfab = &(src[tag.srcIndex]);
                      auto* dfab = &(get(tag.dstIndex));
                if (op == FabArrayBase::COPY)
                {
                    dfab->template copy<RunOn::Host>(*sfab, tag.sbox, scomp, tag.dbox, dcomp, ncomp);
                }
                else
                {
                    dfab->template plus<RunOn::Host>(*sfab, tag.sbox, tag.dbox, scomp, dcomp, ncomp);
                }
            }
        }
    }
    else
    {
        LayoutData<Vector<FabCopyTag<CellFab> > > loc_copy_tags(boxArray(),DistributionMap());
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];
            if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                loc_copy_tags[tag.dstIndex].push_back
                    ({src.fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
            }
        }

/*#ifdef _OPENMP
#pragma omp parallel
#endif*/
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const auto& tags = loc_copy_tags[mfi];
            auto dfab = this->array(mfi);
            if (op == FabArrayBase::COPY)
            {
                for (auto const & tag : tags)
                {
                    auto const sfab = tag.sfab->array();
                    Dim3 offset = tag.offset.dim3();
                    amrex::LoopConcurrentOnCpu (tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,dcomp+n) = sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                    });
                }
            }
            else
            {
                for (auto const & tag : tags)
                {
                    auto const sfab = tag.sfab->array();
                    Dim3 offset = tag.offset.dim3();
                    amrex::LoopConcurrentOnCpu (tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,dcomp+n) += sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                    });
                }
            }
        }
    }
}
