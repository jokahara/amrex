
#ifndef CopyDescriptor_H_
#define CopyDescriptor_H_

#include <AMReX_FACopyDescriptor.H>
#include "CellFabArray.h"

/**
* \brief This class orchestrates filling a destination fab of size destFabBox
* from fabarray on the local processor (myProc).
*/
class CopyDescriptor : public FabArrayCopyDescriptor<CellFab>
{
  public:

    CopyDescriptor () : FabArrayCopyDescriptor<CellFab>() {}

    CopyDescriptor (const CopyDescriptor&) = delete;
    CopyDescriptor& operator= (const CopyDescriptor&) = delete;

    void CollectData ();

    bool DataAvailable () const { return dataAvailable; }
    int CurrentNFabArrays () const { return fabArrays.size(); }
    int nFabComTags () const { return fabComTagList.size(); }
    int nFabCopyDescs () const { return fabCopyDescList.size(); }

    /*
    FabArrayId RegisterFabArray(FabArray<FAB> *fabarray);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp);
                      
    //! Add a box but only from FabArray[fabarrayindex].
    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        fabarrayindex,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp,
                      bool       bUseValidBox = true);


    // copies local data to destFab 
    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab);
    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab,
                  const Box&       destBox);

    void PrintStats () const;
    void clear ();*/

private:
    //!
    //! Helper function for AddBox() routines.
    //!
    /*void AddBoxDoIt (FabArrayId fabarrayid,
                     const Box& destFabBox,
                     BoxList*   returnedUnfilledBoxes,
                     int        faindex,
                     int        srccomp,
                     int        destcomp,
                     int        numcomp,
                     bool       bUseValidBox,
                     BoxDomain& unfilledBoxDomain);*/

    // Some useful typedefs.
    typedef std::multimap<int,FabCopyDescriptor<CellFab>*> FCDMap;
    typedef typename FCDMap::value_type     FCDMapValueType;
    typedef typename FCDMap::iterator       FCDMapIter;
    typedef typename FCDMap::const_iterator FCDMapConstIter;

    typedef std::vector<FabArrayBase::FabComTag> FabComTagContainer;
    typedef std::vector<FabComTagContainer::const_iterator> FabComTagIterContainer;

    // data
    std::vector<FabArray<CellFab>*> fabArrays;
    std::vector<FCDMap>         fabCopyDescList;
    FabComTagContainer          fabComTagList;
    int                         nextFillBoxId;
    bool                        dataAvailable;
};

/*void
InterpAddBox (CopyDescriptor&       fabCopyDesc,
		      BoxList*              returnUnfilledBoxes,
		      Vector<FillBoxId>&    returnedFillBoxIds,
		      const Box&            subbox,
		      FabArrayId            faid1,
		      FabArrayId            faid2,
		      int                   src_comp,
		      int                   dest_comp,
		      int                   num_comp)
{
    if (something)
    {
        returnedFillBoxIds.resize(1);
        returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1,
                                                   subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
    }
    else
    {
        returnedFillBoxIds.resize(2);
        BoxList tempUnfilledBoxes(subbox.ixType());
        returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1,
                                                   subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
        returnedFillBoxIds[1] = fabCopyDesc.AddBox(faid2,
                                                   subbox,
                                                   &tempUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
        //
        // The boxarrays for faid1 and faid2 should be the
        // same so only use returnUnfilledBoxes from one AddBox here.
        //
    }
}

void
InterpFillFab (CopyDescriptor&          fabCopyDesc,
		       const Vector<FillBoxId>& fillBoxIds,
		       FabArrayId               faid1,
		       FabArrayId               faid2,
		       CellFab&                 dest,
		       int                      src_comp,   // these comps need to be removed
		       int                      dest_comp,  // from this routine
		       int                      num_comp)
{
    if (something)
    {
        fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest);
    }
    else
    {
        BL_ASSERT(dest_comp + num_comp <= dest.nComp());

        CellFab dest1(dest.box(), dest.nComp());
        CellFab dest2(dest.box(), dest.nComp());

        fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest1);
        fabCopyDesc.FillFab(faid2, fillBoxIds[1], dest2);

        // dest.linInterp;
    }
}*/

/*void
CopyDescriptor::AddBoxDoIt (FabArrayId fabarrayid,
                            const Box& destFabBox,
                            BoxList*   returnedUnfilledBoxes,
                            int        faindex,
                            int        srccomp,
                            int        destcomp,
                            int        numcomp,
                            bool       bUseValidBox,
                            BoxDomain& unfilledBoxDomain)
{
    const int myProc = ParallelDescriptor::MyProc();

    FabArray<CellFab>* fabArray = fabArrays[fabarrayid.Id()];

    BL_ASSERT(faindex >= 0 && faindex < fabArray->size());

    Box intersect = destFabBox;

    if (bUseValidBox)
    {
        intersect &= fabArray->box(faindex);
    }
    else
    {
        intersect &= fabArray->fabbox(faindex);
    }

    if (intersect.ok())
    {
        FabCopyDescriptor<CellFab>* fcd = new FabCopyDescriptor<CellFab>;

        int remoteProc     = fabArray->DistributionMap()[faindex];
        if(remoteProc >= ParallelDescriptor::NProcs()) {
            std::ostringstream ss;
            ss << ParallelDescriptor::MyProc() << ":: _in AddBoxDoIt:  nProcs remoteProc = "
                << ParallelDescriptor::NProcs() << "  " << remoteProc << "\n";
            amrex::Abort("Bad remoteProc: "+ss.str());
        }
        fcd->fillBoxId     = nextFillBoxId;
        fcd->subBox        = intersect;
        fcd->myProc        = myProc;
        fcd->copyFromProc  = remoteProc;
        fcd->copyFromIndex = faindex;
        fcd->srcComp       = srccomp;
        fcd->destComp      = destcomp;
        fcd->nComp         = numcomp;

        if (ParallelDescriptor::sameTeam(remoteProc))
        {
            //
            // Data is local.
            //
            fcd->fillType       = FillLocally;
            fcd->localFabSource = &(*fabArray)[faindex];
        }
        else
        {
            //
            // Data is remote.
            //
            FabArrayBase::FabComTag fabComTag;

            dataAvailable               = false;
            fcd->fillType               = FillRemotely;
            fcd->localFabSource         = new CellFab(intersect, numcomp);
            fcd->cacheDataAllocated     = true;
            fabComTag.fabArrayId        = fabarrayid.Id();
            fabComTag.fillBoxId         = nextFillBoxId;
            fabComTag.fabIndex          = faindex;
            fabComTag.procThatNeedsData = myProc;
            fabComTag.procThatHasData   = remoteProc;
            fabComTag.box               = intersect;
            fabComTag.srcComp           = srccomp;
            fabComTag.destComp          = destcomp;
            fabComTag.nComp             = numcomp;
            //
            // Do not send the data yet.
            //
            fabComTagList.push_back(fabComTag);
        }

        fabCopyDescList[fabarrayid.Id()].insert(FCDMapValueType(fcd->fillBoxId,fcd));

        if (returnedUnfilledBoxes != 0)
        {
            unfilledBoxDomain.rmBox(intersect);
        }
    }
}*/

inline void
CopyDescriptor::CollectData ()
{
    dataAvailable = true;

    if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI

    BL_PROFILE("FabArrayCopyDescriptor::CollectData()");

    const int MyProc = ParallelDescriptor::MyProc();

    int Total_Rcvs_Size = 0;
    
    // We use this to make finding matching FabComTags more efficient.
    std::map< int, FabComTagIterContainer > RcvTags;

    std::map<int,int> Snds, Rcvs, Npts;
    
    // Set Rcvs[i] to # of blocks needed from CPU i
    for (FabComTagContainer::const_iterator it = fabComTagList.begin(),
             End = fabComTagList.end(); it != End; ++it)
    {
        BL_ASSERT(it->box.ok());
        BL_ASSERT(it->procThatNeedsData == MyProc);
        BL_ASSERT(it->procThatHasData   != MyProc);

        const int Who = it->procThatHasData;
        const int Cnt = (it->box.numPts())*(it->nComp);

        RcvTags[Who].push_back(it);

        Total_Rcvs_Size += Cnt;

        if (Rcvs.count(Who) > 0)
        {
            Rcvs[Who] += 1;
        }
        else
        {
            Rcvs[Who] = 1;
        }

        if (Npts.count(Who) > 0)
        {
            Npts[Who] += Cnt;
        }
        else
        {
            Npts[Who] = Cnt;
        }
    }
    BL_ASSERT(Rcvs.count(MyProc) == 0);

    const int NProcs = ParallelDescriptor::NProcs();
    {
        Vector<int> SndsArray(NProcs,0), RcvsArray(NProcs,0);

        for (std::map<int,int>::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
	    {
            RcvsArray[it->first] = it->second;
	    }

        {
            BL_PROFILE_VAR("CollectData_Alltoall()", blpvCDATA);
            BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
                            BLProfiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoall(RcvsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         SndsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         ParallelDescriptor::Communicator()) );

            BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
                            BLProfiler::AfterCall());

            BL_PROFILE_VAR_STOP(blpvCDATA);
        }
        BL_ASSERT(SndsArray[MyProc] == 0);

        for (int i = 0; i < NProcs; i++) 
            if (SndsArray[i] > 0)
                Snds[i] = SndsArray[i];
    }

    // There are two rounds of send and recv.
    // First, the data receivers need to send the data senders meta-data (e.g., boxes).
    // Then, the senders know what data to send and perform send.
    const int SeqNum_md   = ParallelDescriptor::SeqNum();
    const int SeqNum_data = ParallelDescriptor::SeqNum();

    const int N_snds = Snds.size();
    const int N_rcvs = Rcvs.size();

    if ( N_snds == 0 && N_rcvs == 0 ) return;

    const int Nints = 4 + 3*AMREX_SPACEDIM;  // # of ints in a meta-data

    MPI_Comm comm = ParallelContext::CommunicatorSub();

    // for meta-data
    Vector<int> md_sender, md_offset, md_icnts, md_bcnts;
    int* md_recv_data;
    Vector<int*> md_send_data;
    Vector<MPI_Request> md_recv_reqs, md_send_reqs;

    // for data
    Vector<int> data_sender, data_offset;
    //Vector<CellFab*> recv_data;
    Vector<CellFab*> send_data;
    Vector<MPI_Request> data_recv_reqs, data_send_reqs;

    if (N_snds > 0)
    {
	    // Recv meta-data
	    int N = 0;
        for (std::map<int,int>::const_iterator it = Snds.begin(), End = Snds.end(); it != End; ++it)
        {
            md_sender.push_back(it->first);
            md_bcnts.push_back(it->second);
            int cnt = it->second * Nints;
            md_icnts.push_back(cnt);
            md_offset.push_back(N);
            N += cnt;
        }

        md_recv_data = static_cast<int*>(amrex::The_Arena()->alloc(N*sizeof(int)));

        for (int i = 0; i < N_snds; ++i)
        {
            md_recv_reqs.push_back(ParallelDescriptor::Arecv(&md_recv_data[md_offset[i]],
                                    md_icnts[i], md_sender[i],
                                    SeqNum_md).req());
        }
    }

    if (N_rcvs > 0)
    {
	    // Send meta-data
        for (std::map<int,int>::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
        {
            int rank = it->first;
            int Nmds = it->second;
            int cnt = Nmds * Nints;

            int* p = static_cast<int*>(amrex::The_Arena()->alloc(cnt*sizeof(int)));
            md_send_data.push_back(p);

            const FabComTagIterContainer& tags = RcvTags[rank];

            // initialized the data
            int * md = p;
            for (int i = 0; i < Nmds; ++i, md += Nints)
            {
                md[0] = tags[i]->fabArrayId;
                md[1] = tags[i]->fabIndex;
                md[2] = tags[i]->srcComp;
                md[3] = tags[i]->nComp;
                const int* lo = tags[i]->box.loVect();
                const int* hi = tags[i]->box.hiVect();
                const IntVect& bxtyp = tags[i]->box.type();
                const int* tp = bxtyp.getVect();
                AMREX_D_EXPR(md[4] = lo[0],
                            md[5] = lo[1],
                            md[6] = lo[2]);
                AMREX_D_EXPR(md[4+  AMREX_SPACEDIM] = hi[0],
                            md[5+  AMREX_SPACEDIM] = hi[1],
                            md[6+  AMREX_SPACEDIM] = hi[2]);
                AMREX_D_EXPR(md[4+2*AMREX_SPACEDIM] = tp[0],
                            md[5+2*AMREX_SPACEDIM] = tp[1],
                            md[6+2*AMREX_SPACEDIM] = tp[2]);
            }

            md_send_reqs.push_back(ParallelDescriptor::Asend(p,cnt,rank,SeqNum_md).req());
        }

        //recv_data = static_cast<value_type*>(amrex::The_Arena()->alloc(Total_Rcvs_Size*sizeof(value_type)));
        std::pair<FCDMapIter,FCDMapIter> match;
        std::map< int,FabComTagIterContainer >::const_iterator found;

        for (std::map<int,int>::const_iterator it = Npts.begin(); it != Npts.end(); ++it)
        {
            int Who = it->first;
            int Cnt = it->second;
            data_sender.push_back(Who);
        }

        // Post receives for data
        for (int k = 0; k < N_rcvs; k++)
        {
            const int       Who = data_sender[k];
            //const CellFab* dptr = &recv_data[data_offset[k]];

            found = RcvTags.find(Who);

            BL_ASSERT(found != RcvTags.end());

            const FabComTagIterContainer& tags = found->second;

            for (FabComTagIterContainer::const_iterator it = tags.begin(),
                    End = tags.end(); it != End; ++it)
            {
                const FabArrayBase::FabComTag& tag = **it;

                BL_ASSERT(tag.procThatHasData == Who);

                match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
                {
                    FabCopyDescriptor<CellFab>* fcdp = (*fmi).second;
                    BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                    if (fcdp->subBox == tag.box)
                    {
                        const Box& bx = tag.box;
                        const int number_of_receives = bx.numPts()*tag.nComp;

                        std::vector<void*> addresses(number_of_receives, NULL);
                        std::vector<int> counts(number_of_receives, -1);
                        std::vector<MPI_Datatype> datatypes(number_of_receives, MPI_DATATYPE_NULL);
                        
                        auto dfab = fcdp->localFabSource->array();

                        amrex::Loop( bx, tag.nComp,
                        [&] (int ii, int jj, int kk, int n) noexcept
                        {
                            const IntVect idx = IntVect(ii,jj,kk);
                            const Long cell = bx.index(idx) + bx.numPts()*n;

                            std::tie(
                                addresses[cell], 
                                counts[cell], 
                                datatypes[cell]
                            ) = dfab(idx, n).get_mpi_datatype();
                        });

                        // get displacements in bytes for incoming user data
                        std::vector<MPI_Aint> displacements(number_of_receives, 0);
                        for (size_t i = 0; i < addresses.size(); i++) {
                            displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
                        }

                        MPI_Datatype recv_datatype;
                        MPI_Type_create_struct(
                            number_of_receives,
                            &counts[0],
                            &displacements[0],
                            &datatypes[0],
                            &recv_datatype
                        );
                        
                        MPI_Type_commit(&recv_datatype);

                        MPI_Request recv_reqs;
                        MPI_Irecv(
                            addresses[0],
                            1,
                            recv_datatype,
                            Who,
                            SeqNum_data,
                            comm,
                            &recv_reqs);

                        data_recv_reqs.push_back(recv_reqs);

                        MPI_Type_free(&recv_datatype);

                        //fcdp->localFabSource->template copyFromMem<RunOn::Host>(tag.box,0,tag.nComp,dptr);
	                    /*data_recv_reqs.push_back(ParallelDescriptor::Arecv(&recv_data[Idx],
							                    Cnt,Who,SeqNum_data).req());*/

                        
                        break;  // why break?
                    }
                }
            }
        }
    }

    // Wait on meta-data and do send
    if (N_snds > 0)
    {
	    int send_counter = 0;
        while (send_counter++ < N_snds)
        {
            MPI_Status status;
            int index;
            ParallelDescriptor::Waitany(md_recv_reqs, index, status);

            int rank = status.MPI_SOURCE;
            BL_ASSERT(status.MPI_TAG == SeqNum_md);
            BL_ASSERT(rank == md_sender[index]);

            const int* p = &md_recv_data[md_offset[index]];
            int numboxes = md_bcnts[index];
            Vector<int> faid(numboxes);
            Vector<int> fidx(numboxes);
            Vector<int> scomp(numboxes);
            Vector<int> ncomp(numboxes);
            Vector<int> npts(numboxes);
            Vector<Box> bxs;
            int N = 0;
            const int * md = p;

            for (int i = 0; i < numboxes; ++i, md += Nints)
            {
                faid[i] = md[0]; fidx[i] = md[1];
                scomp[i] = md[2]; ncomp[i] = md[3];
                bxs.push_back(Box(IntVect(&md[4]),
                                  IntVect(&md[4+AMREX_SPACEDIM]),
                                  IntVect(&md[4+AMREX_SPACEDIM*2])));
                npts[i] = bxs.back().numPts()*ncomp[i];
                N += npts[i];
            }

            BL_ASSERT(N < std::numeric_limits<int>::max());

            //CellFab* data = static_cast<CellFab*>(amrex::The_Arena()->alloc(N*sizeof(CellFab)));
            //CellFab* dptr = data;
            //send_data.push_back(data);

            for (int i = 0; i < numboxes; ++i)
            {
                // send data
                CellFab &sfab = (*fabArrays[faid[i]])[fidx[i]]; //template copyToMem<RunOn::Host>(bxs[i],scomp[i],ncomp[i],dptr);

                int number_of_sends = bxs[i].numPts() * ncomp[i];
                std::vector<void*> addresses(number_of_sends, NULL);
                std::vector<int> counts(number_of_sends, -1);
                std::vector<MPI_Datatype> datatypes(number_of_sends, MPI_DATATYPE_NULL);

                amrex::Loop( bxs[i], ncomp[i],
                [&] (int ii, int jj, int kk, int n) noexcept
                {
                    const IntVect idx = IntVect(ii,jj,kk);
                    const Long cell = bxs[i].index(idx) + bxs[i].numPts()*n ;
                    std::tie(
                        addresses[cell], 
                        counts[cell], 
                        datatypes[cell]
                    ) = sfab(idx, n+scomp[i]).get_mpi_datatype();
                });

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

                MPI_Request send_reqs;
                MPI_Isend(
                    addresses[0],
                    1,
                    send_datatype,
                    rank,
                    SeqNum_data,
                    comm,
                    &send_reqs);

                data_send_reqs.push_back(send_reqs); 

                MPI_Type_free(&send_datatype);
            }
        }
        amrex::The_Arena()->free(md_recv_data);
    }

    // Wait
    if (N_rcvs > 0)
    {
	    Vector<MPI_Status> stats(N_rcvs);
        ParallelDescriptor::Waitall(md_send_reqs, stats);
        // free metadata
        for (int i = 0; i < N_rcvs; ++i) {
            amrex::The_Arena()->free(md_send_data[i]);
        }

        ParallelDescriptor::Waitall(data_recv_reqs, stats);
    }

    // Finished send
    if (N_snds > 0)
    {
	    Vector<MPI_Status> stats(N_snds);
        ParallelDescriptor::Waitall(data_send_reqs, stats);
        //MPI_Waitall(data_send_reqs.size(), data_send_reqs.dataPtr(), stats.dataPtr());
    }

#endif /*BL_USE_MPI*/
}

#endif
