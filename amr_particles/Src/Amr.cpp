#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <AMReX_ParmParse.H>

#include "Amr.h"

bool                   Amr::initialized;
Vector<BoxArray>       Amr::initial_ba;
Vector<BoxArray>       Amr::regrid_ba;

void 
Amr::Initialize () // static function
{
    if (initialized) return;
    
    amrex::ExecOnFinalize(Amr::Finalize);
    initialized = true;
}

void
Amr::Finalize () // static function
{
    Amr::regrid_ba.clear();
    Amr::initial_ba.clear();

    initialized = false;
}

Long
Amr::cellCount (int lev) noexcept
{
    return amr_level[lev]->countCells();
}

int
Amr::numGrids (int lev) noexcept
{
    return amr_level[lev]->numGrids();
}

Amr::Amr () : AmrCore()
{
    Initialize();
    InitAmr();
}

Amr::Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord)
    : AmrCore(rb,max_level_in,n_cell_in,coord)
{
    Initialize();
    InitAmr();
}

void
Amr::InitAmr ()
{

    BL_PROFILE("Amr::InitAmr()");
    //
    // Set default values.
    //

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        isPeriodic[i] = false;

    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    pp.query("initial_grid_file", initial_grids_file);
    pp.query("regrid_file"      , regrid_grids_file);

    // Generate internal values from user-supplied values.
    n_comp = 1;
    pp.query("n_comp", n_comp);

    int grow = 1;
    pp.query("n_grow", grow);
    n_grow = {grow, grow, grow};
    
    if (check_input) checkInput();
    
    // Make the default regrid_int = 1 for all levels.
    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (int i = 0; i < max_level; i++)
           regrid_int[i]  = 1;
    }

    // Determine which directions to refine (for 2D simulations)
	Vector<int> ref_dir(AMREX_SPACEDIM, 1);
    pp.getarr("ref_dir",ref_dir,0,AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if(!ref_dir[i]) {
            for (int j = 0; j <= max_level; j++) {
                ref_ratio[j][i] = 1;
            }
            // redefining geometries
            for (int j = 1; j <= max_level; j++) {
                Geometry& geom = Geom(j);
                geom.Domain(Geom(j-1).Domain());
                geom.refine(ref_ratio[j-1]);
            }
        }
    }
    
    // make empty pointers to levels
    int nlev = max_level+1;
    amr_level.resize(nlev);

    // Read in the regrid interval if max_level > 0.
    if (max_level > 0) 
    {
       int numvals = pp.countval("regrid_int");
       if (numvals == 1)
       {
           //
           // Set all values to the single available value.
           //
           int the_regrid_int = 0;
           pp.query("regrid_int",the_regrid_int);
           for (int i = 0; i < max_level; i++)
           {
               regrid_int[i] = the_regrid_int;
           }
       }
       else if (numvals == 0)
       {
           if (verbose) {
               amrex::Print(amrex::ErrorStream()) << "Using default regrid_int = 1 at all levels!\n";
           }
       }
       else if (numvals < max_level)
       {
           amrex::Error("You did not specify enough values of regrid_int");
       }
       else 
       {
           //
           // Otherwise we expect a vector of max_level values
           //
           pp.queryarr("regrid_int",regrid_int,0,max_level);
       }
    }

    if (max_level > 0 && !initial_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' ) {}
        std::ifstream is(initial_grids_file.c_str(),std::ios::in);

        if (!is.good())
            amrex::FileOpenFailed(initial_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        initial_ba.resize(in_finest);

        use_fixed_upto_level = in_finest;
        if (in_finest > max_level)
           amrex::Error("You have fewer levels in your inputs file then in your grids file!");

        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (int i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                bx.refine(ref_ratio[lev-1]);
                bl.push_back(bx);
            }
            initial_ba[lev-1].define(bl);
        }
        is.close();
        if (verbose > 0) {
            amrex::Print() << "Read initial_ba. Size is " << initial_ba.size() << "\n";
        }

#undef STRIP
    }

    if (max_level > 0 && !regrid_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' ) {}
        std::ifstream is(regrid_grids_file.c_str(),std::ios::in);

        if (!is.good())
            amrex::FileOpenFailed(regrid_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        regrid_ba.resize(in_finest);
        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (int i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                 bx.refine(ref_ratio[lev-1]);
                 for (int idim = 0 ; idim < AMREX_SPACEDIM; ++idim)
                 {
                     if (bx.length(idim) > max_grid_size[lev][idim])
                     {
                         std::ostringstream ss;
                         ss << "Grid " << bx << " too large" << '\n';
                         amrex::Error(ss.str());
                     }
                 }
                 bl.push_back(bx);
            }
            regrid_ba[lev-1].define(bl);
        }
        is.close();
#undef STRIP
    }

    loadbalance_with_workestimates = 0;
    pp.query("loadbalance_with_workestimates", loadbalance_with_workestimates);

    loadbalance_level0_int = 2;
    pp.query("loadbalance_level0_int", loadbalance_level0_int);

    loadbalance_max_fac = 1.5;
    pp.query("loadbalance_max_fac", loadbalance_max_fac);
}

Amr::~Amr ()
{
    Amr::Finalize();
}

Long
Amr::cellCount () noexcept
{
    Long cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i]->countCells();
    return cnt;
}

int
Amr::numGrids () noexcept
{
    int cnt = 0;
    for (int i = 0; i <= finest_level; i++)
        cnt += amr_level[i]->numGrids();
    return cnt;
}

void
Amr::checkInput ()
{
    if (max_level < 0)
        amrex::Error("checkInput: max_level not set");
    //
    // Check that blocking_factor is a power of 2.
    //
    for (int i = 0; i <= max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            int k = blocking_factor[i][idim];
            while ( k > 0 && (k%2 == 0) )
                k /= 2;
            if (k != 1)
                amrex::Error("Amr::checkInput: blocking_factor not power of 2");
        }
    }
    //
    // Check level dependent values.
    //
    for (int i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12) {
            amrex::Error("Amr::checkInput: bad ref_ratios");
        }
    }
    const Box& domain = Geom(0).Domain();
    if (!domain.ok()) {
        amrex::Error("level 0 domain bad or not set");
    }
    //
    // Check that domain size is a multiple of blocking_factor[0].
    //
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        int len = domain.length(i);
        if (len%blocking_factor[0][i] != 0)
            amrex::Error("domain size not divisible by blocking_factor");
    }
    //
    // Check that max_grid_size is even.
    //
    for (int i = 0; i <= max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (max_grid_size[i][idim]%2 != 0) {
                amrex::Error("max_grid_size is not even");
            }
        }
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (int i = 0; i <= max_level; i++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (max_grid_size[i][idim]%blocking_factor[i][idim] != 0) {
                amrex::Error("max_grid_size not divisible by blocking_factor");
            }
        }
    }

    if( ! Geom(0).ProbDomain().ok()) {
        amrex::Error("Amr::checkInput: bad physical problem size");
    }

    if(verbose > 0) {
	amrex::Print() << "Successfully read inputs file ... " << '\n';
    }
}

void
Amr::defBaseLevel (const BoxArray* lev0_grids, const Vector<int>* pmap)
{
    BL_PROFILE("defBaseLevel");

    finest_level = 0;
    // Define base level grids.

    // Check that base domain has even number of zones in all directions.
    const Box& domain   = Geom(0).Domain();
    const IntVect& d_len = domain.size();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
        if (d_len[idir]%2 != 0 && !(d_len[idir] == 1))
            amrex::Error("defBaseLevel: must have even number of cells");

    BoxArray lev0;

    if (lev0_grids != 0 && lev0_grids->size() > 0)
    {
        BL_ASSERT(pmap != 0);

        BoxArray domain_ba(domain);
        if (!domain_ba.contains(*lev0_grids))
            amrex::Error("defBaseLevel: domain does not contain lev0_grids!");
        if (!lev0_grids->contains(domain_ba))
            amrex::Error("defBaseLevel: lev0_grids does not contain domain");

        lev0 = *lev0_grids;

        if (refine_grid_layout) {
            ChopGrids(0,lev0,ParallelDescriptor::NProcs());
        }
    }
    else
    {
	    lev0 = MakeBaseGrids();
    }

    this->SetBoxArray(0, lev0);
    this->SetDistributionMap(0, DistributionMapping(lev0));

    // Now build level 0 grids.
    amr_level[0].reset(new AmrLevel(*this, 0, Geom(0), grids[0], dmap[0]));
    amr_level[0]->init();
}

void
Amr::bldFineLevels ()
{
    BL_PROFILE("Amr::bldFineLevels()");
    finest_level = 0;

    Vector<BoxArray> new_grids(max_level+1);
    // Get initial grid placement.
    do
    {
        int new_finest;
        MakeGrids(finest_level, new_finest, new_grids);

        if (new_finest <= finest_level) break;

        // Make new distribution maps
        Vector<DistributionMapping> new_dmap(new_finest+1);
        MakeDistributionMaps(new_finest, new_grids, new_dmap);
        
        InstallNewLevels(new_finest, new_grids, new_dmap);
        finest_level = new_finest;
    }
    while (finest_level < max_level);
    
    // post regridding operations
    for (int lev = 0; lev <= finest_level; ++lev) {
        amr_level[lev]->post_regrid(0, finest_level);
    }
    
    // Iterate grids to ensure fine grids encompass all interesting gunk.
    // but only iterate if we did not provide a grids file.
    if ( regrid_grids_file.empty() || !initial_grids_file.empty() )  
    {
        bool grids_the_same;
        const int MaxCnt = 4;
        int count = 0;

        do
        {
            for (int i = 0; i <= finest_level; i++) {
                new_grids[i] = amr_level[i]->boxArray();
            }

            regrid(0);

            grids_the_same = true;

            for (int i = 0; i <= finest_level && grids_the_same; i++) {
                if (!(new_grids[i] == amr_level[i]->boxArray())) {
                    grids_the_same = false;
                }
            }

            count++;
        }
        while (!grids_the_same && count < MaxCnt);
    }
}

void
Amr::regrid (int  lbase,
             bool initial)
{
    BL_PROFILE("Amr::regrid()");

    if (lbase > std::min(finest_level,max_level-1)) return;
    if (max_level == 0) return;

    if (verbose > 0)
	amrex::Print() << "Now regridding at level lbase = " << lbase << "\n";

    //
    // Compute positions of new grids.
    //
    int                         new_finest;
    Vector<BoxArray>            new_grid_places(max_level+1);
    Vector<DistributionMapping> new_dmap(max_level+1);

    MakeGrids(lbase, new_finest, new_grid_places);

    // copy grids coarser than lbase
    for (int lev = 0; lev < lbase; lev++) {
        new_grid_places[lev] = amr_level[lev]->boxArray();
    }
    // Define the new grids from level start up to new_finest.
    MakeDistributionMaps(new_finest, new_grid_places, new_dmap);

    bool grids_unchanged = finest_level == new_finest;
    for (int lev = lbase, End = std::min(finest_level,new_finest); lev <= End; lev++) {
        if (new_grid_places[lev] != amr_level[lev]->boxArray()) {
            grids_unchanged = false;
        }
    }

    // If use_efficient_regrid flag is set and grids are unchanged, then don't do anything more here.
    if (grids_unchanged)
    {
        if (verbose > 0) {
            amrex::Print() << "Regridding at level lbase = " << lbase 
                << " but grids unchanged\n";
        }
        return;
    }

    // make new levels and load balance
    InstallNewLevels(new_finest, new_grid_places, new_dmap);

    // Reclaim all remaining storage for levels > new_finest.
    for (int lev = new_finest + 1; lev <= finest_level; ++lev) {
        ClearLevel(lev);
    }

    finest_level = new_finest;

    // Check at *all* levels whether we need to do anything special now that the grids
    // at levels lbase+1 and higher may have changed.  
    for (int lev = 0; lev <= finest_level; ++lev) {
        amr_level[lev]->post_regrid(lbase, finest_level);
    }

    // Report creation of new grids.
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << " : REGRID  with lbase = "
                       << lbase
                       << std::endl;

        printGridSummary(amrex::OutStream(),lbase,finest_level);
    }
}

void Amr::MakeDistributionMaps (int new_finest, Vector<BoxArray>& new_grids,
                                Vector<DistributionMapping>& new_dmap)
{   
    const Box& fdomain = geom[new_finest].Domain();
    Vector<BoxList> boxes(new_finest+1);

    // make BoxLists and refine them towards finest level
    boxes[0] = new_grids[0].boxList();
    for (int lev = 1; lev <= new_finest; lev++)
    {
        boxes[lev] = new_grids[lev].boxList();
        for (int i = 0; i < lev; i++) {
            boxes[i].refine(refRatio(lev-1));
        }
    }

    // chop boxes so that levels correctly overlap
    for (int lev = new_finest-1; lev >= 0; --lev)
    {
        BoxList chopped;
        for (const Box &cbox : boxes[lev].data()) {
            // Calculate complement of finer level boxes
            BoxList complement = amrex::complementIn(cbox, boxes[lev+1]);
            chopped.join(complement);
        }
        // Add boxes of finer level to the end the list
        // This conveniently makes all boxes sorted in the order of coarsest to finest level
        chopped.join(boxes[lev+1].data());
        chopped.removeEmpty();
        
        // redefine new grids
        boxes[lev] = std::move(chopped);
        new_grids[lev].define(boxes[lev]);
    }
    // new_grids[0] contains all boxes
    
    // Weigh cells based on the level they belong to.
    // It is also possible to give a weight for each cell separately
    // using MultiFab weights(ba, dm, 1, 0).
    Vector<Real> weights(boxes[0].size(), 0);

    Real volume = 8; //ref_ratio[0][0] * ref_ratio[0][1] * ref_ratio[0][2];
    int weight = 1, start = 0, end = 0;
    for (int lev = 0; lev <= new_finest; lev++)
    {
        start = end;
        end = (lev == new_finest) 
                ? boxes[0].size()
                : boxes[lev].size() - boxes[lev+1].size();
                
        for (int i = start; i < end; i++) {
            weights[i] = weight;
        }
        weight *= volume;
    }

	// calculate distribution mapping based on weights
    DistributionMapping dm = DistributionMapping::makeSFC(weights, new_grids[0]);

    // define distribution of level 0 
    new_dmap[0].define(dm.ProcessorMap());
    
    // make new distribution mappings for fine levels
    for (int lev = 1; lev <= new_finest; lev++) 
    {
        BoxArray& previous_ba = new_grids[lev-1];
        DistributionMapping& previous_dm = new_dmap[lev-1];

        BoxArray& current_ba = new_grids[lev];
        Vector<int> new_map(new_grids[lev].size(), -1);
        
        // We can copy quickly from coarse level because boxes are sorted
        int start = previous_ba.size() - current_ba.size();
        for (int i = 0; i < current_ba.size(); i++) {
            new_map[i] = previous_dm[start++];
        }

        // define distribution map
        new_dmap[lev].define(std::move(new_map));

        // coarsen grids back to correct refinement level
        for (int i = 0; i < lev; i++) {
            new_grids[i].coarsen(ref_ratio[i]);
        }
    }
}

void
Amr::InstallNewLevels (int new_finest, 
                       const Vector<BoxArray>& new_grids,
                       const Vector<DistributionMapping>& new_dmap)
{
    BL_PROFILE("InstallNewLevels()");
    // Before making new levels, 
    // calculate averages to cells which were just unrefined.
    for (int lev = finest_level; lev > 0; lev--)
    {
        AmrLevel& fine_level = *amr_level[lev];
        AmrLevel& crse_level = *amr_level[lev-1];

        const Geometry& fine_geom = geom[lev];
        const Geometry& crse_geom = geom[lev-1];

        IntVect crse_ratio = ref_ratio[lev-1];
        MyInterpolater* mapper = &my_interpolater;

        if (new_finest < lev) {
            // average all data
            Print() << "aver all " << lev << "\n";
            for (AmrIter it(crse_level); it.isValid(); ++it)
            {
                int child = it.childLID();
                //Print() << "    " << it.validbox() << " " <<child << "\n";
                if (child < 0) continue;

                auto& crse = crse_level[it];
                auto& fine = fine_level.getCells().atLocalIdx(child);
                Box fine_box = amrex::refine(it.validbox(), crse_level.fineRatio());
                
                AllPrint() << "    " << fine_box << "\n";
                mapper->average(crse, 0, fine, 0, n_comp, fine_box, 
                                crse_ratio, crse_geom, fine_geom, RunOn::Cpu);
            }
        }
        else
        {
            Box domain = geom[lev].Domain();
            BoxArray coarse_boxes = amrex::complementIn(domain, new_grids[lev]);

            for (AmrIter it(crse_level); it.isValid(); ++it)
            {
                int child = it.childLID();
                if (child < 0) continue;

                auto& crse = crse_level[it];
                auto& fine = fine_level.getCells().atLocalIdx(child);
                Box fine_box = amrex::refine(it.validbox(), crse_level.fineRatio());
                BoxArray fine_region = amrex::intersect(coarse_boxes, fine_box);
                
                for (size_t i = 0; i < fine_region.size(); i++)
                {
                    mapper->average(crse, 0, fine, 0, n_comp, fine_region[i], 
                                    crse_ratio, crse_geom, fine_geom, RunOn::Cpu);
                }
            }
        }
    }
    

    for (int lev = 0; lev <= new_finest; lev++)
    {
        AmrLevel* a = new AmrLevel(*this, lev, Geom(lev), new_grids[lev], new_dmap[lev]);
        
        if (lev <= finest_level)
        {
            // Fill level from old levels
            AmrLevel* old = amr_level[lev].get();
            a->init(old);
        }
        else
        {
            // Fill level only from coarser level.
            a->init();
            a->FillFromCoarsePatch(a->getCells(), 0, n_comp, 0);
        }

        amr_level[lev].reset(a);
        this->SetBoxArray(lev, amr_level[lev]->boxArray());
        this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
    }
    
}

void 
Amr::ClearLevel (int lev) 
{ 
    amr_level[lev].reset();
    this->ClearBoxArray(lev);
    this->ClearDistributionMap(lev);
}

void
Amr::printGridInfo (std::ostream& os,
                    int           min_lev,
                    int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = amr_level[lev]->boxArray();
        int                       numgrid = bs.size();
        Long                      ncells  = amr_level[lev]->countCells();
        double                    ntot    = Geom(lev).Domain().d_numPts();
        Real                      frac    = 100.0_rt*(Real(ncells) / ntot);
        const DistributionMapping& map    = amr_level[lev]->DistributionMap();

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';


        for (int k = 0; k < numgrid; k++)
        {
            const Box& b = bs[k];

            os << ' ' << lev << ": " << b << "   ";
                
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                os << b.length(i) << ' ';

            os << ":: " << map[k] << '\n';
        }
    }

    os << std::endl; // Make sure we flush!
}


void
Amr::FillAllBoundaries() 
{    
    //Print() << "FillFineToCoarse\n";
    // average to coarser levels
    for (int i = finest_level; i >= 0; i--) {
        amr_level[i]->FillFineToCoarse();
    }

    //Print() << "FillBoundary\n";
    if(ParallelDescriptor::NProcs() > 1)
    for (int i = 0; i <= finest_level; i++)
    {
        CellFabArray* cells = &amr_level[i]->getCells();

        // update particle counts between remote neighbors
        Cell::transfer_particles = false;
        cells->FillBoundary(geom[i].periodicity(), true); // cross=true to not fill corners

        // get tags to cells which received data from remote cells
        const auto& receiveTags = cells->get_receive_tags();

        // resize ghost cells 
        for (auto const& kv : receiveTags)
        {
            for (auto const& tag : kv.second)
            {
                const auto& bx = tag.dbox;
                auto dfab = cells->array(tag.dstIndex);

                For( bx, cells->n_comp,
                [&] (int ii, int jj, int kk, int n) noexcept
                {
                    const IntVect idx = IntVect(ii,jj,kk);
                    dfab(idx, n)->resize();
                });
            }
        }
        
        // transmit particle data
        Cell::transfer_particles = true;
        cells->FillBoundary(geom[i].periodicity(), true);
    }
    
    //Print() << "FillCoarseToFine\n";
    // interpolate to fine levels
    for (int i = 0; i <= finest_level; i++) {
        amr_level[i]->FillCoarseToFine();
    }
    
}

void
Amr::MakeGrids (int lbase, int& new_finest,
                Vector<BoxArray>& new_grids)
{
    BL_PROFILE("Amr::MakeGrids()");

    const Real strttime = amrex::second();

    if (lbase == 0) {
	    new_grids[0] = MakeBaseGrids();
    }

    if (!initial_grids_file.empty() && !use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    // Use grids in initial_grids_file as fixed coarse grids.
    if ( !initial_grids_file.empty() && use_fixed_coarse_grids)
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,initial_ba.size());

        for (int lev = lbase+1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);

                if (lev > lbase)
                    bl.push_back(bx);

            }
            if (lev > lbase) new_grids[lev].define(bl);

            new_grids[lev].maxSize(max_grid_size[lev]);
        }
    }
    else if ( !regrid_grids_file.empty() )     // Use grids in regrid_grids_file 
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min<int>(new_finest,regrid_ba.size());
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = regrid_ba[lev-1].size();
            for (int i = 0; i < ngrid; i++)
            {
                Box bx(regrid_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }
    
    // Make new grids based on error estimates.
        Print() << "make new\n";
    MakeNewGrids(lbase, 0, new_finest, new_grids);
    
    if (verbose > 0)
    {
        Real stoptime = amrex::second() - strttime;
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	    amrex::Print() << "MakeGrids() time: " << stoptime << " new finest: " << new_finest<< '\n';
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Amr::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    /* NOTES: 
    * Grids are generated by iterating from current finest level to level 0
    * Be careful that fine areas aren't unrefined do to coarser cells marked as CLEAR!
    * For more see AmrMesh::MakeNewGrids() in AmrCore/AMReX_AmrMesh.cpp
    */
   Print() << "tag " << lev << "\n";
    amr_level[lev]->errorEst(tags, TagBox::CLEAR, TagBox::SET, time, n_error_buf[lev][0], ngrow);
}

BoxArray
Amr::GetAreaNotToTag (int lev)
{
    return BoxArray(amr_level[lev]->getAreaNotToTag());
}

void
Amr::ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev)
{
   Print() << "manual tag " << lev << "\n";
    amr_level[lev]->manual_tags_placement(tags, bf_lev);
}

const Vector<BoxArray>& Amr::getInitialBA() noexcept
{
  return initial_ba;
}
