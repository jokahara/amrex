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
#include "StateData.h"

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

Vector<std::unique_ptr<AmrLevel> >&
Amr::getAmrLevels () noexcept
{
    return amr_level;
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

    int nlev = max_level+1;
    amr_level.resize(nlev);

    // Make the default regrid_int = 1 for all levels.
    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (int i = 0; i < max_level; i++)
           regrid_int[i]  = 1;
    }
    
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
Amr::init (const BoxArray* lev0_grids, const Vector<int>* pmap)
{
    BL_PROFILE("Amr::initialInit()");
    
    ParmParse pp;

    n_comp = 1;
    pp.query("n_comp", n_comp);

    int grow = 1;
    pp.query("n_grow", grow);
    n_grow = {grow, grow, grow};

    if (check_input) checkInput();
    
    // Generate internal values from user-supplied values.
    finest_level = 0;

    // Define base level grids.
    defBaseLevel(0, lev0_grids, pmap);

    // Build fine level grids.
    if (max_level > 0) bldFineLevels(0);

    for (int lev = 0; lev <= finest_level; lev++) {
        amr_level[lev]->post_regrid(0,finest_level);
    }
}

void
Amr::timeStep (int  level,
               Real time,
               int  iteration,
               int  niter,
               Real stop_time)
{
    BL_PROFILE("Amr::timeStep()");
    BL_COMM_PROFILE_NAMETAG("Amr::timeStep TOP");

    // This is used so that the AmrLevel functions can know which level is being advanced 
    //      when regridding is called with possible lbase > level.
    which_level_being_advanced = level;

    int lev_top = std::min(finest_level, max_level-1);
    for (int i(level); i <= lev_top; ++i)
    {
        const int old_finest = finest_level;

        regrid(i,time);

        if (old_finest > finest_level) 
            lev_top = std::min(finest_level, max_level - 1);
    }

    //
    // Load Balance
    //
    if (max_level == 0 && loadbalance_level0_int > 0 && loadbalance_with_workestimates)
    {
        LoadBalanceLevel0(time);
    }

    // Advance grids at this level.
    // amr_level[level]->advance;

    // maybe regrid

    // Advance grids at higher level.
    if (level < finest_level)
    {
        const int lev_fine = level+1;
        timeStep(lev_fine,time,1,1,stop_time);
    }

    //post_timestep

    // Set this back to negative so we know whether we are in fact in this routine
    which_level_being_advanced = -1;
}

void
Amr::defBaseLevel (Real              strt_time, 
                   const BoxArray*   lev0_grids,
                   const Vector<int>* pmap)
{
    BL_PROFILE("Amr::defBaseLevel()");
    // Just initialize this here for the heck of it
    which_level_being_advanced = -1;

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
        Print() << "    MakeBaseGrids\n";
	    lev0 = MakeBaseGrids();
    }

    this->SetBoxArray(0, lev0);
    this->SetDistributionMap(0, DistributionMapping(lev0));

    // Now build level 0 grids.
    amr_level[0].reset(new AmrLevel(*this, 0, Geom(0), grids[0], dmap[0], strt_time));
}

void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{
    BL_PROFILE("Amr::regrid()");

    if (lbase > std::min(finest_level,max_level-1)) return;

    if (verbose > 0)
	amrex::Print() << "Now regridding at level lbase = " << lbase << "\n";

    //
    // Compute positions of new grids.
    //
    int                         new_finest;
    Vector<BoxArray>            new_grid_places(max_level+1);
    Vector<DistributionMapping> new_dmap(max_level+1);

    grid_places(lbase,time,new_finest, new_grid_places);

    bool regrid_level_zero = (!initial) && (lbase == 0)
        && ( loadbalance_with_workestimates || (new_grid_places[0] != amr_level[0]->boxArray()));

    const int start = regrid_level_zero ? 0 : lbase+1;

    bool grids_unchanged = finest_level == new_finest;
    for (int lev = start, End = std::min(finest_level,new_finest); lev <= End; lev++) {
        if (new_grid_places[lev] == amr_level[lev]->boxArray()) {
            new_grid_places[lev] = amr_level[lev]->boxArray();  // to avoid duplicates
            new_dmap[lev] = amr_level[lev]->DistributionMap(); 
        } else {
            grids_unchanged = false;
        }
    }

    //
    // Reclaim all remaining storage for levels > new_finest.
    //
    for(int lev = new_finest + 1; lev <= finest_level; ++lev) {
        amr_level[lev].reset();
        this->ClearBoxArray(lev);
        this->ClearDistributionMap(lev);
    }

    finest_level = new_finest;

    //
    // Define the new grids from level start up to new_finest.
    //
    for(int lev = start; lev <= new_finest; ++lev) {

        // Construct skeleton of new level.
        if (loadbalance_with_workestimates && !initial) {
            new_dmap[lev] = makeLoadBalanceDistributionMap(lev, time, new_grid_places[lev]);
        }
        else if (new_dmap[lev].empty()) {
	        new_dmap[lev].define(new_grid_places[lev]);
	    }

        AmrLevel* a = new AmrLevel(*this, lev, Geom(lev), new_grid_places[lev], new_dmap[lev]);

        amr_level[lev].reset(a);
        this->SetBoxArray(lev, amr_level[lev]->boxArray());
        this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
    }

    // Check at *all* levels whether we need to do anything special now that the grids
    // at levels lbase+1 and higher may have changed.  
    for(int lev(0); lev <= new_finest; ++lev) {
        amr_level[lev]->post_regrid(lbase,new_finest);
    }

    // Report creation of new grids.
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << " : REGRID  with lbase = "
                       << lbase
                       << std::endl;

        printGridSummary(amrex::OutStream(),start,finest_level);
    }
}

DistributionMapping
Amr::makeLoadBalanceDistributionMap (int lev, Real time, const BoxArray& ba) const
{
    BL_PROFILE("makeLoadBalanceDistributionMap()");

    if (verbose) {
        amrex::Print() << "Load balance on level " << lev << " at t = " << time << "\n";
    }

    DistributionMapping newdm;

    if (amr_level[lev])
    {
        DistributionMapping dmtmp;
        if (ba.size() == boxArray(lev).size()) {
            dmtmp = DistributionMap(lev);
        } else {
            dmtmp.define(ba);
        }
        amrex::Print() << "    Calculate work estimates\n";
        // TODO: Fill multifab with weights for dist mapping
        MultiFab workest(ba, dmtmp, 1, 0, MFInfo(), FArrayBoxFactory());
        amr_level[lev]->estimateWork(workest);

        //Real navg = static_cast<Real>(ba.size()) / static_cast<Real>(ParallelDescriptor::NProcs());
        //int nmax = static_cast<int>(std::max(std::round(loadbalance_max_fac*navg), std::ceil(navg)));

        newdm = DistributionMapping::makeKnapSack(workest);
    }
    else
    {
        newdm.define(ba);
    }

    return newdm;
}

void
Amr::LoadBalanceLevel0 (Real time)
{
    BL_PROFILE("LoadBalanceLevel0()");
    const auto& dm = makeLoadBalanceDistributionMap(0, time, boxArray(0));
    InstallNewDistributionMap(0, dm);
    amr_level[0]->post_regrid(0,0);
}

void
Amr::InstallNewDistributionMap (int lev, const DistributionMapping& newdm)
{
    BL_PROFILE("InstallNewDistributionMap()");

    AmrLevel* a = new AmrLevel(*this, lev, Geom(lev), boxArray(lev), newdm);

    AmrLevel* old = amr_level[lev].get();

    CellFabArray& S_new = a->getData(); // should be already defined

    if (lev == 0)
    {
        Print() << "    FillFromLevel0\n";
        // TODO: maybe move this to FillPatch later
        S_new.FillPatchSingleLevel({0,0,0}, old->getData(), 0, 0, nComp(), Geom(lev));
    } else
    {
        AmrLevel::FillPatch(*old, S_new, 0, 0, 0, S_new.nComp()); // Fill S_new with data in old
    }
    
    amr_level[lev].reset(a);

    this->SetBoxArray(lev, amr_level[lev]->boxArray());
    this->SetDistributionMap(lev, amr_level[lev]->DistributionMap());
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
Amr::grid_places (int              lbase,
                  Real             time,
                  int&             new_finest,
                  Vector<BoxArray>& new_grids)
{
    BL_PROFILE("Amr::grid_places()");
    Print() << "    grid_places\n";

    const Real strttime = amrex::second();

    if (lbase == 0) {
        Print() << "    MakeBaseGrids\n";
	    new_grids[0] = MakeBaseGrids();
    }

    if ( time == 0. && !initial_grids_file.empty() && !use_fixed_coarse_grids)
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
    if ( ! initial_grids_file.empty() && use_fixed_coarse_grids)
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
            if (lev > lbase)
                new_grids[lev].define(bl);
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
    
    Print() << "    MakeNewGrids\n";
    MakeNewGrids(lbase, time, new_finest, new_grids);

    if (verbose > 0)
    {
        Real stoptime = amrex::second() - strttime;
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	    amrex::Print() << "grid_places() time: " << stoptime << " new finest: " << new_finest<< '\n';
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Amr::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
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
    amr_level[lev]->manual_tags_placement(tags, bf_lev);
}

void
Amr::bldFineLevels (Real strt_time)
{
    BL_PROFILE("Amr::bldFineLevels()");
    finest_level = 0;

    Vector<BoxArray> new_grids(max_level+1);
    //
    // Get initial grid placement.
    //
    do
    {
        int new_finest;

        grid_places(finest_level,strt_time,new_finest,new_grids);

        if (new_finest <= finest_level) break;
        //
        // Create a new level and link with others.
        //
        finest_level = new_finest;

	    DistributionMapping new_dm {new_grids[new_finest]};

        AmrLevel* level = new AmrLevel(*this, new_finest, Geom(new_finest),
                                       new_grids[new_finest], new_dm, strt_time);

        amr_level[new_finest].reset(level);
        this->SetBoxArray(new_finest, new_grids[new_finest]);
        this->SetDistributionMap(new_finest, new_dm);
    }
    while (finest_level < max_level);
    //
    // Iterate grids to ensure fine grids encompass all interesting gunk.
    // but only iterate if we did not provide a grids file.
    //
    if ( regrid_grids_file.empty() || (strt_time == 0.0 && !initial_grids_file.empty()) )  
    {
        bool grids_the_same;
        const int MaxCnt = 4;
        int count = 0;

        do
        {
            for (int i = 0; i <= finest_level; i++) {
                new_grids[i] = amr_level[i]->boxArray();
            }

            regrid(0,strt_time,true);

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

const Vector<BoxArray>& Amr::getInitialBA() noexcept
{
  return initial_ba;
}
