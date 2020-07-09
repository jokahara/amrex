
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "AmrGrid.hpp"

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrGrid::AmrGrid ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t.resize(nlevs_max, 0.0);
    dt.resize(nlevs_max, 1.e100);

    data.resize(nlevs_max);
    bcs.resize(1);
    
    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    
/*
    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
*/
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }
}

AmrGrid::~AmrGrid ()
{
}

/* advance solution to final time
void
AmrGrid::Evolve ()
{
    Real cur_time = time[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

	// sync up time
	for (lev = 0; lev <= finest_level; ++lev) {
	    t[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
	}

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	WritePlotFile();
    }
}
*/

// initializes multilevel data
void
AmrGrid::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        //AverageDown();

        if (chk_int > 0) {
            //WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        //ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrGrid::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				    const DistributionMapping& dm)
{
    const int ncomp = data[lev-1].nComp();
    const int nghost = data[lev-1].nGrow();

    data[lev].define(ba, dm, ncomp, nghost);
    //phi_old[lev].define(ba, dm, ncomp, nghost);

    t[lev] = time;
    //t_old[lev] = time - 1.e200;

    FillCoarsePatch(lev, time, data[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrGrid::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    const int ncomp = data[lev].nComp();
    const int nghost = data[lev].nGrow();

    CellFabArray new_state(ba, dm, ncomp, nghost);
    //CellFabArray old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, data[lev]);
    //std::swap(old_state, phi_old[lev]);

    t[lev] = time;
    //t_old[lev] = time - 1.e200;
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrGrid::ClearLevel (int lev)
{
    data[lev].clear();
    //phi_old[lev].clear();
    //flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrGrid::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;

    data[lev].define(ba, dm, ncomp, nghost);
    //phi_old[lev].define(ba, dm, ncomp, nghost);

    t[lev] = time;
    //t_old[lev] = time - 1.e200;

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t[lev];

    CellFabArray& state = data[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();
	    //TODO: initdata
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrGrid::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const CellFabArray& state = data[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            auto tagArray  = tags[mfi].array();
            auto dataArray = state[mfi].array();

            //tag cells for refinement
            ParallelFor(box, [&] (int i, int j, int k)
            {
                if (dataArray(i,j,k).number_of_particles >= 4) {
                    tagArray(i,j,k) = TagBox::CLEAR;
                }
                else {
                    tagArray(i,j,k) = TagBox::SET;
                }
            });
        }
    }
}

// read in some parameters from inputs file
void
AmrGrid::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);
	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);
	pp.query("chk_file", chk_file);
	pp.query("chk_int", chk_int);
    pp.query("restart",restart_chkfile);
    }
}
/*
// set covered coarse cells to be the average of overlying fine cells
void
AmrGrid::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(data[lev+1], data[lev],
                            geom[lev+1], geom[lev],
                            0, data[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrGrid::AverageDownTo (int crse_lev)
{
    amrex::average_down(data[crse_lev+1], data[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, data[crse_lev].nComp(), refRatio(crse_lev));
}*/

// compute a new CellFabArray by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrGrid::FillPatch (int lev, Real time, CellFabArray& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        Vector<CellFabArray*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        BndryFuncArray bfunc(phifill); // TODO: lambda
        PhysBCFunct<BndryFuncArray> physbc(geom[lev], bcs, bfunc);
        
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
    }
    else
    {
        Vector<CellFabArray*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                0, icomp, ncomp, geom[lev-1], geom[lev],
                                cphysbc, 0, fphysbc, 0,
                                refRatio(lev-1), mapper, bcs, 0);
    }
}

// fill an entire CellFabArray by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrGrid::FillCoarsePatch (int lev, Real time, CellFabArray& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<CellFabArray*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
	    amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    BndryFuncArray bfunc(phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, 0, fphysbc, 0, refRatio(lev-1),
				 mapper, bcs, 0);
}

// utility to copy in data into another CellFabArray
void
AmrGrid::GetData (int lev, Real time, Vector<CellFabArray*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = 1.e-3;

    if (time > t[lev] - teps && time < t[lev] + teps)
    {
        data.push_back(&this->data[lev]);
        datatime.push_back(t[lev]);
    }
    else
    {
        data.push_back(&this->data[lev]);
        datatime.push_back(t[lev]);
    }
}


// advance a level by dt
// includes a recursive call for finer levels
void
AmrGrid::timeStep (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                        // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                        // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }
	        }
	    }
    }

    if (Verbose()) {
	amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
	amrex::Print() << "ADVANCE with time = " << t[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    //TODO: advance a single level for a single time step

    ++istep[lev];

    if (Verbose())
    {
	amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
            // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }

}

void
AmrGrid::ComputeDt ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt[lev] = 0.1;
    }
}


// get plotfile name
std::string
AmrGrid::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of CellFabArrays for writing
Vector<const CellFabArray*>
AmrGrid::PlotFileMF () const
{
    Vector<const CellFabArray*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(&data[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrGrid::PlotFileVarNames () const
{
    return {"phi"};
}

// write plotfile to disk
void
AmrGrid::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    // TODO
}
