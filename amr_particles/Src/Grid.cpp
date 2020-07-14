
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "Grid.hpp"

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
GridLevel::GridLevel ()
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

    // TODO: Boundary conditions
    /*bcs.resize(1);
    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    
    // walls (Neumann)
    // int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    // int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};

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
    }*/
}

GridLevel::~GridLevel ()
{
}

// initializes multilevel data
void
GridLevel::InitData ()
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
GridLevel::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				                 const DistributionMapping& dm)
{
    Print() << "Making new level: " << lev <<"\n";
    const int ncomp = data[lev-1].nComp();
    const int nghost = data[lev-1].nGrow();
    
    data[lev].define(ba, dm, ncomp, nghost);
    t[lev] = time;

    FillCoarsePatch(lev, data[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
GridLevel::RemakeLevel (int lev, Real time, const BoxArray& ba,
			          const DistributionMapping& dm)
{
    Print() << "Remaking Level: " << lev << "\n";
    const int ncomp = data[lev].nComp();
    const int nghost = data[lev].nGrow();
    
    CellFabArray new_state(ba, dm, ncomp, nghost);

    FillPatch(lev, new_state, 0, ncomp);

    //std::swap(new_state, data[lev]);

    t[lev] = time;
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
GridLevel::ClearLevel (int lev)
{
    Print() << "Clearing Level: " << lev << "\n";
    data[lev].clear();
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void GridLevel::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{
    Print() << "Making new level from scratch: " << lev << "\n";

    const int ncomp = 1;
    const int nghost = 1;

    data[lev].define(ba, dm, ncomp, nghost);
    
    t[lev] = time;
    
    if (lev == 0)
    for (MFIter mfi(data[lev]); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();
        auto a = data[lev][mfi].array();
        
		Loop(box, [=] (int i, int j, int k) 
		{
			Cell* cell_data = &a(i,j,k);
			const unsigned int number_of_particles = i+1;
				//= (unsigned int)ceil(max_particles_per_cell * double(rand()) / RAND_MAX);
			for (unsigned int n = 0; n < number_of_particles; n++) {
				std::array<double, 3> coordinates = {{
					i + double(rand()) / RAND_MAX,
					j + double(rand()) / RAND_MAX,
					k + double(rand()) / RAND_MAX
				}};

				cell_data->particles.push_back(coordinates);
				cell_data->number_of_particles = cell_data->particles.size();
			}
		});
    }
    else FillCoarsePatch(lev, data[lev], 0, ncomp);
}

// 
/**
 * tag all cells for refinement
 * overrides the pure virtual function in AmrCore
 * Tag values: CLEAR, BUF, SET.
 **/
void
GridLevel::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    Print() << "Tagging cells at level: " << lev << "\n";
    static bool first = true;
    
    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const CellFabArray& state = data[lev];
    int tagged = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            auto tagArray  = tags[mfi].array();
            auto dataArray = state[mfi].array();

            //tag cells for refinement
            ParallelFor(box, [&] (int i, int j, int k)
            {
                if (dataArray(i,j,k).number_of_particles >= 4) {
                    tagArray(i,j,k) = TagBox::SET;
                    tagged++;
                }
                else {
                    tagArray(i,j,k) = TagBox::CLEAR;
                }
            });
        }
    }
    Print() << "Tagged cells: " << tagged << "\n";
}

// read in some parameters from inputs file
void
GridLevel::ReadParameters ()
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

// compute a new CellFabArray by coping from valid region and filling ghost cells
void
GridLevel::FillPatch (int lev, CellFabArray& mf, int icomp, int ncomp)
{
    Print() << "FillPatch needed: " << lev << "\n";
    return;
    if (lev == 0)
    {
        CellFabArray& src = data[0];
        mf.FillPatchSingleLevel(src, 0, icomp, ncomp, geom[lev]);
    }
    else
    {
        CellFabArray &coarse = data[lev-1],
                     &fine   = data[lev];

        mf.FillPatchTwoLevels(coarse, fine, 0, icomp, ncomp, 
                              geom[lev-1], geom[lev], refRatio(lev-1));
    }
}

// fill an entire CellFabArray by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
GridLevel::FillCoarsePatch (int lev, CellFabArray& fine, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    CellFabArray& coarse = data[lev-1];

    Print() << "Interpolate to: " << lev << "\n";
    
    const BoxArray& ba = fine.boxArray();
    const DistributionMapping& dm = fine.DistributionMap();

    Geometry &fine_geom = Geom(lev);
    Geometry &coarse_geom = Geom(lev-1);

    IntVect nghost{0,0,0};              // not adding ghost data here
    IntVect ratio = refRatio(lev-1);

    const Box &fine_domain = fine_geom.Domain();

    // make an array of coarsened boxes corresponding to fine boxes
    BoxArray ba_crse_patch(ba.size());
    for (int i = 0, N = ba.size(); i < N; ++i)
    {
        Box bx = amrex::grow(ba[i], nghost);
        bx &= fine_domain;

        ba_crse_patch.set(i, bx.coarsen(ratio));
    }


/*#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif*/
    for (MFIter mfi(ba_crse_patch, dm); mfi.isValid(); ++mfi)
    {
        auto& sfab = coarse[mfi];       // source data
        auto& dfab = fine[mfi];         // destination data
        const Box dbox = dfab.box();    // this box is grown by 1
        const Box sbox = mfi.validbox();

        // copying particles to fine cells
        For(sbox, [&] (int i, int j, int k) 
        {
            const IntVect iv(i,j,k);
            Cell &parent = sfab(iv);
            
            for (uint i = 0; i < parent.number_of_particles; i++)
            {
                auto& particle = parent.particles[i];
                // get the new location corresponding to fine geometry
                IntVect new_iv = fine_geom.CellIndex(particle.data());

                dfab(new_iv).particles.push_back(particle);
                dfab(new_iv).number_of_particles++;
            }

            // clear for now
            parent.particles.clear();
            parent.number_of_particles = 0;
        });
        //interp
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void
GridLevel::TimeStep (int lev, Real time, int iteration)
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
            TimeStep(lev+1, time+(i-1)*dt[lev+1], i);
        }

        //AverageDownTo(lev); // average lev+1 down to lev
    }

}

void
GridLevel::ComputeDt ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt[lev] = 0.1;
    }
}


// get plotfile name
std::string
GridLevel::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// set plotfile variable names
Vector<std::string>
GridLevel::PlotFileVarNames () const
{
    return {"plt"};
}

// write plotfile to disk
void
GridLevel::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    // TODO
}

