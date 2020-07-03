/*
Particle propagator for dccrg.

Copyright 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/

#include "algorithm"
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "iomanip"
#include "iostream"
#include "string"
#include "ctime"
#include "mpi.h"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>

#include "CellFabArray.hpp"

bool Cell::transfer_particles = false;

using namespace std;
using namespace boost;
using namespace amrex;

using PartIter = ParIter<1,1,0,0>;

/*!
Propagates particles between local cells and to/from remote cells into local cells.

The velocity field is constant in space and time and only along x direction.
*/
void propagate_particles(ParticleContainer<1,1,0,0> &grid) {

	const double vx = 0.1;
	// propagate particles in local cells and copies of remote neighbors
	for (PartIter pti(grid, 0); pti.isValid(); ++pti)
	{
		auto &particles = pti.GetArrayOfStructs();
		for (auto &p : particles )
		{
			p.pos(0) += vx;
		}
	}
}


/*!
Saves the local grid and its particles to separate vtk files.

Each process saves its own files.
*/
void save(const int rank, ParticleContainer<1,1,0,0> &grid, unsigned int step)
{
	// write the grid
	const string grid_file_name(
		lexical_cast<string>("simple_")
		+ lexical_cast<string>(rank) + "_"
		+ lexical_cast<string>(step) + "_grid.vtk"
	);

	// write the particles
	const string outname(
		lexical_cast<string>("simple_")
		+ lexical_cast<string>(rank) + "_"
		+ lexical_cast<string>(step) + ".vtk"
	);
	
	ofstream outfile;
	outfile.open(outname.c_str());

	// write header
	outfile << "# vtk DataFile Version 2.0\n"
		"Particle test\nASCII\nDATASET UNSTRUCTURED_GRID\n";

	// calculate the total number local of particles
	uint64_t total_particles = grid.TotalNumberOfParticles(true, true);
	outfile << "POINTS " << total_particles << " float\n";

	// write out coordinates
	uint64_t written_particles = 0;
	
	for (PartIter pti(grid, 0); pti.isValid(); ++pti)
	{
		const auto &particles = pti.GetArrayOfStructs();
		for (auto p : particles )
		{
			outfile << p.pos(0) << "\t"
				<< p.pos(1) << "\t"
				<< p.pos(2) << "\n";
			written_particles++;
		}
	}

	if (written_particles != total_particles) {
		cerr << __FILE__ << ":" << __LINE__
			<< ": Incorrect number of particles written: " << written_particles
			<< ", should be " << total_particles
			<< endl;
		abort();
	}

	outfile << "CELLS " << total_particles << " " << 2 * total_particles << "\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << "1 " << i << "\n";
	}

	outfile << "CELL_TYPES " << total_particles << "\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << "1 ";
	}

	// process numbers of particles
	outfile << "\nCELL_DATA " << total_particles
		<< "\nSCALARS process int 1\nLOOKUP_TABLE default\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << rank << " ";
	}
	outfile << "\n";

	// cell numbers of particles
	outfile << "SCALARS cell int 1\nLOOKUP_TABLE default\n";
	for (PartIter pti(grid, 0); pti.isValid(); ++pti)
	{
		const auto &particles = pti.GetArrayOfStructs();
		for (int i = 0; i < 3; i++)
		{
			outfile << pti.index() << " ";
		}
	}
	outfile << "\n";

	outfile.close();
}


int main_main()
{
	clock_t before = clock();
	amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

	int rank = ParallelDescriptor::MyProc();
	int comm_size = ParallelDescriptor::NProcs();
	
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(0,0,0));
        IntVect dom_hi(AMREX_D_DECL(9,9,0));
        // box containing index space (cell centered by default)
        Box domain(dom_lo, dom_hi);
		
        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
		ba.maxSize(3);
        //ba.maxSize(max_grid_size);
		
       // This defines the physical box.
        RealBox real_box(AMREX_D_DECL(0,0,0),
                         AMREX_D_DECL(10,10,1));

        // periodic in x direction
        Array<int,AMREX_SPACEDIM> is_periodic{{1,1,1}};

        // This defines a Geometry object
        geom.define(domain, real_box, CoordSys::cartesian, is_periodic);
    }
	
    // Nghost = number of ghost cells for each array 
    int Nghost = 1;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
	
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	dm.strategy(DistributionMapping::Strategy::KNAPSACK);
	
	ParticleContainer<1,1,0,0> grid;
	grid.Define(geom, dm, ba);
	ParticleInitType<1,1,0,0> pdata = {};
	grid.InitNRandomPerCell(1, pdata);
	grid.Redistribute();
	
	// value in the ghost box is the same as the real value on the opposite side
	cout << rank << " has particles: " << grid.NumberOfParticlesAtLevel(0, true, true) << endl;
	//BaseFab<Cell> fab;
	
	//FabArray< BaseFab<Cell> > grid(ba, dm, ...);
	
	//grid.FillBoundary(geom.periodicity());
	/*
	AmrCore amr;
	amr.SetGeometry(0, geom);
	amr.SetBoxArray(0, ba);
	amr.SetDistributionMap(0, dm);
	amrex::Print() << amr.CountCells() << "\n";*/

	/*
	Visualize the results for example with visit -o simple_particles.visit
	or visit -o simple_grid.visit or overlay them both from the user interface.
	*/
	const string visit_particles_name("simple_particles.visit"),
		visit_grid_name("simple_grid.visit");

	ofstream visit_particles, visit_grid;

	if (rank == 0) {
		visit_particles.open(visit_particles_name.c_str());
		visit_particles << "!NBLOCKS " << comm_size << "\n";
		visit_grid.open(visit_grid_name.c_str());
		visit_grid << "!NBLOCKS " << comm_size << "\n";
	}

	const unsigned int max_steps = 50;
	for (unsigned int step = 0; step < max_steps; step++) {

		// append current output file names to the visit files
		if (rank == 0) {
			for (int i = 0; i < comm_size; i++) {
				visit_particles << "simple_"
					<< i << "_"
					<< step << ".vtk\n";
				visit_grid << "simple_"
					<< i << "_"
					<< step << "_grid.vtk\n";
			}
		}

		save(rank, grid, step);

		// update particle data between neighboring cells on different processes
		grid.Redistribute(0, -1, 0, 1);
		propagate_particles(grid);
	}

	// append final output file names to the visit files
	if (rank == 0) {
		for (int i = 0; i < comm_size; i++) {
				visit_particles << "simple_"
					<< rank << "_"
					<< max_steps << ".vtk\n";
				visit_grid << "simple_"
					<< rank << "_"
					<< max_steps << "_grid.vtk\n";
		}
		visit_particles.close();
		visit_grid.close();
	}

	save(rank, grid, max_steps);

	clock_t after = clock();
	cout << "Process " << rank << ": " << comm_size << " processes in total"
		<< ": simulation took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds "
		<< endl;
	
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}
	
    amrex::Initialize(argc, argv, false,  MPI_COMM_WORLD);
	
	main_main();

	amrex::Finalize();
	MPI_Finalize();
}