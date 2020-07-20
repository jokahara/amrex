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
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FabArrayBase.H>

#include "CellFabArray.hpp"

bool Cell::transfer_particles = false;

using namespace std;
using namespace boost;
using namespace amrex;

void update_cell_lists(CellFabArray &grid, Geometry &geom) {

	// move particles to the particle list of the cell the particles are currently inside of
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		auto const& a = grid[mfi].array();
		auto box = mfi.validbox();

		LoopOnCpu(box, [&] (int i, int j, int k)
		{
			const IntVect cell(i,j,k);
			auto *previous_data = &a(cell);

			const IntVect neighbor_cell = cell - IntVect(1,0,0);
			auto *neighbor_data = &a(neighbor_cell);

			vector<std::array<double, 3>>::size_type n = 0;
			while (n < neighbor_data->particles.size()) {				
				const IntVect current_cell = 
					geom.CellIndex( neighbor_data->particles[n].data() );
				
				if (current_cell == cell) {
					auto *current_data = &a(current_cell);
					current_data->particles.push_back(neighbor_data->particles[n]);
					current_data->number_of_particles = current_data->particles.size();
				}
				n++;
			}
		});
	}

	// remove particles which don't belong to their current cells
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		auto const& a = grid[mfi].array();
		auto box = mfi.validbox();

		LoopOnCpu(box, [&] (int i, int j, int k)
		{
			const IntVect previous_cell(i,j,k);
			auto *previous_data = &a(previous_cell);

			vector<std::array<double, 3>>::size_type n = 0;
			while (n < previous_data->particles.size()) {

				const IntVect current_cell = 
					geom.CellIndex( previous_data->particles[n].data() );

				if (current_cell == previous_cell) {
					n++;
				}
				else {
					previous_data->particles.erase(previous_data->particles.begin() + n);
					previous_data->number_of_particles = previous_data->particles.size();
				}
			}
		});
	}
}

/*!
Saves the local grid and its particles to separate vtk files.

Each process saves its own files.
*/
void save(const int rank, CellFabArray &grid, unsigned int step)
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
	uint64_t total_particles = 0;
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		auto const& a = grid[mfi].array();
		For(mfi.validbox(), [&] (int i, int j, int k)
		{
			total_particles += a(i,j,k).particles.size();
		});
	}
	outfile << "POINTS " << total_particles << " float\n";
	// write out coordinates
	uint64_t written_particles = 0;
	
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		auto const& a = grid[mfi].array();
		For(mfi.validbox(), [&] (int i, int j, int k)
		{
			auto const particles = a(i,j,k).particles;
			for (int n = 0; n < a(i,j,k).particles.size(); n++)
			{
				outfile << particles[n][0] << "\t"
						<< particles[n][1] << "\t"
						<< particles[n][2] << "\n";
				written_particles++;
			}
		});
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
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		auto const& a = grid[mfi].array();
		For(mfi.validbox(), [&] (int i, int j, int k) {
			for (int n = 0; n < 3; n++)
			{
				outfile << mfi.index() << " ";
			}
		});
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
        IntVect dom_hi(AMREX_D_DECL(3,3,0));
        // box containing index space (cell centered by default)
        Box domain(dom_lo, dom_hi);
		
        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
		ba.maxSize(2);
        //ba.maxSize(max_grid_size);
		
       // This defines the physical box.
        RealBox real_box(AMREX_D_DECL(0,0,0),
                         AMREX_D_DECL(4,4,1));

        // periodic in x direction
        Array<int,AMREX_SPACEDIM> is_periodic{{1,0,0}};

        // This defines a Geometry object
        geom.define(domain, real_box, CoordSys::cartesian, is_periodic);
    }
	
    // Nghost = number of ghost cells for each array 
    IntVect Nghost(1,1,0);
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
	
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	dm.strategy(DistributionMapping::Strategy::KNAPSACK);

	CellFabArray grid;
	grid.define(ba, dm, Ncomp, Nghost);
	
	// initial condition
	const unsigned int max_particles_per_cell = 5;
	for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
		// This is the valid Box of the current FArrayBox.
		// By "valid", we mean the original ungrown Box in BoxArray.
		const Box& box = mfi.validbox();

		// A reference to the current FArrayBox in this loop iteration.
		auto& fab = grid[mfi];

		// Obtain Array4 from CellFab.
		// Array4<Cell> const& a = grid.array(mfi);
		auto const& a = fab.array();
		
		Dim3 lo = lbound(box);
		Dim3 hi = ubound(box);
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
	std::cerr << rank << " initialized " << grid.local_size() << " box(es)" << "\n";

	for (MFIter mfi(grid); mfi.isValid(); ++mfi)
	{
		int non_local = 0;
		auto box = mfi.validbox();

		std::cerr << "Process " << rank << " initialized box: " << box << "\n";
	}
	
	ParallelDescriptor::Barrier();
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

		const double vx = 0.1;

		// move particles in x-direction
		for (MFIter mfi(grid); mfi.isValid(); ++mfi) { // loop over local boxes
			auto const& a = grid[mfi].array();
			auto box = mfi.validbox();
			
			// loop over cells
			LoopOnCpu(box, [&] (int i, int j, int k) 
			{	
				// loop over particles
				for (int n = 0; n < a(i,j,k).particles.size(); n++) 
				{
					auto &particle = a(i,j,k).particles[n];
					particle[0] += vx;

					// hande grid wrap around
					if(particle[0] >= geom.period(0)) {
						particle[0] -= geom.period(0);
					}
				}
			});
		}

		// update particle counts between neighboring cells
		Print() << "update counts: \n";
		Cell::transfer_particles = false;
		grid.FillBoundary(geom.periodicity(), true); // cross=true to not fill corners
		
		// resize ghost cells
    	const auto receiveTags = *grid.getFB(grid.nGrowVect(), geom.periodicity(), true, false).m_RcvTags;
    	for (auto const& kv : receiveTags)
		{
			for (auto const& tag : kv.second)
			{
				const auto& bx = tag.dbox;
				auto dfab = grid.array(tag.dstIndex);
				amrex::Loop( bx, Ncomp,
				[&] (int ii, int jj, int kk, int n) noexcept
				{
					const IntVect idx = IntVect(ii,jj,kk);
					dfab(idx, n).resize(); // resizing
				});
			}
		}

		Print() << "update data: \n";
		// update particle data between neighboring cells
		Cell::transfer_particles = true;
		grid.FillBoundary(geom.periodicity(), true);

		Print() << "done: \n";
		update_cell_lists(grid, geom);
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
	//MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN); 

    amrex::Initialize(argc, argv, false,  MPI_COMM_WORLD);
	
	main_main();

	amrex::Finalize();
	MPI_Finalize();
}