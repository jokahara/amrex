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

#include "Src/Amr.h"

bool Cell::transfer_particles = false;

using namespace std;
using namespace boost;
using namespace amrex;

void update_cell_lists(CellArray &grid, Geometry &geom) {

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
void save(const int rank, AmrGrid &amrGrid, unsigned int step)
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
	for (int lev = 0; lev <= amrGrid.maxLevel(); lev++)
	{
		auto& grid = amrGrid[lev];
		for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k)
			{
				total_particles += a(i,j,k).particles.size();
			});
		}
	}
	
	outfile << "POINTS " << total_particles << " float\n";
	// write out coordinates
	uint64_t written_particles = 0;
	
	for (int lev = 0; lev <= amrGrid.maxLevel(); lev++)
	{
		auto& grid = amrGrid[lev];
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

	outfile << "\nSCALARS level int 1\nLOOKUP_TABLE default\n";
	for (int lev = 0; lev <= amrGrid.maxLevel(); lev++)
	{
		auto& grid = amrGrid[lev];
		for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k) {
				for (int n = 0; n < a(i,j,k).number_of_particles; n++)
				{
					outfile << lev << " ";
				}
			});
		}
	}

	// cell numbers of particles
	outfile << "\nSCALARS cell int 1\nLOOKUP_TABLE default\n";
	for (int lev = 0; lev <= amrGrid.maxLevel(); lev++)
	{
		auto& grid = amrGrid[lev];
		for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k) {
				for (int n = 0; n < 3; n++)
				{
					outfile << mfi.index() << " ";
				}
			});
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
	
	AmrGrid amrGrid;
	amrex::Print() << amrGrid;

	amrGrid.InitData();

	for (int level = 0; level <= amrGrid.maxLevel(); level++)
	{
		std::cout << "Level " << level << ":\n"; 
		auto &grid = amrGrid[level];

		for (MFIter mfi(grid); mfi.isValid(); ++mfi) {
			const Box& box = mfi.validbox();
			std::cout << "  Box: " << box << " ";
			auto a = grid[mfi].array();
			int n = 0;
			LoopOnCpu(box, [&] (int i, int j, int k) 
			{	
				n += a(i,j,k).number_of_particles;
			});

			std::cout << "has " << n << " particles\n";

			Real loc[3];
			amrGrid.Geom(level).CellCenter(IntVect(box.loVect()), loc);
			//std::cout << "  Real low: " << loc[0] << " " << loc[1] << " " << loc[2] << ":\n"; 

			amrGrid.Geom(level).CellCenter(IntVect(box.hiVect()), loc);
			//std::cout << "  Real high: " << loc[0] << " " << loc[1] << " " << loc[2] << ":\n"; 
		}
	}
	
	
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
	
	const unsigned int max_steps = 1;
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

		save(rank, amrGrid, step);

		const double vx = 0.1;

		auto &grid = amrGrid[0];
		int Ncomp = grid.nComp();
		Geometry &geom = amrGrid.Geom(0);
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
		Cell::transfer_particles = false;
		grid.FillBoundary(geom.periodicity(), true); // cross=true to not fill corners
		//ParallelDescriptor::Barrier();

		// resize ghost cells (not working currently with MPI)
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

		// update particle data between neighboring cells
		Cell::transfer_particles = true;
		grid.FillBoundary(geom.periodicity(), true);

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

	save(rank, amrGrid, max_steps);

	clock_t after = clock();
	cout << "Process " << rank << ": " << comm_size << " processes in total"
		<< ": simulation took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds "
		<< endl;
	
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	// argv must include inputs file
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

    amrex::Initialize(argc, argv, true,  MPI_COMM_WORLD);
	
	main_main();

	amrex::Finalize();
	MPI_Finalize();
}