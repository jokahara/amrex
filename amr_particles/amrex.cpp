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
#include "boost/lexical_cast.hpp"
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

void update_cell_lists(AmrLevel& level, Geometry &geom) 
{
	CellFabArray &grid = level.getCells();
	//int added=0;
	// move particles to the particle list of the cell the particles are currently inside of
	for (AmrIter mfi(level); mfi.isValid(); ++mfi) {
		if (!mfi.isPropagated()) continue;

		auto const& a = grid[mfi].array();
		auto box = mfi.validbox();

		For(box.growHi(0,1), [&] (int i, int j, int k)
		{
			const IntVect cell(i,j,k);
			const IntVect neighbor_cell = cell - IntVect(1,0,0);
			auto& neighbor_data = a(neighbor_cell);

			int n = 0;
			while (n < neighbor_data->particles.size()) {				
				const IntVect current_cell = 
					geom.CellIndex( neighbor_data->particles[n].data() );
				
				if (current_cell == cell) {
					auto& current_data = a(current_cell);
					current_data->particles.push_back(neighbor_data->particles[n]);
					current_data->number_of_particles = current_data->particles.size();
					//added++;
					
					neighbor_data->particles.erase(neighbor_data->particles.begin() + n);
					neighbor_data->number_of_particles = neighbor_data->particles.size();
				} 
				else n++;
			}
		});
	}
	//Print() << "moved " << added << "\n";
	
}

void printCounts(Amr& amr) 
{
	sleep(ParallelDescriptor::MyProc());
	int total = 0;
	int nghost = 0;
	int count = 0;
	
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		AllPrint() << "\n\tRANK " << ParallelDescriptor::MyProc() << ", "
				   << "LEVEL " << lev << "\n";

		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) 
		{
			//if (!mfi.isPropagated()) continue;
			const auto& fab = amr[lev].getCells()[mfi];
			int line = 0;
			For(mfi.validbox().grow(IntVect{1,1,0}), amr.nComp(), [&] (int i, int j, int k, int n) 
			{
				if(j != line) {
					line = j;
					AllPrint() << "\n";
				}
				auto& p = fab(IntVect{i,j,k},n);
				if(p->number_of_particles < 10) AllPrint() << " ";
				Print() << p->number_of_particles << " ";
				count += p->number_of_particles;
			});
			AllPrint() << "\n" << mfi.validbox() << "\n";
		}
	}

	Print() << count << " particles\n";
	AllPrint() << "rank " << ParallelDescriptor::MyProc() << " has " << total << " particles in total\n";
	AllPrint() << " + " << nghost << " in ghost cells\n";
}

void count(Amr& amr) 
{
	sleep(ParallelDescriptor::MyProc());
	int total = 0;
	int nghost = 0;
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		AllPrint() << "RANK " << ParallelDescriptor::MyProc() << ", "
				   << "LEVEL " << lev;

		int count = 0;
		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) 
		{
			if (!mfi.isPropagated()) continue;
			const auto& fab = amr[lev].getCells()[mfi];
			int line = 0;
			For(mfi.validbox().grow(IntVect{0,0,0}), amr.nComp(), [&] (int i, int j, int k, int n) 
			{
				auto& p = fab(IntVect{i,j,k},n);
				count += p->number_of_particles;
			});
		}
		AllPrint() << ": " << count << " particles\n";
	}

}

void propagate(const Box& box, CellFab& fab, Geometry& geom) 
{
	const double vx = 0.1;
	//Print() << box << "\n";
	ParallelFor(box, [&] (int i, int j, int k) 
	{	
		// loop over particles
		IntVect iv(i,j,k);
		//Print() << iv << "\n";
		for (uint n = 0; n < fab(iv)->particles.size(); n++) 
		{
			auto &particle = fab(iv)->particles[n];
			particle[0] += vx;

			// hande grid wrap around
			if(particle[0] >= geom.period(0)) {
				particle[0] -= geom.period(0);
			}
		}
	});
}

/*!
Saves the local grid and its particles to separate vtk files.

Each process saves its own files.
*/
void save(const int rank, Amr &amr, unsigned int step)
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
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		auto& grid = amr[lev].getCells();
		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) {
			if (!mfi.isPropagated()) continue;

			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k)
			{
				total_particles += a(i,j,k)->particles.size();
			});
		}
	}
	
	outfile << "POINTS " << total_particles << " float\n";
	// write out coordinates
	uint64_t written_particles = 0;
	
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		auto& grid = amr[lev].getCells();
		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) {
			if (!mfi.isPropagated()) continue;

			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k)
			{
				auto const particles = a(i,j,k)->particles;
				for (uint n = 0; n < a(i,j,k)->particles.size(); n++)
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

	outfile << "\nSCALARS lev int 1\nLOOKUP_TABLE default\n";
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		auto& grid = amr[lev].getCells();
		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) {
			if (!mfi.isPropagated()) continue;

			auto const& a = grid[mfi].array();
			For(mfi.validbox(), [&] (int i, int j, int k) {
				for (int n = 0; n < a(i,j,k)->number_of_particles; n++)
				{
					outfile << lev << " ";
				}
			});
		}
	}

	// cell numbers of particles
	outfile << "\nSCALARS cell int 1\nLOOKUP_TABLE default\n";
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		auto& grid = amr[lev].getCells();
		for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) {
			if (!mfi.isPropagated()) continue;
			
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
	amrex::Print() << "Running " << comm_size << " processes\n";
	
	Amr amr;
	amrex::Print() << "   Domain: " << amr.Geom(0).Domain().size() << "\n";
	amr.defBaseLevel();

	Print() << "Level " << 0 << " has " << amr[0].countCells() << " cells\n"; 
	// Initialize particles
	{
		CellFabArray& grid = amr[0].getCells();
		int y = amr.Geom(0).Domain().length(1);

		for (AmrIter mfi(amr[0]); mfi.isValid(); ++mfi) {
			const Box& box = mfi.validbox();
			auto a = grid[mfi].array();
			int n = 0;
			For(box, [&] (int i, int j, int k) 
			{
				Cell* cell_data = a(i,j,k).get();
				if(!cell_data) Print() << "something went wrong\n";

				const unsigned int number_of_particles = i+y*j;
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
	}

    // Build fine level grids.
    if (amr.maxLevel() > 0) amr.bldFineLevels();

	Print() << "Finest level = " << amr.finestLevel() << "\n";
	
	amr.FillAllBoundaries();
	count(amr);
	//printCounts(amr);

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
		
		//if (step % 5 == 0) amr.LoadBalance();
		//Print() << "step " << step << "\n";
		amr.FillAllBoundaries();

		//count(amr);
	
		save(rank, amr, step);

		int Ncomp = amr.nComp();
		for (int lev = 0; lev <= amr.finestLevel(); lev++)
		{
			//Print() << "	lev " << lev << "\n";
			Geometry &geom = amr.Geom(lev);
			CellFabArray& grid = amr[lev].getCells();

			// move particles in x-direction
			for (AmrIter mfi(amr[lev]); mfi.isValid(); ++mfi) { // loop over local boxes
				if (!mfi.isPropagated()) continue;
				Box box = mfi.validbox();
				propagate(box, grid[mfi], geom);
			}

			auto& fine_bndry = amr[lev].getFineBoudary();
			if (fine_bndry.coarse) {
				BoxList done_boxes;
				
				for (size_t fab = 0; fab < fine_bndry.coarse->size(); fab++)
				{
					CellFab& cfab = fine_bndry.coarse->get(fab);
					BoxList do_boxes(complementIn(cfab.box(), done_boxes));
					for (int i = 0; i < do_boxes.size(); i++)
					{
						propagate(do_boxes.data()[i], cfab, geom);
					}

					done_boxes.join(do_boxes);
				}
			}
			auto& crse_bndry = amr[lev].getCoarseBoudary();
			if (crse_bndry.fine) {
				BoxList done_boxes;
				for (size_t fab = 0; fab < crse_bndry.fine->size(); fab++)
				{
					CellFab& cfab = crse_bndry.fine->get(fab);
					BoxList do_boxes(complementIn(cfab.box(), done_boxes));
					for (int i = 0; i < do_boxes.size(); i++)
					{
						propagate(do_boxes.data()[i], cfab, geom);
					}
					done_boxes.join(do_boxes);
				}
			}
			
			update_cell_lists(amr[lev], amr.Geom(lev));
		}
		
		//count(amr);
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

	save(rank, amr, max_steps);

	clock_t after = clock();
	Print() << "Finalizing"
		<< ": simulation took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds "
		<< endl;
	printCounts(amr);

	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Couldn't initialize MPI." << endl;
		abort();
	}

	// arguments must include inputs file
    amrex::Initialize(argc, argv, true,  MPI_COMM_WORLD);
	
	main_main();

	amrex::Finalize();
	MPI_Finalize();
}