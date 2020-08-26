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
#define NO_CHILD_LID -1

bool Cell::transfer_particles = false;

using namespace std;
using namespace boost;
using namespace amrex;


void propagate(const Box& box, CellFab& fab, const Geometry& geom) 
{
	auto const& a = fab.array();
	const double vx = 0.1;
	ParallelFor(box, [&] AMREX_GPU_DEVICE (int i, int j, int k) 
	{	
		// loop over particles
		IntVect iv(i,j,k);
		for (uint n = 0; n < a(iv)->particles.size(); n++) 
		{
			auto &particle = a(iv)->particles[n];
			particle[0] += vx;

			// hande grid wrap around
			if(particle[0] >= geom.ProbLength(0)) {
				particle[0] -= geom.ProbLength(0);
			}
		}
	});
}

void update_cell_lists(AmrLevel& level, Geometry &geom) 
{
	CellFabArray &grid = level.getCells();
//sleep(ParallelDescriptor::MyProc());
	int count = 0;
	// move particles to the particle list of the cell the particles are currently inside of
	for (AmrIter it(level); it.isValid(); ++it) 
	{	
		bool has_children = (it.childLID() != NO_CHILD_LID);

		auto const& a = grid[it].array();
		auto box = grid[it].box();

		amrex::For(it.validbox(), [&] (int i, int j, int k)
		{
			const IntVect cell(i,j,k);
			auto& cell_data = a(cell);

			int n = 0;
			while (n < cell_data->particles.size()) {				
				IntVect current_cell = 
					geom.CellIndex( cell_data->particles[n].data() );

				if (current_cell != cell) {
					// check for periodicity
					if(!box.contains(current_cell)) {
						current_cell[0] += geom.ProbDomain().length(0);
					}
					if(current_cell[0] == 0) {
						count++;
						//Print() << cell << " adds " << cell_data->particles[n][0] << " to " << current_cell << "\n";
					}
					auto& current_data = a(current_cell);
					current_data->particles.push_back(cell_data->particles[n]);
					current_data->number_of_particles = current_data->particles.size();

					cell_data->particles.erase(cell_data->particles.begin() + n);
					cell_data->number_of_particles = cell_data->particles.size();
				}
				else n++;
			}
			// clear data received from fine level
			if (has_children) cell_data->clear();
		});
	}
	//AllPrint() << ParallelDescriptor::MyProc() << " added first " << count << "\n";
	
	if (ParallelDescriptor::NProcs() > 1)
	for (auto& cont : grid.get_receive_tags())
	{
		for (auto& tag : cont.second)
		{
			//if (level.childLID(tag.dstIndex) != NO_CHILD_LID) continue;
			propagate(tag.dbox, grid[tag.dstIndex], geom);

			auto const& a = grid[tag.dstIndex].array(); 
			const Box& box = grid[tag.dstIndex].box();
			amrex::For(tag.dbox, [&] (int i, int j, int k) 
			{
				const IntVect cell(i,j,k);
				auto& cell_data = a(cell);

				for (int n = 0; n < cell_data->particles.size(); n++) 
				{
					const IntVect current_cell = 
						geom.CellIndex( cell_data->particles[n].data() );
					
					if (current_cell != cell && box.contains(current_cell)) {
						
					if(current_cell[0] == 0) {
						count++;
						//Print() << cell << " remote add " << cell_data->particles[n][0] << " to " << current_cell << "\n";
					}
						auto& current_data = a(current_cell);
						current_data->particles.push_back(cell_data->particles[n]);
						current_data->number_of_particles = current_data->particles.size();

						cell_data->particles.erase(cell_data->particles.begin() + n);
						cell_data->number_of_particles = cell_data->particles.size();
						n--;
					}
				}
				// data not needed anymore
				cell_data->clear();
			});
		}
	}
	//Print() << " added " << count << "\n";
	
	// updating from coarse to fine boundary
	count = 0;
	auto& crse_bndry = level.getCoarseBoudary();
	if (crse_bndry.fine) {
	for (int fab = 0; fab < crse_bndry.fine->size(); fab++)
	{
		Box box = crse_bndry.fine->box(fab);
		auto const& a = crse_bndry.fine->get(fab);
		amrex::For(box, [&] (int i, int j, int k)
		{
			const IntVect cell(i,j,k);
			auto& cell_data = a(cell);
			
			for (int n = 0; n < cell_data->particles.size(); n++) {	
				cell_data->particles[n][0] += 0.1;

				const IntVect current_cell = 
					geom.CellIndex( cell_data->particles[n].data() );

				if (current_cell == cell) continue;

				// find new location
				for (AmrIter it(level); it.isValid(); ++it) {
					if (it.validbox().contains(current_cell)) {
						count++;
						auto& current_data = level[it](current_cell);
						current_data->particles.push_back(cell_data->particles[n]);
						current_data->number_of_particles = current_data->particles.size();
						break;
					}
				}
			}
			cell_data->clear();
		});
	}}
	//ParallelContext::BarrierAll();
}

void printCounts(Amr& amr) 
{
	sleep(ParallelDescriptor::MyProc());
	int total = 0;
	
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		AllPrint() << "\n\tRANK " << ParallelDescriptor::MyProc() << ", "
				   << "LEVEL " << lev << "\n";

		for (AmrIter it(amr[lev]); it.isValid(); ++it) 
		{
			const auto& fab = amr[lev].getCells()[it];
			int line = 0, sum = 0;
			amrex::For(it.validbox().grow(IntVect{1,1,0}), amr.nComp(), [&] (int i, int j, int k, int n) 
			{
				if(j != line) {
					line = j;
					AllPrint() << "\n";
				}
				IntVect iv{i,j,k};
				auto& p = fab(iv,n);
				if(p->number_of_particles < 10) AllPrint() << " ";
				AllPrint() << p->number_of_particles << " ";
				if (it.validbox().contains(iv) ) sum += p->number_of_particles;
				//if (it.validbox().contains(iv)) AllPrint() << "  ";
				//else 
				//AllPrint() << p.use_count() << " ";
			});
			AllPrint() << "\n" << it.validbox() << ": " << sum << "\n";
			
			if (it.childLID() == NO_CHILD_LID) total += sum;
		}
	}

	AllPrint() << total << " particles\n";
}

int count(Amr& amr) 
{
	//sleep(ParallelDescriptor::MyProc());
	int total = 0;
	int nghost = 0;
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		AllPrint() << "RANK " << ParallelDescriptor::MyProc() << ", "
				   << "LEVEL " << lev;

		int count = 0;
		for (AmrIter it(amr[lev]); it.isValid(); ++it) 
		{
			if (it.childLID() != NO_CHILD_LID) continue;

			const auto& fab = amr[lev].getCells()[it];
			int line = 0;
			amrex::For(it.validbox(), amr.nComp(), [&] (int i, int j, int k, int n) 
			{
				auto& p = fab(IntVect{i,j,k},n);
				count += p->number_of_particles;
			});
		}
		AllPrint() << ": " << count << " particles\n";
		total += count;
	}
	return total;
}

/*!
Writes the cells on this process into a vtk file with given name in ASCII format.

The cells are written in ascending order.
Must be called simultaneously on all processes.

Returns true on success, false otherwise.
*/
bool write_grid(Amr& amr, const std::string& file_name)
{
	std::ofstream outfile(file_name);
	if (!outfile.is_open()) {
		std::cerr << "Couldn't open file " << file_name << std::endl;
		return false;
	}

	outfile <<
		"# vtk DataFile Version 2.0\n"
		"Cartesian cell refinable grid\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n";

	int n_points = 0;
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		for (AmrIter it(amr[lev]); it.isValid(); ++it)
		{
			if (it.childLID() == NO_CHILD_LID)
				n_points += it.validbox().numPts();
		}
	}

	// write separate points for every cells corners
	outfile << "POINTS " << n_points * 8 << " float" << std::endl;
	for (int lev = 0; lev <= amr.finestLevel(); lev++)
	{
		for (AmrIter it(amr[lev]); it.isValid(); ++it)
		{
			if (it.childLID() != NO_CHILD_LID) continue;
		
			const Geometry& geom = amr.Geom(lev);
			Box region = it.validbox();
			Vector<Real> x_edge, y_edge, z_edge;
			geom.GetEdgeLoc(x_edge, region, 0);
			geom.GetEdgeLoc(y_edge, region, 1);
			geom.GetEdgeLoc(z_edge, region, 2);
			for (int i = 0; i < x_edge.size() - 1; i++)
			for (int j = 0; j < y_edge.size() - 1; j++)
			for (int k = 0; k < z_edge.size() - 1; k++)
			{
				const std::array<double, 3>
					cell_min = {x_edge[i], y_edge[j], z_edge[k]},
					cell_max = {x_edge[i+1], y_edge[j+1], z_edge[k+1]};

				outfile
					<< cell_min[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
					<< cell_max[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
					<< cell_min[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
					<< cell_max[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
					<< cell_min[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
					<< cell_max[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
					<< cell_min[0] << " " << cell_max[1] << " " << cell_max[2] << "\n"
					<< cell_max[0] << " " << cell_max[1] << " " << cell_max[2] << "\n";
			}
		}
	}
	

	// map cells to written points
	outfile << "CELLS " << n_points << " " << n_points * 9 << std::endl;
	for (unsigned int j = 0; j < n_points; j++) {
		outfile << "8 ";
		for (int i = 0; i < 8; i++) {
				outfile << j * 8 + i << " ";
		}
		outfile << std::endl;
	}

	// cell types
	outfile << "CELL_TYPES " << n_points << std::endl;
	for (unsigned int i = 0; i < n_points; i++) {
		outfile << 11 << std::endl;
	}

	if (!outfile.good()) {
		std::cerr << "Writing of vtk file probably failed" << std::endl;
		// TODO: throw an exception instead
		exit(EXIT_FAILURE);
	}

	outfile.close();

	return true;
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
	write_grid(amr, grid_file_name.c_str());

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
		for (AmrIter it(amr[lev]); it.isValid(); ++it) {
			if (it.childLID() != NO_CHILD_LID) continue;

			auto const& a = grid[it].array();
			amrex::For(it.validbox(), [&] (int i, int j, int k)
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
		for (AmrIter it(amr[lev]); it.isValid(); ++it) {
			if (it.childLID() != NO_CHILD_LID) continue;

			auto const& a = grid[it].array();
			amrex::For(it.validbox(), [&] (int i, int j, int k)
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
		for (AmrIter it(amr[lev]); it.isValid(); ++it) {
			if (it.childLID() != NO_CHILD_LID) continue;

			auto const& a = grid[it].array();
			amrex::For(it.validbox(), [&] (int i, int j, int k) {
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
		for (AmrIter it(amr[lev]); it.isValid(); ++it) {
			if (it.childLID() != NO_CHILD_LID) continue;
			
			auto const& a = grid[it].array();
			amrex::For(it.validbox(), [&] (int i, int j, int k) {
				for (int n = 0; n < 3; n++)
				{
					outfile << it.index() << " ";
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
	Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

	int rank = ParallelDescriptor::MyProc();
	int comm_size = ParallelDescriptor::NProcs();
	
	Amr amr;
	AllPrint() << "Domain: " << amr.Geom(0).Domain().size() << "\n";
	amr.defBaseLevel();
	
	AllPrint() << "Level " << 0 << " has " << amr[0].countCells() << " cells\n"; 
	// Initialize particles
	{
		CellFabArray& grid = amr[0].getCells();
		IntVect length = amr.Geom(0).Domain().size();
		int x = std::min(length[0], 12) / 2;
		int y = length[1] / 2;

		for (AmrIter it(amr[0]); it.isValid(); ++it) 
		{
			const Box& box = it.validbox();
			auto a = grid[it].array();
			int n = 0;
			amrex::For(box, [&] (int i, int j, int k) 
			{
				Cell* cell_data = a(i,j,k).get();
				if(!cell_data) Abort("Cell was not allocated!");

				int nx = (x > i) ? i : abs(i - 2*x);
				int ny = (y > j) ? j : abs(j - 2*y);
				unsigned int number_of_particles = nx*ny;

				if (i>2*x) number_of_particles = 2;
				for (unsigned int n = 0; n < number_of_particles; n++) 
				{
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
    if (amr.maxLevel() > 0) {
    	Print() << "Building fine levels\n";
		amr.bldFineLevels();
	}

	Print() << "Finest level = " << amr.finestLevel() << "\n";
	
	amr.FillAllBoundaries();
	
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
	for (unsigned int step = 0; step < max_steps; step++) 
	{
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

		if (amr.Verbose() > 0) 
			Print() << "step " << step << "\n";
		
		// regrid all levels every 5 steps
		if (step % 5 == 0 && step > 0) {
			amr.regrid(0);
			if(step == 40) printCounts(amr);
		} 
		
		// Fills ghost cells with data from other levels and remote neighbours.
		amr.FillAllBoundaries();
		
		save(rank, amr, step);
		//printCounts(amr);

		int Ncomp = amr.nComp();
		for (int lev = 0; lev <= amr.finestLevel(); lev++)
		{
			Geometry &geom = amr.Geom(lev);
			CellFabArray& grid = amr[lev].getCells();

			// loop over local boxes
			//#pragma omp parallel
			for (AmrIter it(amr[lev], false); it.isValid(); ++it) 
			{ 
				//if (it.childLID() != NO_CHILD_LID) continue;
				Box box = it.validbox();
				// move particles in x-direction
				propagate(box, grid[it], geom);
			}
			
			update_cell_lists(amr[lev], geom);
		}
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

	//printCounts(amr);
	count(amr);
	
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