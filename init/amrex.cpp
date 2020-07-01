/*
Tests the scalability of initializing the grid
*/

#include "cstdlib"
#include "ctime"

#include "boost/program_options.hpp"
#include "zoltan.h"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

struct CellData {

	double variables[3];

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->variables), 3, MPI_DOUBLE);
	}
};

using namespace std;
using namespace amrex;

int main(int argc, char* argv[])
{   
	clock_t before = clock();

    amrex::Initialize(argc,argv);
    
	amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    
	int rank = -1;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	/*
	Options
	*/
	uint64_t x_length, y_length, z_length;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("x_length",
			boost::program_options::value<uint64_t>(&x_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the x direction")
		("y_length",
			boost::program_options::value<uint64_t>(&y_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the y direction")
		("z_length",
			boost::program_options::value<uint64_t>(&z_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the z direction");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
    	amrex::Finalize();
		return EXIT_SUCCESS;
	}

    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(x_length-1, y_length-1, z_length-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        //ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0_rt,-1.0_rt,-1.0_rt)},
                         {AMREX_D_DECL( 1.0_rt, 1.0_rt, 1.0_rt)});

        // periodic in all direction
        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

        // This defines a Geometry object
        geom.define(domain,real_box,CoordSys::cartesian,is_periodic);
    }
	BaseFab<CellData> data;
	
	
	int myproc = ParallelDescriptor::MyProc();  // Return the rank

	int nprocs = ParallelDescriptor::NProcs();  // Return the number of processes

	if (ParallelDescriptor::IOProcessor()) {
		// Only the I/O process executes this
	}

	int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

	if (nprocs > 1)
	{
	
	// Broadcast 100 ints from the I/O Processor
	Vector<int> a(100), b(100);
	a.assign(a.size(), myproc+10);
	
    const int SeqNum = ParallelDescriptor::SeqNum();

	cout << "Process " << myproc << " to " << (myproc+1) % nprocs << "\n";
	ParallelDescriptor::Send(a, (myproc+1) % nprocs, SeqNum);

	MPI_Request req = ParallelDescriptor::Recv(b, (myproc+1) % nprocs, SeqNum).req();
	MPI_Status stat;
	cout << "Process " << myproc << " waits" << "\n";
	ParallelDescriptor::Wait(req, stat);
	
	cout << "Received b[5]=" << b[5] << "\n";

	// See AMReX_ParallelDescriptor.H for many other Reduce functions
	//ParallelDescriptor::ReduceRealSum(x);

	}

	clock_t after = clock();

	cout << "Process " << myproc << " " << nprocs << " " << ioproc
		<< ": grid initialization took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds (total grid size " << x_length * y_length * z_length << ")"
		<< endl;

	//MPI_Finalize();

    amrex::Finalize();

	return EXIT_SUCCESS;
}

