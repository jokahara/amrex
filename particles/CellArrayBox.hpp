
#include <AMReX_Box.H>
#include <AMReX_BaseFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

using namespace amrex;

class Cell {
public:
	unsigned int number_of_particles;

	// coordinates of particles in this cell
	std::vector<std::array<double, 3> > particles;

	/*
	The number of particles is transferred over MPI if
	this is false, otherwise the particle coordinates
	are transferred.
	*/
	static bool transfer_particles;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		void* address = NULL;
		int count = -1;
		MPI_Datatype datatype = MPI_DATATYPE_NULL;
	
		if (Cell::transfer_particles) {
			
			if (number_of_particles > 0) {
				if(particles.size() != number_of_particles) {
					std::cerr << "particles not resized correctly " << particles.size() 
							<< " != " << number_of_particles << "\n";
					abort();
				}
				address = &(particles[0]);
			} else {
				// return a sane address just in case
				address = &(number_of_particles);
			}

			count = number_of_particles * 3;
			datatype = MPI_DOUBLE;

		} else {

			address = &(number_of_particles);
			count = 1;
			datatype = MPI_UNSIGNED;

		}

		return std::make_tuple(address, count, datatype);
	}

	Cell() {
		number_of_particles = 0;
	};

	Cell(const Cell &c) {
		*this = c;
	}

	Cell& operator=(const Cell &c) {
		number_of_particles = c.number_of_particles;
		resize();

		std::memcpy(&particles[0], &c.particles[0], 3 * number_of_particles * sizeof(double) );
		
		return *this;
	}

	// reserves space for particle data coming over MPI.
	void resize()
	{
		this->particles.resize(this->number_of_particles);
	}
};


class cArrayBox : public BaseFab<Cell>
{
public:
    //! Construct an invalid FAB with no memory.
    cArrayBox () noexcept;

    //explicit cArrayBox (Arena* ar) noexcept;
    //cArrayBox (const Box& b, int ncomp, Arena* ar);

    /**
    * \brief Construct an initial FAB with the data space allocated but
    * not inititialized. ncomp is the number of components
    * (variables) at each data point in the Box.
    */
    explicit cArrayBox (const Box& b,
                        int ncomp=1,
						bool alloc=true,
						bool shared=false,
            			Arena* ar = nullptr)
    : BaseFab<Cell>(b,ncomp,alloc,shared,ar) {}

    cArrayBox (const cArrayBox& rhs, MakeType make_type, int scomp, int ncomp)
	: BaseFab<Cell>(rhs, make_type, scomp, ncomp) {}

    // explicit = cannot convert or copy
    explicit cArrayBox (Array4<Cell> const& a) noexcept : BaseFab<Cell>(a) {}

    explicit cArrayBox (Array4<Cell> const& a, IndexType t) noexcept : BaseFab<Cell>(a,t) {}

    explicit cArrayBox (Array4<Cell const> const& a) noexcept : BaseFab<Cell>(a) {}

    explicit cArrayBox (Array4<Cell const> const& a, IndexType t) noexcept : BaseFab<Cell>(a,t) {}

    //!  The destructor.
    ~cArrayBox () noexcept {}

    cArrayBox (cArrayBox&& rhs) noexcept = default;

    cArrayBox (const cArrayBox&) = delete;  	// copying not allowed here
    cArrayBox& operator= (const cArrayBox&) = delete;
    cArrayBox& operator= (cArrayBox&&) = delete;

private:

    static bool do_initval;

};

/*!
\brief Returns the error string of given MPI error.
*/
class Error_String
{
public:

	/*!
	mpi_return_value is the value returned by a failed MPI function.
	*/
	std::string operator()(int mpi_return_value)
	{
		char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
		int string_length;
		MPI_Error_string(mpi_return_value, mpi_error_string, &string_length);
		mpi_error_string[string_length + 1] = '\0';

		const std::string result(mpi_error_string);
		return result;
	}
};