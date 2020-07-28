
#ifndef Cell_H
#define Cell_H

#include <vector>
#include <tuple>

#include "mpi.h"

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
		particles = c.particles;
		
		return *this;
	}

	Cell& operator+=(const Cell &c) {
		for (size_t i = 0; i < c.number_of_particles; i++)
		{
			particles.push_back(c.particles[i]);
		}
		number_of_particles = particles.size();
		
		return *this;
	}

	// efficiently move src data to this cell 
	void swap(Cell& src) {
		particles.swap(src.particles);
		number_of_particles = particles.size();
		src.number_of_particles = src.particles.size();
	}

	// reserves space for particle data coming over MPI.
	void resize()
	{
		this->particles.resize(this->number_of_particles);
	}
};

#endif