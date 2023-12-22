#ifndef SUBDOMAIN_HPP
#define SUBDOMAIN_HPP

#include "defines.hpp"
#include "world.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

class SubDomain : public World_LC
{
public:
    SubDomain();
    ~SubDomain(){};

    int local_particle_num = 0;
    int local_cell_num = 0;
    // number of cells per direction (without the ghost) 
    int local_cell_N[DIM];
    // number of cells in the ghost
    int ghost_cell_num;
    // rank of process
    int myrank;
    // number of processors
    int myprocessors;
    // total number of processes
    int numprocs;
    // multiindex of proc
    int ip[DIM];
    // Multiindex of process hyperrectangle
    int N_p[DIM];
    // Neighboring processes
    int ip_lower[DIM];
    int ip_upper[DIM];

    // vector of breadth of border_cells - at same time multiindex of first
    // interior cell
    int ic_start[DIM];
    // multiindex of bordercell just outside interior - diagonally opposed to
    // ic_start
    int ic_stop[DIM];
    // multiindex of cell hyperrectangle on proc
    int ic_number[DIM];
    // global index of first interior cell.
    int ic_lower_global[DIM];

    MPI_Datatype MPI_Particle;
    std::vector<Particle> received_particles;

    // local variables for reduction
    real local_e_kin;
    real local_e_pot;
	
    void construct_MPI_Particle();

    void read_Parameter(const std::string& filename);
	// compute single indices if neighboring processes
    void compute_neighboring_procs();
    /// Init only cells appropriate for subdomain.
    void init_local_cells_in_subdomain();

    std::vector<unsigned> compute_cell_neighbors(const int linear_cell_index);
    /// Overwrites World, because want to decide immediately whether to keep
    /// particle or not.
    void read_Particles(const std::string& filename);
    /// Sorting particles at init stage
    void init_particle_in_cells();

    /// Sorting particles during simulation
    void sort_particles_in_cells();


	// index computation
    int compute_local_cell_single_index(const int (&cell_index)[DIM]);
    void compute_local_cell_multi_index(const int local_cell_single_index,
                                        int (&cell_index)[DIM]);
	int local_linear_cell_index(Particle *p); 

    int compute_processor_single_index(const int (&processor_index)[DIM]);
    void compute_processor_multi_index(const int processor_single_index,
                                       int (&processor_index)[DIM]);
	
	
    void collect_global_data();


protected:
    void print(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, SubDomain& S);
};

#endif  // SUBDOMAIN_HPP
