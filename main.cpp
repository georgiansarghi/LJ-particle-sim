#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <mpi.h>

#include "ljpotential.hpp"
#include "observer.hpp"
#include "subdomain.hpp"
#include "thermostat.hpp"
#include "velocityverlet.hpp"
#include "world.hpp"

int main(int argc, char* argv[]) {
    double wall_time;
    int my_rank;

    // check arguments
    if (argc < 2) {
        std::cerr << "Error: missing arguments." << std::endl;
        std::cerr << "Usage: " << std::endl
                  << "\t" << argv[0] << " <parameter_file> <particle_data_file>"
                  << std::endl;
        return EXIT_FAILURE;
    }

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
        wall_time = MPI_Wtime();
    
    // instanciate Potential
    LJPotential Pot;

    // create World
    SubDomain W;

    // read Parameters
    W.read_Parameter(argv[1]);

    // read Particles
    W.read_Particles(argv[2]);

    std::cout << W << std::endl;

    W.T.set_start_velocities(W.particles);

    // create the Observer
    Observer_Sub O(W);
    
    O.open_files();


    // instanciate timediscretization
    // remark: & is used to get the address of Pot
    VelocityVerlet_Sub Verlet(W, &Pot, O);

    // run the simulation
    Verlet.simulate();

    MPI_Finalize();

    if (my_rank == 0) {
        wall_time = MPI_Wtime() - wall_time;
        std::cout << std::endl
                  << "Elapsed wall clock time = " << wall_time << " seconds."
                  << std::endl
                  << std::endl;
    }
    return EXIT_SUCCESS;
}
