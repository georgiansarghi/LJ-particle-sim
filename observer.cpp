#include "observer.hpp"

Observer::Observer(World& W) : W(W) {}

Observer::~Observer()
{
    // close the statistics file
    if (statistics.is_open())
        statistics.close();
    // and the coordinates file
    if (coordinates.is_open())
        coordinates.close();
    // and the xyz file
    if (xyz.is_open())
        xyz.close();
}

void Observer::open_files() {

    // open statistics file
    std::string statistics_filename = W.name + ".statistics";
    // open file, overwrite existing files, take no prisioners
    statistics.open(statistics_filename.c_str());
    // and tell the world
    std::cout << "Opened " << statistics_filename << " for writing."
              << std::endl;

    std::string coordinates_filename = W.name + ".coordinates";
    coordinates.open(coordinates_filename.c_str());
    std::cout << "Opened " << coordinates_filename << " for writing."
              << std::endl;

    // open xyz file
    std::string xyz_filename = W.name + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyz_filename.c_str());
    // and tell the world
    std::cout << "Opened " << xyz_filename << " for writing."
              << std::endl;


}

void Observer::output_statistics() {
    // write statistics into the filestream, seperated with tabulars
    statistics << W.t << "\t" << W.e_pot << "\t" << W.e_kin << "\t" << W.e_tot << "\t"
               << std::endl;
}

void Observer::output_coordinates() {
    // write position of particles into the filestream, seperated with tabulars
    coordinates << W.t << "\t";
    
    for(unsigned i = 0; i < W.particles.size(); i++)
        for(unsigned j = 0; j < DIM; j++)
            coordinates << W.particles[i].x[j] << "\t";

    coordinates << std::endl;
}

void Observer::output_xyz() {
    // number of particles
    xyz << W.particles.size() << std::endl;
    // time stamp
    xyz << "Zeit: " << W.t << std::endl;
    
    // write position of particles into the filestream, seperated with tabulars
    for(unsigned i = 0; i < W.particles.size(); i++) {
        // particle type
        xyz << "H\t";
        // particle position
        for(unsigned j = 0; j < DIM; j++)
            xyz << W.particles[i].x[j] << "\t";
        xyz << std::endl;
    }
}

void Observer::notify() {
    // call output functions
    output_statistics();
    output_coordinates();
    output_xyz();
}

///////////////////////// Observer_LC starting here.


void Observer_LC::output_statistics() {
    // write statistics into the filestream, seperated with tabulars
	// For the first 99 timesteps, the average energies are outputted as the actual energies
    statistics << W.t << "\t" 
                << W.e_pot << "\t" 
                << W.e_kin << "\t" 
                << W.e_tot << "\t" 
                << W.T.current_temperature;
	if (W.timestep < W.out_interval)
		statistics << "\t\t" << W.e_pot << "\t" << W.e_kin << "\t" << W.e_tot;
	else
		statistics << "\t\t" << W.e_pot_avg << "\t" << W.e_kin_avg << "\t" << W.e_tot_avg;
    statistics << std::endl;
}



void Observer_LC::output_coordinates() {
    // write position of particles into the filestream, seperated with tabulars
    coordinates << W.t << "\t";

    for(unsigned i = 0; i < W.particles.size(); i++)
    {
        for(unsigned d = 0; d < DIM; d++)
            coordinates << W.particles[i].x[d] << "\t";
        coordinates << std::endl;
    }
}



void Observer_LC::output_xyz() {
    // number of particles

    xyz << W.number_of_particles << std::endl;
    
    // time stamp
    xyz << "Zeit: " << W.t << std::endl;
    
    // write position of particles into the filestream, seperated with tabulars
    for(unsigned i = 0; i < W.particles.size(); i++) {
        if(W.particles[i].active) {
            // particle type
            xyz << "H\t";
            // particle position
            for(unsigned j = 0; j < DIM; j++)
                xyz << W.particles[i].x[j] << "\t";
            xyz << std::endl;
        }
    }
}

//////////////// Observer for subdomain
Observer_Sub::Observer_Sub(SubDomain& W) : Observer_LC(W), W(W) {}


void Observer_Sub::open_files() {


    // For multiprocs add proc number to output files.
    char rankstring[6];
    sprintf(rankstring, ".%04d", W.myrank);

    //  Close Observer stuff
    if (statistics.is_open())
        statistics.close();
    if (W.myrank == 0) {
        // open statistics file
        std::string statistics_filename = "./results/" + W.name + ".statistics";
        // open file, overwrite existing files, take no prisioners
        statistics.open(statistics_filename.c_str());
        // and tell the world
        std::cout << "Opened " << statistics_filename << " for writing."
                  << std::endl;
    }

    if (coordinates.is_open())
        coordinates.close();
    // open coordinates file
    std::string coordinates_filename = "./results/" + W.name + rankstring + ".coordinates";
    coordinates.open(coordinates_filename.c_str());
    std::cout << "Opened " << coordinates_filename << " for writing."
              << std::endl;

    if(xyz.is_open())
        xyz.close();
    // open xyz file
    std::string xyz_filename = "./results/" + W.name + rankstring + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyz_filename.c_str());
    // and tell the world
    std::cout << "Opened " << xyz_filename << " for writing." << std::endl;

}



void Observer_Sub::notify() {
    // only one statistic output.
    if (W.myrank == 0)
        output_statistics();
    //  But all write output of coordinates.
    Observer_Sub::output_xyz();
}



void Observer_Sub::output_xyz() {

    if ((out % W.out_interval) == 0) {
        int count = 0;
        for(unsigned i = 0; i < W.cells.size(); i++) {
            if(!W.cells[i].ghost_flag)
                for(unsigned j = 0; j < W.cells[i].particles.size(); j++) {
                    count++;
                }
        }
        xyz << count << std::endl;
        // time stamp
        xyz << "Zeit: " << W.t << std::endl;
        
        // write position of particles into the filestream, seperated with tabulars
        for(unsigned i = 0; i < W.cells.size(); i++) {
            if(!W.cells[i].ghost_flag)
                for(unsigned j = 0; j < W.cells[i].particles.size(); j++) {
                    // particle type
                    xyz << "H\t";
                    // particle position
                    for(unsigned n = 0; n < DIM; n++)
                        xyz << W.cells[i].particles[j].x[n] << "\t";
                    xyz << std::endl;
                }
        }
    } 
    out++;
}



void Observer_Sub::output_statistics() {
    // write statistics into the filestream, seperated with tabulars
    statistics << W.t << "\t" 
                << W.e_pot << "\t" 
                << W.e_kin << "\t" 
                << W.e_tot << "\t" 
                << std::endl;
}
