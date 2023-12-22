
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <unordered_set>

#include "world.hpp"



World::World()
    : name("unknown"), t(0), delta_t(0), t_end(0), timestep(0), e_kin(0), e_pot(0), e_tot(0), e_kin_avg(0), e_pot_avg(0), e_tot_avg(0) {
    // empty constructor
}



BorderType World::string_to_BorderType(std::string s) {
    if (s == "unknown" ) return unknown;
    else if (s == "periodic" ) return periodic;
    else return leaving;
}



void World::read_Parameter(const std::string& filename) {
    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("read_Parameter(): Can't open file '" +
                                 filename + "' for reading.");

    // helper strings
    std::string line, option, s;

    // read file till eof
    while (parfile.good()) {
        // init empty option
        option = "";
        // read line from file
        getline(parfile, line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option and value from stringstream
        strstr >> option;
        // check options and read values
        if (option == "delta_t") {
            strstr >> delta_t;
        } else if (option == "t_end") {
            strstr >> t_end;
        } else if (option == "name") {
            strstr >> name;
        } else if (option == "length") {
            for(int i = 0 ; i < DIM ; i++) {
                strstr >> length_x[i];
            }
        } else if (option == "upper_border") {
            for(int i = 0 ; i < DIM ; i++) {
                    strstr >> s;
                    upper[i] = string_to_BorderType(s);
                }
        } else if (option == "lower_border") {
            for(int i = 0 ; i < DIM ; i++) {
                    strstr >> s;
                    lower[i] = string_to_BorderType(s);
                }
        } else if (option == "epsilon") {
            strstr >> eps;
        } else if (option == "sigma") {
            strstr >> sigma;
        } else if (option == "output_interval") {
            strstr >> out_interval;
        }
    }

    parfile.close();
}



void World::read_Particles(const std::string& filename) {
    std::ifstream parfile(filename.c_str());
    if(!parfile.is_open())
        throw std::runtime_error("read_Particles(): Can't open file '" +
                                 filename + "' for reading.");

	// helper string
    std::string line, option;
	
	// helper particle
	Particle p;

     // read file till eof
    while(parfile.good()) {
		// init empty option
        option = "";
        // read line from file
        getline(parfile, line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
		// test for empty line
		strstr >> option;
        
		if( option == "" ) 
			continue;
        // read values from stringstream and write them in particles

		// mass
        strstr >> p.m;
		// position
		for(int k = 0; k < DIM; k++)
			strstr >> (p.x)[k];
		// velocity
		for(int k = 0; k < DIM; k++)
			strstr >> (p.v)[k];
		// add entry p to vector particles
		particles.push_back(p);
    }
    number_of_particles = particles.size();
    // close file
    parfile.close();
}



void World::print(std::ostream& os) const {
    os << "t=" << t << " delta_t=" << delta_t << " t_end=" << t_end
       << " Number of Particles=" << particles.size() << std::endl;
    for(int d = 0; d < DIM; d++)
        os << "length_x[" << d << "] = " << length_x[d] << std::endl;
    for(int d = 0; d < DIM; d++)
        os << "upper[" << d << "] = " << upper[d] << std::endl;
    for(int d = 0; d < DIM; d++)
        os << "lower[" << d << "] = " << lower[d] << std::endl;
}



void World_LC::print(std::ostream& os) const {
    World::print(os);
    os << "Cells per dim: ";
    for(int d = 0; d < DIM; d++)
        os << cell_N[d] << "\t";
    os << std::endl;

    os << "Length per cell: ";
    for(int d = 0; d < DIM; d++)
        os << cell_length[d] << "\t";
    os << std::endl;

    os << "Total number of cells: " << global_cell_number << std::endl;
}



std::ostream& operator<<(std::ostream& os, World_LC& W) {
    W.print(os);
    return os;
}



std::ostream& operator<<(std::ostream& os, World& W) {
    W.print(os);
    return os;
}



World_LC::World_LC() : World(), T() 
{
    // empty constructor
}



void World_LC::read_Parameter(const std::string& filename) {

    // Call function from baseclass for all member parameters of world
    World::read_Parameter(filename);
    T.gamma = delta_t;
    // continue with parameters only known to world_lc

    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("read_Parameter(): Can't open file '" +
                                 filename + "' for reading.");
    // helper strings
    std::string line, option;

    // read file till eof
    while (parfile.good()) {
        // init empty option
        option = "";
        // read line from file
        getline(parfile, line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option and value from stringstream
        strstr >> option;
        // check options and read values
        if (option == "cell_r_cut") { // Currently only one option.
            strstr >> cell_r_cut;
        } else if (option == "set_start_temperature") {
            strstr >> T.start_temperature ;
        } else if (option == "thermostat_step_interval") {
            strstr >> T.thermostat_interval;
            T.thermostat_interval *= delta_t;
        } else if (option == "thermostat_target_temperature") {
            strstr >> T.target_temperature;
        } else if (option == "random_seed") {
            strstr >> T.maxwell_boltzmann_seed;
        }
    }

    parfile.close();



    // computing:
    // - total number of cells: global_cell_number
    // - number of cells in each dimension: cell_N[d]
    // - size of each cell: cell_length[d]

    global_cell_number = 1;

    for(int d = 0; d < DIM; d++) {
        // number of possible cells with given r_cut
        cell_N[d] = std::max(1, int(std::floor(length_x[d] / cell_r_cut)));
        // length of cell per dim.
        cell_length[d] = length_x[d] / (double(cell_N[d]));

        global_cell_number *= cell_N[d];
    }


    // initialize appropriate number of cells in vector.
    for(int i = 0; i < global_cell_number; i++) {
        cells.push_back(Cell());
    }
}



void World_LC::read_Particles(const std::string& filename) {
    // Get all particles from file and store them in World
    World::read_Particles(filename);

    for(unsigned i = 0; i < particles.size(); i++) {
        particles[i].active = true;
    }

    initialize_Cells();
}



void World_LC::initialize_Cells(void) {
    for(unsigned i = 0; i < particles.size(); i++) {
        cells[linear_cell_index(&particles[i])].particles.push_back(particles[i]);
    }

    for(unsigned i = 0; i < cells.size(); i++) {
        cells[i].neighbors = cell_neighbors(i);
    }
}



void World_LC::sort_particles_in_cells(void) {
    unsigned global_cell_index;

    for(unsigned n = 0; n < cells.size(); n++) {
        for(auto i = cells[n].particles.begin(); i != cells[n].particles.end();) {

            if((*i).active) {
                // we compute the index of the cell in the vector cells
                global_cell_index = linear_cell_index(&(*i));

                if(global_cell_index != n) {
                    cells[global_cell_index].particles.push_back(*i);
                    i = cells[n].particles.erase(i);
                } else {
                    i++;
                }

            } else {
                i = cells[n].particles.erase(i);
                number_of_particles -= 1;
            }
        }
    }
}



std::vector<unsigned> World_LC::cell_neighbors(const int linear_cell_index) {
    std::vector<unsigned> neighbor_indices;
    std::unordered_set<unsigned> indices;
    
    int cell_index[DIM];
    int neighbor[DIM];
    int neighbor_index;

    compute_cell_index(linear_cell_index, cell_index);

    for(int i = 0; i < int(pow(3, DIM)); i++) {
        for(int d = 0; d < DIM; d++) {
            neighbor[d] = cell_index[d] - 1 + int( floor(i / pow(3, d)) ) % 3;
            // positive modulo
            neighbor[d] = (neighbor[d] % cell_N[d] + cell_N[d]) % cell_N[d];
        }

        neighbor_index = compute_global_cell_index(neighbor);

        if(neighbor_index != linear_cell_index) {
            indices.insert(neighbor_index);
        }
    }

    neighbor_indices.insert(neighbor_indices.end(), indices.begin(), indices.end());
    return neighbor_indices;
}



void World_LC::compute_cell_index(int linear_cell_index, int (&cell_index)[DIM]) {
    int lcsi = linear_cell_index;
    cell_index[DIM - 1] = lcsi % (cell_N[DIM - 1]);

    for(int d = DIM - 2; d >= 0; d--) {
        lcsi -= cell_index[d + 1];
        lcsi = (int) lcsi / cell_N[d + 1];
        cell_index[d] = lcsi % cell_N[d];
    }
}



int World_LC::compute_global_cell_index(const int (&cell_index)[DIM]) {
    int global_cell_index = cell_index[0];

    for(int j = 1; j <  DIM; j++) {
        global_cell_index *= cell_N[j]; 
        global_cell_index += cell_index[j];
    }

    return global_cell_index;
}



int World_LC::linear_cell_index(Particle *p) {
    int dim_index[DIM];
    int lin_index;

    for(int k = 0; k < DIM; k++) {
        dim_index[k] = int(floor(p->x[k] * cell_N[k] / length_x[k]));
    }

    lin_index = dim_index[0];

    for(int j = 1; j <  DIM; j++) {
        lin_index *= cell_N[j]; 
        lin_index += dim_index[j];
    }

    return lin_index;
}
