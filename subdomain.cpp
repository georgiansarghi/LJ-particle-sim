#include "subdomain.hpp"
#include "particle.cpp"
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

SubDomain::SubDomain() : World_LC(), local_particle_num(0)
{
    construct_MPI_Particle();
}


void SubDomain::read_Parameter(const std::string& filename)
{
    // Call function from baseclass for all member parameters of world
    World_LC::read_Parameter(filename);

    // continue with parameters only known to subdomain

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
        if (option == "num_procs") {
            for (int d = 0; d < DIM; d++)
                strstr >> N_p[d];
        }
    }

	// compute number of local cells in each dimension
    for(int d = 0; d < DIM; d++) {
        local_cell_N[d] = (int) cell_N[d] / N_p[d];
    }

    // Close file.
    parfile.close();

    //// set derived values
    MPI_Comm_size(MPI_COMM_WORLD, &myprocessors);
	
    // set number of processors
    numprocs = 1;
    for (int d = 0; d < DIM; d++)
        numprocs *= N_p[d];

    if (numprocs != myprocessors) {
        std::cerr << "Processor distribution does not correspond to number of "
                     "processors!"
                  << std::endl;
        std::cerr << "You have " << myprocessors
                  << " processors, but " << numprocs << " are expected!"
                  << std::endl;
        exit(1);
    }
    // Get myrank
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Get my multiindex
    compute_processor_multi_index(myrank, ip);

    // Compute neighboring indices, set temp
    compute_neighboring_procs();

    // Compute Cells to be used
    init_local_cells_in_subdomain();
}


void SubDomain::print(std::ostream& os) const
{
    if (myrank == 0)
        World_LC::print(os);
    os << std::endl;

    os << "local_particle_num: ";
    os << local_particle_num << std::endl;
    os << "My rank:  ";
    os << myrank << std::endl;

    os << "My proc multiindex:  ";
    for (int d = 0; d < DIM; d++)
        os << ip[d] << "\t";
    os << std::endl;

    os << "Number of procs per dim:  ";
    for (int d = 0; d < DIM; d++)
        os << N_p[d] << "\t";
    os << std::endl;

    os << "Neighboring processes below:  ";
    for (int d = 0; d < DIM; d++)
        os << ip_lower[d] << "\t";
    os << std::endl;

    os << "Neighboring processes above:  ";
    for (int d = 0; d < DIM; d++)
        os << ip_upper[d] << "\t";
    os << std::endl;

    os << "ic_start:  ";
    for (int d = 0; d < DIM; d++)
        os << ic_start[d] << "\t";
    os << std::endl;

    os << "ic_stop:  ";
    for (int d = 0; d < DIM; d++)
        os << ic_stop[d] << "\t";
    os << std::endl;

    os << "ic_number:  ";
    for (int d = 0; d < DIM; d++)
        os << ic_number[d] << "\t";
    os << std::endl;

    os << "ic_lower_global:  ";
    for (int d = 0; d < DIM; d++)
        os << ic_lower_global[d] << "\t";
    os << std::endl;
}


std::ostream& operator<<(std::ostream& os, SubDomain& S)
{
    S.print(os);
    return os;
}


// compute the processor single index according to the formula
// proc_single_index = proc_index[0] + N_p[0] * ( proc_index[1] + N_p[1] * proc_index[2])
int SubDomain::compute_processor_single_index(const int (&proc_index)[DIM]) {
    int proc_single_index = proc_index[DIM - 1];

    for(int j = DIM - 2; j >= 0; j--) {
        proc_single_index *= N_p[j]; 
        proc_single_index += proc_index[j];
    }

    return proc_single_index;
}


// compute the processor multi index according to the formula
// processor_single_index = processor_index[0] + N_p[0] * ( processor_index[1] + N_p[1] * processor_index[2])
void SubDomain::compute_processor_multi_index(const int processor_single_index,
                                              int (&processor_index)[DIM]) {
    
    int psi = processor_single_index;
    processor_index[0] = psi % N_p[0];

    for(int d = 1; d < DIM; d++) {
        psi -= processor_index[d - 1];
        psi = (int) psi / N_p[d - 1];
        processor_index[d] = psi % N_p[d];
    }   
    
}


// modification due to enlarged domain
void SubDomain::compute_local_cell_multi_index(const int local_cell_single_index, int (&cell_index)[DIM]) 
{
    int lcsi = local_cell_single_index;
    cell_index[DIM - 1] = lcsi % (local_cell_N[DIM - 1] + 2);

    for(int d = DIM - 2; d >= 0; d--) {
        lcsi -= cell_index[d + 1];
        lcsi = (int) lcsi / (local_cell_N[d + 1] + 2);
        cell_index[d] = lcsi % (local_cell_N[d] + 2);
    }

}


// modification due to enlarged domain
int SubDomain::compute_local_cell_single_index(const int (&cell_index)[DIM]) {
    
    int cell_single_index = cell_index[0];

    for(int j = 1; j <  DIM; j++) {
        cell_single_index *= (local_cell_N[j] + 2); 
        cell_single_index += cell_index[j];
    }

    return cell_single_index;
    
}


// gets called in read_Particles
void SubDomain::init_particle_in_cells() { 

    for(unsigned i = 0; i < particles.size(); i++) {
        cells[local_linear_cell_index(&particles[i])].particles.push_back(particles[i]);
    }

    for(unsigned i = 0; i < cells.size(); i++) {
        if (!cells[i].ghost_flag)
            cells[i].neighbors = compute_cell_neighbors(i);
    }
}



// will handle also particles in the ghost during iteration zero
// it is called only when we already know if the particle is in the subdomain (or ghost)
int SubDomain::local_linear_cell_index(Particle *p) {

    int dim_index[DIM];
    int lin_index;

    for(int d = 0; d < DIM; d++) {
		
        if(ip[d] == N_p[d] - 1 && p->x[d] >= length_x[d]) {

            dim_index[d] = ic_number[d] + 1;
        
        } else if(ip[d] == 0 && p->x[d] < 0) {
        
            dim_index[d] = 0;
        
        } else {

            dim_index[d] = int(floor(p->x[d] * (double) cell_N[d] / length_x[d]));

            if (dim_index[d] == ic_lower_global[d] - 1) {
                
                dim_index[d] = 0;

            } else if (dim_index[d] == ic_lower_global[d] + ic_number[d]) {
                
                dim_index[d] = ic_number[d] + 1;

            } else {

                dim_index[d] %= ic_number[d];
                // we take count of the outer boundary
                dim_index[d] += 1;
            }
        }
    }

    lin_index = dim_index[0];

    for(int j = 1; j <  DIM; j++) {
        lin_index *= (ic_number[j] + 2);
        lin_index += dim_index[j];
    }

    return lin_index;
}


void SubDomain::compute_neighboring_procs() {
    int proc_index[DIM];

	// proc_index = ip
    compute_processor_multi_index(myrank, proc_index);

    for(int d = 0; d < DIM; d++) {
        if (lower[d] == leaving && ip[d] == 0) {
			
            ip_lower[d] = - 1;
			
        } else if (lower[d] == periodic && ip[d] == 0) {
			
			proc_index[d] = N_p[d] - 1;
            ip_lower[d] = compute_processor_single_index(proc_index);
			proc_index[d] = 0;
			
        } else {
			
			proc_index[d]--;
			ip_lower[d] = compute_processor_single_index(proc_index);
			proc_index[d]++;
			
		}

        if (upper[d] == leaving && ip[d] == N_p[d] - 1) {
			
            ip_upper[d] = - 1;    
			
        } else if(upper[d] == periodic && ip[d] == N_p[d] - 1) {
			
			proc_index[d] = 0;
            ip_upper[d] = compute_processor_single_index(proc_index);
			proc_index[d] = N_p[d] - 1;
			
        } else {
			
			proc_index[d]++;
			ip_upper[d] = compute_processor_single_index(proc_index);
			proc_index[d]--;
			
		}
    }
}


std::vector<unsigned> SubDomain::compute_cell_neighbors(const int linear_cell_index) {

    std::vector<unsigned> neighbor_indices;
    std::unordered_set<unsigned> indices;
    
    int cell_index[DIM];
    int neighbor[DIM];
    int neighbor_index;

    compute_local_cell_multi_index(linear_cell_index, cell_index);

    for(int i = 0; i < int(pow(3, DIM)); i++) {
        for(int d = 0; d < DIM; d++) {
            neighbor[d] = cell_index[d] - 1 + int( floor(i / pow(3, d)) ) % 3;
        }

        neighbor_index = compute_local_cell_single_index(neighbor);

        if(neighbor_index != linear_cell_index) {
            indices.insert(neighbor_index);
        }
    }

    neighbor_indices.insert(neighbor_indices.end(), indices.begin(), indices.end());
	
    return neighbor_indices;
}


void SubDomain::init_local_cells_in_subdomain()
{
	local_cell_num = 1;
    int total_cell_num = 1;
	// set initial data
	for(int d = 0; d < DIM; d++) {
		ic_start[d] = 1;
		ic_stop[d] = cell_N[d]/N_p[d] + ic_start[d];
		ic_number[d] = ic_stop[d] - ic_start[d];

        ic_lower_global[d] = ip[d] * cell_N[d] / N_p[d];

		local_cell_num *= ic_number[d];
	}
    
    // we compute the number of cells in ghost (ideally)make
    // some of them will be empty in the leaving bounary type case

    for(unsigned d = 0; d < DIM; d++)
        total_cell_num *= (ic_number[d] + 2*ic_start[d]);

    ghost_cell_num = total_cell_num - local_cell_num;

    cells.clear();
	
    // initialize appropriate number of cells in vector
    for(int i = 0; i < local_cell_num + ghost_cell_num; i++) {

        cells.push_back(Cell());

    }
    

    // we set a ghost_flag == true in every outer border cell
    
    int         index_lower = 0, index_lower_d[DIM] = {0, 0, 0},
                index_upper = 0, index_upper_d[DIM] = {0, 0, 0},
                d1, d2, increment[DIM - 1];
    for (int i = 0; i < DIM - 1; i++)
        increment[i] = 0;
    for (unsigned d = 0; d < DIM; d++) {
        /*
        d1 = (d + 1) % DIM;
        d2 = (d + 2) % DIM;
        */
        if (d == 0) {d1 = 1; d2 = 2;}
        else if (d == 1) {d1 = 0; d2 = 2;}
        else {d1 = 0; d2 = 1;}
        // loop over outer boundary
        index_lower_d[d] = 0;
        index_upper_d[d] = ic_number[d] + 1;
        for (int i = 1 - increment[0]; i < ic_number[d1] + 1 + increment[0]; i++)
            for (int j = 1 - increment[1]; j < ic_number[d2] + 1 + increment[1]; j++) {

                index_lower_d[d1] = index_upper_d[d1] = i;
                index_lower_d[d2] = index_upper_d[d2] = j;

                index_lower = compute_local_cell_single_index(index_lower_d);
                index_upper = compute_local_cell_single_index(index_upper_d);

                // set the flag on lower outer boundary
                cells[index_lower].ghost_flag = true;
                // set the flag on lower outer boundary
                cells[index_upper].ghost_flag = true;
            }
        if (increment[0] == 0) increment[0]++;
        else if (increment[1] == 0) increment[1]++;
    }


}


// reads the particles belonging to the subdomain (only inner for cells)
void SubDomain::read_Particles(const std::string& filename)
{

    std::ifstream parfile(filename.c_str());
    if (!parfile.is_open())
        throw std::runtime_error("read_Particles(): Can't open file '" +
                                 filename + "' for reading.");
	
	// helper string
    std::string line;
	
	// helper particle
	Particle p;

     // read file till eof
    while(parfile.good()) {

        // read line from file
        getline(parfile, line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;

        // read values from stringstream and write them in particles
        strstr >> p.id;
		// mass
		strstr >> p.m;
		// position
		for(int d = 0; d < DIM; d++)
			strstr >> p.x[d];
		// velocity
		for(int d = 0; d < DIM; d++) {
				strstr >> p.v[d];
        }
		
		// add particle p to vector particles if p is in subdomain
		bool flag = true;
		
		for (int d = 0; d <  DIM; d++) {
			if(     (p.x[d] <    ip[d]      * ((double)length_x[d] / N_p[d])) //- (double)length_x[d] / cell_N[d])
                ||   (p.x[d] >=  (ip[d] + 1) * ((double)length_x[d] / N_p[d]))) { //+ (double)length_x[d] / cell_N[d])) {
                flag = false;
            }
        }
        // init Force
        for(int d = 0; d < DIM; d++) {
                p.F[d] = 0;
        }
        
        // init Force_old
        for(int d = 0; d < DIM; d++) {
                p.F_old[d] = 0;
        }
		
        if (flag) {
            particles.push_back(p);
        }
    }

    local_particle_num = particles.size();

    // close file
    parfile.close();
	
	// set active flags for particles (only relevenat for leaving boundary conditions)
    for(unsigned i = 0; i < particles.size(); i++) {
        particles[i].active = true;
    }
	
	// initialize particles in cells 
	init_particle_in_cells();
}


void SubDomain::construct_MPI_Particle(void) {
    Particle p[2];

    // build new mpi datatype
    MPI_Datatype XD;

    MPI_Datatype type[5] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int blocklen[5]      = {1, 1, 1, DIM, DIM};

    MPI_Aint disp[5];  // displacements
    MPI_Aint base;
    MPI_Aint sizeofentry;

    /* compute displacements */
    MPI_Get_address(&p[0],      &disp[0]);
    MPI_Get_address(&p[0].id,   &disp[1]);
    MPI_Get_address(&p[0].m,    &disp[2]);
    MPI_Get_address(&p[0].x[0], &disp[3]);
    MPI_Get_address(&p[0].v[0], &disp[4]);

    base = disp[0];
    for (unsigned i = 0; i < 5; i++)
        disp[i] -= base;

    MPI_Type_create_struct(5, blocklen, disp, type, &XD);
    MPI_Get_address(p+1, &sizeofentry);
    sizeofentry = MPI_Aint_diff(sizeofentry, base); 
    MPI_Type_create_resized(XD, 0, sizeofentry, &MPI_Particle); 
    MPI_Type_commit(&MPI_Particle);
}


void SubDomain::sort_particles_in_cells() {

    unsigned subdomain_cell_index;

    for(unsigned n = 0; n < cells.size(); n++) {

        for(auto i = cells[n].particles.begin(); i != cells[n].particles.end();) {


            if(i->active) {
                // we compute the index of the cell in the vector cells
                subdomain_cell_index = local_linear_cell_index(&(*i));

                if(subdomain_cell_index != n) {
                    cells[subdomain_cell_index].particles.push_back(*i);
                    i = cells[n].particles.erase(i);
                } else {
                    i++;
                }

            } else {
                i = cells[n].particles.erase(i);
                local_particle_num -= 1;
            }
        }
    }
}


// compute global kinetic and potential energies
void SubDomain::collect_global_data() {

    MPI_Reduce(&local_e_kin, &e_kin, 1, MPI_DOUBLE,
                MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_e_pot, &e_pot, 1, MPI_DOUBLE,
                MPI_SUM, 0, MPI_COMM_WORLD);
}
