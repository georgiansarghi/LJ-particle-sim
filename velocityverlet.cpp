#include "velocityverlet.hpp"

#include <algorithm>
#include <cstring>



VelocityVerlet::VelocityVerlet(World& W, Potential& Pot, Observer& O)
    : TimeDiscretization(W, Pot, O) {
    // empty constructor
}



VelocityVerlet::VelocityVerlet(World& W, Potential* Pot, Observer& O)
    : TimeDiscretization(W, (*Pot), O) {
    // empty constructor
}



void VelocityVerlet::simulate() {
    std::cout << "Start." << std::endl;
    O.notify();
    compute_Force();

    while (W.t < W.t_end && W.particles.size() > 0) {
        timestep(W.delta_t);
    }
    
    std::cout << std::endl << "End." << std::endl;
}



void VelocityVerlet::handle_borders() {
    for(unsigned i = 0; i < W.particles.size(); i++) {
        for(int j = 0; j < DIM; j++) {
            if( (W.lower[j] == leaving && W.particles[i].x[j] < 0) 
             || (W.upper[j] == leaving && W.particles[i].x[j] > W.length_x[j]) ) {

                W.particles[i].active = false;
            
            } else if( W.lower[j] == periodic || W.upper[j] == periodic ) {
            
                W.particles[i].x[j] = std::fmod(W.particles[i].x[j] + W.length_x[j], W.length_x[j]);
         
            }
        }
    }
}



void VelocityVerlet::timestep(real delta_t) {
    std::cout << "\rTime: " << W.t << "        ";
	
    W.e_kin = 0.0;
    W.e_pot = 0.0;

    compute_Force();
    update_V();
    update_X();
    // compute total energy
    W.e_tot = W.e_kin + W.e_pot;
	
    // update particles, deleting from the vector all the particles who left the box
    // set active=false for particles out of the domain of simulation
    handle_borders();

    // increase time
    W.t += delta_t;
    // notify observer
    O.notify();
}



void VelocityVerlet::compute_Force() {
    real diff[DIM];

	for(unsigned i = 0; i < W.particles.size(); i++)	{
        memset(W.particles[i].F, 0, DIM * sizeof(real));
		// compute the potential energy
		for(unsigned j = 0; j < i; j++) {
            compute_diff(&W.particles[i], &W.particles[j], diff);
			W.e_pot += Pot.force(W.particles[i], W.particles[j], diff, 1, 1);
        }
		// compute remaining forces
		for(unsigned j = i+1; j < W.particles.size(); j++) {
            compute_diff(&W.particles[i], &W.particles[j], diff);
			Pot.force(W.particles[i], W.particles[j], diff, 1, 1);
        }
	}
}



void VelocityVerlet::update_V() {
    Particle *p;

	for(unsigned i = 0; i < W.particles.size(); i++) {
		for(unsigned j = 0; j < DIM; j++) {
            p = &W.particles[i];

			p->v[j] += W.delta_t  * .5/p->m * (p->F[j] + p->F_old[j]);
			// compute the kinetic energy
			W.e_kin += 0.5 * p->m * sqr(p->v[j]);
		}
    }
}



void VelocityVerlet::update_X() {
    Particle *p;

	for(unsigned i = 0; i < W.particles.size(); i++)	{
		for(unsigned j = 0; j < DIM; j++) {
            p = &W.particles[i];

			p->x[j] += W.delta_t * (p->v[j] + .5/p->m * (p->F[j] * W.delta_t)); 
			// update F_old
			p->F_old[j] = p->F[j];
		}
	}

}



real VelocityVerlet::compute_distance(Particle *p, Particle *q) {
    real d_sqr = 0;
    real xx[DIM];
    compute_diff(p, q, xx);

    for(unsigned k = 0; k < DIM; k++)
        d_sqr += sqr(xx[k]);

    return pow(d_sqr, 0.5);
}




// calculates (vectorial) difference of two particles and stores it in the array diff
void VelocityVerlet::compute_diff(Particle *p, Particle *q, real (&diff)[DIM]) {
    int s;

    for(unsigned k = 0; k < DIM; k++) {
        diff[k] = q->x[k] - p->x[k];
        if((W.lower[k] == periodic || W.upper[k] == periodic) && W.length_x[k] - abs(diff[k]) < abs(diff[k])) {
            // the sign of diff[k]
            s = ((diff[k] > 0) ? 1 : (diff[k] < 0) ? -1 : 0);
            diff[k] = s * (abs(diff[k]) - W.length_x[k]);
        }
    }
}



////////////////// Functions of Velocity_Verlet_LC


VelocityVerlet_LC::VelocityVerlet_LC(World_LC& W, Potential& Pot, Observer& O)
    : VelocityVerlet(W, Pot, O), W(W) {
    // empty constructor
}



VelocityVerlet_LC::VelocityVerlet_LC(World_LC& W, Potential* Pot, Observer& O)
    : VelocityVerlet(W, (*Pot), O), W(W) {
    // empty constructor
}



void VelocityVerlet_LC::simulate() {
    std::cout << "Start LC." << std::endl;
    O.notify();

    // for all cells itialize compute force
    
    for(unsigned cell_index = 0; cell_index < W.cells.size(); cell_index++) {
        compute_Force(cell_index);
    }

    // run simulation 
    while (W.t < W.t_end) {
        timestep(W.delta_t);
    }

    std::cout << std::endl << "End." << std::endl;
}



void VelocityVerlet_LC::timestep(real delta_t) {
    std::cout << "\rTime: " << W.t << "        ";

    W.e_kin = 0.0;
    W.e_pot = 0.0;
	
    // For all cells do integration steps.
    // Need several loops for correct order.
    // (Forces depend on all x updated!)

    // This updates only coordinates.

for(unsigned cell_index = 0; cell_index < W.cells.size(); cell_index++) {
        update_X(cell_index);
    }

	// update particles, deleting from the vector all the particles who left the box
    // set active=false for particles out of the domain of simulation
    handle_borders();

    // This updates cells of moved particles.
    W.sort_particles_in_cells();

    
    // Not called with iterator because we actually need the index!
    for(unsigned cell_index = 0; cell_index < W.cells.size(); cell_index++) {
            compute_Force(cell_index);
        }


    for(auto cell_it = W.cells.begin(); cell_it != W.cells.end(); cell_it++) {
        if(!(cell_it->particles.empty())) {
            update_V(cell_it->particles);
        }
    }
    

    W.T.execute_thermostat(W.t, W.cells, W.number_of_particles, W.e_kin);

	
	// increase number of timesteps
	W.timestep++;
	
	// compute average potential, kinetic and total energy over the last "value of W.out_interval" timesteps
	if (W.timestep < W.out_interval) {
		W.e_pot_list.push_back(W.e_pot);
		W.e_kin_list.push_back(W.e_kin);
		W.e_pot_avg += W.e_pot;
		W.e_kin_avg += W.e_kin;
	}
	else if (W.timestep == W.out_interval) {
		W.e_pot_list.push_back(W.e_pot);
		W.e_kin_list.push_back(W.e_kin);
		W.e_pot_avg = (W.e_pot_avg + W.e_pot)/W.out_interval;
		W.e_kin_avg = (W.e_kin_avg + W.e_kin)/W.out_interval;
	}
	else {
		W.e_pot_list.push_back(W.e_pot);
		W.e_kin_list.push_back(W.e_kin);
		W.e_pot_avg = (W.out_interval*W.e_pot_avg - W.e_pot_list.front() + W.e_pot)/W.out_interval;
		W.e_kin_avg = (W.out_interval*W.e_kin_avg - W.e_kin_list.front() + W.e_kin)/W.out_interval;
		W.e_pot_list.pop_front();
		W.e_kin_list.pop_front();
	}

    // compute total energy
    W.e_tot = W.e_kin + W.e_pot;
	
	// compute average total energy
	W.e_tot_avg = W.e_pot_avg + W.e_kin_avg;

    // compute average total energy
	W.e_tot_avg = W.e_pot_avg + W.e_kin_avg;
    
    // increase time
    W.t += delta_t;

    // notify observer
    O.notify();
}



void VelocityVerlet_LC::compute_Force(const unsigned int cell_index) {

    compute_Force_from_same_cell(cell_index);
    
    for(unsigned n = 0; n < W.cells[cell_index].neighbors.size(); n++) {
        compute_Force_from_neighboring_cell(cell_index, W.cells[cell_index].neighbors[n]);
    }    
}


// exact copy of update_V for the whole box, restricted to a cell
void VelocityVerlet_LC::update_V(std::vector<Particle>& particles) {
    Particle p;

    for(unsigned i = 0; i < particles.size(); i++)
        for(unsigned j = 0; j < DIM; j++) {
            p = particles[i];

            p.v[j] += W.delta_t  * .5/p.m * (p.F[j] + p.F_old[j]);
            // compute the kinetic energy
            W.e_kin += 0.5 * p.m * sqr(p.v[j]);
            // update velocity in global vector of particles
        }

}


// exact copy of update_X for  the whole box, resticted to a cell, we update the particles who moved in sort_particles_in cells
void VelocityVerlet_LC::update_X(const unsigned int cell_index) {
    Particle p;

    for(unsigned i = 0; i < W.cells[cell_index].particles.size(); i++)	{
        for(unsigned j = 0; j < DIM; j++) {
            p = W.cells[cell_index].particles[i];

            // update position
            p.x[j] += W.delta_t * (p.v[j] + .5/p.m * (p.F[j] * W.delta_t));

            // update F_old
            p.F_old[j] = p.F[j];
        }
    }
} 



void VelocityVerlet_LC::compute_Force_from_neighboring_cell(const unsigned int cell_index, int other_cell_index) {
    real diff[DIM];
    Particle p, q;
    real d;

   for(unsigned i = 0; i < W.cells[cell_index].particles.size(); i++) {
        for(unsigned j = 0; j < W.cells[other_cell_index].particles.size(); j++) {

            p = W.cells[other_cell_index].particles[j];
            q = W.cells[cell_index].particles[i];

            d = compute_distance(&p, &q);

            if(d < W.cell_r_cut) {
                compute_diff(&q, &p, diff);
                W.e_pot += 0.5 * Pot.force(q, p, diff, 1, 1);
            }
        }
    }
}

void VelocityVerlet_LC::compute_Force_from_same_cell(const unsigned int cell_index) {
    real diff[DIM];
    Particle p, q;
    real d;

    for( unsigned i = 0; i < W.cells[cell_index].particles.size(); i++)	{

        // initialize force vector for every particle
        memset(W.cells[cell_index].particles[i].F, 0, DIM * sizeof(real));

        // compute the potential energy and forces
        for( unsigned j = 0; j < W.cells[cell_index].particles.size(); j++) {
            // compute distance between particles
            if(i != j) {
                p = W.cells[cell_index].particles[j];
                q = W.cells[cell_index].particles[i];

                d = compute_distance(&p, &q);
                // checks if the distance between particles is sfficiently small
                if (d < W.cell_r_cut) {
                    compute_diff(&q, &p, diff);
                    W.e_pot += 0.5 * Pot.force(q, p, diff, 1, 1);
                }
            }
        }
    }
}



////////////////// Functions of Velocity_Verlet_Sub


VelocityVerlet_Sub::VelocityVerlet_Sub(SubDomain& W, Potential& Pot,
                                       Observer& O)
    : VelocityVerlet_LC(W, Pot, O), W(W) {
    // empty constructor
}



VelocityVerlet_Sub::VelocityVerlet_Sub(SubDomain& W, Potential* Pot,
                                       Observer& O)
    : VelocityVerlet_LC(W, (*Pot), O), W(W) {
    // empty constructor
}


// set active flag of particles which leave the domain to false
// only relevant for leaving boundaries
void VelocityVerlet_Sub::handle_borders() {
    for(unsigned i = 0; i < W.cells.size(); i++) {
        for(unsigned j = 0; j < W.cells[i].particles.size(); j++) {
            for(int d = 0; d < DIM; d++) {
                if( (W.lower[d] == leaving && W.cells[i].particles[j].x[d] < 0)
                 || (W.upper[d] == leaving && W.cells[i].particles[j].x[d] > W.length_x[d]) ) {

                    W.cells[i].particles[j].active = false;

                }
            }
        }
    }
}



void VelocityVerlet_Sub::simulate() {
    std::cout << "Start LC." << std::endl;
    O.notify();

    // helper iterator for cells
    unsigned int cell_index_iterator;
    // for all cells itialize compute force
    // no need to create and clear ghost since we already initialized an enlarged subdomain

    send_inner_receive_outer();


    for (cell_index_iterator = 0; cell_index_iterator < W.cells.size(); cell_index_iterator++) {
        for(unsigned i = 0; i < W.cells[cell_index_iterator].particles.size(); i++)  {
            // initialize force vector for every particle
            memset(W.cells[cell_index_iterator].particles[i].F, 0, DIM * sizeof(real));
        }
    }


    for (cell_index_iterator = 0; cell_index_iterator < W.cells.size(); cell_index_iterator++) {
        if(not W.cells[cell_index_iterator].ghost_flag) {
            compute_Force(cell_index_iterator);
        }
    }

    clear_ghosts();

    // while simulation end time not reached
    
    while (W.t < W.t_end) {

        if (W.myrank == 0)
            std::cout << "\rTime: " << W.t << "        ";

        W.local_e_kin = 0.0;
        W.local_e_pot = 0.0;


        timestep(W.delta_t);
		
        // updates W.local_e_kin and W.local_e_pot
        W.collect_global_data();

        W.T.execute_thermostat(W.t, W.cells, W.number_of_particles, W.local_e_kin);

        // compute total energy
        if (W.myrank == 0){
            W.e_tot = W.e_kin + W.e_pot;
        }

        //increase time
        W.t += W.delta_t;

        //notify observer
        O.notify();
    }

    std::cout << std::endl << "End." << std::endl;
}

void VelocityVerlet_Sub::timestep(real delta_t) {
	
    // Required in multiple loops.
    std::vector<Cell>::iterator cell_iterator;
    unsigned int cell_index_iterator;

    // For all cells do integration steps.
    // Need several loops for correct order.
    // (Forces depend on all x updated!)

    // This updates only coordinates.
    for (cell_index_iterator = 0; cell_index_iterator < W.cells.size();
         cell_index_iterator++) {
        W.cells[cell_index_iterator].force_flag = true;


        if (!(W.cells[cell_index_iterator].ghost_flag)) {
            update_X(cell_index_iterator);
        }
    }

    // deal with leaving case
    handle_borders();

    // This updates cells of moved particles prior to sending.
    W.sort_particles_in_cells();

    send_outer_receive_inner();
	
	// clear ghosts before sending innner particles
    clear_ghosts();
	
    send_inner_receive_outer();
	
    // Not called with iterator because we actually need the index!
    // Create ghosts only for compute_Force
    for (cell_index_iterator = 0; cell_index_iterator < W.cells.size();
         cell_index_iterator++) {
        if (not (W.cells[cell_index_iterator].ghost_flag)) {
            compute_Force(cell_index_iterator);
        }
    }
    
    // Clear ghosts immediately after compute force.
    clear_ghosts();

    // Now update velocities. No extra handling.
    for (cell_iterator = W.cells.begin(); cell_iterator < W.cells.end(); cell_iterator++) {
        if (not cell_iterator->particles.empty() && not cell_iterator->ghost_flag)
            update_V(cell_iterator->particles);
    }
}

void VelocityVerlet_Sub::clear_ghosts() {
    
    int         d1, d2, bound[DIM-1] = {0, 0},
                index_lower = 0, index_lower_d[DIM] = {0, 0, 0},
                index_upper = 0, index_upper_d[DIM] = {0, 0, 0};
				
	// loop over outer boundary in dimension d			
    for (unsigned d = 0; d < DIM; d++) {
        /*
        d1 = (d + 1) % DIM;
        d2 = (d + 2) % DIM;
        */
        if (d == 0) {d1 = 1; d2 = 2;}
        else if (d == 1) {d1 = 0; d2 = 2;}
        else {d1 = 0; d2 = 1;}
		
        index_lower_d[d] = 0;
        index_upper_d[d] = W.ic_number[d] + 1;
		
		// loop over outer boundary in dimension d1
        for (int i = 1 - bound[0]; i < W.ic_number[d1] + 1 + bound[0]; i++)
			// loop over outer boundary in dimension d2
            for (int j = 1 - bound[1]; j < W.ic_number[d2] + 1 + bound[1]; j++) {

                index_lower_d[d1] = index_upper_d[d1] = i;
                index_lower_d[d2] = index_upper_d[d2] = j;

                index_lower = W.compute_local_cell_single_index(index_lower_d);
                index_upper = W.compute_local_cell_single_index(index_upper_d);

                if (!W.cells[index_lower].ghost_flag) std::cout << "wrong: the cell should be on outer boundary" <<std::endl;
                if (!W.cells[index_upper].ghost_flag) std::cout << "wrong: the cell should be on outer boundary" <<std::endl;

                // clear lower d component of ghost (checking if there is a lower neighbor)
                if (W.ip_lower[d] >= 0) {
                    W.cells[index_lower].particles.clear(); 
                }
                // clear upper d component of ghost (checking if there is an upper neighbor)
                if (W.ip_upper[d] >= 0) {
                    W.cells[index_upper].particles.clear();
                }
            }
        if (bound[0] == 0) bound[0]++;
        else if (bound[1] == 0) bound[1]++;
    }
}

// sort particles in the vector received_particles in cells
void VelocityVerlet_Sub::create_ghosts(int d) {
    Particle p;
    int index;
    real c;	// displacement in periodic case

    for (unsigned i = 0; i < W.received_particles.size(); i++) {
        p = W.received_particles[i];

		if(W.lower[d] == periodic && W.upper[d] == periodic) {
			c = W.ic_start[d] * W.length_x[d] / W.cell_N[d];

			if(p.x[d] > W.length_x[d] - c && W.ip[d] == 0) {

				p.x[d] -= W.length_x[d];

			} else if(p.x[d] < c && W.ip[d] == W.N_p[d] - 1) {

				p.x[d] += W.length_x[d];
			}
		}

        index = W.local_linear_cell_index(&p);

        if(not W.cells[index].ghost_flag) {
            W.local_particle_num++;
        }

        p.active = true;
        W.cells[index].particles.push_back(p);

    }

    W.received_particles.clear();
}



void VelocityVerlet_Sub::send_outer_receive_inner() {

	//////// variables needed for sending
    // sending buffers: store particles which are to be sent to neighboring processes
    std::vector <Particle> lower_buffer;
    std::vector <Particle> upper_buffer;
	// sizes of sending buffers
	unsigned vsu, vsl;
	
	//// variables needed for receiving
	// receiving buffers: store particles which are received from neighboring processes
    std::vector <Particle> upper_data;
    std::vector <Particle> lower_data;
	// sizes of receiving buffers
	unsigned su, sl;

	// helper particle
    Particle p;
	
	// MPI request variables
    MPI_Request lower_request, upper_request, 
                lower_s_request, upper_s_request;    

	// variables needed in loop over outer cells
	int     d1, d2, bound[DIM-1] = {1, 1},
			index_lower = 0, index_lower_d[DIM] = {0, 0, 0},
			index_upper = 0, index_upper_d[DIM] = {0, 0, 0};


    //////// assemble upper buffer and lower buffer for each dimension
    //////// moreover update the vector of received particles for each dimension
	// loop over inner boundary in dimension d
    for (int d = DIM - 1; d >= 0; d--) {
		// sending buffers
        upper_buffer.clear();
        lower_buffer.clear();

        // compute other indices
        /*
        d1 = (d + 1) % DIM;
        d2 = (d + 2) % DIM;
        */
        if (d == 0) {d2 = 1; d1 = 2;}
        else if (d == 1) {d2 = 0; d1 = 2;}
        else {d2 = 0; d1 = 1;}
		
        index_lower_d[d] = 0;
        index_upper_d[d] = W.ic_number[d] + W.ic_start[d];
	
		// loop over inner boundary in dimension d1
        for (int i = 1 - bound[0]; i < W.ic_number[d1] + 1 + bound[0]; i++) {
			// loop over inner boundary in dimension d2
            for (int j = 1 - bound[1]; j < W.ic_number[d2] + 1 + bound[1]; j++) {

                index_lower_d[d1] = index_upper_d[d1] = i;
                index_lower_d[d2] = index_upper_d[d2] = j;

                index_lower = W.compute_local_cell_single_index(index_lower_d);
                index_upper = W.compute_local_cell_single_index(index_upper_d);

                // assemble lower buffer (checking if there is a lower neighbor)
                if (W.ip_lower[d] >= 0)
                    for (unsigned k = 0; k < W.cells[index_lower].particles.size(); k++) {
                        p = W.cells[index_lower].particles[k];
                        lower_buffer.push_back(p);
                    }                    
                // assemble upper buffer (checking if there is an upper neighbor)
                if (W.ip_upper[d] >= 0)
                    for (unsigned k = 0; k < W.cells[index_upper].particles.size(); k++) {
                        p = W.cells[index_upper].particles[k];
                        upper_buffer.push_back(p);
                    }
            }
        }


		//// non-blocking receiving of buffer sizes
        if (W.ip_lower[d] >= 0) {
            MPI_Irecv(&sl,  1, MPI_UNSIGNED, W.ip_lower[d], 
                            1, MPI_COMM_WORLD, &lower_s_request);
        } else {
            sl = 0; 
            lower_data.empty();
            lower_s_request = MPI_REQUEST_NULL;
            lower_request = MPI_REQUEST_NULL;
        }
    
        if (W.ip_upper[d] >= 0) {
            MPI_Irecv(&su,                   1, MPI_UNSIGNED, W.ip_upper[d], 
                                                1, MPI_COMM_WORLD, &upper_s_request);
        } else {
            su = 0; 
            upper_data.empty();
            upper_s_request = MPI_REQUEST_NULL;
            upper_request = MPI_REQUEST_NULL;
        }


		//// blocking sending of buffer sizes
        if (W.ip_lower[d] >= 0) {

            vsl = lower_buffer.size();

            MPI_Send(&vsl,                1, MPI_UNSIGNED, W.ip_lower[d], 
                                        1, MPI_COMM_WORLD);
        }

        if (W.ip_upper[d] >= 0) {

            vsu = upper_buffer.size();
            MPI_Send(&vsu,                1, MPI_UNSIGNED, W.ip_upper[d], 
                                        1, MPI_COMM_WORLD);
        }
		

		//// wait until Irecv of buffer sizes is finished
        MPI_Wait(&lower_s_request, MPI_STATUS_IGNORE);
        MPI_Wait(&upper_s_request, MPI_STATUS_IGNORE);


		//// non-blocking receiving of data
        if (W.ip_lower[d] >= 0) {
            // resizing
            lower_data.resize(sl);
            MPI_Irecv(&lower_data.front(),   sl, W.MPI_Particle, W.ip_lower[d],
                                            1, MPI_COMM_WORLD, &lower_request);
        }
        
        if (W.ip_upper[d] >= 0) {
            // resizing
            upper_data.resize(su);
            MPI_Irecv(&upper_data.front(),   su, W.MPI_Particle, W.ip_upper[d],
                                            1, MPI_COMM_WORLD, &upper_request);
        }


		//// blocking sending data
        if (W.ip_lower[d] >= 0)
            // resizing
            MPI_Send(&lower_buffer.front(),     vsl, W.MPI_Particle, W.ip_lower[d],
                                    1, MPI_COMM_WORLD);
        

        if (W.ip_upper[d] >= 0)
            // resizing
            MPI_Send(&upper_buffer.front(),     vsu, W.MPI_Particle, W.ip_upper[d],
                                    1, MPI_COMM_WORLD);


		//// waiting until Irecv of buffers is finished
        if (MPI_Wait(&lower_request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            std::cout << "Waitall lower: failed" << std::flush << std::endl;

        if (MPI_Wait(&upper_request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            std::cout << "Waitall upper: failed" << std::flush << std::endl;


		//// sort received particles in vector received_particles
        if (W.ip_lower[d] >= 0)
            for (unsigned i = 0; i < sl; i++) {
                W.received_particles.push_back(lower_data[i]);
            }
        
        if (W.ip_upper[d] >= 0)
            for (unsigned i = 0; i < su; i++) {
                W.received_particles.push_back(upper_data[i]);
            }

		// clear receiving buffers
        upper_data.clear();
        lower_data.clear();

		// sort received particles in inner boundary cells
        create_ghosts(d);

        // increase size of information that will be send on next iteration
        // this is to do once we fixed the (out of range) cell indices in the subdomain
        // on the ghost (outer boundary of subdomain)
        if (bound[1] == 1) bound[1]--;
        else if (bound[0] == 1) bound[0]--;
    }
}


// inverse to send_outer_receive_inner()
void VelocityVerlet_Sub::send_inner_receive_outer() {

    //////// variables needed for sending
    // sending buffers: store particles which are to be sent to neighboring processes
    std::vector <Particle> lower_buffer;
    std::vector <Particle> upper_buffer;
	// sizes of sending buffers
	unsigned vsu, vsl;
	
	//// variables needed for receiving
	// receiving buffers: store particles which are received from neighboring processes
    std::vector <Particle> upper_data;
    std::vector <Particle> lower_data;
	// sizes of receiving buffers
	unsigned su, sl;

	// helper particle
    Particle p;
	
	// MPI request variables
    MPI_Request lower_request, upper_request, 
                lower_s_request, upper_s_request;    

	// variables needed in loop over inner cells
	int     d1, d2, bound[DIM-1] = {1, 1},
			index_lower = 0, index_lower_d[DIM] = {0, 0, 0},
			index_upper = 0, index_upper_d[DIM] = {0, 0, 0};


    //////// assemble upper buffer and lower buffer for each dimension
    //////// moreover update the vector of received particles for each dimension
	// loop over inner boundary cells in dimension d
    for (int d = 0; d < DIM; d++) {
        upper_buffer.clear();
        lower_buffer.clear();

        // compute other indices
        /*
        d1 = (d + 1) % DIM;
        d2 = (d + 2) % DIM;
        */
        if (d == 0) {d1 = 1; d2 = 2;}
        else if (d == 1) {d1 = 0; d2 = 2;}
        else {d1 = 0; d2 = 1;}

        index_lower_d[d] = 1;
        index_upper_d[d] = W.ic_number[d];
		// loop over inner boundary cells in dimension d1
        for (int i = 1 - bound[0]; i < W.ic_number[d1] + 1 + bound[0]; i++) {
			// loop over inner boundary cells in dimension d1
            for (int j = 1 - bound[1]; j < W.ic_number[d2] + 1 + bound[1]; j++) {

                index_lower_d[d1] = index_upper_d[d1] = i;
                index_lower_d[d2] = index_upper_d[d2] = j;

                index_lower = W.compute_local_cell_single_index(index_lower_d);
                index_upper = W.compute_local_cell_single_index(index_upper_d);

                // assemble lower buffer (checking if there is a lower neighbor)
                if (W.ip_lower[d] >= 0)
                    for (unsigned k = 0; k < W.cells[index_lower].particles.size(); k++) {
                        p = W.cells[index_lower].particles[k];
                        lower_buffer.push_back(p);
                    }                    
                // assemble upper buffer (checking if there is an upper neighbor)
                if (W.ip_upper[d] >= 0)
                    for (unsigned k = 0; k < W.cells[index_upper].particles.size(); k++) {
                        p = W.cells[index_upper].particles[k];
                        upper_buffer.push_back(p);
                    }
            }
        }

		//// non-blocking receiving of buffer sizes
        if (W.ip_lower[d] >= 0) {
            MPI_Irecv(&sl,  1, MPI_UNSIGNED, W.ip_lower[d], 
                            1, MPI_COMM_WORLD, &lower_s_request);
        } else {
            sl = 0; 
            lower_data.empty();
            lower_s_request = MPI_REQUEST_NULL;
            lower_request = MPI_REQUEST_NULL;
        }
    
        if (W.ip_upper[d] >= 0) {
            MPI_Irecv(&su,                   1, MPI_UNSIGNED, W.ip_upper[d], 
                                                1, MPI_COMM_WORLD, &upper_s_request);
        } else {
            su = 0; 
            upper_data.empty();
            upper_s_request = MPI_REQUEST_NULL;
            upper_request = MPI_REQUEST_NULL;
        }

		//// blocking sending of buffer sizes
        if (W.ip_lower[d] >= 0) {
            vsl = lower_buffer.size();

            MPI_Send(&vsl,                1, MPI_UNSIGNED, W.ip_lower[d], 
                                        1, MPI_COMM_WORLD);
        }

        if (W.ip_upper[d] >= 0) {

            vsu = upper_buffer.size();
            MPI_Send(&vsu,                1, MPI_UNSIGNED, W.ip_upper[d], 
                                        1, MPI_COMM_WORLD);
        }


		//// wait until Irecv of buffer sizes is finished
        MPI_Wait(&lower_s_request, MPI_STATUS_IGNORE);
        MPI_Wait(&upper_s_request, MPI_STATUS_IGNORE);


		//// non-blocking receiving of data
        if (W.ip_lower[d] >= 0) {
            // resizing
            lower_data.resize(sl);
            MPI_Irecv(&lower_data.front(),   sl, W.MPI_Particle, W.ip_lower[d],
                                            1, MPI_COMM_WORLD, &lower_request);
        }
        
        if (W.ip_upper[d] >= 0) {
            // resizing
            upper_data.resize(su);
            MPI_Irecv(&upper_data.front(),   su, W.MPI_Particle, W.ip_upper[d],
                                            1, MPI_COMM_WORLD, &upper_request);
        }

		//// blocking sending of data
        if (W.ip_lower[d] >= 0)
            // resizing
            MPI_Send(&lower_buffer.front(),     vsl, W.MPI_Particle, W.ip_lower[d],
                                    1, MPI_COMM_WORLD);
        

        if (W.ip_upper[d] >= 0)
            // resizing
            MPI_Send(&upper_buffer.front(),     vsu, W.MPI_Particle, W.ip_upper[d],
                                    1, MPI_COMM_WORLD);


		//// wait until Irecv of data is finished
        if (MPI_Wait(&lower_request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            std::cout << "Waitall lower: failed" << std::flush << std::endl;

        if (MPI_Wait(&upper_request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
            std::cout << "Waitall upper: failed" << std::flush << std::endl;

        //// sort received particles in vector received_particles
        if (W.ip_lower[d] >= 0)
            for (unsigned i = 0; i < sl; i++) {
                W.received_particles.push_back(lower_data[i]);
            }
        
        if (W.ip_upper[d] >= 0)
            for (unsigned i = 0; i < su; i++) {
                W.received_particles.push_back(upper_data[i]);
            }

		// clear receiving buffers
        upper_data.clear();
        lower_data.clear();

		// sort received particles in ghost cells (= outer boundary cells)
        create_ghosts(d);

        // increase size of information that will be send on next iteration
        // this is to do once we fixed the (out of range) cell indices in the subdomain
        // on the ghost (outer boundary of subdomain)
        if (bound[0] == 0) bound[0]++;
        else if (bound[1] == 0) bound[1]++;
    }
}



// exact copy of update_V for the whole box, restricted to a cell
void VelocityVerlet_Sub::update_V(std::vector<Particle>& particles) {
    Particle * p;

    for(unsigned i = 0; i < particles.size(); i++)
        for(unsigned j = 0; j < DIM; j++) {
            p = &particles[i];
            
            p->v[j] += W.delta_t  * .5/p->m * (p->F[j] + p->F_old[j]);

            //compute the kinetic energy
            W.local_e_kin += 0.5 * p->m * sqr(p->v[j]);
            // update velocity in global vector of particles
        }

}


// exact copy of update_X for  the whole box, resticted to a cell, we update the particles who moved in sort_particles_in cells
void VelocityVerlet_Sub::update_X(const unsigned int cell_index) {
    Particle * p;

    for(unsigned i = 0; i < W.cells[cell_index].particles.size(); i++)  {
        if(not W.cells[cell_index].ghost_flag) {
            for(unsigned j = 0; j < DIM; j++) {
                p = &W.cells[cell_index].particles[i];

                // update position
                p->x[j] += W.delta_t * (p->v[j] + .5/p->m * (p->F[j] * W.delta_t));

                // update F_old
                p->F_old[j] = p->F[j];
            }
        }

        // initialize force vector for every particle
        memset(W.cells[cell_index].particles[i].F, 0, DIM * sizeof(real));
    }
} 


// force computation
void VelocityVerlet_Sub::compute_Force(const unsigned int cell_index) {

    compute_Force_from_same_cell(cell_index);
    
    for(unsigned n = 0; n < W.cells[cell_index].neighbors.size(); n++) {
        if(W.cells[W.cells[cell_index].neighbors[n]].force_flag) {
            compute_Force_from_neighboring_cell(cell_index, W.cells[cell_index].neighbors[n]);
        }
    }

    W.cells[cell_index].force_flag = false;
}

// force computation with particles from neighboring cell
void VelocityVerlet_Sub::compute_Force_from_neighboring_cell(const unsigned int cell_index, int other_cell_index) {
    real diff[DIM];
    Particle *p, *q;
    real dist;

   for(unsigned i = 0; i < W.cells[cell_index].particles.size(); i++) {
        for(unsigned j = 0; j < W.cells[other_cell_index].particles.size(); j++) {

            p = &W.cells[other_cell_index].particles[j];
            q = &W.cells[cell_index].particles[i];

            dist = compute_distance(p, q);
			// checks if the distance between particles is sufficiently small
            if(dist < W.cell_r_cut) {
                compute_diff(q, p, diff);

                if(W.cells[other_cell_index].ghost_flag)
                    W.local_e_pot += 0.5 * Pot.force(*q, *p, diff, W.sigma, W.eps);
                else
                    W.local_e_pot += Pot.force(*q, *p, diff, W.sigma, W.eps);
            }
        }
    }
}

// force computation with particles from same cell
void VelocityVerlet_Sub::compute_Force_from_same_cell(const unsigned int cell_index) {
    real diff[DIM];
    Particle *p, *q;
    real dist;

    for( unsigned i = 0; i < W.cells[cell_index].particles.size(); i++) {



        // compute the potential energy and forces
        for( unsigned j = 0; j < i; j++) {
            // compute distance between particles
            p = &W.cells[cell_index].particles[j];
            q = &W.cells[cell_index].particles[i];

            dist = compute_distance(p, q);   
            // checks if the distance between particles is sufficiently small
            if (dist < W.cell_r_cut) {
                compute_diff(q, p, diff);
                W.local_e_pot += Pot.force(*q, *p, diff, W.sigma, W.eps);
            }   
        }
    }
}


// computes (scalar) distance between two particles
real VelocityVerlet_Sub::compute_distance(Particle *p, Particle *q) {
    real d_sqr = 0;
    real diff[DIM];
    compute_diff(p, q, diff);

    for(unsigned k = 0; k < DIM; k++)
        d_sqr += sqr(diff[k]);

    return pow(d_sqr, 0.5);
}


// calculates (vectorial) difference of two particles and stores it in the array diff
void VelocityVerlet_Sub::compute_diff(Particle *p, Particle *q, real (&diff)[DIM]) {

    for(unsigned k = 0; k < DIM; k++) {
        diff[k] = q->x[k] - p->x[k];
    }
}
