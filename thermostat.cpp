#include "thermostat.hpp"

Thermostat::Thermostat()
    : timestamp(0.0), target_temperature(-1.0), start_temperature(-1.0)
{
}

/// Empty destructor
Thermostat::~Thermostat() {}

void Thermostat::execute_thermostat(real time, std::vector<Cell>& cells,
                                    const int number_of_particles,
                                    const real e_kin)
{
    if (target_temperature < 0.0)
        return;
    if (time > timestamp) {
        scale_particle_velocities(cells, number_of_particles, e_kin);
        timestamp += thermostat_interval;
    }
}

void Thermostat::scale_particle_velocities(std::vector<Cell>& cells,
                                           const int number_of_particles,
                                           const real e_kin)
{
	real beta;
	// compute current temperature
	current_temperature	= 2 * e_kin / (3 * number_of_particles);
	// compute scaling factor beta
    
    // beta = sqrt(1 + gamma * (target_temperature / current_temperature - 1));
    beta = sqrt(target_temperature / current_temperature);

    //std::cout << "\nbeta'>\t" << beta
    //            << "\nbeta>\t" << sqrt(target_temperature / current_temperature)
    //            << "\nt-now>\t" << current_temperature
    //            << "\nt-goal>\t" << target_temperature << "\n" << std::flush;

	// update the velocity of all particles by multiplying them with beta

	for (auto cell_it = cells.begin(); cell_it < cells.end(); cell_it++) {
		for (auto particle_it = cell_it->particles.begin(); particle_it < cell_it->particles.end(); particle_it++) {
			for (unsigned d = 0; d < DIM; d++) {
				(*particle_it).v[d] *= beta;
			}
		}
	}
}



void Thermostat::normal_sample(real v)
{
    real a, b, r;
    do {
        // we generate two real values in the interval [0,1]
        a = real(rand()) / RAND_MAX;
        b = real(rand()) / RAND_MAX;
        // we shift them so they lie in [-1,1)
        a = 2. * a - 1.;
        b = 2. * b - 1.;
        // we compute the squared radius
        r = sqr(a) + sqr(b);
    } while (r >= 1.0 || r == 0);

    sample = sqrt(v) * a * sqrt(-2 * log(r) / r);
}



void Thermostat::set_start_velocities(std::vector<Particle>& particles)
{
    if (target_temperature < 0.0)
        return;

    if (maxwell_boltzmann_seed < 1)
		srand(time(NULL));
	else 
		srand(maxwell_boltzmann_seed);

    real var;
    for (auto particle_it = particles.begin(); particle_it < particles.end(); particle_it++) {
        // we compute the variance, depending on the mass of the particle
        var = start_temperature / particle_it->m;

        // we initialize velocity components according to random normal distribution
        for (int d = 0; d < DIM; d++) {
            normal_sample(var);
            particle_it->v[d] = sample;
        }
    }
}
