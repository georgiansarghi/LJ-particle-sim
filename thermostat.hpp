#ifndef THERMOSTAT_HPP
#define THERMOSTAT_HPP

#include "cell.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>

/**
 * @brief Implementation of thermostat
 */

class Thermostat
{
public:
    Thermostat();
    ~Thermostat();

    void execute_thermostat(real time, std::vector<Cell>& cells,
                            const int number_of_particles, const real e_kin);
    void scale_particle_velocities(std::vector<Cell>& cells,
                                   const int number_of_particles,
                                   const real e_kin);
    void normal_sample(real v);

    void set_start_velocities(std::vector<Particle>& particles);

    /// interval of time between scaling of velocities
    real thermostat_interval;
    /// timestamp of last thermostat execution
    real timestamp;
    /// needed to compute beta
    real gamma;
    /// current temperature of the system
    real current_temperature;
    /// desired temperature to which velocities are scaled.
    real target_temperature;
    /// temperature to which the ensemble is set in the beginning
    real start_temperature;
    /// seed for random number generator to control maxwell-boltzmann
    /// distribution
    int maxwell_boltzmann_seed;
    /// value of sample from normal distribution
    /// it is assigned every time we call the normal_sample function
    real sample;
};

#endif  // THERMOSTAT_HPP
