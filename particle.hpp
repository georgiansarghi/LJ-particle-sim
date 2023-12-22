#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "defines.hpp"
#include <fstream>
#include <iostream>

/**
 * @brief Particle data
 *
 * Contains the particle data.
 */
class Particle
{
public:
    int id;
    /// Mass
    real m;
    /// Position
    real x[DIM];
    /// Velocity
    real v[DIM];
    /// Force
    real F[DIM];
    /// Force (previous step)
    real F_old[DIM];

    // only used for BorderType leaving to indicate
    // inactive (out of simulation domain) particles
	bool active;

protected:
    void print(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, Particle& p);
};

#endif  // PARTICLE_HPP
