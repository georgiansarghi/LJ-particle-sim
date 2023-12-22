#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include "particle.hpp"
#include <iostream>
#include <math.h>



/**
 * @brief Abstract Potential class.
 */
class Potential {
public:
    /**
     * @brief Calculate the force between the two particles and add it to p.
     *
     * @param p particle p
     * @param q particl q
     *
     * @return potential energy
     */
    virtual real force(Particle& p, Particle& q, real *difference, real sigma, real eps) = 0;
};

#endif  // POTENTIAL_HPP
