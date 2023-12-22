#ifndef LJPOTENTIAL_HPP
#define LJPOTENTIAL_HPP

#include "potential.hpp"

class LJPotential : public Potential
{
public:
    /**
     * @brief Calculate the force between the two particles and add it to p.
     *
     * @param p particle p
     * @param q particle q
	 * @param sigma root parameter of LJ potential
	 * @param eps depth parameter of LJ potential
     *
     * @return potential energy
     */
    virtual real force(Particle& p, Particle& q, real *difference, real sigma, real eps);
};

#endif  // LJPOTENTIAL_HPP
