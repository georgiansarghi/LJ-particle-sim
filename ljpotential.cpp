#include "ljpotential.hpp"

real LJPotential::force(Particle& p, Particle& q, real *difference, real sigma, real eps) {
    real U, sigma_term, dist_sqr = 0.0;

    for (unsigned i = 0; i < DIM; i++) {
        dist_sqr += sqr(difference[i]);
    }


    sigma_term = sqr(sqr(sigma) * sigma) / (sqr(dist_sqr) * dist_sqr);

    for (unsigned i = 0; i < DIM; i++) {
        p.F[i] += 24 * eps * sigma_term * ( 1 - 2 * sigma_term ) * difference[i] / dist_sqr;

        // we use simmetry of the force
        q.F[i] -= 24 * eps * sigma_term * ( 1 - 2 * sigma_term ) * difference[i] / dist_sqr;
    }

    U = 4 * eps * sigma_term * (sigma_term - 1);


    return U;
}
