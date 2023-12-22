#ifndef TIMEDISCRETIZATION_HPP
#define TIMEDISCRETIZATION_HPP

#include "observer.hpp"
#include "potential.hpp"
#include "world.hpp"
#include <iostream>

/**
 * @brief Base Class for TimeDiscretization Algorithms
 */
class TimeDiscretization
{
public:
    /**
     * @brief constructor
     *
     * @param W world configuration
     * @param Pot potential used for force calculation
     * @param O Observer of the simulation
     */
    TimeDiscretization(World& W, Potential& Pot, Observer& O);

    /**
     * @brief run a single timestep
     *
     * @param delta_t length of the timestep
     */
    virtual void timestep(real delta_t) = 0;

    /**
     * @brief run the simulation
     */
    virtual void simulate() = 0;

    /**
     * @brief calculares the forces affecting the particles at the current time
     */
    virtual void compute_Force() = 0;

    /**
     * @brief calculates the new velocity of the particles
     */
    virtual void update_V() = 0;

    /**
     * @brief calculate the new position of all particles according to their
     * velocity
     */
    virtual void update_X() = 0;

protected:
    // data structures
    /// the world where the particles live in
    World& W;
    /// the potential used for force calculation
    Potential& Pot;
    /// the observer of the simulation
    Observer& O;

private:
    TimeDiscretization();
};

#endif  // TIMEDISCRETIZATION_HPP
