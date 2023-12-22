/** \mainpage LJ-particle-sim
 *
 * \section main_sec Overview
 *
 * Simulating particle interactions using the 
 * Lennard-Jones potential.
 *
 * Key components:
 *
 * - \b World: This class represents the entire simulation environment.
 *
 * - \b Observer: The Observer class is designed to monitor the simulation's 
 *   progress. Specifically, it observes the time discretization of the simulation, 
 *   tracking and recording relevant metrics at each step.
 *
 * - \b Potential: This class is responsible for calculating the 
 *   potential energy between particles.
 *
 * - \b Particle: The Particle class represents individual particles within the 
 *   simulation. Each Particle instance holds its own state, such as position, 
 *   velocity, and any other relevant attributes.
 *
 * - \b TimeDiscretization: This is a base class for different time discretization 
 *   strategies, like the Velocity Verlet algorithm.
 *
 * - \b Thermostat: The Thermostat class is responsible for regulating the temperature 
 *   of the simulation environment.
 *
 */
#ifndef DEFINES_HPP
#define DEFINES_HPP

// define the dimension of the particles
#define DIM 3
// reals in double precision
typedef double real;
// square define
#define sqr(_x) ((_x) * (_x))
// border types
enum BorderType { unknown = 0, leaving = 1, periodic = 2 };

#endif  // DEFINES_HPP
