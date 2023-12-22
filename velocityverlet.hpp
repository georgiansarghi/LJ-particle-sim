#ifndef VELOCITYVERLET_HPP
#define VELOCITYVERLET_HPP

#include "subdomain.hpp"
#include "timediscretization.hpp"

/**
 * @brief Implementation of the Velocity Verlet Algorithm
 */
class VelocityVerlet : public TimeDiscretization
{
public:
    /**
     * @brief constructor
     *
     * @param W world configuration
     * @param Pot potential used for force calculation
     * @param O Observer of the simulation
     */
    VelocityVerlet(World& W, Potential& Pot, Observer& O);

    /**
     * @brief constructor
     *
     * This is an example for Constructor overloading. If you have read until
     * here you can use the other constructor and change the blatt1 main
     * function.
     *
     * @param W world configuration
     * @param Pot potential used for force calculation
     * @param O Observer of the simulation
     */
    VelocityVerlet(World& W, Potential* Pot, Observer& O);

    /**
     * @brief run a single timestep
     *
     * @param delta_t length of the timestep
     */
    virtual void timestep(real delta_t);

    /**
     * @brief run the simulation
     */
    virtual void simulate();

    /**
     * @brief calculares the forces affecting the particles at the current time
     */
    virtual void compute_Force();

    /**
     * @brief calculates the new velocity of the particles
     */
    virtual void update_V();

    /**
     * @brief calculate the new position of all particles according to their
     * velocity
     */
    virtual void update_X();

    /**
     * @brief check for each particle if has left the box
     */
    virtual void handle_borders();

    // data structures inherited from TimeDiscretization

    virtual void compute_diff(Particle *p, Particle *q, real (&diff)[DIM]);
    virtual real compute_distance(Particle *p, Particle *q);
protected:



private:
    VelocityVerlet();
};

/**
 * @brief Implementation of the Velocity Verlet Algorithm
 */
class VelocityVerlet_LC : public VelocityVerlet
{
public:
    /**
     * @brief constructor
     *
     * @param W world configuration
     * @param Pot potential used for force calculation
     * @param O Observer of the simulation
     */
    VelocityVerlet_LC(World_LC& W, Potential& Pot, Observer& O);

    /**
     * @brief constructor
     *
     * This is an example for Constructor overloading. If you have read until
     * here you can use the other constructor and change the blatt1 main
     * function.
     *
     * @param W world configuration
     * @param Pot potential used for force calculation
     * @param O Observer of the simulation
     */
    VelocityVerlet_LC(World_LC& W, Potential* Pot, Observer& O);

    /**
     * @brief run a single timestep
     *
     * @param delta_t length of the timestep
     */
    virtual void timestep(real delta_t);

    /**
     * @brief run the simulation
     */
    virtual void simulate();

    /**
     * @brief calculates the forces affecting the particles at the current time
     */
    virtual void compute_Force(const unsigned int cell_index);

    /**
     * @brief calculates the new velocity of the particles
     */
    virtual void update_V(std::vector<Particle>& particles);

    /**
     * @brief calculate the new position of all particles according to their
     * velocity
     */
    virtual void update_X(const unsigned int cell_index);
	
	//virtual void compute_diff(Particle *p, Particle *q, real (&diff)[DIM]);

    // data structures inherited from TimeDiscretization
protected:
    /**
     *  @brief computation of forces for neighboring cells
     */
    void compute_Force_from_neighboring_cell(
        const unsigned int cell_index, const int other_cell_index);
    /**
     *  @brief computation of forces for same cell
     */
    void compute_Force_from_same_cell(const unsigned int cell_index);

    // data structures inherited from TimeDiscretization

    /// Overwriting world in timediscretization with linked_cell world
    World_LC& W;

private:
    VelocityVerlet_LC();
};

class VelocityVerlet_Sub : public VelocityVerlet_LC
{
public:
    VelocityVerlet_Sub(SubDomain& W, Potential& Pot, Observer& O);
    VelocityVerlet_Sub(SubDomain& W, Potential* Pot, Observer& O);

    virtual void timestep(real delta_t);
    virtual void simulate();

    virtual void compute_Force(const unsigned int cell_index);
    
    void compute_Force_from_neighboring_cell(
        const unsigned int cell_index, const int other_cell_index);

    void compute_Force_from_same_cell(const unsigned int cell_index);

    virtual void update_X(const unsigned int current_cell_index);

    virtual void update_V(std::vector<Particle>& particles);

    /// Clearing border cells from particles
    void clear_ghosts();
    /// Sending particles to boundaries to serve as ghosts.
    void create_ghosts(int d);

    void send_inner_receive_outer();
    void send_outer_receive_inner();

	// set active flag of particles which leave the domain to false
	// only relevant for leaving boundaries
    void handle_borders();

	
    real compute_distance(Particle *p, Particle *q);
    // calculates (vectorial) difference of two particles and stores it in the array diff
    void compute_diff(Particle *p, Particle *q, real (&diff)[DIM]);
    
protected:
    SubDomain& W;

private:
    VelocityVerlet_Sub();
};

#endif  // VELOCITYVERLET_HPP
