#ifndef OBSERVER_HPP
#define OBSERVER_HPP

#include "subdomain.hpp"
#include "world.hpp"
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <vector>



/**
 * @brief Observer for the time discretization algorithm.
 */
class Observer {
public:
    /**
     * @brief constructor
     *
     * opens and creates the files written during observation
     *
     * @param W
     */
    Observer(World& W);

    /**
     * @brief destructor
     *
     * closes the files written during the obervation
     */
    ~Observer();

    /**
     * @brief notify the observer that the world has changed
     */
    virtual void notify();

    /**
     * @brief output the particel configuration in the xyz format
     */
    virtual void output_xyz();

    /**
     * @brief output statistics like kinetic, potential and total energy
     */
    virtual void output_statistics();

    virtual void open_files();
    /**
     * @brief output the coordinates of the particles
     */
    virtual void output_coordinates();

protected:
    /// The world we are observing
    World& W;
    /// Statistics filestream
    std::ofstream statistics;
    /// coordiantes filestream
    std::ofstream coordinates;
    /// xyz filestream
    std::ofstream xyz;

private:
    /// Disabled Constructor
    Observer();
};



class Observer_LC : public Observer {
public:
    Observer_LC(World_LC& W) : Observer(W), W(W){};

	/**
     * @brief output the particel configuration in the xyz format
     */
    virtual void output_xyz();
	
	/**
     * @brief output statistics like kinetic, potential and total energy
     */
    virtual void output_statistics();
	
	/**
     * @brief output the coordinates of the particles
     */
    virtual void output_coordinates();

protected:
    World_LC& W;
};

class Observer_Sub : public Observer_LC
{
public:
    Observer_Sub(SubDomain& W);
	
	/**
     * @brief counts time steps to control when we output data
     */
    int out = 0;

    virtual void open_files();
	
	/**
     * @brief output statistics like kinetic, potential and total energy
     */
    virtual void output_statistics();
	
	/**
     * @brief output the particel configuration in the xyz format
     */
    virtual void output_xyz();
	
    virtual void notify();

protected:
    SubDomain& W;
};

#endif  // OBSERVER_HPP
