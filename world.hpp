#ifndef WORLD_HPP
#define WORLD_HPP

#include "cell.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include "thermostat.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <list>

/**
 * @brief Holds all information of the simulation environment.
 */
class World
{
public:
    World();

    /**
     * @brief convert given string to border type, return leaving if string is
     * unknown
     *
     * @param s
     *
     * @return border type
     */
    BorderType string_to_BorderType(std::string s);

    /**
     * @brief read the world parameters from the given parameter file
     *
     * parameter file example
     * \code
     * delta_t 0.1
     * t_end 1.0
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    virtual void read_Parameter(const std::string& filename);

    /**
     * @brief read the particles from the given data file
     *
     * @param filename filename of the particle data file
     */
    virtual void read_Particles(const std::string& filename);

    // data structures
    /// Name of the simulated world
    std::string name;
    /// Current time
    real t;
    /// Timestep
    real delta_t;
    /// epsilon
    real eps;
    /// sigma
    real sigma;
    // output interval (we average the energies over the last out_interval timesteps)
    int out_interval;
    /// End of simulation
    real t_end;
	/// number of timesteps in the simulation
	int timestep;
    /// kinetic energy
    real e_kin;
    /// potential energy
    real e_pot;
    /// total energy
    real e_tot;
	/// list of kinetic energies
	std::list<real> e_kin_list;
	/// list of potential energies
	std::list<real> e_pot_list;
	/// average kinetic energy
	real e_kin_avg;
    /// average potential energy
    real e_pot_avg;
    /// average total energy
    real e_tot_avg;
    /// Vector of particles
    std::vector<Particle> particles;
    /// Border conditions
    BorderType upper[DIM], lower[DIM];
    // domain
    real length_x[DIM];
    // number of particles
    int number_of_particles;

protected:
    virtual void print(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, World& W);
};

class World_LC : public World
{
public:
    /// Constructor of World_LC
    World_LC();
    void read_Parameter(const std::string& filename);
    void read_Particles(const std::string& filename);
    void initialize_Cells();


    /// Number of cells in every dimension
    int cell_N[DIM];
    /// length of cells
    real cell_length[DIM];
    /// Cut of radius used for calculation of the cell length
    real cell_r_cut;
    // Total number of cells in world
    int global_cell_number;
    /// Vector of cells
    std::vector<Cell> cells;
    /// 2d array of neighboring cells, will be of size 3^DIM - 1 rows and DIM coloumns
    //int** neighboring_indices; 
    std::vector<int> neighboring_indices;

    /// updates cells of moved particles
    void sort_particles_in_cells();

    std::vector<unsigned> cell_neighbors(const int linear_cell_index);

    int compute_global_cell_index(const int (&cell_index)[DIM]);

    Thermostat T;

    /// computes global cell index
    int linear_cell_index(Particle *p);

    void compute_cell_index(int linear_cell_index, int (&cell_index)[DIM]);


protected:
    virtual void print(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, World_LC& W);
};

/**
 * @brief a ostream operator for the World class
 *
 * @param os stream object
 * @param W the world
 *
 * @return resulting stream object
 */
std::ostream& operator<<(std::ostream& os, World& W);
std::ostream& operator<<(std::ostream& os, World_LC& W);

#endif  // WORLD_HPP
