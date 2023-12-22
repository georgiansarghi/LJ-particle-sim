#ifndef CELL_HPP
#define CELL_HPP

#include "particle.hpp"
#include <vector>

class Cell
{
public:
    Cell() : particles(0){};
    std::vector<Particle> particles;

	/**
	* @brief vector of global cell_indices of the neighboring cells,
	* these will depend on the BoundaryType
	*/
    std::vector<unsigned> neighbors;
    
	/**
	* @brief shows if a particle is in a ghost cell or not
	*/
    bool ghost_flag = false;

	/**
	* @brief true if the forces in the particles in the cell have to be updated
	*/
    bool force_flag = true;
};

#endif  // CELL_HPP
