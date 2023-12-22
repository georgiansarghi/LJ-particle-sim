#include "particle.hpp"
#include <sstream>
#include <stdexcept>

void Particle::print(std::ostream& os) const
{
    //os << "Particle: " << id << std::endl;
    //os << "Mass: " << m << std::endl;
    os << " Postion: ";
    for (int d = 0; d < DIM; d++)
        os << x[d] << "\t";
    os << std::endl;

    os << "Velocity: ";
    for (int d = 0; d < DIM; d++)
        os << v[d] << "\t";
    os << std::endl;

    os << "Force: ";
    for (int d = 0; d < DIM; d++)
        os << F[d] << "\t";
    os << std::endl;

    os << "Force_old: ";
    for (int d = 0; d < DIM; d++)
        os << F_old[d] << "\t";
    os << std::endl;

    os << " Id: ";
    os << id << "\t";
    os << std::endl;
}

std::ostream& operator<<(std::ostream& os, Particle& p)
{
    p.print(os);
    return os;
}
