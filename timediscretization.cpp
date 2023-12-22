#include "timediscretization.hpp"

TimeDiscretization::TimeDiscretization(World& W, Potential& Pot, Observer& O)
    : W(W), Pot(Pot), O(O)
{
}
