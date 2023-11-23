#include "RHS.hpp"
#include "Grid.hpp"

RHS::RHS(Grid& grid)
{
	rhs_vector.resize(grid.num_points);
}