#include "RHS.hpp"
#include "Grid.hpp"

RHS::RHS(Grid& grid)
{
	rhs_vector.resize(grid.num_points);
}

void RHS::forcingFunction(const Grid& grid, const Connectivity connect)
{
	for (auto& index : connect.ID_vector)
	{
		double x = grid.X[index];
		double y = grid.Y[index];

		double value = 0.02 * std::exp(-((x - 0.5) * (x - 0.5) / 0.09 + (y - 0.5) * (y - 0.5) / 0.25));
		rhs_vector[index] = value;
	}
}