#include "AssmbleLHS.hpp"

std::vector<double> laplaceLHS(
	Grid& grid,
	Connectivity& connect,
	std::vector<BoundaryCondition>& BC_vector,
	GenericOperator& op,
	std::vector<double>& x)
{
	return op.laplace.calculateOperator(grid, connect, BC_vector, x);
}

std::vector<double> R_LHS(
	Grid& grid,
	Connectivity& connect,
	std::vector<BoundaryCondition>& BC_vector,
	GenericOperator& op,
	std::vector<double> x)
{
	std::vector<double> result(grid.num_points, 0.0);
	std::vector<double> laplace = op.laplace.calculateOperator(grid, connect, BC_vector, x);

	for (int i = 0; i < grid.num_points; ++i)
	{
		laplace[i] = 0.5 * grid.dt * grid.nu * laplace[i];
		result[i] = x[i] - laplace[i];
	}

	return result;
}