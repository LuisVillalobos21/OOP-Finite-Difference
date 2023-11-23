#include "LHS.hpp"

LHS::LHS(Grid& grid)
{
	lhs_vector.resize(grid.num_points);
}

std::vector<double> LHS::assembleLHS(Grid& grid, LaplacianOperator& laplace)
{
	laplace.calculateOperator(grid, connect, const std::vector<int>&nodeIDs, std::vector<double>&function_values);
}