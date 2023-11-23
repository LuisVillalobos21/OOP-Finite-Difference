#include "LHS.hpp"

LHS::LHS(Grid& grid)
{
	lhs_vector.resize(grid.num_points);
}

void LHS::assembleLHS(Grid& grid, Connectivity& connect, Field& field, LaplacianOperator& laplace)
{
	std::fill(lhs_vector.begin(), lhs_vector.end(), 0.0);

	laplace.calculateOperator(grid, connect, field.vector);

	lhs_vector = laplace.laplace_vector;
}