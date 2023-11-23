#include "Grid.hpp"
#include "Connectivity.hpp"
#include "LaplacianOperator.hpp"
#include "RHS.hpp"
#include "LHS.hpp"
#include "Field.hpp"
#include "ConjugateGradient.hpp"

int main()
{
	double dx = 0.25;
	double dy = 0.25;
	double start_x = 0.0;
	double end_x = 1.0;
	double start_y = 0.0;
	double end_y = 1.0;

	Grid grid(start_x, end_x, start_y, end_y, dx, dy);

	Connectivity connect(grid);

	LaplacianOperator laplace(grid);

	RHS rhs(grid);

	LHS lhs(grid);

	Field u_velocity(grid);

	lhs.assembleLHS(grid, connect, u_velocity, laplace);

	rhs.forcingFunction(grid, connect);

	return 0;
}