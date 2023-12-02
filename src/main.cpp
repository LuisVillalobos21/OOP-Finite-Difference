#include "Grid.hpp"
#include "Connectivity.hpp"
#include "LaplacianOperator.hpp"
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

	//ConjugateGradient(, , 1e-3, [&]() { ; });

	return 0;
}