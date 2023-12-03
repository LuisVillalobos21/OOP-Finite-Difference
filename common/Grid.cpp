#include "Grid.hpp"

Grid::Grid(double x1, double x2, double y1, double y2, double dx, double dy, double dt, double nu)
	: dx(dx), dy(dy), dt(dt), nu(nu)
{
	nx = static_cast<int>((x2 - x1) / dx) - 1;
	ny = static_cast<int>((y2 - y1) / dy) - 1;
	num_points = nx * ny;
	Eigen::VectorXd x_spacing = Eigen::VectorXd::LinSpaced(nx, x1 + dx, x2 - dx);
	Eigen::VectorXd y_spacing = Eigen::VectorXd::LinSpaced(ny, y1 + dy, y2 - dx);

	X.resize(num_points);
	Y.resize(num_points);

	createMeshgrid(x_spacing, y_spacing);
	createMeshdata();
}

void Grid::createMeshgrid(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
	X_grid.resize(y.size(), x.size());
	Y_grid.resize(y.size(), x.size());

	for (int i = 0; i < y.size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			X_grid(i, j) = x(j);
			Y_grid(i, j) = y(i);
		}
	}
}


void Grid::createMeshdata()
{
	int k = 0;

	for (int i = 0; i < Y_grid.rows(); ++i)
	{
		for (int j = 0; j < Y_grid.cols(); ++j)
		{
			X[k] = X_grid(i, j);
			Y[k] = Y_grid(i, j);
			k += 1;
		}
	}
}

void Grid::printMeshData(std::vector<double> x)
{
	std::cout << "Mesh values are: " << '\n';
	for (auto& value : x)
	{
		std::cout << value << '\n';
	}
}
