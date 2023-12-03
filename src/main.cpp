#include "Grid.hpp"
#include "Connectivity.hpp"
#include "GenericOperator.hpp"
#include "LaplacianOperator.hpp"
#include "ConjugateGradient.hpp"
#include "BoundaryCondition.hpp"

std::vector<double> rhsForcing(const Grid& grid, const Connectivity& connect)
{
	std::vector<double> rhs(grid.num_points, 0.0);

	for (auto& currentID : connect.ID_vector)
	{
		double x = grid.X[currentID];
		double y = grid.Y[currentID];
		double value = -std::exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / (2 * 0.1 * 0.1));
		rhs[currentID] = value;
	}

	return rhs;
}

void csvOut(const std::vector<std::vector<double>>& results, const std::string& filename)
{
	std::ofstream outFile(filename);

	if (!outFile.is_open())
	{
		std::cerr << "Error: Could not open file for writing." << std::endl;
		return;
	}

	for (size_t i = 0; i < results[0].size(); ++i) 
	{
		outFile << results[0][i] << ", " << results[1][i] << ", " << results[2][i] << "\n";
	}

	outFile.close();
}

std::vector<std::vector<double>> collectResults(std::vector<double>& x, std::vector<double>& y, std::vector<double>& data)
{
	std::vector<std::vector<double>> results;
	results.resize(3);
	results[0] = x;
	results[1] = y;
	results[2] = data;

	return results;
}

int main()
{
	double tolerance = 1e-4;

	double dx = 0.005;
	double dy = 0.005;
	double start_x = 0.0;
	double end_x = 1.0;
	double start_y = 0.0;
	double end_y = 1.0;

	Grid grid(start_x, end_x, start_y, end_y, dx, dy);

	std::vector<double> u_velocity;
	std::vector<double> v_velocity;
	std::vector<double> pressure;

	u_velocity.resize(grid.num_points, 0.0);
	v_velocity.resize(grid.num_points, 0.0);
	pressure.resize(grid.num_points, 0.0);

	Connectivity connect(grid);

	std::vector<std::vector<int>> boundary_nodes;
	boundary_nodes.resize(1);
	boundary_nodes[0] = connect.boundary_ID;

	LaplacianOperator laplace(grid);

	std::vector<double> RHS = rhsForcing(grid, connect);

	// Create a Dirichlet boundary condition
	double dirichletValue = -5.0; 
	int dirichletFlag = -2; 
	DirichletCondition dirichletCondition(dirichletValue, dirichletFlag);
	BoundaryCondition boundaryCondition(&dirichletCondition, dirichletFlag);

	// Create a Neumann boundary condition
	//double neumannValue = 3.0; 
	//int neumannFlag = 2.0; 
	//NeumannCondition neumannCondition(neumannValue, neumannFlag);
	//BoundaryCondition boundaryCondition(&neumannCondition, neumannFlag);

	laplace.applyBoundaryCondition(grid, connect, connect.boundary_ID, pressure, boundaryCondition);

	for (int i = 0; i < RHS.size(); ++i)
	{
		RHS[i] += laplace.rhs_vector[i];
	}

	std::vector<double> solution = ConjugateGradient(
		laplace,
		RHS,
		tolerance,
		grid,
		connect,
		boundary_nodes,
		pressure
	);

	std::vector<std::vector<double>> results = collectResults(grid.X, grid.Y, solution);

	csvOut(results, "results.csv");

	return 0;
}