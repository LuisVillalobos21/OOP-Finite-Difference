#include "Grid.hpp"
#include "Connectivity.hpp"
#include "GenericOperator.hpp"
#include "ConjugateGradient.hpp"
#include "BoundaryCondition.hpp"
#include "AssmbleLHS.hpp"
#include "ArrayOperations.hpp"

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
		outFile << results[0][i] << ", " << results[1][i] << ", " << results[2][i] << ", " << results[3][i] << ", " << results[4][i] << "\n";
	}

	outFile.close();

	std::cout << "Results have been written to file: " << filename << '\n';
}

std::vector<std::vector<double>> collectResults(std::vector<double>& x, std::vector<double>& y, std::vector<double>& u, std::vector<double>& v, std::vector<double>& p)
{
	std::vector<std::vector<double>> results;
	results.resize(5);
	results[0] = x;
	results[1] = y;
	results[2] = u;
	results[3] = v;
	results[4] = p;

	return results;
}

int main()
{
	double dx = 0.0078125;
	double dy = 0.0078125;
	double start_x = 0.0;
	double end_x = 1.0;
	double start_y = 0.0;
	double end_y = 1.0;
	 
	double dt = 0.001; 
	double nu = .01; 

	double u_lid = 1;

	Grid grid(start_x, end_x, start_y, end_y, dx, dy, dt, nu);

	Connectivity connect(grid);

	//connect.printIDMatrix();

	// Creation of boundary condition struct vectors for u,v,pressure
	BoundaryCondition u_wall_boundary = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -1);
	BoundaryCondition u_move_wall_boundary = makeBoundaryCondition(BoundaryConditionType::Dirichlet, u_lid, -2); 
	BoundaryCondition u_wall_boundary3 = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -3);
	std::vector<BoundaryCondition> u_bc_struct_vector;
	u_bc_struct_vector.push_back(u_wall_boundary);
	u_bc_struct_vector.push_back(u_move_wall_boundary);
	u_bc_struct_vector.push_back(u_wall_boundary3);

	BoundaryCondition v_wall_boundary1 = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -1);
	BoundaryCondition v_wall_boundary2 = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -2);
	BoundaryCondition v_wall_boundary3 = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -3);
	std::vector<BoundaryCondition> v_bc_struct_vector;
	v_bc_struct_vector.push_back(v_wall_boundary1);
	v_bc_struct_vector.push_back(v_wall_boundary2);
	v_bc_struct_vector.push_back(v_wall_boundary3);

	BoundaryCondition p_wall_boundary_1 = makeBoundaryCondition(BoundaryConditionType::Neumann, 0.0, -1);
	BoundaryCondition p_wall_boundary_2 = makeBoundaryCondition(BoundaryConditionType::Neumann, 0.0, -2);
	BoundaryCondition p_wall_boundary_3 = makeBoundaryCondition(BoundaryConditionType::Dirichlet, 0.0, -3);
	std::vector<BoundaryCondition> p_bc_struct_vector;
	p_bc_struct_vector.push_back(p_wall_boundary_1);
	p_bc_struct_vector.push_back(p_wall_boundary_2);
	p_bc_struct_vector.push_back(p_wall_boundary_3);

	GenericOperator op(grid);
	GenericOperator op2(grid);
	GenericOperator op3(grid);

	std::vector<double> u_velocity;
	std::vector<double> v_velocity;
	std::vector<double> pressure;

	u_velocity.resize(grid.num_points, 0.0);
	v_velocity.resize(grid.num_points, 0.0);
	pressure.resize(grid.num_points, 0.0);

	AssembleLHS assemblyFunction_Laplace = laplaceLHS;
	AssembleLHS assemblyFunction_Roperator = R_LHS;

	std::vector<double> u_fractional(grid.num_points, 0.0);
	std::vector<double> v_fractional(grid.num_points, 0.0);
	std::vector<double> p_corrected(grid.num_points, 0.0);

	std::vector<double> u_old(grid.num_points, 0.0);
	std::vector<double> v_old(grid.num_points, 0.0);

	std::vector<double> RHS(grid.num_points, 0.0);
	std::vector<double> RHS2(grid.num_points, 0.0);

	double tolerance = 1e-4;
	int max_num_iter = 500'000;

	// START OF NAVIER STOKES SOLVE LOOP
	for (int k = 0; k < max_num_iter; ++k)
	{
		// U FRACTIONAL VELOCITY CALCULATION 

		op.laplace.calculateOperator(grid, connect, u_bc_struct_vector, u_velocity);
		op2.laplace.calculateOperator(grid, connect, u_bc_struct_vector, u_velocity);
		op.advec.calculateOperator(grid, connect, u_bc_struct_vector, u_velocity, v_velocity, u_velocity);
		op2.advec.calculateOperator(grid, connect, u_bc_struct_vector, u_velocity, v_velocity, u_old);

		for (int i = 1; i < grid.num_points; ++i)
		{
			RHS2[i] = u_velocity[i] + grid.dt * -(1.5 * op.advec.advection_vector[i] - 0.5 * op2.advec.advection_vector[i]) +
				grid.dt * grid.nu * (op.laplace.laplace_vector[i] - op.laplace.rhs_vector[i]);
			u_fractional[i] = RHS2[i];

			RHS[i] = u_velocity[i] + grid.dt * -(1.5 * op.advec.advection_vector[i] - 0.5 * op2.advec.advection_vector[i]) +
				0.5 * grid.dt * grid.nu * (op.laplace.laplace_vector[i] - 2 * op.laplace.rhs_vector[i]/* - op2.laplace.rhs_vector[i]*/);
		}

		//u_fractional = ConjugateGradient(
		//	assemblyFunction_Roperator,
		//	RHS,
		//	tolerance,
		//	grid,
		//	connect,
		//	op3,
		//	u_bc_struct_vector,
		//	u_velocity, 2
		//);

		// V FRACTIONAL VELOCITY CALCULATION 

		op.laplace.calculateOperator(grid, connect, v_bc_struct_vector, v_velocity);
		op.advec.calculateOperator(grid, connect, v_bc_struct_vector, u_velocity, v_velocity, v_velocity);
		op2.advec.calculateOperator(grid, connect, v_bc_struct_vector, u_velocity, v_velocity, v_old);

		for (int i = 1; i < grid.num_points; ++i)
		{
			RHS2[i] = v_velocity[i] + grid.dt * -(1.5 * op.advec.advection_vector[i] - 0.5 * op2.advec.advection_vector[i]) +
				grid.dt * grid.nu * (op.laplace.laplace_vector[i] - op.laplace.rhs_vector[i]);
			v_fractional[i] = RHS2[i];

			RHS[i] = v_velocity[i] + grid.dt * -(1.5 * op.advec.advection_vector[i] - 0.5 * op2.advec.advection_vector[i]) +
				0.5 * grid.dt * grid.nu * (op.laplace.laplace_vector[i] - 2 * op.laplace.rhs_vector[i]/* - op2.laplace.rhs_vector[i]*/);
		}

		//v_fractional = ConjugateGradient(
		//	assemblyFunction_Roperator,
		//	RHS,
		//	tolerance,
		//	grid,
		//	connect,
		//	op3,
		//	v_bc_struct_vector,
		//	v_velocity, 2
		//);

		// PRESSURE CALCULATION

		op.laplace.calculateOperator(grid, connect, p_bc_struct_vector, pressure);
		op.div.calculateOperator(grid, connect, u_bc_struct_vector, v_bc_struct_vector, u_fractional, v_fractional);
		op2.div.calculateOperator(grid, connect, u_bc_struct_vector, v_bc_struct_vector, u_fractional, v_fractional);

		for (int i = 0; i < grid.num_points; ++i)
		{
			op.laplace.rhs_vector[i] = (1 / grid.dt) * op.laplace.rhs_vector[i];  
			op.div.divergence_vector[i] = (1 / grid.dt) * op.div.divergence_vector[i];
			RHS[i] = op.div.divergence_vector[i] + op.laplace.rhs_vector[i];
		}

		//std::cout << "Solving Pressure!" << '\n';
		p_corrected = ConjugateGradient(
			assemblyFunction_Laplace,
			RHS,
			tolerance,
			grid,
			connect,
			op,
			p_bc_struct_vector,
			pressure, 2
		);

		// CORRECT VELOCITY CALCULATION

		u_old = u_velocity;
		v_old = v_velocity;

		op.grad2.calculateOperator(grid, connect, p_bc_struct_vector, p_corrected); // bc vector is dummy placeholder

		for (int i = 0; i < grid.num_points; ++i)
		{
			pressure[i] = p_corrected[i];
			u_velocity[i] = u_fractional[i] - grid.dt * op.grad2.gradient_vector_x[i];
			v_velocity[i] = v_fractional[i] - grid.dt * op.grad2.gradient_vector_y[i];
		}

		// CHECK FOR TIME INVARIANT SOLUTION

		double l2Norm = 0.0;
		double normUold = 0.0;

		for (size_t i = 0; i < u_velocity.size(); ++i) {
			l2Norm += (u_velocity[i] - u_old[i]) * (u_velocity[i] - u_old[i]);
			normUold += u_old[i] * u_old[i];
		}

		l2Norm = std::sqrt(l2Norm);
		normUold = std::sqrt(normUold);
		double relative_error_tolerance = 1e-5;

		double relativeError = normUold > 0 ? l2Norm / normUold : l2Norm;

		if (relativeError < relative_error_tolerance) {
			std::cout << "Relative error below tolerance: " << relativeError << std::endl;
			break;
		}
		else {
			std::cout << "Time step: "<< k + 1 << " Relative error between steps: " << relativeError << std::endl;
		}
	}

	std::vector<std::vector<double>> results = collectResults(grid.X, grid.Y, u_velocity, v_velocity, pressure);

	csvOut(results, "results.csv");

	return 0;
}