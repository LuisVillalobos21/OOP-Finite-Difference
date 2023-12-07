#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "BoundaryCondition.hpp"

struct Grid;
struct Connectivity;
struct BoundaryCondition;

struct Divergence2Order
{
	    std::vector<double> divergence_vector;
	
	    Divergence2Order(Grid& grid);
	
	    std::vector<double> calculateOperator(Grid& grid,
	        Connectivity& connect,
	        std::vector<BoundaryCondition>& u_BC_vector,
			std::vector<BoundaryCondition>& v_BC_vector,
	        std::vector<double>& u_values,
	        std::vector<double>& v_values);
	
	    void calculateInnerPoints(Grid& grid,
	        Connectivity& connect,
	        std::vector<double>& u_values,
	        std::vector<double>& v_values);
	
	    void calculateBoundaryPoints(Grid& grid,
	        Connectivity& connect,
	        std::vector<double>& u_values,
	        std::vector<double>& v_values);
	
		void applyBoundaryCondition(Grid& grid,
			Connectivity& connect,
			const std::vector<int>& nodeIDs,
			std::vector<double>& u_values,
			std::vector<double>& v_values,
			BoundaryCondition& u_BC_vector,
			BoundaryCondition& v_BC_vector);
};
