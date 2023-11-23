#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "LaplacianOperator.hpp"

struct Grid;
struct LaplacianOperator;

struct LHS 
{
	std::vector<double> lhs_vector;
	LHS(Grid& grid);
	std::vector<double> assembleLHS(Grid& grid, LaplacianOperator& laplace);
};


