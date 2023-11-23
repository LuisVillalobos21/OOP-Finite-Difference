#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "LaplacianOperator.hpp"
#include "Field.hpp"

struct Grid;
class LaplacianOperator;
struct Field;

struct LHS 
{
	std::vector<double> lhs_vector;
	LHS(Grid& grid);
	void assembleLHS(Grid& grid, Connectivity& connect, Field& field, LaplacianOperator& laplace);
};


