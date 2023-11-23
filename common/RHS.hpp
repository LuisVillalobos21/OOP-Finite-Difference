#pragma once
#include "ProjectIncludes.hpp"

struct Grid;

struct RHS
{
	std::vector<double> rhs_vector;

	RHS(Grid& grid);
};

