#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"

struct Grid;
struct Connectivity;

struct RHS
{
	std::vector<double> rhs_vector;

	RHS(Grid& grid);

	void forcingFunction(const Grid& grid, const Connectivity connect);
};

