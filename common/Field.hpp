#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"

struct Field
{
	std::vector<double> vector;

	Field(Grid& grid);
};


