#pragma once
#include <ProjectIncludes.hpp>

struct Grid;

struct Connectivity {

	std::vector<int> ID_vector;				// the IDs for all points in the domain
	std::vector<std::vector<int>> ID_matrix;

	std::vector<int> inner_ID;
	std::vector<int> left_ID;
	std::vector<int> right_ID;
	std::vector<int> bot_ID;
	std::vector<int> top_ID;
	std::vector<int> bot_left_ID;
	std::vector<int> bot_right_ID;
	std::vector<int> top_left_ID;
	std::vector<int> top_right_ID;

	std::vector<int> neighbor_ID;			// neighbor arrays 

	std::vector<int> neighbor_index;		// neighbor indexing arrays

	Connectivity(Grid& grid);

	void constructIDMatrix(int nx, int ny);

	void construct_Solution_PointIDs(Grid& grid);

	void constructNeighborIDs();

	void printIDData(std::vector<int> x);
};