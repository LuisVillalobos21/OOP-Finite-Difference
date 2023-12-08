#pragma once
#include <ProjectIncludes.hpp>

struct Grid;

struct Connectivity {

	std::vector<int> ID_vector;				// the IDs for all points in the domain
	std::vector<std::vector<int>> ID_matrix;

	std::vector<int> inner_ID;
	std::vector<int> boundary_ID;
	std::vector<int> left_ID;
	std::vector<int> right_ID;
	std::vector<int> bot_ID;
	std::vector<int> top_ID;
	std::vector<int> bot_left_ID;
	std::vector<int> bot_right_ID;
	std::vector<int> top_left_ID;
	std::vector<int> top_right_ID;

	std::vector<std::vector<int>> neighbor_IDs;			// neighbor arrays 

	Connectivity(Grid& grid);

	void constructIDMatrix(int nx, int ny);

	void construct_Solution_PointIDs(Grid& grid);

	void constructNeighborIDs(Grid& grid);

	void Connectivity::printIDMatrix() const;

	void printIDData(std::vector<int> x);

	void Connectivity::printNeighborIDs();
};