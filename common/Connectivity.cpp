#include <Connectivity.hpp>
#include "Grid.hpp"

Connectivity::Connectivity(Grid& grid)
{
	ID_vector.resize(grid.nx * grid.ny);
	std::iota(ID_vector.begin(), ID_vector.end(), 0);

	constructIDMatrix(grid.nx, grid.ny);

	construct_Solution_PointIDs(grid);
	constructNeighborIDs();
}


void Connectivity::constructIDMatrix(int nx, int ny)
{
	ID_matrix.resize(ny);
	for (auto& row : ID_matrix)
	{
		row.resize(nx);
	}

	int k = 0;
	for (int i = 0; i < ny; ++i)
	{
		for (int j = 0; j < nx; ++j)
		{
			int shuffled_index = ID_vector[k];
			ID_matrix[i][j] = shuffled_index;
			k += 1;
		}
	}
}


void Connectivity::construct_Solution_PointIDs(Grid& grid)
{
	int i;
	int j;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		 // all solution points (excluding 
	{													 // the ACTUAl boundary condition points)
		for (int j = 1; j < ID_matrix[i].size() - 1; ++j)
		{
			soln_ID.push_back(ID_matrix[i][j]);
		}
	}

	for (int i = 2; i < ID_matrix.size() - 2; ++i)		// inner solution points 
	{
		for (int j = 2; j < ID_matrix[i].size() - 2; ++j)
		{
			inner_ID.push_back(ID_matrix[i][j]);
		}
	}

	j = 1;
	for (int i = 2; i < ID_matrix.size() - 2; ++i)		// left solution points 
	{
		left_ID.push_back(ID_matrix[i][j]);
	}

	j = ID_matrix[1].size() - 2;
	for (int i = 2; i < ID_matrix.size() - 2; ++i)		// right solution points 
	{
		right_ID.push_back(ID_matrix[i][j]);
	}

	i = 1;
	for (int j = 2; j < ID_matrix[i].size() - 2; ++j)	// bot solution points 
	{
		bot_ID.push_back(ID_matrix[i][j]);
	}

	i = ID_matrix.size() - 2;
	for (int j = 2; j < ID_matrix[i].size() - 2; ++j)	// top solution points
	{
		top_ID.push_back(ID_matrix[i][j]);
	}

	i = 1;
	j = 1;
	bot_left_ID.push_back(ID_matrix[i][j]);				// bot left corner

	i = 1;
	j = ID_matrix[1].size() - 2;
	bot_right_ID.push_back(ID_matrix[i][j]);			// bot right corner

	i = ID_matrix.size() - 2;
	j = 1;
	top_left_ID.push_back(ID_matrix[i][j]);				// top left corner

	i = ID_matrix.size() - 2;
	j = ID_matrix[1].size() - 2;
	top_right_ID.push_back(ID_matrix[i][j]);			// top right corner
}

void Connectivity::constructNeighborIDs()
{
	int i;
	int j;

	int k = 0;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// all solution points, edges will have neighbors which point to 
	{													// fixed boundary bounds
		for (int j = 1; j < ID_matrix[i].size() - 1; ++j)
		{
			neighbor_index.push_back(k);				// this k is the index used to access the neighbors

			neighbor_ID.push_back(ID_matrix[i][j - 1]);	// west neighbor
			k += 1;

			neighbor_ID.push_back(ID_matrix[i][j + 1]);	// east neighbor
			k += 1;

			neighbor_ID.push_back(ID_matrix[i - 1][j]);; // south neightbor
			k += 1;

			neighbor_ID.push_back(ID_matrix[i + 1][j]);	// north neighbor
			k += 1;
		}
	}
}
