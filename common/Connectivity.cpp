#include <Connectivity.hpp>
#include "Grid.hpp"

Connectivity::Connectivity(Grid& grid)
{
	ID_vector.resize(grid.nx * grid.ny);
	std::iota(ID_vector.begin(), ID_vector.end(), 0);

	constructIDMatrix(grid.nx, grid.ny);

	construct_Solution_PointIDs(grid);
	constructNeighborIDs(grid);
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

	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// inner solution points 
	{
		for (int j = 1; j < ID_matrix[i].size() - 1; ++j)
		{
			inner_ID.push_back(ID_matrix[i][j]);
		}
	}

	j = 0;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// left solution points 
	{
		left_ID.push_back(ID_matrix[i][j]);
	}

	j = ID_matrix[1].size() - 1;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// right solution points 
	{
		right_ID.push_back(ID_matrix[i][j]);
	}

	i = 0;
	for (int j = 1; j < ID_matrix[i].size() - 1; ++j)	// bot solution points 
	{
		bot_ID.push_back(ID_matrix[i][j]);
	}

	i = ID_matrix.size() - 1;
	for (int j = 1; j < ID_matrix[i].size() - 1; ++j)	// top solution points
	{
		top_ID.push_back(ID_matrix[i][j]);
	}

	i = 0;
	j = 0;
	{
		bot_left_ID.push_back(ID_matrix[i][j]);				// bot left corner
	}

	i = 0;
	j = ID_matrix[1].size() - 1;
	{
		bot_right_ID.push_back(ID_matrix[i][j]);			// bot right corner
	}

	i = ID_matrix.size() - 1;
	j = 0;
	{
		top_left_ID.push_back(ID_matrix[i][j]);				// top left corner
	}

	i = ID_matrix.size() - 1;
	j = ID_matrix[1].size() - 1;
	{
		top_right_ID.push_back(ID_matrix[i][j]);			// top right corner
	}
}

void Connectivity::constructNeighborIDs(Grid& grid)
{
	int i;
	int j;
	int k = 0;
	neighbor_IDs.resize(grid.num_points);

	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// inner solution points 
	{
		for (int j = 1; j < ID_matrix[i].size() - 1; ++j)
		{
			std::vector<int> neighbors;
			neighbors.push_back(ID_matrix[i][j - 1]);
			neighbors.push_back(ID_matrix[i][j + 1]);
			neighbors.push_back(ID_matrix[i - 1][j]);
			neighbors.push_back(ID_matrix[i + 1][j]);

			int currentNode = ID_matrix[i][j];

			neighbor_IDs[currentNode] = neighbors;
		}
	}

	j = 0;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// left solution points 
	{
		std::vector<int> neighbors;
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i][j + 1]);
		neighbors.push_back(ID_matrix[i - 1][j]);
		neighbors.push_back(ID_matrix[i + 1][j]);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	j = ID_matrix[1].size() - 1;
	for (int i = 1; i < ID_matrix.size() - 1; ++i)		// right solution points 
	{
		std::vector<int> neighbors;
		neighbors.push_back(ID_matrix[i][j - 1]);
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i - 1][j]);
		neighbors.push_back(ID_matrix[i + 1][j]);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	i = 0;
	for (int j = 1; j < ID_matrix[i].size() - 1; ++j)	// bot solution points 
	{
		std::vector<int> neighbors;
		neighbors.push_back(ID_matrix[i][j - 1]);
		neighbors.push_back(ID_matrix[i][j + 1]);
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i + 1][j]);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	i = ID_matrix.size() - 1;
	for (int j = 1; j < ID_matrix[i].size() - 1; ++j)	// top solution points
	{
		std::vector<int> neighbors;
		neighbors.push_back(ID_matrix[i][j - 1]);
		neighbors.push_back(ID_matrix[i][j + 1]);
		neighbors.push_back(ID_matrix[i - 1][j]);
		neighbors.push_back(-1);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	i = 0;
	j = 0;
	{
		std::vector<int> neighbors;							// bot left corner
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i][j + 1]);
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i + 1][j]);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;			
	}

	i = 0;
	j = ID_matrix[1].size() - 1;
	{
		std::vector<int> neighbors;							// bot right corner
		neighbors.push_back(ID_matrix[i][j - 1]);
		neighbors.push_back(-1);
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i + 1][j]);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	i = ID_matrix.size() - 1;
	j = 0;
	{
		std::vector<int> neighbors;							// top left corner
		neighbors.push_back(-1);	
		neighbors.push_back(ID_matrix[i][j + 1]);
		neighbors.push_back(ID_matrix[i - 1][j]);
		neighbors.push_back(-1);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}

	i = ID_matrix.size() - 1;
	j = ID_matrix[1].size() - 1;
	{
		std::vector<int> neighbors;							// top right corner
		neighbors.push_back(ID_matrix[i][j - 1]);
		neighbors.push_back(-1);
		neighbors.push_back(ID_matrix[i - 1][j]);
		neighbors.push_back(-1);

		int currentNode = ID_matrix[i][j];

		neighbor_IDs[currentNode] = neighbors;
	}
}

void Connectivity::printIDData(std::vector<int> x)
{
	std::cout << "ID values are: " << '\n';
	for (auto& value : x)
	{
		std::cout << value << '\n';
	}
}

void Connectivity::printNeighborIDs()
{
	std::cout << "Neighbor IDs are: ";
	for (auto& ID : neighbor_IDs)
	{
		for (auto& neighbor : ID)
		{
			std::cout << neighbor << ", ";
		}
		std::cout << '\n';
	}
}
