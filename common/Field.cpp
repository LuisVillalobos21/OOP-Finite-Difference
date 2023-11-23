#include "Field.hpp"

Field::Field(Grid& grid)
{
	vector.resize(grid.num_points, 0.0);
}