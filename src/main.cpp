#include <Grid.hpp>
#include <Connectivity.hpp>

int main()
{
	double dx = 0.25;
	double dy = 0.25;
	double start_x = 0.0;
	double end_x = 1.0;
	double start_y = 0.0;
	double end_y = 1.0;

	Grid grid(start_x, end_x, start_y, end_y, dx, dy);

	//grid.printMeshData(grid.X);

	//grid.printMeshData(grid.Y);

	Connectivity connect(grid);

	connect.printIDData(connect.top_right_ID);

	connect.printNeighborIDs();

	return 0;
}