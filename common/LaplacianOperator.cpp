#include "LaplacianOperator.hpp"

LaplacianOperator::LaplacianOperator(Grid& grid)
{
    laplace_vector.resize(grid.num_points, 0.0);
}

std::vector<double> LaplacianOperator::calculateOperator(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values)
{
    std::fill(laplace_vector.begin(), laplace_vector.end(), 0.0);
    calculateInnerPoints(grid, connect, function_values);
    calculateBoundaryPoints(grid, connect, connect.boundary_ID, function_values);

    return laplace_vector;
}

void LaplacianOperator::calculateInnerPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values)
{
    double invdx2 = 1.0 / (grid.dx * grid.dx); 
    double invdy2 = 1.0 / (grid.dy * grid.dy); 
    std::vector<double> invD = { invdx2, invdx2, invdy2, invdy2 };

    for (const auto& nodeID : connect.inner_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double laplace_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            laplace_value += invD[i] * function_values[neighborID];
            ++i;
        }

        laplace_value -= 2 * (invdx2 + invdy2) * function_values[nodeID];

        laplace_vector[nodeID] = laplace_value;
    }
}

void LaplacianOperator::calculateBoundaryPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values)
{
    double invdx2 = 1.0 / (grid.dx * grid.dx);
    double invdy2 = 1.0 / (grid.dy * grid.dy);
    std::vector<double> invD = { invdx2, invdx2, invdy2, invdy2 };

    for (const auto& nodeID : nodeIDs)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double laplace_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors) 
        {
            if (neighborID == -1) {
                ++i;
                continue; 
            }

            laplace_value += invD[i] * function_values[neighborID];
            ++i;
        }

        laplace_value -= 2 * (invdx2 + invdy2) * function_values[nodeID];

        laplace_vector[nodeID] = laplace_value;
    }
}
