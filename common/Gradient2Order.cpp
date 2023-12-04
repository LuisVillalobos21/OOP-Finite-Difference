#include "Gradient2Order.hpp"

Gradient2Order::Gradient2Order(Grid& grid)
{
    gradient_vector_x.resize(grid.num_points, 0.0);
    gradient_vector_y.resize(grid.num_points, 0.0);
}

std::pair<std::vector<double>, std::vector<double>> Gradient2Order::calculateOperator(Grid& grid, Connectivity& connect, std::vector<BoundaryCondition>& BC_vector, std::vector<double>& function_values)
{
    std::fill(gradient_vector_x.begin(), gradient_vector_x.end(), 0.0);
    std::fill(gradient_vector_y.begin(), gradient_vector_y.end(), 0.0);
    calculateInnerPoints(grid, connect, function_values);
    calculateBoundaryPoints(grid, connect, function_values);

    for (int i = 0; i < BC_vector.size(); ++i)
    {
        applyBoundaryCondition(grid, connect, connect.boundary_ID, function_values, BC_vector[i]);
    }
    return {gradient_vector_x, gradient_vector_y};
}

void Gradient2Order::calculateInnerPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.inner_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double grad_value_x = 0.0;
        double grad_value_y = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            grad_value_x += invD_x[i] * function_values[neighborID];
            grad_value_y += invD_y[i] * function_values[neighborID];

            ++i;
        }

        gradient_vector_x[nodeID] = grad_value_x;
        gradient_vector_y[nodeID] = grad_value_y;
    }
}

void Gradient2Order::calculateBoundaryPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.boundary_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double grad_value_x = 0.0;
        double grad_value_y = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID < 0) {
                ++i;
                continue;
            }

            grad_value_x += invD_x[i] * function_values[neighborID];
            grad_value_y += invD_y[i] * function_values[neighborID];
            ++i;
        }

        gradient_vector_x[nodeID] = grad_value_x;
        gradient_vector_y[nodeID] = grad_value_y;
    }
}

void Gradient2Order::applyBoundaryCondition(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values, BoundaryCondition& BC)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };
    std::vector<int> opposite = {1, 0, 3, 2};

    for (const auto& nodeID : nodeIDs)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double rhs_value, operator_value;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID == BC.boundary_flag)
            {
                //std::tie(rhs_value, operator_value) = BC.returnBCValue(grid, nodeID, function_values, i);

                gradient_vector_x[nodeID] += invD_x[i] * ((2 * function_values[nodeID]) - function_values[neighbors[opposite[i]]]);
                gradient_vector_y[nodeID] += invD_y[i] * ((2 * function_values[nodeID]) - function_values[neighbors[opposite[i]]]);
                ++i;
                continue;
            }
            ++i;
        }
    }
}