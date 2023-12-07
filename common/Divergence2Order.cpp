#include "Divergence2Order.hpp"

Divergence2Order::Divergence2Order(Grid& grid)
{
	divergence_vector.resize(grid.num_points, 0.0);
}

std::vector<double> Divergence2Order::calculateOperator(Grid& grid,
    Connectivity& connect,
    std::vector<BoundaryCondition>& u_BC_vector,
    std::vector<BoundaryCondition>& v_BC_vector,
    std::vector<double>& u_values,
    std::vector<double>& v_values)
{
    std::fill(divergence_vector.begin(), divergence_vector.end(), 0.0);
    calculateInnerPoints(grid, connect, u_values, v_values);
    calculateBoundaryPoints(grid, connect, u_values, v_values);

    for (int i = 0; i < u_BC_vector.size(); ++i)
    {
        applyBoundaryCondition(grid, connect, connect.boundary_ID, u_values, v_values, u_BC_vector[i], v_BC_vector[i]);
    }
    return divergence_vector;
}

void Divergence2Order::calculateInnerPoints(Grid& grid,
    Connectivity& connect,
    std::vector<double>& u_values,
    std::vector<double>& v_values)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.inner_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double div_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            div_value += invD_x[i] * u_values[neighborID];
            div_value += invD_y[i] * v_values[neighborID];

            ++i;
        }

        divergence_vector[nodeID] = div_value;
    }
}

void Divergence2Order::calculateBoundaryPoints(Grid& grid,
    Connectivity& connect,
    std::vector<double>& u_values,
    std::vector<double>& v_values)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.boundary_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double div_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID < 0)
            {
                ++i;
                continue;
            }

            div_value += invD_x[i] * u_values[neighborID];
            div_value += invD_y[i] * v_values[neighborID];

            ++i;
        }
        divergence_vector[nodeID] = div_value;
    }
}

void Divergence2Order::applyBoundaryCondition(Grid& grid,
    Connectivity& connect,
    const std::vector<int>& nodeIDs,
    std::vector<double>& u_values,
    std::vector<double>& v_values,
    BoundaryCondition& u_BC_vector,
    BoundaryCondition& v_BC_vector)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };
    std::vector<int> opposite = { 1, 0, 3, 2 };

    for (const auto& nodeID : nodeIDs)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double rhs_value_u, operator_value_u, rhs_value_v, operator_value_v;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID == u_BC_vector.boundary_flag)
            {
                std::tie(rhs_value_u, operator_value_u) = u_BC_vector.returnBCValue(grid, nodeID, u_values, i);
                std::tie(rhs_value_v, operator_value_v) = v_BC_vector.returnBCValue(grid, nodeID, v_values, i);

                divergence_vector[nodeID] += invD_x[i] * -rhs_value_u;
                divergence_vector[nodeID] += invD_y[i] * -rhs_value_v;

                //divergence_vector[nodeID] += invD_x[i] * ((2 * u_values[nodeID]) - u_values[neighbors[opposite[i]]]);
                //divergence_vector[nodeID] += invD_y[i] * ((2 * v_values[nodeID]) - v_values[neighbors[opposite[i]]]);
                ++i;
                continue;
            }
            ++i;
        }
    }
}