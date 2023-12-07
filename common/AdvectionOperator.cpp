#include "AdvectionOperator.hpp"


AdvectionOperator::AdvectionOperator(Grid& grid)
{
	advection_vector.resize(grid.num_points, 0.0);
}

std::vector<double> AdvectionOperator::calculateOperator(Grid& grid,
    Connectivity& connect,
    std::vector<BoundaryCondition>& BC_vector,
    std::vector<double>& u_values,
    std::vector<double>& v_values,
    std::vector<double>& quantity)
{
    std::fill(advection_vector.begin(), advection_vector.end(), 0.0);
    calculateInnerPoints(grid, connect, u_values, v_values, quantity);
    calculateBoundaryPoints(grid, connect, u_values, v_values, quantity);

    for (int i = 0; i < BC_vector.size(); ++i)
    {
        applyBoundaryCondition(grid, connect, connect.boundary_ID, u_values, v_values, quantity, BC_vector[i]);
    }
    return advection_vector;
}

void AdvectionOperator::calculateInnerPoints(Grid& grid,
    Connectivity& connect,
    std::vector<double>& u_values,
    std::vector<double>& v_values,
    std::vector<double>& quantity)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.inner_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double advection_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            advection_value += invD_x[i] * u_values[nodeID] * quantity[neighborID];
            advection_value += invD_y[i] * v_values[nodeID] * quantity[neighborID];

            ++i;
        }

        advection_vector[nodeID] = advection_value;
    }
}

void AdvectionOperator::calculateBoundaryPoints(Grid& grid,
    Connectivity& connect,
    std::vector<double>& u_values,
    std::vector<double>& v_values,
    std::vector<double>& quantity)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : connect.boundary_ID)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double advection_value = 0.0;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID < 0)
            {
                ++i;
                continue;
            }

            advection_value += invD_x[i] * u_values[nodeID] * quantity[neighborID];
            advection_value += invD_y[i] * v_values[nodeID] * quantity[neighborID];

            ++i;
        }
        advection_vector[nodeID] = advection_value;
    }
}

void AdvectionOperator::applyBoundaryCondition(Grid& grid,
    Connectivity& connect,
    const std::vector<int>& nodeIDs,
    std::vector<double>& u_values,
    std::vector<double>& v_values,
    std::vector<double>& quantity,
    BoundaryCondition& BC)
{
    double invdx = 0.5 / grid.dx;
    double invdy = 0.5 / grid.dy;
    std::vector<double> invD_x = { -invdx, invdx, 0.0, 0.0 };
    std::vector<double> invD_y = { 0.0, 0.0, -invdy, invdy };

    for (const auto& nodeID : nodeIDs)
    {
        const std::vector<int>& neighbors = connect.neighbor_IDs[nodeID];

        double rhs_value, operator_value;
        int i = 0;
        for (const auto& neighborID : neighbors)
        {
            if (neighborID == BC.boundary_flag)
            {
                std::tie(rhs_value, operator_value) = BC.returnBCValue(grid, nodeID, quantity, i);

                advection_vector[nodeID] += invD_x[i] * u_values[nodeID] * -rhs_value;
                advection_vector[nodeID] += invD_x[i] * u_values[nodeID] * operator_value;
                advection_vector[nodeID] += invD_y[i] * v_values[nodeID] * -rhs_value;
                advection_vector[nodeID] += invD_y[i] * v_values[nodeID] * operator_value;
                ++i;
                continue;
            }
            ++i;
        }
    }
}