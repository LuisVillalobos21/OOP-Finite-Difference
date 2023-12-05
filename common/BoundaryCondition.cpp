#include "BoundaryCondition.hpp"

DirichletCondition::DirichletCondition(double value, int flag) : boundary_value(value), boundary_flag(flag) {}

std::pair<double, double> DirichletCondition::returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor)
{
    return { -boundary_value, 0.0 };
}

NeumannCondition::NeumannCondition(double value, int flag) : boundary_value(value), boundary_flag(flag) {}

std::pair<double, double> NeumannCondition::returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor)
{
    //std::vector<double> signNeumann = { 1.0, 1.0, 1.0, -1.0 };
    //std::vector<double> invhNeumann = { grid.dx, grid.dx, grid.dy, grid.dy };
    //double rhs_value;
    //rhs_value = boundary_value * invhNeumann[neighbor] * signNeumann[neighbor];

    return { 0.0, function_values[nodeID]}; // the ONLY neumann condition is zero for wall boundary condition for pressure, fuck it
}

BoundaryCondition::BoundaryCondition(IBoundaryCondition* strategy, int flag) : strategy(strategy), boundary_flag(flag) {}

std::pair<double, double> BoundaryCondition::returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor)
{
    return strategy->returnBCValue(grid, nodeID, function_values, neighbor);
}

BoundaryCondition makeBoundaryCondition(BoundaryConditionType type, double value, int flag) 
{
    IBoundaryCondition* condition = nullptr;

    switch (type) 
    {
    case BoundaryConditionType::Dirichlet:
        condition = new DirichletCondition(value, flag);
        break;
    case BoundaryConditionType::Neumann:
        condition = new NeumannCondition(value, flag);
        break;
    default:
        throw std::invalid_argument("Unknown Boundary Condition Type");
    }

    return BoundaryCondition(condition, flag);
}
