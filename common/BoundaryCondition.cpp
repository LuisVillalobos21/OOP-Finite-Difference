#include "BoundaryCondition.hpp"

DirichletCondition::DirichletCondition(double value, double flag) : boundary_value(value), boundary_flag(flag) {}

std::pair<double, double> DirichletCondition::returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors)
{
    return { boundary_value, 0.0 };
}

NeumannCondition::NeumannCondition(double value, double flag) : boundary_value(value), boundary_flag(flag) {}

std::pair<double, double> NeumannCondition::returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) 
{
    return { boundary_value, function_values[nodeID]};
}

BoundaryCondition::BoundaryCondition(IBoundaryCondition* strategy, double flag) : strategy(strategy), boundary_flag(flag) {}

std::pair<double, double> BoundaryCondition::returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors)
{
    return strategy->returnBCValue(nodeID, function_values, neighbors);
}
