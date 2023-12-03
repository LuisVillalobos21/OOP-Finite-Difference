#pragma once
#include "ProjectIncludes.hpp"

struct IBoundaryCondition {
    virtual ~IBoundaryCondition() {}
    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) = 0;
};

struct DirichletCondition : IBoundaryCondition {
    double boundary_value;
    int boundary_flag;

    DirichletCondition(double value, double flag);
    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) override;
};

struct NeumannCondition : IBoundaryCondition {
    double boundary_value;
    int boundary_flag;

    NeumannCondition(double value, double flag);
    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) override;
};

struct BoundaryCondition {
    IBoundaryCondition* strategy;
    int boundary_flag;

    BoundaryCondition(IBoundaryCondition* strategy, double flag);
    std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors);
};
