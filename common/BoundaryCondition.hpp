#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"

struct Grid;

struct IBoundaryCondition {
    virtual ~IBoundaryCondition() {}
    virtual std::pair<double, double> returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor) = 0;
};

struct DirichletCondition : IBoundaryCondition {
    double boundary_value;
    int boundary_flag;

    DirichletCondition(double value, int flag);
    virtual std::pair<double, double> returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor) override;
};

struct NeumannCondition : IBoundaryCondition {
    double boundary_value;
    int boundary_flag;

    NeumannCondition(double value, int flag);
    virtual std::pair<double, double> returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor) override;
};

struct BoundaryCondition {
    IBoundaryCondition* strategy;
    int boundary_flag;

    BoundaryCondition(IBoundaryCondition* strategy, int flag);
    std::pair<double, double> returnBCValue(Grid& grid, int nodeID, const std::vector<double>& function_values, int neighbor);
};
