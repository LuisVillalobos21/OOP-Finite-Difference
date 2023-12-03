#pragma once
#include "ProjectIncludes.hpp"

// Boundary condition interface
struct IBoundaryCondition {
    virtual ~IBoundaryCondition() {}
    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) = 0;
};

// Dirichlet boundary condition
struct DirichletCondition : IBoundaryCondition {
    double fixedValue;

    DirichletCondition(double value) : fixedValue(value) {}

    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) override {
        return { fixedValue, 0.0 }; // Return fixed value for rhs_value and 0.0 for operator value
    }
};

// Neumann boundary condition
struct NeumannCondition : IBoundaryCondition {
    // Specific properties and methods for NeumannCondition
    double boundaryValue;

    NeumannCondition(double value) : boundaryValue(value) {}

    virtual std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) override {
        double rhs_value = /* logic for calculating rhs_value based on Neumann condition */;
        double operator_value = /* logic for calculating operator_value based on Neumann condition */;

        return { rhs_value, operator_value };
    }
};

// BoundaryCondition class with a strategy pattern
struct BoundaryCondition {
    IBoundaryCondition* strategy;

    BoundaryCondition(IBoundaryCondition* strategy) : strategy(strategy) {}

    std::pair<double, double> returnBCValue(int nodeID, const std::vector<double>& function_values, const std::vector<int>& neighbors) {
        return strategy->returnBCValue(nodeID, function_values, neighbors);
    }
};
