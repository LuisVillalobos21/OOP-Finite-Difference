#pragma once
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "BoundaryCondition.hpp"

struct Grid;
struct Connectivity;
struct BoundaryCondition;

struct LaplacianOperator
{
    std::vector<double> laplace_vector;
    std::vector<double> rhs_vector;

    LaplacianOperator(Grid& grid);

    std::vector<double> LaplacianOperator::calculateOperator(Grid& grid, Connectivity& connect, std::vector<BoundaryCondition>& BC_vector, std::vector<double>& function_values);
    void calculateInnerPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values);
    void calculateBoundaryPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values);
    void LaplacianOperator::applyBoundaryCondition(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values, BoundaryCondition& boundaryCondition);
};
