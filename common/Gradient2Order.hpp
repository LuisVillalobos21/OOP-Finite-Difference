#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "BoundaryCondition.hpp"

struct Grid;
struct Connectivity;
struct BoundaryCondition;

struct Gradient2Order
{
    std::vector<double> gradient_vector;
    //std::vector<double> rhs_vector;

    Gradient2Order(Grid& grid);

    std::vector<double> calculateOperator(Grid& grid, Connectivity& connect, std::vector<BoundaryCondition>& BC_vector, std::vector<double>& function_values);
    void calculateInnerPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values);
    void calculateBoundaryPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values);
    //void applyBoundaryCondition(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values, BoundaryCondition& boundaryCondition);
};
