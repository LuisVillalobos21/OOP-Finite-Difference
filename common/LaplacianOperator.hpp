#pragma once
#include "Grid.hpp"
#include "Connectivity.hpp"

struct Grid;
struct Connectivity;

struct LaplacianOperator
{
    std::vector<double> laplace_vector;

    LaplacianOperator(Grid& grid);

    std::vector<double> LaplacianOperator::calculateOperator(Grid& grid, Connectivity& connect, std::vector<std::vector<int>>& boundaryIDs, std::vector<double>& function_values);
    void calculateInnerPoints(Grid& grid, Connectivity& connect, std::vector<double>& function_values);
    void calculateBoundaryPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values);
};
