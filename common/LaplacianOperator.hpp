#pragma once
#include "Operator.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"

struct Grid;
struct Connectivity;

class LaplacianOperator : public Operator
{
public:
    std::vector<double> laplace_vector;

    LaplacianOperator(Grid& grid);

    void calculateOperator(Grid& grid, Connectivity& connect, std::vector<double>& function_values) override;
    void calculateInnerPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values) override;
    void calculateBoundaryPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values) override;
};
