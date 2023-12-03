#pragma once
#include "ProjectIncludes.hpp"
#include "LaplacianOperator.hpp"  


double dotProduct(std::vector<double>& x, std::vector<double>& y);

double dotProduct(std::vector<double>& x, std::vector<double>& y);

std::vector<double> ConjugateGradient(
    LaplacianOperator& laplaceOperator,
    std::vector<double>& b,
    double tolerance,
    Grid& grid,
    Connectivity& connect,
    std::vector<BoundaryCondition>& BC_struct_vector,
    std::vector<double>& function_values);

