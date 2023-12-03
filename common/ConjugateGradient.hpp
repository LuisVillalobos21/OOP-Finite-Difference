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
    std::vector<std::vector<int>>& boundaryIDs,
    std::vector<double>& function_values);

