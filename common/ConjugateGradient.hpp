#pragma once
#include "ProjectIncludes.hpp"
#include "GenericOperator.hpp"  
#include "AssmbleLHS.hpp"
#include "ArrayOperations.hpp"

//double dotProduct(std::vector<double>& x, std::vector<double>& y);
//
//void subtractArrays(std::vector<double>& result, const std::vector<double>& lhs, const std::vector<double>& rhs);

std::vector<double> ConjugateGradient(
    AssembleLHS& assemblyFunction,
    std::vector<double>& b,
    double tolerance,
    Grid& grid,
    Connectivity& connect,
    GenericOperator& op,
    std::vector<BoundaryCondition>& BC_struct_vector,
    std::vector<double>& function_values);

