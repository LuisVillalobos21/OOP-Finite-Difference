#pragma once
#include "ProjectIncludes.hpp"
#include "GenericOperator.hpp"  
#include "AssmbleLHS.hpp"
#include "ArrayOperations.hpp"

std::vector<double> ConjugateGradient(
    AssembleLHS& assemblyFunction,
    std::vector<double>& b,
    double tolerance,
    Grid& grid,
    Connectivity& connect,
    GenericOperator& op,
    std::vector<BoundaryCondition>& BC_struct_vector,
    std::vector<double>& function_values,
    double pin_option);

