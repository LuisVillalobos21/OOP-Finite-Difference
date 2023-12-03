#pragma once
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "LaplacianOperator.hpp"

struct GradientOperator;
struct DivergenceOperator;
struct AdvectionOperator;

struct GenericOperator
{
    LaplacianOperator laplacianOp;
    GradientOperator* gradientOp;      
    DivergenceOperator* divergenceOp; 
    AdvectionOperator* advectionOp;   

    GenericOperator(Grid& grid) : laplacianOp(grid), gradientOp(nullptr), divergenceOp(nullptr), advectionOp(nullptr) {}

    // Destructor - make sure to properly manage memory if you allocate it for the placeholders
    ~GenericOperator() {
        delete gradientOp;
        delete divergenceOp;
        delete advectionOp;
    }

    // Methods for LaplacianOperator
    std::vector<double> applyLaplacian(Grid& grid, Connectivity& connect, std::vector<std::vector<int>>& boundaryIDs, std::vector<double>& function_values) {
        return laplacianOp.calculateOperator(grid, connect, boundaryIDs, function_values);
    }

    // Placeholder methods for other operators
    // void applyGradient(...); 
    // void applyDivergence(...);
    // void applyAdvection(...);

};
