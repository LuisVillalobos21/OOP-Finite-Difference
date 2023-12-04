#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "BoundaryCondition.hpp"

struct Grid;
struct Connectivity;
struct BoundaryCondition;

struct AdvectionOperator
{
    std::vector<double> advection_vector;

    AdvectionOperator(Grid& grid);

    std::vector<double> calculateOperator(Grid& grid, 
        Connectivity& connect, 
        std::vector<BoundaryCondition>& BC_vector, 
        std::vector<double>& u_values, 
        std::vector<double>& v_values,
        std::vector<double>& quantity);

    void calculateInnerPoints(Grid& grid, 
        Connectivity& connect, 
        std::vector<double>& u_values,
        std::vector<double>& v_values,
        std::vector<double>& quantity);

    void calculateBoundaryPoints(Grid& grid, 
        Connectivity& connect, 
        std::vector<double>& u_values,
        std::vector<double>& v_values,
        std::vector<double>& quantity);

    void applyBoundaryCondition(Grid& grid, 
        Connectivity& connect, 
        const std::vector<int>& nodeIDs, 
        std::vector<double>& u_values, 
        std::vector<double>& v_values, 
        std::vector<double>& quantity,
        BoundaryCondition& BC);
};
