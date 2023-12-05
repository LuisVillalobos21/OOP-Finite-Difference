#pragma once
#include "ProjectIncludes.hpp"
#include "LaplacianOperator.hpp"
#include "Gradient2Order.hpp"
#include "AdvectionOperator.hpp"
#include "Divergence2Order.hpp"

struct Grid;

struct GenericOperator 
{
    LaplacianOperator laplace;
    Gradient2Order grad2;
    AdvectionOperator advec;
    Divergence2Order div;

    GenericOperator(Grid& grid) : laplace(grid), grad2(grid), advec(grid), div(grid) 
    {
    }
};