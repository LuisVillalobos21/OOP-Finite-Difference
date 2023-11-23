#pragma once
#include "ProjectIncludes.hpp"

struct Grid;
struct Connectivity;

class Operator {
public:
    virtual ~Operator() {}

    virtual void calculateOperator() = 0;
    virtual void calculateInnerPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values) = 0;
    virtual void calculateBoundaryPoints(Grid& grid, Connectivity& connect, const std::vector<int>& nodeIDs, std::vector<double>& function_values) = 0;
};