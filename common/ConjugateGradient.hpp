#pragma once
#include "ProjectIncludes.hpp"
#include "LHS.hpp"
#include "RHS.hpp"

double dotProduct(std::vector<double>& x, std::vector<double>& y);

std::vector<double> ConjugateGradient(LHS& lhs, RHS& rhs, double tolerance);

