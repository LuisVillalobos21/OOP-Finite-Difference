#pragma once
#include "ProjectIncludes.hpp"

double dotProduct(std::vector<double>& x, std::vector<double>& y);

std::vector<double> ConjugateGradient(std::vector<double> b, double tolerance, std::function<void()> assembleLHSFunc);

