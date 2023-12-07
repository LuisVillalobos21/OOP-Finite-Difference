#pragma once
#include "ProjectIncludes.hpp"
#include "Grid.hpp"
#include "Connectivity.hpp"
#include "BoundaryCondition.hpp"
#include "GenericOperator.hpp"
#include "ArrayOperations.hpp"

struct Grid;
struct Connectivity;
struct BoundaryCondition;
struct GenericOperator;

using AssembleLHS = std::function<std::vector<double>(
	Grid&,
	Connectivity&,
	std::vector<BoundaryCondition>&,
	GenericOperator&,
	std::vector<double>&)>;

std::vector<double> laplaceLHS(
	Grid& grid,
	Connectivity& connect,
	std::vector<BoundaryCondition>& BC_vector,
	GenericOperator& op,
	std::vector<double>& x);

std::vector<double> R_LHS(
	Grid& grid,
	Connectivity& connect,
	std::vector<BoundaryCondition>& BC_vector,
	GenericOperator& op,
	std::vector<double> x);