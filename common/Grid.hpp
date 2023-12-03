#pragma once
#include "ProjectIncludes.hpp"

struct Grid
{
    int nx;
    int ny;
    double dx;
    double dy;
    double dt;
    double nu;
    int num_points;
    Eigen::MatrixXd X_grid;
    Eigen::MatrixXd Y_grid;
    std::vector<double> X;
    std::vector<double> Y;

    Grid(double x1, double x2, double y2, double y1, double dx, double dy, double dt, double nu);

    void createMeshgrid(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
    void createMeshdata();
    void printMeshData(std::vector<double> x);
};
