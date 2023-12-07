#include "ConjugateGradient.hpp"

std::vector<double> ConjugateGradient(
    AssembleLHS& assembleLHSFunction,
    std::vector<double>& b,
    double tolerance,
    Grid& grid,
    Connectivity& connect,
    GenericOperator& op,
    std::vector<BoundaryCondition>& BC_struct_vector,
    std::vector<double> x,
    double pin_option) 
{
    std::vector<double> r;
    std::vector<double> q;
    r.resize(grid.num_points);
    q.resize(grid.num_points);
    subtractArrays(r, b, assembleLHSFunction(grid, connect, BC_struct_vector, op, x));
    std::vector<double> d = r;
    double delta_new = dotProduct(r, r);
    double delta_0 = delta_new;
    double delta_old;

    size_t i = 0; 
    while (i < b.size())
    {
        q = assembleLHSFunction(grid, connect, BC_struct_vector, op, d);

        if (dotProduct(d, q) == 0.0) {
            std::cout << "Early termination: dot product of d and q is zero.\n";
            return x;
        }

        double alpha = delta_new / dotProduct(d, q);

        for (size_t j = 0; j < x.size(); ++j)
        {
            x[j] += alpha * d[j];
        }

        //if (i % 50 == 0) 
        //{
            subtractArrays(r, b, assembleLHSFunction(grid, connect, BC_struct_vector, op, x));
        //}
        //else
        //{
        //    for (size_t j = 0; j < r.size(); ++j)
        //    {
        //        r[j] -= alpha * q[j];
        //    }
        //}

        delta_old = delta_new;
        delta_new = dotProduct(r, r);

        //std::cout << "Iteration " << i + 1 << ", Residual norm: " << sqrt(delta_new) << std::endl;

        if (sqrt(delta_new) < tolerance)
        {
            //std::cout << "CONVERGED! Tolerance is below: " << tolerance << '\n';
            return x;
        }

        double beta = delta_new / delta_old;

        for (size_t j = 0; j < d.size(); ++j)
        {
            d[j] = r[j] + beta * d[j];
        }

        ++i; 
    }
    std::cout << "CONJ GRAD FAILED TO CONVERGE IN TOTAL # OF CELL ITERATIONS " << "Residual = " << sqrt(delta_new) << '\n';

    return x;
}
