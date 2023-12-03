#include "ConjugateGradient.hpp"

double dotProduct(std::vector<double>& x, std::vector<double>& y)
{
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    double dotProduct = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        dotProduct += x[i] * y[i];
    }

    return dotProduct;
}

void subtractArrays(std::vector<double>& result, const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vectors must be of the same size to subtract");
    }

    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] - rhs[i];
    }
}

std::vector<double> ConjugateGradient(
    LaplacianOperator& laplaceOperator,
    std::vector<double>& b,
    double tolerance,
    Grid& grid,
    Connectivity& connect,
    std::vector<std::vector<int>>& boundaryIDs,
    std::vector<double>& x) 
{
    std::vector<double> r;
    std::vector<double> q;
    r.resize(grid.num_points);
    q.resize(grid.num_points);
    subtractArrays(r, b, laplaceOperator.calculateOperator(grid, connect, boundaryIDs, x));
    std::vector<double> d = r;
    double delta_new = dotProduct(r, r);
    double delta_0 = delta_new;
    double delta_old;

    size_t i = 0; 
    while (i < b.size())
    {
        q = laplaceOperator.calculateOperator(grid, connect, boundaryIDs, d);

        double alpha = delta_new / dotProduct(d, q);

        for (size_t j = 0; j < x.size(); ++j)
        {
            x[j] += alpha * d[j];
        }

        if (i % 50 == 0) 
        {
            subtractArrays(r, b, laplaceOperator.calculateOperator(grid, connect, boundaryIDs, x));
        }
        else
        {
            for (size_t j = 0; j < r.size(); ++j)
            {
                r[j] -= alpha * q[j];
            }
        }

        delta_old = delta_new;
        delta_new = dotProduct(r, r);

        std::cout << "Iteration " << i + 1 << ", Residual norm: " << sqrt(delta_new) << std::endl;

        if (sqrt(delta_new) < tolerance)
        {
            std::cout << "Converged! Tolerance is below: " << tolerance << '\n';
            break;
        }

        double beta = delta_new / delta_old;

        for (size_t j = 0; j < d.size(); ++j)
        {
            d[j] = r[j] + beta * d[j];
        }

        ++i; 
    }

    return x;
}
