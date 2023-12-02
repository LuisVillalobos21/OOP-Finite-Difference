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

std::vector<double> ConjugateGradient(const std::vector<double>& b, double tolerance, std::function<void(std::vector<double>&)> assembleLHSFunc)
{
    std::vector<double> x(b.size(), 0.0);
    std::vector<double> r = b;
    std::vector<double> d = r;
    std::vector<double> lhs_vector(b.size());

    double delta_new = dotProduct(r, r);
    double delta_0 = delta_new;

    for (size_t i = 0; i < b.size(); ++i)
    {
        assembleLHSFunc(d);

        double denom = dotProduct(d, lhs_vector);
        if (denom == 0.0) {
            break;
        }
        double alpha = delta_new / denom;

        for (size_t j = 0; j < x.size(); ++j)
        {
            x[j] += alpha * d[j];
            r[j] -= alpha * lhs_vector[j];
        }

        double delta_old = delta_new;
        delta_new = dotProduct(r, r);

        std::cout << "Iteration " << i + 1 << ", Residual norm: " << sqrt(delta_new) << std::endl;

        if (delta_new < tolerance * tolerance * delta_0)
        {
            break;
        }

        double beta = delta_new / delta_old;
        for (size_t j = 0; j < d.size(); ++j)
        {
            d[j] = r[j] + beta * d[j];
        }
    }

    return x;
}
