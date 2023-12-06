#include "ArrayOperations.hpp"

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

void addArrays(std::vector<double>& result, const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vectors must be of the same size to subtract");
    }

    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] + rhs[i];
    }
}