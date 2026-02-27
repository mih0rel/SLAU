#pragma once
#include "matrices/TridiagonalMatrix.hpp"
#include <vector>

std::vector<double> solve_tridiagonal_matrix(const TridiagonalMatrix& matrix,
                                             const std::vector<double>& d) {
    const std::size_t n = matrix.size();
    std::vector<double> x(n);
    std::vector<double> p(n);
    std::vector<double> q(n);

    p[1] = -matrix.c(0) / matrix.b(0);
    q[1] = d[0] / matrix.b(0);
    for (std::size_t i = 1; i + 1 < n; i++) {
        const double denom = matrix.a(i - 1) * p[i] + matrix.b(i);
        p[i + 1] = -matrix.c(i) / denom;
        q[i + 1] = (d[i] - matrix.a(i - 1) * q[i]) / denom;
    }

    x[n - 1] = (d[n - 1] - matrix.a(n - 2) * q[n - 1]) /
               (matrix.a(n - 2) * p[n - 1] + matrix.b(n - 1));
    for (std::size_t i = n - 2; i >= 1; i--) {
        x[i] = p[i + 1] * x[i + 1] + q[i + 1];
    }
    x[0] = p[1] * x[1] + q[1];

    return x;
}
