#pragma once

#include "matrices/CSR.hpp"

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

template<typename T>
std::vector<T> gseidel(CSR<T>& matrix,
                       const std::vector<T>& b,
                       std::vector<T> x,
                       double break_condition)
{
    double r_cur = std::numeric_limits<double>::infinity();
    std::vector<double> res(x.size());

    while (r_cur > break_condition) {
        std::vector<T> x_old = x;

        for (std::size_t j = 0; j < x.size(); ++j) {
            T L = 0;
            T U = 0;
            T d = 0;

            for (std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j + 1); ++k) {
                std::size_t col = matrix.get_cols(k);
                T a = matrix.get_values(k);

                if (col < j) {
                    L += a * x[col];
                }
                else if (col > j) {
                    U += a * x_old[col];
                }
                else {
                    d = a;
                }
            }

            if (d == 0) {
                throw std::runtime_error("0");
            }

            x[j] = (b[j] - L - U) / d;
        }

        res = (matrix * x) - b;
        r_cur = std::sqrt(dot(res, res));
    }

    return x;
}