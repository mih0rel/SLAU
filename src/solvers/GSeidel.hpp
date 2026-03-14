#pragma once

#include "matrices/CSR.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

template<typename T>
std::vector<T> gseidel( CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0, double break_condition){
    double r_cur = 1e8;
    std::vector<double> res(x_0.size()); 
    std::vector<T> x_old(x_0.size());
    for(; r_cur > break_condition;){
        std::copy(x_0.begin(), x_0.end(), x_old.begin());
        for(std::size_t j = 0; j < x_0.size(); j++){
            T L = 0;
            T U = 0;
            T d = 0;
            for(std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j + 1); k++){
                if(matrix.get_cols(k) < j)
                    L += matrix.get_values(k) * x_0[matrix.get_cols(k)];
                else if(matrix.get_cols(k) > j)
                    U += matrix.get_values(k) * x_old[matrix.get_cols(k)];
                else
                    d = matrix.get_values(k);
            }
            if(d == 0) throw std::runtime_error("0");
            x_0[j] = (b[j] - L - U) / d;
        }
         res = (matrix * x_0) - b;
         r_cur = std::sqrt(dot(res, res));
    }
    return x_0;
}
