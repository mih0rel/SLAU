#pragma once

#include "matrices/CSR.hpp"

#include <algorithm>
#include <cmath>

template<typename T>
std::vector<double> jacobi( CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0, double break_condition){
    std::vector<double> x(x_0.size());
    double tmp = 0;
    double r_cur = 1e8;
    std::vector<double> res(x_0.size()); 
    T d = 0;
    for(; r_cur > break_condition;){
        for(std::size_t j = 0; j < x_0.size(); j++){
            tmp = 0;
            for(std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j+1); k++){
                if(matrix.get_cols(k) != j)
                    tmp += matrix.get_values(k) * x_0[matrix.get_cols(k)];
                else
                    d = matrix.get_values(k);
              }
            x[j] = (b[j] - tmp) / d;
        }
         std::copy(x.begin(), x.end(), x_0.begin());
         res = (matrix * x) - b;
         r_cur = std::sqrt(dot(res, res));
    }
    return x_0;
}
