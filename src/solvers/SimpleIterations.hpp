#pragma once

#include "matrices/CSR.hpp"

#include <cmath>

template<typename T>
std::vector<double> FPI(CSR<T>& matrix, const std::vector<T>& b,  std::vector<T> x_0, double tau, double break_condition){
    double r_cur = 1e8;
    std::vector<double> res(x_0.size()); 
    for(; r_cur > break_condition;){
        res = (matrix * x_0 - b);
        x_0 = x_0 - tau * (matrix * x_0 - b);
        r_cur = std::sqrt(dot(res, res));
    }
    return x_0;
}


