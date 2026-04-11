#pragma once

#include "matrices/CSR.hpp"

template<typename T>
std::vector<double> CG(CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0, double break_condition){
    std::vector<double> r_cur = matrix * x_0 - b, r_prev;
    std::vector<double> d = r_cur;

    for (;norm(r_cur) > break_condition;) {
        x_0 = x_0 - dot(r_cur, r_cur) / (dot(d, (matrix * d))) * d;
        
        r_prev = r_cur;
        r_cur = matrix * x_0 - b;

       
        d = r_cur + ((dot(r_cur,r_cur)) / (dot(r_prev,r_prev))) * d;
    }

    return x_0;
}
