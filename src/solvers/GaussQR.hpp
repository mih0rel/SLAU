#pragma once

#include "decompositions/Householder.hpp"

template<typename T>
std::vector<double> gauss_qr(const DenseMatrix<T>& A, const std::vector<T>& b){
    std::vector<double> x(b.size());
    const auto& [Q, R] = householder(A);
    std::vector<double> y = Q.transpose() * b;
    double tmp;
    for(std::size_t i = b.size() - 1; i > 0; --i){
        tmp = 0;
        for(std::size_t  j = i + 1; j < b.size(); j++){
            tmp += R(i, j) * x[j];
        }
        x[i] = (y[i] - tmp) / R(i,i);
    }
    tmp = 0;
    for(std::size_t j =  1; j < b.size(); j++){
            tmp += R(0, j) * x[j];
        }
    x[0] = (y[0] - tmp) / R(0,0);
    return x;
}
