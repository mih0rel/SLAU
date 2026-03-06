#pragma once

#include "matrices/DenseMatrix.hpp"

#include <cmath>
#include <utility>

template<typename T>
std::pair<DenseMatrix<double>, DenseMatrix<double>> householder(const DenseMatrix<T>& matrix_in){
    std::size_t width =  matrix_in.get_width();
    std::size_t height = matrix_in.get_height();
    std::vector<double> v_i;
    std::vector<double> R = matrix_in.get_values();
    std::vector<double> Q(height * height, 0.0);
    for(std::size_t i = 0; i < height; i++) Q[i * width + i] = 1;

    double a = 0;
    for(std::size_t i = 0; i < width; i++){
        std::vector<T> x_i(height);
        for(std::size_t j = i; j < height; j++) x_i[j - i] = R[j * width + i];
        v_i = x_i;
        v_i[0] = (std::abs(x_i[0]) > 1e-10) ? v_i[0] + x_i[0] / std::abs(x_i[0]) * std::sqrt(dot(x_i, x_i)) : v_i[0] + std::sqrt(dot(x_i, x_i));
  
        if (std::sqrt(dot(x_i, x_i)) < 1e-10) continue;
        for(std::size_t j = i; j < width; j ++ ){
            a = 0.0;
            for(std::size_t k = i; k < height; k ++) a += R[k * width + j] * v_i[k - i];
            for(std::size_t k = i; k < height; k ++) R[k * width + j] -= 2.0 / (dot(v_i, v_i)) * a * v_i[k - i];
        }

        for(std::size_t j = 0; j < height; j++){
            a = 0.0;
            for(std::size_t k = i; k < height; k ++) a += Q[k * width + j] * v_i[k - i];
            for(std::size_t k = i; k < height; k ++) Q[k * width + j] -= 2.0 / (dot(v_i, v_i)) * a * v_i[k - i];
        }
        x_i.pop_back();
    }
    double tmp;
    for(std::size_t i = 0; i < height; i ++) {
        for(std::size_t j = i; j < height; j ++){
            tmp = Q[i * height + j];
            Q[i * height + j] = Q[j * height + i];
            Q[j * height + i] = tmp;
        }
    }
    DenseMatrix<double> Q_ans(Q, height, height);
    DenseMatrix<double> R_ans(R, height, width);
    std::pair<DenseMatrix<double>, DenseMatrix<double>> QR(Q_ans, R_ans);
    return QR;

}
