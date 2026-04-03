#include "matrices/CSR.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

template<typename T>
std::vector<double> gseidel_iter(CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0){
    std::vector<T> x_old(x_0.size());
    std::copy(x_0.begin(), x_0.end(), x_old.begin());
    for(std::size_t i = 0; i < x_0.size(); i++){
        std::size_t j = x_0.size() - 1 - i ;
        T L = 0;
        T U = 0;
        T d = 0;
        for(std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j + 1); k++){
            if(matrix.get_cols(k) < j)
                L += matrix.get_values(k) * x_old[matrix.get_cols(k)];
            else if(matrix.get_cols(k) > j)
                U += matrix.get_values(k) * x_0[matrix.get_cols(k)];
            else
                d = matrix.get_values(k);
        }
        x_0[j] = (b[j] - L - U) / d;
    }

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
        x_0[j] = (b[j] - L - U) / d;
    }
    return x_0;
}


template<typename T>
std::vector<double> symmetric_gseidel( CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0, double break_condition){
    double r_cur = 1e8;
    std::vector<double> res(x_0.size()); 
    for(; r_cur > break_condition;){
         x_0 = gseidel_iter(matrix, b, x_0);
         res = (matrix * x_0) - b;
         r_cur = std::sqrt(dot(res, res));
    }
    return x_0;
}


template<typename T>
std::vector<double> accelerated_gseidel( CSR<T>& matrix, const std::vector<T>& b, std::vector<T> y_0, double rho, double break_condition){
    double r_cur = 1e8;
    std::vector<double> res(y_0.size()); 

    double mu_0 = 1, mu_1 = 1 / rho;

    std::vector<double> y_1 = gseidel_iter(matrix, b, y_0);
    res = (matrix * y_1) - b;
    r_cur = std::sqrt(dot(res, res));

    std::vector<T> y(y_0.size());
    for(; r_cur > break_condition;){
        y_1 = gseidel_iter(matrix, b, y_1);
        y = 2 * mu_1 / rho * y_1 - mu_0 * y_0;
        mu_0 = 2 / rho * mu_1 - mu_0;
        y /= mu_0;

        y_0 = y_1;
        y_1 = y;

        std::swap(mu_0, mu_1);


        res = (matrix * y) - b;
        r_cur = std::sqrt(dot(res, res));

    }

    return y_1;
}
