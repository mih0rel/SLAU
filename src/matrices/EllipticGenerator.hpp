#pragma once

#include "CSR.hpp"
#include <stdexcept>
template<typename T>
struct Elliptic{
    CSR<T> matrix;
    std::vector<T> rhs;
    T h;
};
template<typename T, typename F, typename G>
Elliptic<T> dirichlet(std::size_t n, T L, F f, G g){
    if(n == 0 || L <= 0) throw std::runtime_error("0");

    T h = L / (n + 1);
    T inv_h2 = T(1) / (h * h);
    std::size_t size = n * n;

    std::vector<T> values;
    std::vector<long unsigned int> cols;
    std::vector<long unsigned int> rows(size + 1);
    std::vector<T> rhs(size);

    long unsigned int amount = 0;
    rows[0] = 0;

    for(std::size_t j = 0; j < n; j++){
        for(std::size_t i = 0; i < n; i++){
            std::size_t id = j * n + i;
            T x = (i + 1) * h;
            T y = (j + 1) * h;

            rhs[id] = f(x, y);

            if(j > 0){
                values.push_back(-inv_h2);
                cols.push_back(id - n);
                amount++;
            }
            else rhs[id] += inv_h2 * g(x, 0);

            if(i > 0){
                values.push_back(-inv_h2);
                cols.push_back(id - 1);
                amount++;
            }
            else rhs[id] += inv_h2 * g(0, y);

            values.push_back(4 * inv_h2);
            cols.push_back(id);
            amount++;

            if(i + 1 < n){
                values.push_back(-inv_h2);
                cols.push_back(id + 1);
                amount++;
            }
            else rhs[id] += inv_h2 * g(L, y);

            if(j + 1 < n){
                values.push_back(-inv_h2);
                cols.push_back(id + n);
                amount++;
            }
            else rhs[id] += inv_h2 * g(x, L);

            rows[id + 1] = amount;
        }
    }

    return {CSR<T>(values, cols, rows, size, size), rhs, h};
}
