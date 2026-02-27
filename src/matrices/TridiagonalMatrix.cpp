// TridiagonalMatrix.cpp
#include "TridiagonalMatrix.hpp"

#include <cmath>
#include <stdexcept>

TridiagonalMatrix::TridiagonalMatrix(std::size_t n)
    : a_(n == 0 ? 0 : n - 1, 0.0), b_(n, 0.0), c_(n == 0 ? 0 : n - 1, 0.0) {}

TridiagonalMatrix::TridiagonalMatrix(std::vector<double> a,
                                     std::vector<double> b,
                                     std::vector<double> c)
    : a_(std::move(a)), b_(std::move(b)), c_(std::move(c)) {
    check_size_consistency_();
}

std::size_t TridiagonalMatrix::size() const noexcept {
    return b_.size();
}

void TridiagonalMatrix::check_size_consistency_() const {
    if (b_.empty()) {
        if (!a_.empty() || !c_.empty()) {
            throw std::invalid_argument("TridiagonalMatrix: for n=0, a and c must be empty.");
        }
        return;
    }

    if (a_.size() + 1 != b_.size() || c_.size() + 1 != b_.size()) {
        throw std::invalid_argument("TridiagonalMatrix: expected |a|=|c|=|b|-1.");
    }
}

void TridiagonalMatrix::check_index_(std::size_t i) const {
    if (i >= size()) {
        throw std::out_of_range("TridiagonalMatrix: index i out of range.");
    }
}

void TridiagonalMatrix::check_offdiag_index_(std::size_t i) const {
    if (i >= a_.size()) {
        throw std::out_of_range("TridiagonalMatrix: off-diagonal index out of range.");
    }
}

double TridiagonalMatrix::a(std::size_t i) const {
    check_offdiag_index_(i);
    return a_[i];
}
double TridiagonalMatrix::b(std::size_t i) const {
    check_index_(i);
    return b_[i];
}
double TridiagonalMatrix::c(std::size_t i) const {
    check_offdiag_index_(i);
    return c_[i];
}

void TridiagonalMatrix::set_a(std::size_t i, double v) {
    check_offdiag_index_(i);
    a_[i] = v;
}
void TridiagonalMatrix::set_b(std::size_t i, double v) {
    check_index_(i);
    b_[i] = v;
}
void TridiagonalMatrix::set_c(std::size_t i, double v) {
    check_offdiag_index_(i);
    c_[i] = v;
}

double TridiagonalMatrix::at(std::size_t i, std::size_t j) const {
    check_index_(i);
    check_index_(j);

    if (i == j) return b_[i];
    if (j + 1 == i) return a_[i - 1]; // A[i][i-1]
    if (i + 1 == j) return c_[i];     // A[i][i+1]
    return 0.0;
}

std::vector<double> TridiagonalMatrix::mul(const std::vector<double>& x) const {
    const std::size_t n = size();
    if (x.size() != n) {
        throw std::invalid_argument("TridiagonalMatrix::mul: x size mismatch.");
    }
    if (n == 0) return {};

    std::vector<double> y(n, 0.0);

    if (n == 1) {
        y[0] = b_[0] * x[0];
        return y;
    }

    // i=0
    y[0] = b_[0] * x[0] + c_[0] * x[1];

    // 1..n-2
    for (std::size_t i = 1; i + 1 < n; ++i) {
        y[i] = a_[i - 1] * x[i - 1] + b_[i] * x[i] + c_[i] * x[i + 1];
    }

    // i=n-1
    y[n - 1] = a_[n - 2] * x[n - 2] + b_[n - 1] * x[n - 1];

    return y;
}
