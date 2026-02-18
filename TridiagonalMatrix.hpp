#pragma once
#include <vector>

class TridiagonalMatrix {
public:
    explicit TridiagonalMatrix(std::size_t n);
    TridiagonalMatrix(std::vector<double> a,
                      std::vector<double> b,
                      std::vector<double> c);

    std::size_t size() const noexcept;

    double a(std::size_t i) const;
    double b(std::size_t i) const;
    double c(std::size_t i) const;

    void set_a(std::size_t i, double v);
    void set_b(std::size_t i, double v);
    void set_c(std::size_t i, double v);

    double at(std::size_t i, std::size_t j) const;
    std::vector<double> mul(const std::vector<double>& x) const;

private:
    void check_size_consistency_() const;
    void check_index_(std::size_t i) const;

    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
};
