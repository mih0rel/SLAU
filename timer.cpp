#include "src/matrices/CSR.hpp"
#include "src/matrices/DenseMatrix.hpp"

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using steady_clock = std::chrono::steady_clock;

namespace {

constexpr std::size_t N = 3000;
constexpr std::size_t FILL_TO = 9000000;
constexpr std::size_t FILL_STEP = 2500;
constexpr int WARMUP_RUNS = 2;
constexpr int MEASURE_RUNS = 8;

volatile long long benchmark_sink = 0;

std::vector<int> make_vector_values(std::size_t n, int zero_percent, std::uint32_t seed) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> prob_dist(0, 99);
    std::uniform_int_distribution<int> value_dist(1, 100);

    std::vector<int> x;
    x.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        x.push_back(prob_dist(gen) < zero_percent ? 0 : value_dist(gen));
    }
    return x;
}

template <typename MatrixT>
long long measure_us(MatrixT& matrix, const std::vector<int>& x) {
    for (int i = 0; i < WARMUP_RUNS; ++i) {
        std::vector<int> y = matrix * x;
        if (!y.empty()) benchmark_sink += y[0];
    }

    auto start = steady_clock::now();
    for (int i = 0; i < MEASURE_RUNS; ++i) {
        std::vector<int> y = matrix * x;
        if (!y.empty()) benchmark_sink += y[0];
    }
    auto finish = steady_clock::now();

    return std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count() / MEASURE_RUNS;
}

}  // namespace

int main() {
    std::ofstream dense_file("Dense_time_fill.csv");
    std::ofstream csr_file("CSR_time_fill.csv");
    dense_file << "size,time_us\n";
    csr_file << "size,time_us\n";

    std::mt19937 gen(424242u);
    std::uniform_int_distribution<std::size_t> pos_dist(0, N * N - 1);
    std::uniform_int_distribution<int> value_dist(1, 100);

    std::vector<int> dense_values(N * N, 0);
    std::vector<int> x = make_vector_values(N, 20, 123456u);

    std::size_t previous_fill = 0;
    for (std::size_t fill_count = 0; fill_count <= FILL_TO; fill_count += FILL_STEP) {
        for (std::size_t k = previous_fill; k < fill_count; ++k) {
            dense_values[pos_dist(gen)] = value_dist(gen);
        }
        previous_fill = fill_count;

        DenseMatrix<int> dense_matrix(dense_values, N, N);
        CSR<int> csr_matrix(dense_values, N, N);

        long long dense_time = measure_us(dense_matrix, x);
        long long csr_time = measure_us(csr_matrix, x);

        dense_file << fill_count << "," << dense_time << "\n";
        csr_file << fill_count << "," << csr_time << "\n";

        std::cout << "size=" << fill_count
                  << " dense_us=" << dense_time
                  << " csr_us=" << csr_time << '\n';
    }

    std::cout << "benchmark_sink=" << benchmark_sink << '\n';
    return 0;
}
