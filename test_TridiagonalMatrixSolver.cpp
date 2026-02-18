#include "gtest/gtest.h"
#include "TridiagonalMatrixSolver.hpp"


TEST(TRIDIAGONAL, TEST_1) {
    std::vector<double> vect1{5, 1};
    std::vector<double> vect2{2 , 4, -3};
    std::vector<double> vect3{-1, 2};
    std::vector<double> d{3, 6, 2};
    TridiagonalMatrix matrix(vect1, vect2, vect3);
    std::vector<double> x = solve_tridiagonal_matrix(matrix, d);
    ASSERT_NEAR(x[0], 1.48837209302326, 1e-14);
    ASSERT_NEAR(x[1], -0.023255813953488, 1e-14);
    ASSERT_NEAR(x[2], -0.674418604651163, 1e-14);
}

TEST(TRIDIAGONAL, TEST_2) {
    std::vector<double> vect1{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    std::vector<double> vect2{2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9};
    std::vector<double> vect3{1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8};
    std::vector<double> d{5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9};
    TridiagonalMatrix matrix(vect1, vect2, vect3);
    std::vector<double> x = solve_tridiagonal_matrix(matrix, d);
    ASSERT_NEAR(x[0], 1.6936085168050307, 1e-12);
    ASSERT_NEAR(x[1], 1.6127829663899385, 1e-12);
    ASSERT_NEAR(x[2], 1.4034499262732964, 1e-12);
    ASSERT_NEAR(x[3], 1.4915446407673001, 1e-12);
    ASSERT_NEAR(x[4], 1.1141633448870929, 1e-12);
    ASSERT_NEAR(x[5], 1.5209929399743264, 1e-12);
    ASSERT_NEAR(x[6], 0.7602906517470917, 1e-12);
    ASSERT_NEAR(x[7], 1.6941553384206036, 1e-12);
    ASSERT_NEAR(x[8], 0.3491630176714153, 1e-12);
    ASSERT_NEAR(x[9], 1.9261218221019746, 1e-12);
}
