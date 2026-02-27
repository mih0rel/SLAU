#include "gtest/gtest.h"
#include "matrices/DenseMatrix.hpp"


TEST(DENSE, TEST_1) {
    std::vector<std::vector<double>> a{{1, 0, 2, 0}, {0, 0, 3, 4}, {5, 6,0, 8}};
    std::vector<std::vector<double>> b{{3, 4, 8, 0}, {1, 0, 0, 0}, {0, 3, 7, 1}};
    DenseMatrix matrix1(a);
    DenseMatrix matrix2(b);
    int c = 2;
    int d = 3;
    std::vector<double> vect2{2 , 4, -3, 1};
    DenseMatrix a1  = matrix1 * c;
    DenseMatrix a2 = d * matrix2;
    DenseMatrix a3 = a1 + a2;
    DenseMatrix a4 = a3 - matrix1;
    std::vector<double> ans_ =  a4 * vect2;
    std::vector<double> ans{-10,1,18};
    ASSERT_EQ(ans_, ans);
}
