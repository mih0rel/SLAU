cmake_minimum_required(VERSION 3.20)
project(SLE)

enable_testing()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(sle_lib
    src/matrices/TridiagonalMatrix.cpp
)
target_sources(sle_lib PUBLIC
    src/decompositions/Householder.hpp
    src/matrices/CSR.hpp
    src/matrices/DenseMatrix.hpp
    src/matrices/TridiagonalMatrix.hpp
    src/matrices/Vector.hpp
    src/solvers/GSeidel.hpp
    src/solvers/GaussQR.hpp
    src/solvers/Jacobi.hpp
    src/solvers/SimpleIterations.hpp
    src/solvers/TridiagonalMatrixSolver.hpp
)
target_include_directories(sle_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(GTest QUIET)
if (GTest_FOUND)
    add_executable(test_DenseMatrix
        tests/test_DenseMatrix.cpp
    )
    target_link_libraries(test_DenseMatrix PRIVATE GTest::gtest_main)
    target_include_directories(test_DenseMatrix PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_DenseMatrix COMMAND test_DenseMatrix)

    add_executable(test_CSR
        tests/test_CSR.cpp
    )
    target_link_libraries(test_CSR PRIVATE GTest::gtest_main)
    target_include_directories(test_CSR PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_CSR COMMAND test_CSR)

    add_executable(test_TridiagonalMatrixSolver
        tests/test_TridiagonalMatrixSolver.cpp
    )
    target_link_libraries(test_TridiagonalMatrixSolver PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_TridiagonalMatrixSolver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_TridiagonalMatrixSolver COMMAND test_TridiagonalMatrixSolver)

    add_executable(test_Jacobi
        tests/test_Jacobi.cpp
    )
    target_link_libraries(test_Jacobi PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_Jacobi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_Jacobi COMMAND test_Jacobi)

    add_executable(test_SimpleIterations
        tests/test_SimpleIterations.cpp
    )
    target_link_libraries(test_SimpleIterations PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_SimpleIterations PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_SimpleIterations COMMAND test_SimpleIterations)

    add_executable(test_GSeidel
        tests/test_GSeidel.cpp
    )
    target_link_libraries(test_GSeidel PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_GSeidel PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_GSeidel COMMAND test_GSeidel)
endif()
