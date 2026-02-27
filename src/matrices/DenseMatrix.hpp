#pragma once
#include "Vector.hpp"

template <typename T>
class DenseMatrix{
private:
    std::vector<T> matrix_;
    std::size_t height_, width_;   
public:
    DenseMatrix(const std::vector<T>& matrix, const std::size_t height, const std::size_t width):
          matrix_(matrix), height_(height), width_(width){}
    DenseMatrix(const std::vector<std::vector<T>>& data){
          height_ = data.size();
          width_ = data[0].size();
          for(size_t i = 0; i < height_; i++){
            matrix_.insert(matrix_.end(), data[i].begin(), data[i].end());
          }
    }
    DenseMatrix(const DenseMatrix& matrix): matrix_(matrix.get_values()), height_(matrix.get_height()), width_(matrix.get_width()) {}
    DenseMatrix(const T value, const std::size_t height, const std::size_t width):
          matrix_(std::vector(height * width, value)), height_(height), width_(width){}

    size_t get_width() const  {return width_;}
    size_t get_height() const {return height_;}
    const std::vector<T>& get_values() const {return matrix_;}

    std::vector<T> get_col(std::size_t i){
        std::vector<T> ans(height_);
        for(std::size_t j = 0; j < height_;j++){
            ans[j] = matrix_[width_ * j + i];
        }
        return ans;
    }

    void write_col(const std::vector<T>& vector, std::size_t i){
        for(std::size_t j = 0; j < height_;j++){
            matrix_[width_ * j + i] = vector[j];
        }
    }

    std::vector<T> operator*(const std::vector<T> &vector) const {
        std::vector<T> result(height_, 0);
        for (std::size_t i = 0; i < height_; i++) {
            for (std::size_t j = 0; j < width_; j++) {
                result[i] += matrix_[width_ * i + j] * vector[j];
                }
            }
        return result;
    }

    template<typename U>
    DenseMatrix<T> operator*(U a) const{
        std::vector<T> res(height_* width_, 0);
        for (std::size_t i = 0; i < matrix_.size(); i++) {
            res[i] = (matrix_[i] * a);
        }
        DenseMatrix<T> ans(res,height_,width_);
        return ans;
    }

    DenseMatrix<T>& operator=(const DenseMatrix<T>& a){
        matrix_ = a.get_values();
        width_ = a.get_width();
        height_ = a.get_height();
        return *this;
    }

    template<typename U>
    DenseMatrix<T> operator*(const DenseMatrix<U>& other) const {
        std::vector<T> res_matrix(height_ * other.get_width());
        for (std::size_t i = 0; i < height_; i++) {
            for (std::size_t k = 0; k < width_; k++) {
                for (std::size_t j = 0; j < other.get_width(); ++j) {
                    res_matrix[i * other.get_width() + j] += (matrix_[i * width_ + k] * other(k, j));
                }
            }
        }
        DenseMatrix ans(res_matrix, height_, other.get_width());
        return ans;
    } 


    T operator()(long unsigned int i,  long unsigned int j) const{
        return(matrix_[i * width_ + j]);
    }

    T& change_elem(long unsigned int i,  long unsigned int j){
        return(matrix_[i * width_ + j]);
    }

    DenseMatrix<T> transpose() const{
        std::vector<T> matrix(width_ * height_);
        for (std::size_t i = 0; i < height_; i++) {
            for (std::size_t j = 0; j < width_; j++) {
                matrix[j * height_ + i] = matrix_[i * width_ + j];
            }
        }
        DenseMatrix ans(matrix, width_, height_);
        return ans; 

    }
    
    
};



template <typename T>
DenseMatrix<T> operator+(const DenseMatrix<T> & left, const DenseMatrix<T> & right){
    std::vector<T> res(left.get_height() * left.get_width(), 0);
    for (std::size_t i = 0; i < left.get_height(); i ++){
        for (std::size_t j = 0; j < left.get_width(); j ++){
                res[i * left.get_width() + j] = (left(i,j) + right(i, j));
        }
    }
    DenseMatrix<T> ans(res, left.get_height(), left.get_width());
    return ans; 
}

template <typename T>
DenseMatrix<T> operator-(const DenseMatrix<T> & left, const DenseMatrix<T> & right){
    std::vector<T> res(left.get_height() * left.get_width(), 0);
    for (std::size_t i = 0; i < left.get_height(); i ++){
        for (std::size_t j = 0; j < left.get_width(); j ++){
                res[i * left.get_width() + j] = (left(i,j) - right(i, j));
        }
    }
    DenseMatrix<T> ans(res, left.get_height(), left.get_width());
    return ans; 
}

template<typename T, typename U>
DenseMatrix<T> operator*(U a, DenseMatrix<T> matrix){
    return matrix * a;
}

template <typename T>   
std::ostream& operator<<(std::ostream& os, const DenseMatrix<T>& dense_matrix){
    for (size_t i = 0; i < dense_matrix.get_height(); i++) {
        for(size_t j = 0; j < dense_matrix.get_width(); j++) os << dense_matrix(i,j) << " ";
        std::cout << std::endl;
    }
    return os;
}
