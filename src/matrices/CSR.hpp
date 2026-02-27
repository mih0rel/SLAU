#pragma once
#include "Vector.hpp"
#include<cmath>

template<typename T>
class CSR {
private:
    std::vector<T> values_;
    std::vector<long unsigned int> cols_;
    std::vector<long unsigned int> rows_;
    std::size_t width_;
    std::size_t height_;

public:
    CSR (const std::vector<T>& values, const std::vector<long unsigned int>& cols,const std::vector<long unsigned int>& rows,std::size_t height, std::size_t width):
        values_(values), cols_(cols), rows_(rows), width_(width), height_(height) {}
    CSR (const std::vector<std::vector<T>>& data){
        long unsigned int amount = 0;
        width_ = data[0].size();
        height_= data.size();
        rows_.resize(height_ + 1);
        rows_[0] = 0;
        for(std::size_t i = 0; i < height_; i++){
            for(std::size_t j = 0; j< width_; j++){
                if (data[i][j] != 0) {
                    values_.push_back(data[i][j]);
                    cols_.push_back(j);
                    amount ++;
                   }
            }
            rows_[i+1] = (amount);
        }
    }
    CSR (const std::vector<T>& data, size_t height, size_t width){
        long unsigned int amount = 0;
        height_ = height;
        width_ = width;
        rows_.resize(height_ + 1);
        rows_[0] = 0;
        for(std::size_t i = 0; i < height; i++){
            for(std::size_t j = 0; j< width; j++){
                if (data[i * width + j] != 0) {
                    values_.push_back(data[i * width + j]);
                    cols_.push_back(j);
                    amount ++;
                   }
            }
            rows_[i+1] = (amount);
        }
    }

    const std::vector<long unsigned int>& get_cols() const {return cols_;}
    long unsigned int get_cols(long unsigned int i) const {return cols_[i];}
    const std::vector<long unsigned int>& get_rows() const {return rows_;}
    long unsigned int get_rows(long unsigned int i) const {return rows_[i];}
    const std::vector<T>& get_values() const {return values_;}
    T  get_values(long unsigned int i) const {return values_[i];}


    size_t get_width() const  {return width_;}
    size_t get_height() const {return height_;}

    std::vector<T> operator*(const std::vector<T> &vec){
        T item;
        std::vector<T> res(height_);
        for (std::size_t i = 0; i + 1 < rows_.size(); i++) {
            item = 0;
            for (std::size_t j = rows_[i]; j < rows_[i + 1]; j++) {
                item += values_[j] * vec[cols_[j]];
            }
            res[i] = (item);
        }
        return res;
    }

    template<typename U>
    CSR<T> operator*(U a){
        std::vector<T> res(values_.size());
        for (std::size_t i = 0; i < values_.size(); i++) {
            res[i] = (values_[i] * a);
        }
        CSR<T> ans(res,cols_,rows_, height_, width_);
        return ans;
    }

    T operator()(size_t i, size_t j) const{
        for (long unsigned int k = rows_[i]; k < rows_[i + 1]; k++) {
            if (cols_[k] == j) {
                return values_[k];
            }
        }
        return 0;
    }

    void print_CSR(){
        std::cout << "values: ";
        for(std::size_t i = 0; i < values_.size(); i++) std::cout << values_[i] << " ";
        std::cout << std::endl;
        std::cout << "cols: ";
        for(std::size_t i = 0; i < cols_.size(); i++) std::cout << cols_[i] << " ";
        std::cout << std::endl;
        std::cout << "rows: ";
        for(std::size_t i = 0; i < rows_.size(); i++) std::cout << rows_[i] << " ";
        std::cout << std::endl;    
    }

};

template<typename T>
CSR<T> operator+ (CSR<T>& left, CSR<T>& right){
    long unsigned int n = left.get_height();
    long unsigned int m = left.get_width();
    std::vector<T> res_values;
    std::vector<long unsigned int> res_cols;
    std::vector<long unsigned int> res_rows(n + 1);
    long unsigned int amount = 0;
    res_rows[0]=amount;
    for (unsigned int i = 0; i  < n; i ++){
        for (unsigned int j = 0; j <= m; j++){
            if ((right(i,j) + left(i,j))!=0){
                amount++;
                res_values.push_back(right(i,j) + left(i,j));
                res_cols.push_back(j);
            }
        }
    res_rows[i + 1] = (amount);
    }
    return CSR(res_values, res_cols, res_rows, n, m);
}


template<typename T>
CSR<T> operator- (CSR<T> left, CSR<T> right){
    long unsigned int n = left.get_height();
    long unsigned int m = left.get_width();
    std::vector<T> res_values;
    std::vector<long unsigned int> res_cols;
    std::vector<long unsigned int> res_rows(n + 1);
    long unsigned int amount = 0;
    res_rows[0] = (amount);
    for (unsigned int i = 0; i < n; i ++){
        for (unsigned int j = 0; j <= m; j++){
            if ((right(i,j) - left(i,j))!=0){
                amount++;
                res_values.push_back(left(i,j) - right(i,j));
                res_cols.push_back(j);
            }
        }
    res_rows[i+1] = (amount);
    }
    return CSR(res_values, res_cols, res_rows, n, m);
} 

template<typename T, typename U>
CSR<T> operator*(U a, CSR<T> matrix){
    return matrix * a;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const CSR<T>& csr_matrix){
    for (unsigned int i = 0; i < (csr_matrix.get_height()); i ++){
        for (unsigned int j = 0; j <= *std::max_element(csr_matrix.get_cols().begin(), csr_matrix.get_cols().end()); j++) std::cout << csr_matrix(i, j) << ' ';
    std::cout << std::endl;
    }
    return os;
}
