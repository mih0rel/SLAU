#pragma once
#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>


template<typename T>
void operator+=(std::vector<T>& rhs, const std::vector<T>& lhs){
    for (unsigned int i = 0; i < rhs.size(); i++) rhs[i] += lhs[i];
}
template<typename T>
std::vector<T> operator+(const std::vector<T>& rhs, const std::vector<T>& lhs){
    std::vector<T> res;
    for (unsigned int i = 0; i < rhs.size(); i++) res.push_back(rhs[i] + lhs[i]);
    return res;
}
template<typename T>
void operator-=(std::vector<T>& rhs, const std::vector<T>& lhs){
    for (unsigned int i = 0; i < rhs.size(); i++) rhs[i] -= lhs[i];
}
template<typename T>
std::vector<T> operator-(const std::vector<T>& rhs, const std::vector<T>& lhs){
    std::vector<T> res;
    for (unsigned int i = 0; i < rhs.size(); i++) res.push_back(rhs[i] - lhs[i]);
    return res;
}
template<typename T>
std::vector<T> operator*(const std::vector<T> &rhs, T lhs){
	std::vector<T> res;
	for(size_t i = 0; i < rhs.size(); i++) res.push_back(rhs[i] * lhs);
	return res;
}
template<typename T>
std::vector<T> operator*(T rhs, const std::vector<T> &lhs){
	return lhs*rhs;
}

template<typename T>
void operator*=(std::vector<T> &rhs, T lhs){
	for(size_t i = 0; i < rhs.size(); i++) rhs[i] *= lhs;
}

template<typename T>
void operator/=(std::vector<T> &rhs, T lhs){
	for(size_t i = 0; i < rhs.size(); i++) rhs[i] /= lhs;
}

template<typename T>
T dot(const std::vector<T>& rhs, const std::vector<T>& lhs){
    T res = 0;
    for(size_t i = 0; i < rhs.size(); i++) res += (rhs[i] * lhs[i]);
    return res;
}

template<typename T>
double norm(const std::vector<T> rhs){
    double res = 0;
    for(size_t i = 0; i < rhs.size(); i++) res += (rhs[i] * rhs[i]);
    return std::sqrt(res);
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& lhs){
    for (unsigned int i = 0; i < lhs.size(); i ++) std::cout << lhs[i] << ' ';
    std::cout << std::endl;
    return os;
}