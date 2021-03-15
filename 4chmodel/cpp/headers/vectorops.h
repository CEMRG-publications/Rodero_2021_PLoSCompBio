#ifndef __VECTOROPS_H_INCLUDED__
#define __VECTOROPS_H_INCLUDED__

#include <vector>
#include <math.h>
#include "ReadFiles.h"

std::vector<double> operator+(const std::vector<double> &,const std::vector<double> &);
std::vector<double> operator-(const std::vector<double> &,const std::vector<double> &);
std::vector<double> operator*(const double &, const std::vector<double> &);
std::vector<double> operator*( const std::vector<double> &, const double &);
std::vector<double> operator/( const std::vector<double> &, const double &);
double norm(const std::vector<double> &);
std::vector<double> pts2vec(const points &, const int &);
std::vector<double> abs(const std::vector<double> &);
int min_pos(const std::vector<double> &);
int factorial(const int&);
double max_value_vec(const std::vector<double> &);
double min_positive(const std::vector<double> &);
std::vector<double> min_vec(const std::vector<double> &, const std::vector<double> &);
std::vector<double> min_positive(const std::vector<double> &, const std::vector<double> &);
std::vector<double> cross(const std::vector<double> &,const std::vector<double> &);
double sum(const std::vector<double> &);
double dot(const std::vector<double> &,const std::vector<double> &);
double average(const std::vector<double> &);

#endif