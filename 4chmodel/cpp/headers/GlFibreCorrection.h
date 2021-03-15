#ifndef __GLFIBRECORRECTION_H_INCLUDED__
#define __GLFIBRECORRECTION_H_INCLUDED__

#include "./ReadFiles.h"

void normalizeLon(const std::string&,const std::string&);
std::vector<double> matVec_prod(const std::vector<std::vector<double>>&,const std::vector<double>&);
double dot_prod(const std::vector<double> &,const std::vector<double>&);
std::vector<double> cross_prod(const std::vector<double>&,const std::vector<double>&);
std::vector<double> scalar_prod(const std::vector<double>&,const double&);
std::vector<double> orthogonalise(const std::vector<double>&,const std::vector<double>&);
void orthogonalise_lon(const std::string&,const std::string&);
void GlFibreCorrection(const std::string&,const std::string&);

#endif