#ifndef __HEART_METRICS_H_INCLUDED__
#define __HEART_METRICS_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <limits>
#include <fstream>

#include "./ReadFiles.h"
#include "./vectorops.h"
#include "./convert.h"

struct AT_vol
{
	double el_ActTime;
    double el_Volume;
};

struct matrix_vec
{
    std::vector<double> vec;
};

std::vector<double> surfaceVolume(const std::string &,const std::string &,const std::string &);
double dist2pts(const std::vector<double> &, const points&);
double volume_mesh(const std::string &,const std::string &,const std::string &,const std::string &,const bool &,const std::string &,const std::string &);
double aortaHeight(const std::string &path, const std::string &, const std::string &);
std::vector<double> AT_metrics(const std::string &,const std::string &,const std::string &,const std::string &,const std::string &,const std::string &);
std::vector<double> AT_metrics_BiV(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void BuildHClusterTable(const std::string&);
void AT_metrics_midseptum(const std::string &, const std::string &, std::string &);
std::vector<double> AT_metrics_BiV_optimised(const elem&,std::vector<double>,const std::string&, const std::string&);
void BuildHClusterTable_optimised(const std::string&,const std::string&);
void BuildHClusterTable_RV(const std::string&,const std::string&,const std::string&,const std::string&);
void get_RVapex_UVC();
void get_midseptum_AHA(const std::string &,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
#endif
