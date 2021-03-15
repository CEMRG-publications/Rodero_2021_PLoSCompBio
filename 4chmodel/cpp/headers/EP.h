#ifndef __EP_H_INCLUDED__
#define __EP_H_INCLUDED__

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "ReadFiles.h"
#include "vectorops.h"

void electrodeVtx(const std::string&,const std::string&,const int&,const std::string&);
void electrodeEndoVtx(const std::string&,const std::string&,const std::string&,const int&,const std::string&,const double&);
void vtkActivationTime(const std::string&,const std::string&,const std::string&,const std::string&);
void endoLayer(const std::string&,const std::string&,const std::string&,const std::string&,const int&,const std::string&,const std::string&,const double&,const std::string&);
void createElectrodesThread(const std::string &, const std::string &, const std::string &x, const int &, const int &, const double &, const int &, const std::string &, const std::string &);
void createQuadpoles(const std::string&,const std::string&,const int&);
void removeFEC(const std::string&,const std::string&,const std::string&,const int&,const int&);
void createFEC(const std::string &,const std::string &,const std::string &,const std::string &, const std::string&, const std::string&,const double &,const double&);
std::vector<double> scale_UVC( UVC, std::string);
void scale_UVC_file(const std::string&,const std::string&);
int single_electrode_UVC(const UVC&,const std::vector<double>&,const double&,const double&,const double&,const double&,const double&,const double&);
void mergeTwoElectrodes(const std::vector<std::string> &paths, const std::vector<std::string> &electrodes, const std::string &, const std::string&);
void lead2electrodes(const std::string&,const std::string&,const std::string&,const std::string&);
int findRVapex(const UVC&);
int findLVapex(const UVC&);

#endif
