#ifndef __CONVERT_H_INCLUDED__
#define __CONVERT_H_INCLUDED__

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>

#include "./ReadFiles.h"

void Vtx2pts(const std::string&,const std::string&,const std::string&,const std::string&);
void vtx2vtk(const std::string&,const std::string&,const std::string&,const std::string&);
void surf2vtk(const std::string&,const std::string&,const std::string&,const std::string&);
void neubc2vtk(const std::string&,const std::string&,const std::string&,const std::string&);
void mesh2vtk(const std::string&,const std::string&,const std::string&,const std::string&,const double&);
void vtx2surf(const std::string&,const std::string&,const std::string&,const std::string&);
bool str2bool(const std::string&);
void elem2vtx(const std::string&,const std::string&,const std::vector<int>&,const std::string&, const std::string&);
int UVC2vtx(const UVC&,const double&,const double&,const double&,const int&);
int UVC2vtx(const std::string&,const std::string&,const double&,const double&,const double&,const int&);

#endif