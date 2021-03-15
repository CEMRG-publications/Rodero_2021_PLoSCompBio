#ifndef __COMMON_H_INCLUDED__
#define __COMMON_H_INCLUDED__

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
#include <iterator>
#include "./ReadFiles.h"

void commonVertices(const std::string&,const std::string&,const std::string&,const std::string&,const int&,const std::string&);
void rmCommonTriangles(const std::string&,const std::string&,const std::string&,const std::string&);
void commonTriangles(const std::string&,const std::string&,const std::string&,const std::string&);
void rmVtxTriangles(std::string path2vtx, std::string vtxName, std::string path2surf2, std::string surfName2);
void rmSeptum(const std::string& path_elem, const std::string& elem_file, const std::string& path_septum, const std::string& septum_file, const std::string& out_path, const std::string& out_file);

#endif