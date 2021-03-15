#ifndef __MAP_H_INCLUDED__
#define __MAP_H_INCLUDED__

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

struct mapNodes;
std::vector<int> mapPoints(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void mapSurface(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
std::vector<int> mapElem(const std::string&,const std::string&,const std::string&,const std::string&);
void mapBiV(const std::string&,const std::string&,const std::string&,const std::string&);
void mapExtendedBase(const std::string&,const std::string&,const std::string&,const std::string&);
void mapFibres(const std::string&,const std::string&,const std::string&,const std::string&,const std::vector<int>&,const std::string&);
void mapNeubcBack(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void mapNeubcBack_fast(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void mapSurfBack(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void mapVtxBack(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void mapUVCBack(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
mapNodes mapDuplicate(const std::string&,const std::string&,const std::string&,const std::string&);
void newActSeq(const mapNodes&,const std::string&,const std::string&,const std::string&);
void vtx2neubc(const std::string&,const std::string&,const std::string&,const std::string&);
void map2bivptsdatfile(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
void map_error_varifold(const std::string &case_num);
#endif