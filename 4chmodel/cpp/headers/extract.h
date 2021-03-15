#ifndef __EXTRACT_H_INCLUDED__
#define __EXTRACT_H_INCLUDED__

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
#include <sys/stat.h>

#include "./ReadFiles.h"
#include "./convert.h"

void extract(const std::string&,const std::string&,const std::vector<int>&,const std::string&,const std::string&,const int&);
void midseptum_point(const std::string&,const std::string&,const std::string&);

#endif