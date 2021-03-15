#ifndef __READFILES_H_INCLUDED__
#define __READFILES_H_INCLUDED__

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
#include <sys/types.h>
#include <dirent.h>

#include <numeric>      // std::iota
#include <functional>   // std::bind
using namespace std::placeholders;


// #include <boost/numeric/ublas/matrix_sparse.hpp>
// #include <boost/numeric/ublas/io.hpp>
// #include <eigen3/Eigen/Sparse>
/**
 * @brief Points with the form <VAR>(x,y,z)</VAR>.
 * 
 */
struct points
{
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
};
/**
 * @brief Surface - Indices of the nodes forming the triangles. 
*/

struct surf
{
	std::vector<int> v1;
	std::vector<int> v2;
	std::vector<int> v3;
};
/**
 * @brief Surface for Neumann boundary conditions in CARP - indices of the nodes 
   forming the triangles, the index of the fourth node of the tetrahedron
   and the index of the tetrahedron containing the surface triangle.
 */
struct neubc
{	
	int nElem;
	int nNodes;
	std::vector<int> v1;
	std::vector<int> v2;
	std::vector<int> v3;
	std::vector<int> v4;
	std::vector<int> tet;
};
/**
 * @brief Elements - indices of the fourth nodes forming the tetrahedra and the 
   element tags defining different regions in the mesh. 
 */
struct elem
{	
	int nElem;
	std::vector<int> v1;
	std::vector<int> v2;
	std::vector<int> v3;
	std::vector<int> v4;
	std::vector<int> tag;
};

/**
 * @brief Fibres direction file - a unitary vector defining the fiber direction
   <VAR>(x1,y1,z1)</VAR> and a unitary vector defining the sheet direction <VAR>(x2,y2,z2)</VAR>. 
 */
struct lon
{	
	int nDir;
	std::vector<double> x1;
	std::vector<double> y1;
	std::vector<double> z1;
	std::vector<double> x2;
	std::vector<double> y2;
	std::vector<double> z2;
};

struct edges
{
	std::vector<int> n1;
	std::vector<int> n2;
};

struct UVC
{
	std::vector<double> Z; //apicobasal
	std::vector<double> RHO; //endoepi
	std::vector<double> PHI; //rotational
	std::vector<double> V; // LV-RV
};

struct EP_metrics_table
{
    std::vector<double> TAT;
    std::vector<double> TAT_LV;
    std::vector<double> AT_10_90;
    std::vector<std::string> HF; // 0 or 1
    std::vector<std::string> heart; // 01 to 20
    std::vector<std::string> lead_position; // Anterior, Antero-Lateral, Lateral, Latero-Posterior, Posterior
    std::vector<std::string> comb_elec; // 1, 2, ..., 8
};


struct HD_metrics_table
{

    std::vector<std::string> heart; // 01 to 20
	std::vector<double> EDV_LV;
	std::vector<int> t_EDV_LV;
	std::vector<int> t_EDV_RV;
	std::vector<double> EDV_RV;
	std::vector<double> ESV_LV;
	std::vector<int> t_ESV_LV;
	std::vector<int> t_ESV_RV;
	std::vector<double> ESV_RV;
	std::vector<double> SV;
	std::vector<double> EF;

};



std::vector<int> ReadVtx(const std::string &,const std::string &);
points ReadPts(const std::string&,const std::string&);
surf ReadSurfElem(const std::string&,const std::string&);
surf ReadSurf(const std::string &,const std::string&);
elem ReadElem(const std::string&,const std::string&,bool outflag=1);
lon ReadLon(const std::string&,const std::string&);
neubc ReadNeubc(const std::string&,const std::string&);
std::vector<double> ReadActTime(const std::string &,const std::string&,bool outflag=1);
std::vector<double> ReadVolFile(const std::string &,const std::string&);
void mergeVtx(const std::string&,const std::vector<std::string>&,const std::string&,const std::string&);
void mergeNeubc(const std::string&,const std::vector<std::string>&,const std::string&,const std::string&);
void mergeSurf(const std::string&,const std::vector<std::string>&,const std::string&,const std::string&);
void mergeLon(const std::string&,const std::string&,const std::string&,const std::string&,const std::vector<int>&,const std::string&,const std::string&,const std::vector<int>&);
void mergeElem(const std::string&,const std::vector<std::string>&,const std::string&,const std::string&);
void mergePts(const std::string&,const std::vector<std::string>&,const std::string&,const std::string&);
void defaultLon(const std::string&,const std::string&);
void write2file(std::string&,std::fstream&);
std::vector<double> extract_point(const points &, const int &);
void change_tags(const std::string &, const std::string &, const std::vector<int> &, const std::vector<int> &);
edges getEdges(const std::string &, const std::string &, const std::string &, const std::string &, const bool &);
bool compare(double,double,std::vector<double>);
std::string GotoLine(const std::string&,const int&);
UVC ReadUVC(const std::string&,const std::string&);
void correctmeshtoolfuckup(const std::string&,const std::string&,const std::string&,const std::string&,const std::vector<int>&,const std::string&,const std::string&);
EP_metrics_table ReadEP_metrics_table(const std::string&,const std::string&);
std::vector<std::string> read_directory(const std::string&);
void dist_betw_pts(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,const std::string&);
std::vector<std::vector<double>> ReadErrorTable(const std::string &path,const std::string &filename);
void FindFirstWave(std::string& , std::string& , std::string&, std::string&,std::string&, std::string&, std::string& , std::string&);

#endif