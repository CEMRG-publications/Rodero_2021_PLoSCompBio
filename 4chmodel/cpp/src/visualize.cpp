/**
 * @brief Set of functions to facilitate the visualization of the data produced.
 * 
 * @file visualize.cpp
 * @author Marina Strocchi
 * @date 2018
 */

#include "../headers/visualize.h"

/**
 * @brief Scales the coordinates of a <CODE>.pts</CODE> file of a given scale factor.
 * 
 * @param path2pts (<CODE>string</CODE>) path to the <CODE>.pts</CODE> file. 
 * @param ptsName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file to scale.
 * @param scale (<CODE>double</CODE>) scale factor to apply.
 */
void scalePts(const std::string &path2pts,const std::string &ptsName,const double &scale)
{
	std::string outPath;
	points pts = ReadPts( path2pts, ptsName );

	std::fstream fOutPts;
	outPath = path2pts + "/" + ptsName + "_scaled.pts";
	write2file(outPath, fOutPts );
	fOutPts << pts.x.size() << std::endl;

	for( int i = 0; i < pts.x.size(); i++ )
	{

		fOutPts << pts.x[i]*scale << " " << pts.y[i]*scale << " " << pts.z[i]*scale << std::endl;

	}

	fOutPts.close();

}