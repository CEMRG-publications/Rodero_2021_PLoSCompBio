#include "../headers/convert.h"

// cd /home/crg17/Desktop/scripts/4chmodel/cpp; g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./src/convert.cpp -o ./bin/convert_formats.o;


void Vtx2pts(const std::string &path2vtx,const std::string &vtxFile, const std::string &path2pts, const std::string &ptsFile)
{
	std::string outPath;
	std::vector<int> vtx = ReadVtx(path2vtx,vtxFile);
	points points = ReadPts(path2pts,ptsFile);
	std::fstream fOutPts;

	outPath = path2vtx + "/" + vtxFile + ".pts";
	write2file(outPath, fOutPts );

	fOutPts << vtx.size() << std::endl;

	int index;

	for( int i = 0; i < vtx.size(); i++ ) 
    {			

    	index = vtx[i];
		// std::cout << i << " " << vtx[i] << std::endl;

   		fOutPts << points.x[index] << " " << points.y[index] << " " << points.z[index] << std::endl;

    } 

    fOutPts.close();

}

/**
 * @brief Converts a <CODE>.surf</CODE> file into a <CODE>.vtk</CODE> surface file.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.surf</CODE> file.
 * @param mshName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file.
 * @param surfName (<CODE>string</CODE>) name of the <CODE>.surf</CODE> file that will be converted to <CODE>.vtk</CODE>.
 */
void surf2vtk(const std::string &path2pts,const std::string &ptsFile,const std::string &path2surf, const std::string &surfFile)
{	
	std::string outPath;
	/* Read points of the original mesh */
	points meshPts = ReadPts(path2pts,ptsFile);

	/* Read .surf file */
	surf surf = ReadSurf(path2surf,surfFile);

	/* Open .vtk file output */
	std::fstream fOutVtk;
	outPath = path2surf + "/" + surfFile + ".vtk";
	write2file(outPath, fOutVtk );

	/* Header */
    fOutVtk << "# vtk DataFile Version 3.0" << std::endl;
    fOutVtk << "vtk output" << std::endl;
	fOutVtk << "ASCII" << std::endl;
	fOutVtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fOutVtk << " " << std::endl;

	fOutVtk << "POINTS " << meshPts.x.size() << " float" << std::endl;
	for( int i = 0; i < meshPts.x.size(); i++ )
	{
		fOutVtk << meshPts.x[i] << " " << meshPts.y[i] << " " << meshPts.z[i] << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELL_TYPES " << surf.v1.size() << std::endl;
	for( int i = 0; i < surf.v1.size(); i++ )
	{
		fOutVtk << 5 << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELLS " << surf.v1.size() << " " << surf.v1.size()*4 << std::endl;
	for( int i = 0; i < surf.v1.size(); i++ )
	{
		fOutVtk << 3 << " " << surf.v1[i] << " " << surf.v2[i] << " " << surf.v3[i] << std::endl;
	}	

	fOutVtk.close();
}

void neubc2vtk(const std::string &path2pts,const std::string &ptsFile,const std::string &path2neubc, const std::string &neubcFile)
{	
	std::string outPath;
	/* Read points of the original mesh */
	points meshPts = ReadPts(path2pts,ptsFile);

	/* Read .surf file */
	neubc surf = ReadNeubc(path2neubc,neubcFile);

	/* Open .vtk file output */
	std::fstream fOutVtk;
	outPath = path2neubc + "/" + neubcFile + ".vtk";
	write2file(outPath, fOutVtk );

	/* Header */
    fOutVtk << "# vtk DataFile Version 3.0" << std::endl;
    fOutVtk << "vtk output" << std::endl;
	fOutVtk << "ASCII" << std::endl;
	fOutVtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fOutVtk << " " << std::endl;

	fOutVtk << "POINTS " << meshPts.x.size() << " float" << std::endl;
	for( int i = 0; i < meshPts.x.size(); i++ )
	{
		fOutVtk << meshPts.x[i] << " " << meshPts.y[i] << " " << meshPts.z[i] << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELL_TYPES " << surf.v1.size() << std::endl;
	for( int i = 0; i < surf.v1.size(); i++ )
	{
		fOutVtk << 5 << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELLS " << surf.v1.size() << " " << surf.v1.size()*4 << std::endl;
	for( int i = 0; i < surf.v1.size(); i++ )
	{
		fOutVtk << 3 << " " << surf.v1[i] << " " << surf.v2[i] << " " << surf.v3[i] << std::endl;
	}	

	fOutVtk.close();
}

/**
 * @brief  writes in output a new vtk file with the new fiber field giving in input a new fiber field 
 * 
 * @param path (<CODE>string</CODE>) path to the CARP mesh.
 * @param meshName (<CODE>string</CODE>) name of the CARP mesh.
 * @param pathFibre (<CODE>string</CODE>) path to the new fiber field.
 * @param fibreField (<CODE>string</CODE>) name of the <CODE>.lon</CODE> file containing the new fiber field.
 * @param scale (<CODE>double</CODE>) scale factor for the <CODE>.pts</CODE> file.
 */
void mesh2vtk(const std::string &path,const std::string &meshName,const std::string &pathFibre,const std::string &fibreField,const double &scale)
{	
	std::string outPath;
	/* Read .pts file */
	points points  = ReadPts(path,meshName);

	/* Read .elem file */
	elem  elem = ReadElem(path,meshName);

	/* Read .lon file */
	lon lon = ReadLon(pathFibre,fibreField);

	/* Open .vtk file output */
	std::fstream fOutVtk;
	outPath = path + "/" + meshName + ".vtk";
	write2file(outPath, fOutVtk );

	/* Header */
    fOutVtk << "# vtk DataFile Version 3.0" << std::endl;
    fOutVtk << "vtk output" << std::endl;
	fOutVtk << "ASCII" << std::endl;
	fOutVtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fOutVtk << " " << std::endl;

	fOutVtk << "POINTS " << points.x.size() << " float" << std::endl;

	for( int i = 0; i < points.x.size(); i++ )
	{
		fOutVtk << points.x[i]*scale << " " << points.y[i]*scale << " " << points.z[i]*scale << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELL_TYPES " << elem.v1.size() << std::endl;
	for( int i = 0; i < elem.v1.size(); i++ )
	{
		fOutVtk << 10 << std::endl;
	}

	fOutVtk << " " << std::endl;
	fOutVtk << "CELLS " << elem.v1.size() << " " << elem.v1.size()*5 << std::endl;
	for( int i = 0; i < elem.v1.size(); i++ )
	{
		fOutVtk << 4 << " " << elem.v1[i] << " " << elem.v2[i] << " " << elem.v3[i] << " " << elem.v4[i] << std::endl;
	}	

	fOutVtk << " " << std::endl;
	fOutVtk << "CELL_DATA " << elem.v1.size() << std::endl;
	fOutVtk << "SCALARS tags int 1" << std::endl;
	fOutVtk << "LOOKUP_TABLE default" << std::endl;
	for( int i = 0; i < elem.v1.size(); i++ )
	{
		fOutVtk << elem.tag[i] << std::endl;
	}	

	fOutVtk << " " << std::endl;
	fOutVtk << "VECTORS fibres float" << std::endl;
	for( int i = 0; i < lon.x1.size(); i++ )
	{
		fOutVtk << lon.x1[i] << " " << lon.y1[i] << " " << lon.z1[i] << std::endl;
	}	

	fOutVtk.close();
}

void vtx2surf(const std::string &path2msh,const std::string &mshName,const std::string &path2vtx,const std::string &vtxName)
{  

    /* Read input mesh */
    points mshPoints = ReadPts(path2msh,mshName);
    elem mshElem = ReadElem(path2msh,mshName);

    /* Read input vtx */
    std::vector<int> vtxSurf = ReadVtx(path2vtx,vtxName);
    std::sort(vtxSurf.begin(),vtxSurf.end());

    /* Initialize surf elements */
    std::vector<int> v1(mshElem.v1.size()), v2(mshElem.v1.size()), v3(mshElem.v1.size());
    int count = 0;

    #pragma omp parallel for
    for( int i = 0; i < mshElem.v1.size(); i++)
    {   

        std::vector<int> diff;
        std::vector<int> tet = {mshElem.v1[i],mshElem.v2[i],mshElem.v3[i],mshElem.v4[i]};
        std::sort(tet.begin(),tet.end());

        std::set_difference(tet.begin(), tet.end(), vtxSurf.begin(), vtxSurf.end(), std::inserter(diff,diff.begin()));

        if( diff.size() < 2 )
        {

            std::vector<int> tr;
            
            std::set_intersection(tet.begin(),tet.end(),vtxSurf.begin(),vtxSurf.end(),std::back_inserter(tr));
            tr.resize(3);
            diff.resize(1);

            v1[count] = tr[0];
			v2[count] = tr[1]; 
			v3[count] = tr[2];

            count = count + 1;

        }

    }

    v1.resize(count); v2.resize(count); v3.resize(count);

    surf surfOut;
    surfOut.v1 = v1; surfOut.v2 = v2; surfOut.v3 = v3;


    /* Open output .surf file */
    std::fstream fOutSurf;
    std::string outName = path2vtx + "/" + vtxName + ".surf";
    write2file(outName, fOutSurf);
    fOutSurf << count << std::endl;

    for( int i = 0; i < v1.size(); i++ )
    {   
        fOutSurf << surfOut.v1[i] << " " << surfOut.v2[i] << " " << surfOut.v3[i] << std::endl;
    }

    fOutSurf.close();

}

bool str2bool(const std::string &s){
	if(s == "1" || s == "t" || s == "T" || s == "true" || s == "True")
		return true;
	else if(s == "0" || s == "f" || s == "F" || s == "false" || s == "False")
		return false;
}

void elem2vtx(const std::string &path2elem, const std::string &elemName, const std::vector<int> &tags, const std::string &outpath, const std::string &vtxName){

	elem elemfile = ReadElem(path2elem,elemName);
	int n = elemfile.v1.size();
	n *= 4;
	std::vector<int> vtxvec(n,-1);
	
	#pragma omp parallel
	for(int i=0; i< elemfile.v1.size(); i++){
		if( std::find(tags.begin(), tags.end(), elemfile.tag[i]) != tags.end() ){
			vtxvec[i] = elemfile.v1[i];
			vtxvec[i+1] = elemfile.v2[i];
			vtxvec[i+2] = elemfile.v3[i];
			vtxvec[i+3] = elemfile.v4[i];
		}
	}

    std::sort(vtxvec.begin(),vtxvec.end());
    std::vector<int>::iterator ip = std::unique(vtxvec.begin(), vtxvec.end());
    vtxvec.resize(std::distance(vtxvec.begin(), ip));

	vtxvec.erase(std::remove(vtxvec.begin(),vtxvec.end(),-1),vtxvec.end());

	n = vtxvec.size();
	
	std::string fulloutpath = outpath + "/" + vtxName + ".vtx";
	std::fstream fOut(fulloutpath,std::fstream::out);
	fOut << n << "\nintra\n";

	for(int i = 0; i < vtxvec.size(); i++)
		fOut << vtxvec[i] << std::endl;

	fOut.close();


}

int UVC2vtx(const UVC &BiV_UVC, const double& Z, const double& RHO, const double& PHI, const int& V){

	int vertex = -1;
	int close_enough;

	for(int i=0; i < BiV_UVC.Z.size(); i++){
		if(BiV_UVC.V[i] == V){
				if(std::abs(BiV_UVC.RHO[i] - RHO) < 1e-2){
					close_enough = i;
					if(std::abs(BiV_UVC.Z[i] - Z) < 1e-2){
						close_enough = i;
						if(std::abs(BiV_UVC.PHI[i] - PHI) < 1e-2){
						vertex = i;
						break;
					}
				}
			}
		}
	}

	if(vertex == -1)
		vertex = close_enough;
	
	return vertex;

}

int UVC2vtx(const std::string& path2UVC, const std::string& UVCname, const double& Z, const double& RHO, const double& PHI, const int& V){

	int vertex;

	UVC BiV_UVC = ReadUVC(path2UVC,UVCname);
	
	vertex = UVC2vtx(BiV_UVC,Z,RHO,PHI,V);
	
	return vertex;

}

void vtx2vtk(const std::string &path2vtx,const std::string &vtxFile, const std::string &path2pts, const std::string &ptsFile)
{
	std::string outPath;
	std::vector<int> vtx = ReadVtx(path2vtx,vtxFile);
	points points = ReadPts(path2pts,ptsFile);
	std::fstream fOutPts;

	outPath = path2vtx + "/" + vtxFile + ".vtk";
	write2file(outPath, fOutPts );

	fOutPts << "# vtk DataFile Version 3.0 \nvtk output \nASCII \nDATASET UNSTRUCTURED_GRID \n\nPOINTS " << vtx.size() << " float\n";

	int index;

	for( int i = 0; i < vtx.size(); i++ ) 
    {			

    	index = vtx[i];
   		fOutPts << points.x[index] << " " << points.y[index] << " " << points.z[index] << std::endl;

    } 

    fOutPts.close();

}

/*
int main(int argc,char* argv[]){
	
	std::string ifmt, ofmt, path2vtx, vtxFile, path2pts, ptsFile, path2surf, surfFile, path2neubc, 
		    neubcFile, path2elem, elemName, outpath, vtxname;

	std::cout << "\nSelect the file format to convert from (vtx, surf, neubc, elem):\n";
	std::cin >> ifmt;

	
	if(ifmt == "elem")
		ofmt = "vtx";
	else
		ofmt = "vtk";

	if(ifmt == "vtx" && ofmt == "vtk"){
		std::cout << "\nType the path to the vtx file:\n";
		std::cin >> path2vtx;

		std::cout << "\nType the name of the vtx file:\n";
		std::cin >> vtxFile;

		std::cout << "\nType the path of the pts file where that vtx comes from:\n";
		std::cin >> path2pts;

		std::cout << "\nType the name of the pts file where that vtx comes from:\n";
		std::cin >> ptsFile;

		vtx2vtk(path2vtx, vtxFile, path2pts, ptsFile);

		std::cout << "\nThe output file is " << path2vtx << "/" << vtxFile << ".vtk\n";
		
	}

	else if(ifmt == "surf" && ofmt == "vtk"){
		std::cout << "\nType the path to the surf file:\n";
		std::cin >> path2surf;

		std::cout << "\nType the name of the surf file:\n";
		std::cin >> surfFile;

		std::cout << "\nType the path of the pts file where that surf comes from:\n";
		std::cin >> path2pts;

		std::cout << "\nType the name of the pts file where that surf comes from:\n";
		std::cin >> ptsFile;

		surf2vtk(path2pts, ptsFile, path2surf, surfFile);

		std::cout << "\nThe output file is " << path2surf << "/" << surfFile << ".vtk\n";
		
	}

	else if(ifmt == "neubc" && ofmt == "vtk"){
		std::cout << "\nType the path to the neubc file:\n";
		std::cin >> path2neubc;

		std::cout << "\nType the name of the neubc file:\n";
		std::cin >> neubcFile;

		std::cout << "\nType the path of the pts file where that neubc comes from:\n";
		std::cin >> path2pts;

		std::cout << "\nType the name of the pts file where that neubc comes from:\n";
		std::cin >> ptsFile;

		neubc2vtk(path2pts, ptsFile, path2neubc, neubcFile);

		std::cout << "\nThe output file is " << path2neubc << "/" << neubcFile << ".vtk\n";
		
	}

	else if(ifmt == "elem" && ofmt == "vtx"){
		std::cout << "\nType the path to the elem file:\n";
		std::cin >> path2elem;

		std::cout << "\nType the name of the elem name file:\n";
		std::cin >> elemName;

		std::cout << "\nType the output path:\n";
		std::cin >> outpath;

		std::cout << "\n Type the output name:\n";
		std::cin >> vtxname;


		elem2vtx(path2elem, elemName, {1,2,25,26,27,28,29},outpath,vtxname);

		std::cout << "\n The output file is " << outpath << "/" << vtxname << ".vtx\n";
	}

	else
		std::cout << "\nCombination of input-output formats not recognised or implemented.";
	
	
	return 0;
}*/

