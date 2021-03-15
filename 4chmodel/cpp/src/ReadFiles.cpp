/**
 * @brief Routines to deal with the input of files.
 * 
 * @file ReadFiles.cpp
 * @author Marina Strocchi
 * @author Cristobal Rodero
 * @date 2018-08-31
 */


#include "../headers/ReadFiles.h"

/**
 * @brief Reads a <CODE>.vtx</CODE> file.
 * 
 *  The function reads in a <CODE> .vtx </CODE> file expecting two header lines containing the number of vertices to read and
    a string specifying if the nodes are in the intracellular space ("intra") or in the extracellular space
    ("extra").
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE> .vtx </CODE> file to read.
 * @param filename (<CODE>string</CODE>) name of the <CODE> .vtx </CODE> to read.
 * @return <CODE>std::vector<int></CODE> Vector with the indices.
 */

std::vector<int> ReadVtx(const std::string & path,const std::string & filename)
{

	std::cout << "Reading " << filename << ".vtx file...\n";

	std::ifstream inputFile(path + "/" + filename + ".vtx");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path << "/" << filename + ".vtx" << std::endl;
		exit(0);
	}

	std::string line;	

	/* Read first line with the number of nodes */
	
	getline(inputFile,line);
	std::istringstream iss1(line);	
	int nPoints;	
	iss1 >> nPoints;

	getline(inputFile,line);
	std::istringstream iss2(line);	
	std::string intra;
	iss2 >> intra; 

	std::vector<int> vtx(nPoints);	

	int i = 0;

	while(getline(inputFile,line))
	{

		std::istringstream iss(line);
		iss >> vtx[i];
		i = i+1;

	}

	return vtx;

}

/**
 * @brief Reads a <CODE>.pts</CODE> file.
 * @param path (<CODE>string</CODE>) path to the <CODE>.pts</CODE> file to read.
 * @param filename (<CODE>string</CODE>) name of the <CODE>.pts</CODE> to read
 * @return <CODE>points</CODE> Structure containing the list of points.
 */

points ReadPts(const std::string &path,const std::string &filename)
{	

	points points;

	std::cout << "Reading " << filename << ".pts file...\n";

	std::ifstream inputFile(path + "/" + filename + ".pts");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << filename + ".pts" << std::endl;
		exit(0);
	}

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	int nPoints;

	iss >> nPoints;

	std::vector<double> x(nPoints);
	std::vector<double> y(nPoints);
	std::vector<double> z(nPoints);

	int i = 0;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		     

        lineStream >> std::setprecision(1) >> x[i] >> y[i] >> z[i];

        i = i+1;

    }
    points.x = x;
    points.y = y;
    points.z = z;

	inputFile.close();

	return points;

}

/**
 * @brief Reads a <CODE>.elem</CODE> file.
 * 
 * The function reads in a <CODE>.elem</CODE> file expecting one header line containing the number of surface triangles
    to read.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.elem</CODE> file to read.
 * @param filename fileName (<CODE>string</CODE>) name of the <CODE>.elem</CODE> to read.
 * @return <CODE>surf</CODE> Structure containing the surface. 
 */

surf ReadSurfElem(const std::string &path,const std::string &filename)
{	

	surf elem;

	std::cout << "Reading " << filename << ".elem file...\n";

	std::ifstream inputFile(path + "/" + filename + ".elem");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << filename + ".elem" << std::endl;
		exit(0);
	}

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	int nElem;

	iss >> nElem;

	std::vector<int> v1(nElem);
	std::vector<int> v2(nElem);
	std::vector<int> v3(nElem);

	int i = 0;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

		std::string Tr;

        lineStream >> Tr >> v1[i] >> v2[i] >> v3[i];

        i = i+1;

    }

    elem.v1 = v1;
    elem.v2 = v2;
    elem.v3 = v3;

	inputFile.close();

	return elem;

}
/**
 * @brief Reads a <CODE>.surf</CODE> file.
 * 
 * The function reads in a <CODE>.surf</CODE> file expecting one header line containing the number of triangles
    to read.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.surf</CODE> file to read.
 * @param filename fileName (<CODE>string</CODE>) name of the <CODE>.surf</CODE> to read.
 * @return <CODE>surf</CODE> Structure containing the surface. 
 */

surf ReadSurf(const std::string &path,const std::string &filename)
{	

	surf elem;

	std::cout << "Reading " << filename << ".surf file...\n";

	std::ifstream inputFile(path + "/" + filename + ".surf");

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path << "/" << filename + ".surf" << std::endl;
		exit(0);
	}	

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	int nElem;

	iss >> nElem;

	std::vector<int> v1(nElem);
	std::vector<int> v2(nElem);
	std::vector<int> v3(nElem);

	int i = 0;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

		std::string Tr;

        lineStream >> Tr >> v1[i] >> v2[i] >> v3[i];

        i = i+1;

    }

    elem.v1 = v1;
    elem.v2 = v2;
    elem.v3 = v3;

	inputFile.close();

	return elem;

}
/**
 * @brief Reads a <CODE>.elem</CODE> file.
 * 
 * The function reads in a <CODE>.elem</CODE> file expecting one header line containing the number of tetrahedra
    to read.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.elem</CODE> file to read.
 * @param filename fileName (<CODE>string</CODE>) name of the <CODE>.elem</CODE> to read.
 * @return <CODE>elem</CODE> Structure containing the elements. 
*/

elem ReadElem(const std::string &path1,const std::string &filename,bool outflag)
{	

	elem elem;
	std::string path = path1;

	if(outflag)
	std::cout << "Reading " << path + "/" + filename << ".elem file...\n";

	std::ifstream inputFile(path + "/" + filename + ".elem");
	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	iss >> elem.nElem;

	std::vector<int> v1(elem.nElem);
	std::vector<int> v2(elem.nElem);
	std::vector<int> v3(elem.nElem);
	std::vector<int> v4(elem.nElem);
	std::vector<int> tag(elem.nElem);

	int i = 0;
    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

		std::string Tr;

        lineStream >> Tr >> v1[i] >> v2[i] >> v3[i] >> v4[i] >> tag[i];

        i = i+1;

    }
    elem.v1 = v1;
    elem.v2 = v2;
    elem.v3 = v3;
    elem.v4 = v4;
    elem.tag = tag;

	inputFile.close();

	return elem;

}

/**
 * @brief Reads a <CODE>.lon</CODE> file.
 * 
 * The function reads in a <CODE>.lon</CODE> file expecting one header line containing the number of of directions given 
	for each element
 * @param path (<CODE>string</CODE>) path to the <CODE>.lon</CODE> file to read.
 * @param filename (<CODE>string</CODE>) name of the </CODE>.lon</CODE> to read.
 * @return <CODE>lon</CODE> structure with the list of fibres.
 * @attention Both fiber and sheet directions are given.
 */

lon ReadLon(const std::string &path,const std::string &filename)
{	

	lon lon;

	std::cout << "Reading " << filename << ".lon file...\n";

	std::ifstream inputFile(path + "/" + filename + ".lon");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << filename + ".lon" << std::endl;
		exit(0);
	}

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	iss >> lon.nDir;

	std::vector<double> x1;
	std::vector<double> y1;
	std::vector<double> z1;
	std::vector<double> x2;
	std::vector<double> y2;
	std::vector<double> z2;

	double x1Temp, y1Temp, z1Temp, x2Temp, y2Temp, z2Temp;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

        lineStream >> x1Temp >> y1Temp >> z1Temp >> x2Temp >> y2Temp >> z2Temp;

        x1.push_back(x1Temp);
        y1.push_back(y1Temp);
        z1.push_back(z1Temp);
        x2.push_back(x2Temp);
        y2.push_back(y2Temp);
        z2.push_back(z2Temp);

    }


    lon.x1 = x1;
    lon.y1 = y1;
    lon.z1 = z1;
    lon.x2 = x2;
    lon.y2 = y2;
    lon.z2 = z2;

	inputFile.close();

	return lon;

}
/**
 * @brief Reads a <CODE>.neubc</CODE> file.
 * 
 * The function reads in a <CODE>.neubc</CODE> file expecting one header line of the following format:<BR>
 * <VAR>nTriangles # mesh elast nElements nNodes </VAR> :<BR>
 * where <VAR>nTriangles</VAR>, <VAR>nElements</VAR> and <VAR>nNodes</VAR> are the number of surface elements, the number of tetrahedra in the 
	initial three-dimensional mesh and the number of nodes in the initial three-dimensional mesh.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.neubc</CODE> file to read.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.neubc</CODE> tor read.
 * @return <CODE>neubc</CODE> Structur with the <CODE>.neubc</CODE> file content.
 */
neubc ReadNeubc(const std::string &path,const std::string &filename)
{	

	neubc neubcOut;

	std::cout << "Reading " << filename << ".neubc file..\n.";

	std::ifstream inputFile(path + "/" + filename + ".neubc");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << filename + ".neubc" << std::endl;
		exit(0);
	}	

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	int nElem;
	int nElemMsh;
	int nNodesMsh;
	std::string Temp;

	iss >> nElem >> Temp >> Temp >> Temp >> nElemMsh >> nNodesMsh;

	std::vector<int> v1(nElem);
	std::vector<int> v2(nElem);
	std::vector<int> v3(nElem);
	std::vector<int> v4(nElem);
	std::vector<int> tet(nElem);

	int i = 0;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

		int Temp;

        lineStream >> v1[i] >> v2[i] >> v3[i] >> Temp >> v4[i] >> tet[i];

        i = i+1;

    }

    neubcOut.nElem = nElemMsh;
    neubcOut.nNodes = nNodesMsh;
    neubcOut.v1 = v1;
    neubcOut.v2 = v2;
    neubcOut.v3 = v3;
    neubcOut.v4 = v4;
    neubcOut.tet = tet;

	inputFile.close();

	return neubcOut;

}
/**
 * @brief Reads a <CODE>.dat</CODE> file.
 * 
 * The function reads in a <CODE>.dat</CODE> file expecting no header line.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.dat</CODE> file to read.
 * @param filename (<CODE>string</CODE>) name of the <CODE>.dat</CODE> to read.
 * @return <CODE>std::vector<double></CODE> structure containing the <CODE>.dat</CODE> file content.
 */

std::vector<double> ReadActTime(const std::string &path,const std::string &filename, bool outflag)
{
	
	if(outflag)
	std::cout << "Reading " << path << "/" << filename << ".dat file...\n";

	std::ifstream inputFile(path + "/" + filename + ".dat");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path + "/" + filename + ".dat" << std::endl;
		exit(0);
	}	

	std::string line;	

	double tAct;
	std::vector<double> actTime;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

        lineStream >> tAct;

        actTime.push_back(tAct);

    }

	// std::cout << "\n File read successfully.\n";
	return actTime;
}
std::vector<double> ReadVolFile(const std::string &path,const std::string &filename )
{

	std::cout << "Reading " << filename << ".dat file...\n";

	std::ifstream inputFile(path + "/" + filename + ".dat");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path + "/" + filename + ".dat" << std::endl;
		exit(0);
	}	

	std::string line;	

	double tAct;
	std::vector<double> actTime;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		

        lineStream >> tAct;

        actTime.push_back(tAct*1e-9);

    }


	return actTime;
}
/**
 * @brief Merges <CODE>.vtx</CODE> files.
 * 
 * The function reads in <CODE>.vtx</CODE> files, merges them and writes the output <CODE>.vtx</CODE> file to the given path.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.vtx</CODE> files to merge.
 * @param fileNames (<CODE>vector<string></CODE>) vector containing the names of the <CODE>.vtx</CODE> files to merge.
 * @param pathOut (<CODE>string</CODE>) path to the output file.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file.
 */
UVC ReadUVC(const std::string &path,const std::string &filename)
{	

	UVC UVC_file;

	std::cout << "Reading " << filename << ".dat file...\n";

	std::ifstream inputFile(path + "/" + filename + ".pts");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path << "/" + filename + ".dat" << std::endl;
		exit(0);
	}

	std::string line;	

	/* Read header line */

	std::getline(inputFile, line);

	std::istringstream iss(line);	

	int nPoints;

	iss >> nPoints;

	std::vector<double> Z(nPoints);
	std::vector<double> RHO(nPoints);
	std::vector<double> PHI(nPoints);
	std::vector<double> V(nPoints);

	int i = 0;

    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);		     

        lineStream >> std::setprecision(6) >> Z[i] >> RHO[i] >> PHI[i] >> std::setprecision(1) >> V[i];

        i = i+1;

    }
    
	UVC_file.Z = Z;
	UVC_file.RHO = RHO;
	UVC_file.PHI = PHI;
	UVC_file.V = V;
	
	inputFile.close();

	return UVC_file;
}

EP_metrics_table ReadEP_metrics_table(const std::string &path2table, const std::string &tableName){

	std::ifstream inputFile(path2table + "/" + tableName + ".txt");

	EP_metrics_table tab;
	std::string line, header;
	std::vector<double> TAT, TAT_LV, AT_10_90;
	std::vector<std::string> HF, heart, lead_position, comb_elec;
	double aux1, aux2, aux3;
	std::string aux4, aux5, aux6, aux7;

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << path2table + "/" + tableName + ".txt" << std::endl;
		exit(0);
	}

	/* Read header line */

	std::getline(inputFile, line);
	std::istringstream iss(line);	

	iss >> header;

    while(std::getline(inputFile, line)){
		std::istringstream lineStream(line);		     
        lineStream >> aux1 >> aux2 >> aux3 >> aux4 >> aux5 >> aux6 >> aux7;
        TAT.push_back(aux1);
		TAT_LV.push_back(aux2);
		AT_10_90.push_back(aux3);
		HF.push_back(aux4);
		heart.push_back(aux5);
		lead_position.push_back(aux6);
		comb_elec.push_back(aux7);
	}
	inputFile.close();

	tab.AT_10_90 = AT_10_90;
	tab.comb_elec = comb_elec;
	tab.heart = heart;
	tab.HF = HF;
	tab.lead_position = lead_position;
	tab.TAT = TAT;
	tab.TAT_LV = TAT_LV;

	// for(int j = 0; j < tab.TAT.size(); j++)
	// 	std::cout << tab.TAT[j] << std::endl << tab.TAT_LV[j] << std::endl << tab.AT_10_90[j] << std::endl << tab.HF[j] << std::endl << tab.heart[j] << std::endl << tab.lead_position[j] << std::endl << tab.comb_elec[j] << std::endl;

	return tab;
}

void mergeVtx(const std::string &path,const std::vector<std::string> &fileNames,const std::string &pathOut,const std::string &fileNameOut )
{

	std::vector<int> vtxOut;

	for(int i = 0; i < fileNames.size(); i++)
	{
		std::vector<int> vtxTemp = ReadVtx(path,fileNames.at(i));

		for(int j = 0; j < vtxTemp.size(); j++)
		{
			vtxOut.push_back(vtxTemp[j]);
		}
	}

	std::fstream fOut(pathOut + "/" + fileNameOut + ".vtx",std::fstream::out);
	fOut << vtxOut.size() << std::endl << "intra" << std::endl;

	for(int i = 0; i < vtxOut.size(); i++)
	{
		fOut << vtxOut[i] << std::endl;
	}

	fOut.close();

}
/**
 * @brief Merges <CODE>.nebc</CODE> files.
 * 
 * The function reads in <CODE>.neubc</CODE> files, merges them and writes the output <CODE>.neubc</CODE> file to the given path.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.neubc</CODE> files to merge.
 * @param fileNames (<CODE>vector<string></CODE>) vector containing the names of the <CODE>.neubc</CODE> files to merge.
 * @param pathOut (<CODE>string</CODE>) path to the output file.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file.
 */
void mergeNeubc(const std::string &path,const std::vector<std::string> &fileNames,const std::string &pathOut,const std::string &fileNameOut )
{

	neubc neubcOut;

	int nElem;
	int nNodes;

	for(int i = 0; i < fileNames.size(); i++)
	{
		neubc neubcTemp = ReadNeubc(path,fileNames.at(i));

		nElem = neubcTemp.nElem;
		nNodes = neubcTemp.nNodes;

		for(int j = 0; j < neubcTemp.v1.size(); j++)
		{
			neubcOut.v1.push_back(neubcTemp.v1[j]);
			neubcOut.v2.push_back(neubcTemp.v2[j]);
			neubcOut.v3.push_back(neubcTemp.v3[j]);
			neubcOut.v4.push_back(neubcTemp.v4[j]);
			neubcOut.tet.push_back(neubcTemp.tet[j]);
		}
	}

	std::fstream fOut(pathOut + "/" + fileNameOut + ".neubc",std::fstream::out);
	fOut << neubcOut.v1.size() << " # " << "mesh elast " << nElem << " " << nNodes << " :" << std::endl;

	for(int i = 0; i < neubcOut.v1.size(); i++)
	{
		fOut << neubcOut.v1[i] << " " << neubcOut.v2[i] << " " << neubcOut.v3[i] << " " << 1 << " " << neubcOut.v4[i] << " " << neubcOut.tet[i] << std::endl;
	}

	fOut.close();

}

/**
 * @brief Merges <CODE>.surf</CODE> files.
 * 
 * The function reads in <CODE>.surf</CODE> files, merges them and writes the output <CODE>.surf</CODE> file to the given path.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.surf</CODE> files to merge.
 * @param fileNames (<CODE>vector<string></CODE>) vector containing the names of the <CODE>.surf</CODE> files to merge.
 * @param pathOut (<CODE>string</CODE>) path to the output file.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file.
 */
void mergeSurf(const std::string &path,const std::vector<std::string> &fileNames,const std::string &pathOut,const std::string &fileNameOut )
{

	surf surfOut;

	for(int i = 0; i < fileNames.size(); i++)
	{
		surf surfTemp = ReadSurf(path,fileNames.at(i));

		for(int j = 0; j < surfTemp.v1.size(); j++)
		{
			surfOut.v1.push_back(surfTemp.v1[j]);
			surfOut.v2.push_back(surfTemp.v2[j]);
			surfOut.v3.push_back(surfTemp.v3[j]);
		}
	}

	std::fstream fOut(pathOut + "/" + fileNameOut + ".surf",std::fstream::out);
	fOut << surfOut.v1.size() << std::endl;

	for(int i = 0; i < surfOut.v1.size(); i++)
	{
		fOut << "Tr " << surfOut.v1[i] << " " << surfOut.v2[i] << " " << surfOut.v3[i] << std::endl;
	}

	fOut.close();

}

/**
 * @brief Merges two fibres files.
 * 
 * The function reads in two <CODE>.lon</CODE> files, merges them and writes the output <CODE>.lon</CODE> file to the given path.
 * For the elements out of the fibres tags, it sets default fibres.
 * 
 * @param path (<CODE>string</CODE>) path to the original mesh.
 * @param mshName (CODE>string</CODE>) name of the original mesh.
 * @param path2Lon1 (<CODE>string</CODE>) path to the first <CODE>.lon</CODE> file. 
 * @param nameLon1 (<CODE>string</CODE>) name of the first <CODE>.lon</CODE> file.
 * @param tags1 (<CODE>std::vector<int></CODE>) Element tags corresponding to the first set of fibres.
 * @param path2Lon2  path2Lon2 (<CODE>string</CODE>) path to the second <CODE>.lon</CODE> file.
 * @param nameLon2 (<CODE>string</CODE>) name of the second <CODE>.lon</CODE> file.
 * @param tags2 (<CODE>std::vector<int></CODE>) Element tags corresponding to the second set of fibres.
 */

void mergeLon(const std::string &path,const std::string &mshName,const std::string &path2Lon1,const std::string &nameLon1,const std::vector<int> &tags1,const std::string &path2Lon2,const std::string &nameLon2,const std::vector<int> &tags2 )
{

	/* Read .lon file */
    lon lon1 = ReadLon(path2Lon1,nameLon1);
    lon lon2 = ReadLon(path2Lon2,nameLon2);
    lon lonOut;

    /* Read heart.elem for element tag */
    elem mshElem = ReadElem(path,mshName);

    std::vector<double> x1(mshElem.v1.size());
    std::vector<double> y1(mshElem.v1.size());
    std::vector<double> z1(mshElem.v1.size());
    std::vector<double> x2(mshElem.v1.size());
    std::vector<double> y2(mshElem.v1.size());
    std::vector<double> z2(mshElem.v1.size()); 

    int k1 = 0;
    int k2 = 0;  

    for(int i = 0; i < mshElem.v1.size(); i++)
    {
    	if( std::find(tags1.begin(), tags1.end(), mshElem.tag[i]) != tags1.end() )
    	{
    		x1[i] = lon1.x1[k1];
    		y1[i] = lon1.y1[k1];
    		z1[i] = lon1.z1[k1];
    		x2[i] = lon1.x2[k1];
    		y2[i] = lon1.y2[k1];
    		z2[i] = lon1.z2[k1];    	
    		k1 = k1+1;	
    	} 
    	else if( std::find(tags2.begin(), tags2.end(), mshElem.tag[i]) != tags2.end() )
    	{
    		x1[i] = lon2.x1[k2];
    		y1[i] = lon2.y1[k2];
    		z1[i] = lon2.z1[k2];
    		x2[i] = lon2.x2[k2];
    		y2[i] = lon2.y2[k2];
    		z2[i] = lon2.z2[k2];      		
    		k2 = k2+1;
    	}
    	else
    	{
    		x1[i] = 1.0;
    		y1[i] = 0.0;
    		z1[i] = 0.0;
    		x2[i] = 0.0;
    		y2[i] = 1.0;
    		z2[i] = 0.0;   
    	}
    }

    lonOut.x1 = x1;
    lonOut.y1 = y1;
    lonOut.z1 = z1;
    lonOut.x2 = x2;
    lonOut.y2 = y2;
    lonOut.z2 = z2;

    /* Write output */ 
    std::fstream fOutFibres(path + "/" + mshName + ".lon",std::fstream::out);   
    fOutFibres << 2 << std::endl;

    for( int i = 0; i < lonOut.x1.size(); i++ )
    {
        fOutFibres << lonOut.x1[i] << " " << lonOut.y1[i] << " " << lonOut.z1[i] << " " << lonOut.x2[i] << " " << lonOut.y2[i] << " " << lonOut.z2[i] << std::endl;
    }

    fOutFibres.close();

}

/**
 * @brief Creates a default <CODE>.lon</CODE> file for a given mesh.
 * 
 * @param path (<CODE>string</CODE>) path to the mesh.
 * @param mshName (<CODE>string</CODE>) name of the mesh.
 */
/**
 * @brief Merges <CODE>.elem</CODE> files.
 * 
 * The function reads in <CODE>.elem</CODE> files, merges them and writes the output <CODE>.elem</CODE> file to the given path.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.elem</CODE> files to merge.
 * @param fileNames (<CODE>vector<string></CODE>) vector containing the names of the <CODE>.elem</CODE> files to merge.
 * @param pathOut (<CODE>string</CODE>) path to the output file.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file.
 */
void mergeElem(const std::string &path,const std::vector<std::string> &fileNames,const std::string &pathOut,const std::string &fileNameOut )
{

	elem elemOut;

	for(int i = 0; i < fileNames.size(); i++)
	{
		elem elemTemp = ReadElem(path,fileNames.at(i));

		for(int j = 0; j < elemTemp.v1.size(); j++)
		{
			elemOut.v1.push_back(elemTemp.v1[j]);
			elemOut.v2.push_back(elemTemp.v2[j]);
			elemOut.v3.push_back(elemTemp.v3[j]);
			elemOut.tag.push_back(elemTemp.tag[j]);
		}
	}

	std::fstream fOut(pathOut + "/" + fileNameOut + ".elem",std::fstream::out);
	fOut << elemOut.v1.size() << std::endl;

	for(int i = 0; i < elemOut.v1.size(); i++)
		fOut << "Tt " << elemOut.v1[i] << " " << elemOut.v2[i] << " " << elemOut.v3[i] << " " << elemOut.tag[i] << std::endl;

	fOut.close();

}




/**
 * @brief Merges <CODE>.pts</CODE> files.
 * 
 * The function reads in <CODE>.pts</CODE> files, merges them and writes the output <CODE>.pts</CODE> file to the given path.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.pts</CODE> files to merge.
 * @param fileNames (<CODE>vector<string></CODE>) vector containing the names of the <CODE>.pts</CODE> files to merge.
 * @param pathOut (<CODE>string</CODE>) path to the output file.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file.
 */
void mergePts(const std::string &path,const std::vector<std::string> &fileNames,const std::string &pathOut,const std::string &fileNameOut )
{

	points ptsOut;

	for(int i = 0; i < fileNames.size(); i++)
	{
		points ptsTemp = ReadPts(path,fileNames.at(i));

		for(int j = 0; j < ptsTemp.x.size(); j++)
		{
			ptsOut.x.push_back(ptsTemp.x[j]);
			ptsOut.y.push_back(ptsTemp.y[j]);
			ptsOut.z.push_back(ptsTemp.z[j]);
		}
	}

	std::fstream fOut(pathOut + "/" + fileNameOut + ".pts",std::fstream::out);
	fOut << ptsOut.x.size() << std::endl;

	for(int i = 0; i < ptsOut.x.size(); i++)
		fOut << ptsOut.x[i] << " " << ptsOut.y[i] << " " << ptsOut.z[i] << std::endl;

	fOut.close();

}

/**
 * @brief Writes a file if it doesn't exist.
 * 
 *  The function checks if the given file exists. If it does, it gives the option to proceed and rewrite on it
 *  or not. If it does not exist, the file is opened 
 * 
 * @param fileName (<CODE>string</CODE>) name of the output file.
 * @param file (<CODE>std::fstream</CODE>) file to write to.
 * 
 * @attention Uncoment to get the message of overwriting warning.
 */

void defaultLon(const std::string &path,const std::string &mshName)
{

    /* Read elem for element tag */
    elem mshElem = ReadElem(path,mshName);

    /* Write output */ 
    std::fstream fOutFibres(path + "/" + mshName + "_default.lon",std::fstream::out);   
    fOutFibres << 2 << std::endl;

    for( int i = 0; i < mshElem.v1.size(); i++ )
    {
        fOutFibres << 0.0 << " " << 0.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << 0.0 << std::endl;
    }

    fOutFibres.close();

}

/**
 * @brief Function to write a file. It looks so simple because it's a modification from an earlier version so to be backwards compatible.
 * 
 * @param fileName Path+filename of the file to write
 * @param file fstream where to write the file.
 */
void write2file(std::string &fileName, std::fstream &file )
{	
	file.open(fileName,std::ios::out);
}

/**
 * @brief Extracts the <VAR>i</VAR>-th point from a <CODE>.pts</CODE> file.
 * 
 * @param p (<CODE>points</CODE>) Structure with the <CODE>.pts</CODE> file content.
 * @param i (<CODE>int</CODE>) Component to extract from <CODE>points</CODE> structure.
 * @return std::vector<double> 
 */

std::vector<double> extract_point(const points &p, const int &i){
	
	std::vector<double> v(3); /**< Results vector */

	v.at(0) = (p.x).at(i);
	v.at(1) = (p.y).at(i);
	v.at(2) = (p.z).at(i);

	return v;
}
/**
 * @brief Changes the numbering of the specified tags.
 * 
 * @param path (<CODE>string</CODE>) Path to the <CODE>.elem</CODE> file to change the tags.
 * @param elem_file (<CODE>string</CODE>) Elements file to change the tags.
 * @param tags_old (<CODE>std::vector<int></CODE>) Set of tags to change.
 * @param tags_new (<CODE>std::vector<int></CODE>) New tags to get set.
 * 
 * @attention <CODE>tags_old</CODE> and <CODE>tags_new</CODE> must have the same length.
 * @attention The change is done <B>sequentially</B>, so exchange of tags is allowed.
 */
void change_tags(const std::string &path, const std::string &elem_file, const std::vector<int> &tags_old, const std::vector<int> &tags_new){

	elem el_old = ReadElem(path,elem_file);
	elem el_new = el_old;

	#pragma omp parallel for
	for(long int i = 0; i < el_old.v1.size(); i++){
		for(int t = 0; t < tags_old.size(); t++){
			if(el_old.tag[i] == tags_old.at(t)){
				el_new.tag[i] = tags_new.at(t);}
		}
	}

	std::fstream fOut(path + "/" + elem_file + "_retagged.elem",std::fstream::out);
	fOut << el_new.v1.size() << std::endl;

	for(long int i = 0; i < el_new.v1.size(); i++){
		fOut << "Tt " << el_new.v1[i] << " " << el_new.v2[i] << " " << el_new.v3[i] << " " << el_new.v4[i] << " " << el_new.tag[i] << std::endl;
	}
	fOut.close();
	
}
/**
 * @brief Obtains the edges of a mesh as pairs of the connected nodes.
 * 
 * @bug NOT FUNCTIONAL. Sparses matrices give problems when parallelized (boost and eigen) and without sparse the sizes are way too big.
 * 
 * @param path2elem 
 * @param elemName 
 * @param path2pts 
 * @param ptsName 
 * @param wanna_write 
 * @return edges 
 */
edges getEdges(const std::string &path2elem, const std::string &elemName, const std::string &path2pts, const std::string &ptsName, const bool &wanna_write){

	// using namespace boost::numeric::ublas; // For the sparse matrix

	elem el_file = ReadElem(path2elem,elemName);
	points pts_file = ReadPts(path2pts,ptsName);

	int n_cols = pts_file.x.size();
	int n_rows = n_cols;

	long int size_tril = n_cols*(n_cols-1)/2;
	std::cout << size_tril;
	std::vector<int> n1(5,0);
/*	std::vector<int> n2(size_tril,0);
	// std::vector<int> nnz_per_col(n_cols);
	
	std::cout << "cool\n";
	// mapped_matrix<bool> adj_matrix(pts_file.x.size(), pts_file.x.size());
	std::vector<std::vector<bool>> adj_matrix(n_rows,std::vector<bool>(n_cols));
*/
	edges edge_struct;
/*	
	
	std::cout << "Creating adjacency matrix...\n";
	// Create a sparse adjacency matrix

	#pragma omp parallel
	for(long int e = 0; e < el_file.nElem; e++){
		// if(100*double(e)/double(el_file.nElem) > cout_int){
		// 	std::cout << cout_int << "\%\n";
		// 	cout_int++; 
		// }
		adj_matrix[(el_file.v1)[e]][(el_file.v2)[e]] = 1;
		adj_matrix[(el_file.v1)[e]][(el_file.v3)[e]] = 1;
		adj_matrix[(el_file.v1)[e]][(el_file.v4)[e]] = 1;
		adj_matrix[(el_file.v2)[e]][(el_file.v3)[e]] = 1;
		adj_matrix[(el_file.v2)[e]][(el_file.v4)[e]] = 1;
		adj_matrix[(el_file.v3)[e]][(el_file.v4)[e]] = 1;
		adj_matrix[(el_file.v2)[e]][(el_file.v1)[e]] = 1;
		adj_matrix[(el_file.v3)[e]][(el_file.v1)[e]] = 1;
		adj_matrix[(el_file.v4)[e]][(el_file.v1)[e]] = 1;
		adj_matrix[(el_file.v3)[e]][(el_file.v2)[e]] = 1;
		adj_matrix[(el_file.v4)[e]][(el_file.v2)[e]] = 1;
		adj_matrix[(el_file.v4)[e]][(el_file.v3)[e]] = 1;
	}
	
	std::cout<< "Creating edges....\n";
	// Create the edge structure


	#pragma omp parallel
	for(long int i = 1; i < n_rows; i++){
		#pragma omp parallel
		for(long int j = 0; j < i; j++){
			n1[(i*(i-1))/2+j] = adj_matrix[i][j]*i;
			n2[(i*(i-1))/2+j] = adj_matrix[i][j]*j;
		}
	}
	std::cout<< "Done.\n";
	edge_struct.n1 = n1;
	edge_struct.n2 = n2;
	std::cout << "Writing...\n";
	if(wanna_write){
		std::fstream fOut(path2elem + "/" + elemName + ".edges", std::fstream::out);
		for(long int i = 0; i < n1.size(); i++)
			fOut << n1[i] << " " << n2[i] << " " << std::endl;
		fOut.close();
	}
*/
	
	return edge_struct;
}

	bool compare(double a, double b, std::vector<double> data){
    	return data[int(a)]<data[int(b)];
	}



std::string GotoLine(const std::string &file, const int &num){

	std::ifstream myfile(file);
	std::string line, result;
	
	for(int lineno = 0; getline(myfile,line) && lineno < num; lineno++)
		if(lineno == num - 1)
			result = line;

	return result;
	
}




void correctmeshtoolfuckup(const std::string &path2elemnoFEC, const std::string &nameelemnoFEC, const std::string &path2elemFEC, const std::string &nameelemFEC, const std::vector<int> &FEC_tags, const std::string &path2output, const std::string &nameoutput){

	elem elemnoFEC = ReadElem(path2elemnoFEC,nameelemnoFEC);
	elem elemFEC = ReadElem(path2elemFEC,nameelemFEC);
	elem elemcorrected = elemFEC;

	// Go through elements
	// if the tag is not in the FEC_tags
	// if they are different 
	// keep the old

	#pragma omp parallel
	for(int i = 0; i < elemcorrected.tag.size(); i++)
		if(elemnoFEC.tag[i] != elemFEC.tag[i]) // If tags are different...
			if( std::find(FEC_tags.begin(), FEC_tags.end(), elemFEC.tag[i]) == FEC_tags.end() ) // ..but they are not the FEC tags
				elemcorrected.tag[i] = elemnoFEC.tag[i];

	std::fstream fOut(path2output + "/" + nameoutput + ".elem",std::fstream::out);
	fOut << elemcorrected.v1.size() << std::endl;

	for(long int i = 0; i < elemcorrected.v1.size(); i++)
		fOut << "Tt " << elemcorrected.v1[i] << " " << elemcorrected.v2[i] << " " << elemcorrected.v3[i] << " " << elemcorrected.v4[i] << " " << elemcorrected.tag[i] << std::endl;
	
	fOut.close();
}

// Returns a vector with all the files in the directory name.

std::vector<std::string> read_directory(const std::string& name){

	std::vector<std::string> v;
	std::vector<std::string> v_nohidden;

     DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }

    closedir(dirp);

	for(int i=0; i < v.size(); i++){
		if((v[i])[0] != '.')
		v_nohidden.push_back(v[i]);
	}

	return v_nohidden;
}


void dist_betw_pts(const std::string &path2pts1, const std::string &pts1name, const std::string &path2pts2, const std::string &pts2name, const std::string &path2out, const std::string &outname){
	
	points pts1 = ReadPts(path2pts1,pts1name);
	points pts2 = ReadPts(path2pts2,pts2name);

	std::vector<double> distvec(pts1.x.size());

	#pragma omp parallel
	for(long int i = 0; i < distvec.size(); i++)
		distvec[i] = sqrt(pow(pts1.x[i] - pts2.x[i],2)+pow(pts1.y[i] - pts2.y[i],2)+pow(pts1.z[i] - pts2.z[i],2));
	
	
	std::fstream fOut(path2out + "/" + outname + ".dat",std::fstream::out);

	for(long int i = 0; i < distvec.size(); i++)
		fOut << distvec[i] << std::endl;
	
	fOut.close();


}

std::vector<std::vector<double>> ReadErrorTable(const std::string &path,const std::string &filename){


	std::cout << "Reading " << filename << ".txt file...\n";

	std::ifstream inputFile(path + "/" + filename + ".txt");	

	if (!inputFile.is_open()) 
	{	
		std::cout << "Error opening file " << filename + ".txt" << std::endl;
		exit(0);
	}

	std::string line;	


	std::vector<std::vector<double>> mat(19,std::vector<double>(24));

	int i = 0;
    while(std::getline(inputFile, line))
    {

		std::istringstream lineStream(line);

		for(int j = 0; j < 24; j++){
			lineStream >> std::setprecision(10) >> mat[i][j];
		}
	i++;
	}

	inputFile.close();

	return mat;

}

void FindFirstWave(std::string& path1, std::string& file1, std::string& path2, std::string& file2, std::string& path3, std::string& file3, std::string& outPath, std::string& outName){

	std::vector<double> AT1 = ReadActTime(path1,file1);
	std::vector<double> AT2 = ReadActTime(path2,file2);
	std::vector<double> AT3 = ReadActTime(path3,file3);
	std::vector<int> vec_tag(AT1.size());

	for(int i = 0; i < AT1.size(); i++){
		if(AT1[i] < AT2[i]){
			if(AT1[i] < AT3[i]){
				vec_tag[i] = 1;
			}
			else{
				vec_tag[i] = 3;
			}
		}
		else if (AT2[i] < AT3[i])
		{
			vec_tag[i] = 2;
		}
		else
			vec_tag[i] = 3;
	}

	std::fstream fOut(outPath + "/" + outName + ".dat",std::fstream::out);

	for(long int i = 0; i < vec_tag.size(); i++)
		fOut << vec_tag[i] << std::endl;
	
	fOut.close();
}