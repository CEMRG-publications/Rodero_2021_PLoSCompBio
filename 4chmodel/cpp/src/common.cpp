/**
 * @brief Script with the functions to remove common vertices and triangles.
 * 
 * @file common.cpp
 * @author Marina Strocchi
 * @date 2017-03-21
 */

#include "../headers/common.h"
/**
 * @brief Finds common vertices between two <CODE>.vtx</CODE> files.
 * 
 * The function finds common vertices from the two <CODE>.vtx</CODE> files given in input. <VAR>flag</VAR> gives the option to save 
 * the <CODE>.vtx</CODE> file with the common vertices (<VAR>flag</VAR> = 1) or to save the <CODE>.vtx</CODE> files in input without
 * the common vertices (<VAR>flag</VAR> = 0).
 * 
 * @param pathFile1 (<CODE>string</CODE>) path to the first <CODE>.vtx</CODE> file.
 * @param fileName1 (<CODE>string</CODE>) name of the first <CODE>.vtx</CODE> file.
 * @param pathFile2 (<CODE>string</CODE>) path to the second <CODE>.vtx</CODE> file.
 * @param fileName2 (<CODE>string</CODE>) name of the second <CODE>.vtx</CODE> file.
 * @param flag (<CODE>int</CODE>) whether to save the <CODE>.vtx</CODE> file of common vertices (1) or to save the vtx 
 * files without the common vertices.
 * @param fileNameOut (<CODE>string</CODE>) name of the output file for <VAR>flag</VAR> = 1.
 */

void commonVertices(const std::string &pathFile1,const std::string &fileName1,const std::string &pathFile2,const std::string &fileName2,const int &flag,const std::string &fileNameOut)
{	
	std::string outPath;
	/* Read in the files */
	std::vector<int> vtx1 = ReadVtx(pathFile1,fileName1);
	std::vector<int> vtx2 = ReadVtx(pathFile2,fileName2);

	/* Initialize vectors without common verticies */
	std::vector<int> Vtx1;
	std::vector<int> Vtx2;
	std::vector<int> commonVtx(vtx1.size());

	/* Find common vertices */
	
	for( int i = 0; i < vtx1.size(); i++ )
	{	
		// int flag = 1;
		for( int j = 0; j < vtx2.size(); j++ )
		{
			if(vtx1[i]==vtx2[j])
			{
				commonVtx.push_back(vtx1[i]);
			}
		}
	}

	if(flag==0)
	{
		/* Remove common vertices */
		std::set_difference(vtx1.begin(), vtx1.end(), commonVtx.begin(), commonVtx.end(), std::inserter(Vtx1, Vtx1.begin()));
		std::set_difference(vtx2.begin(), vtx2.end(), commonVtx.begin(), commonVtx.end(), std::inserter(Vtx2, Vtx2.begin()));	

		/* Write .vtx files without common nodes */
		std::fstream fOutVtx1;
		outPath = pathFile1 + "/" + fileName1 + "_rm" + ".vtx";
		write2file(outPath, fOutVtx1 );
		fOutVtx1 << Vtx1.size() << std::endl;
		fOutVtx1 << "intra" << std::endl;

		if(Vtx1.size() > 0){
		for( int i = 0; i < Vtx1.size(); i++ )
			fOutVtx1 << Vtx1.at(i) << std::endl;
		}	

		fOutVtx1.close();	

		std::fstream fOutVtx2;
		outPath = pathFile2 + "/" + fileName2 + "_rm" + ".vtx";
		write2file(outPath, fOutVtx2 );
		fOutVtx2 << Vtx2.size() << std::endl;
		fOutVtx2 << "intra" << std::endl;	

		if(Vtx2.size() > 0){
		for( int i = 0; i < Vtx2.size(); i++ )
			fOutVtx2 << Vtx2.at(i) << std::endl;
		}	

		fOutVtx2.close();
	}
	else if(flag==1)
	{
		/* Write common vertices to an output file */ 
		std::fstream fOutVtx;
		outPath = pathFile1 + "/" + fileNameOut + ".vtx";
		write2file(outPath, fOutVtx );
		fOutVtx << commonVtx.size() << std::endl;
		fOutVtx << "intra" << std::endl;	

		for( int i = 0; i < commonVtx.size(); i++ )
		{
			fOutVtx << commonVtx.at(i) << std::endl;
		}	

		fOutVtx.close();	

	}

}
/**
 * @brief Finds and removes common triangles to two surfaces defined by both <CODE>.surf</CODE> and <CODE>.neubc</CODE> files.
 * 
 * @param pathFile1 (<CODE>string</CODE>) path to the first <CODE>.surf</CODE> and <CODE>.neubc</CODE> file.
 * @param fileName1 (<CODE>string</CODE>) name of the first <CODE>.surf</CODE> and <CODE>.neubc</CODE> file.
 * @param pathFile2 (<CODE>string</CODE>) path to the second <CODE>.surf</CODE> and <CODE>.neubc</CODE> file.
 * @param fileName2 (<CODE>string</CODE>) name of the second <CODE>.surf</CODE> and <CODE>.neubc</CODE> file. 
 */

void rmCommonTriangles(const std::string &pathFile1,const std::string &fileName1,const std::string &pathFile2,const std::string &fileName2)
{	
	std::string outPath;
	/* Read in the files */
	surf surf1 = ReadSurf(pathFile1,fileName1);
	surf surf2 = ReadSurf(pathFile2,fileName2);
	// surf surf1 = ReadSurfElem(pathFile1,fileName1);
	// surf surf2 = ReadSurfElem(pathFile2,fileName2);

	neubc neubc1 = ReadNeubc(pathFile1,fileName1);
	neubc neubc2 = ReadNeubc(pathFile2,fileName2);


	/* Find common triangles */
	surf commonTr;

	for( int i = 0; i < surf1.v1.size(); i++ )
	{
		for( int j = 0; j < surf2.v1.size(); j++ )
		{
			if( (neubc1.v1[i]==neubc2.v1[j]) && (neubc1.v2[i]==neubc2.v2[j]) && (neubc1.v3[i]==neubc2.v3[j]) )
			{
				commonTr.v1.push_back(neubc1.v1[i]);
				commonTr.v2.push_back(neubc1.v2[i]);
				commonTr.v3.push_back(neubc1.v3[i]);
			}
		}
	}
		std::cout << commonTr.v1.size() <<" common triangles found" << std::endl;


	/* Write .surf and .neubc files without common nodes */
	std::fstream fOutSurf1;
	outPath = pathFile1 + "/" + fileName1 + "_rm" + ".surf";
	write2file(outPath, fOutSurf1 );
	fOutSurf1 << surf1.v1.size() - commonTr.v1.size() << std::endl;
	
	std::fstream fOutNeubc1(pathFile1 + "/" + fileName1 + "_rm" + ".neubc",std::fstream::out);
	fOutNeubc1 << neubc1.v1.size() - commonTr.v1.size() << " # " << "mesh elast " << neubc1.nElem << " " << neubc1.nNodes << " :" << std::endl;
	
	for(int i = 0; i < neubc1.v1.size(); i++)
	{
		int flag = 1;
		for(int j = 0; j < commonTr.v1.size(); j++)
		{
			if( (neubc1.v1[i] == commonTr.v1[j]) && (neubc1.v2[i] == commonTr.v2[j]) && (neubc1.v3[i] == commonTr.v3[j]))
			{
				flag = 0;
			}
		}
		if(flag == 1)
		{
			fOutSurf1 << "Tr " << surf1.v1.at(i) << " " << surf1.v2.at(i) << " "  << surf1.v3.at(i) << std::endl;
			fOutNeubc1 << neubc1.v1.at(i) << " " << neubc1.v2.at(i) << " "  << neubc1.v3.at(i) << " " << 1 << " " << neubc1.v4.at(i) << " " << neubc1.tet.at(i) << std::endl;
		}
	}

	fOutSurf1.close();
	fOutNeubc1.close();

	std::fstream fOutSurf2;
	outPath = pathFile2 + "/" + fileName2 + "_rm" + ".surf";
	write2file(outPath, fOutSurf2 );
	fOutSurf2 << surf2.v1.size() - commonTr.v1.size() << std::endl;

	std::fstream fOutNeubc2(pathFile2 + "/" + fileName2 + "_rm" + ".neubc",std::fstream::out);
	fOutNeubc2 << neubc2.v1.size() - commonTr.v1.size() << " # " << "mesh elast " << neubc2.nElem << " " << neubc2.nNodes << " :" << std::endl;

	for( int i = 0; i < neubc2.v1.size(); i = i + 1)
	{
		int flag = 1;
		for(int j = 0; j < commonTr.v1.size(); j = j + 1)
		{
			if( (neubc2.v1[i] == commonTr.v1[j]) && (neubc2.v2[i] == commonTr.v2[j]) && (neubc2.v3[i] == commonTr.v3[j]))
			{
				flag = 0;
			}
		}
		if(flag == 1)
		{
			fOutSurf2 << "Tr " << surf2.v1.at(i) << " " << surf2.v2.at(i) << " "  << surf2.v3.at(i) << std::endl;
			fOutNeubc2 << neubc2.v1.at(i) << " " << neubc2.v2.at(i) << " "  << neubc2.v3.at(i) << " " << 1 << " " << neubc2.v4.at(i) << " " << neubc2.tet.at(i) << std::endl;
		}
	}

	fOutSurf2.close();
	fOutNeubc2.close();

}

void commonTriangles(const std::string &pathFile1,const std::string &fileName1,const std::string &pathFile2,const std::string &fileName2, const std::string &outName)
{	
	std::string outPath;
	/* Read in the files */
	// surf surf1 = ReadSurfElem(pathFile1,fileName1);
	// surf surf2 = ReadSurfElem(pathFile2,fileName2);

	surf surf1 = ReadSurf(pathFile1,fileName1);
	surf surf2 = ReadSurf(pathFile2,fileName2);

	neubc neubc1 = ReadNeubc(pathFile1,fileName1);
	neubc neubc2 = ReadNeubc(pathFile2,fileName2);


	/* Find common triangles */
	surf commonTr;

	for( int i = 0; i < surf1.v1.size(); i++ )
	{
		for( int j = 0; j < surf2.v1.size(); j++ )
		{
			if( (neubc1.v1[i]==neubc2.v1[j]) && (neubc1.v2[i]==neubc2.v2[j]) && (neubc1.v3[i]==neubc2.v3[j]) )
			{
				commonTr.v1.push_back(neubc1.v1[i]);
				commonTr.v2.push_back(neubc1.v2[i]);
				commonTr.v3.push_back(neubc1.v3[i]);
			}
		}
	}

	std::cout << commonTr.v1.size() <<" common triangles found" << std::endl;

	/* Write .surf and .neubc files without common nodes */
	std::fstream fOutSurf1;
	outPath = pathFile1 + "/" + outName + ".surf";
	write2file(outPath, fOutSurf1 );
	fOutSurf1 << commonTr.v1.size() << std::endl;
	
	for(int i = 0; i < neubc1.v1.size(); i++)
	{
		int flag = 1;
		for(int j = 0; j < commonTr.v1.size(); j++)
		{
			if( (neubc1.v1[i] == commonTr.v1[j]) && (neubc1.v2[i] == commonTr.v2[j]) && (neubc1.v3[i] == commonTr.v3[j]))
			{
				flag = 0;
			}
		}
		if(flag == 1)
		{
			fOutSurf1 << "Tr " << surf1.v1.at(i) << " " << surf1.v2.at(i) << " "  << surf1.v3.at(i) << std::endl;
		}
	}

	fOutSurf1.close();

}


void rmVtxTriangles(std::string path2vtx, std::string vtxName, std::string path2surf2, std::string surfName2)
{	

	/* Read in the files */
	surf surf2 = ReadSurf(path2surf2,surfName2);

    std::fstream vtx2File;
    vtx2File.open(path2surf2 + "/" + surfName2 + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtxSurf2; 
    if (vtx2File.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtxSurf2 = ReadVtx(path2surf2,surfName2);
    } else
    {
	    /* Generate a .vtx from the input surface */ 
		vtxSurf2.reserve( surf2.v1.size()*3  ); // preallocate memory
		vtxSurf2.insert( vtxSurf2.end(), surf2.v1.begin(), surf2.v1.end() );
		vtxSurf2.insert( vtxSurf2.end(), surf2.v2.begin(), surf2.v2.end() );
		vtxSurf2.insert( vtxSurf2.end(), surf2.v3.begin(), surf2.v3.end() );
		std::sort( vtxSurf2.begin(), vtxSurf2.end() );
	    std::vector<int>::iterator ip_surf2;
	    ip_surf2 = std::unique(vtxSurf2.begin(), vtxSurf2.end());
	 
	    vtxSurf2.resize(std::distance(vtxSurf2.begin(), ip_surf2)); 
    }
	
    std::vector<int> vtxSurf1 = ReadVtx(path2vtx,vtxName);

    /* Correct the vertices forming a triangle on the second surface. This 
       ensures that the two surfaces touch only with a line */
    int count_vtx = 0;
	std::vector<int> v2(3);
	std::vector<int> intersection2(3);
	std::vector<int> vtxRemove;
	std::vector<int>::iterator it2;

	std::vector<int> vtxSurf2_rep; 
	std::vector<int> score(3);
	vtxSurf2_rep.reserve( surf2.v1.size()*3  ); // preallocate memory
	vtxSurf2_rep.insert( vtxSurf2_rep.end(), surf2.v1.begin(), surf2.v1.end() );
	vtxSurf2_rep.insert( vtxSurf2_rep.end(), surf2.v2.begin(), surf2.v2.end() );
	vtxSurf2_rep.insert( vtxSurf2_rep.end(), surf2.v3.begin(), surf2.v3.end() );
	std::sort( vtxSurf2_rep.begin(), vtxSurf2_rep.end() );


    for(int i = 0; i < surf2.v1.size(); i++)
    {	
    	v2 = {surf2.v1[i],surf2.v2[i],surf2.v3[i]};
		std::sort(v2.begin(),v2.end());
		it2 = std::set_intersection(v2.begin(),v2.end(),vtxSurf1.begin(),vtxSurf1.end(), intersection2.begin());
		intersection2.resize(it2-intersection2.begin());  
		if(intersection2.size() > 2)
		{	
			score[0] = std::count(vtxSurf2_rep.begin(), vtxSurf2_rep.end(), v2[0]);
			score[1] = std::count(vtxSurf2_rep.begin(), vtxSurf2_rep.end(), v2[1]);
			score[2] = std::count(vtxSurf2_rep.begin(), vtxSurf2_rep.end(), v2[2]);

			std::vector<int>::iterator max = std::max_element(score.begin(), score.end());
    		int max_ind = std::distance(score.begin(),max);
			vtxRemove.push_back(v2[max_ind]);
			count_vtx = count_vtx+1;
		}
		v2.clear();
		intersection2.clear();
    }

    std::vector<int> vtxFinal;
    std::sort(vtxRemove.begin(),vtxRemove.end());
    std::vector<int>::iterator ip_r;
    ip_r = std::unique(vtxRemove.begin(), vtxRemove.end());
 
    vtxRemove.resize(std::distance(vtxRemove.begin(), ip_r));    

    std::cout << "Found " << vtxRemove.size() <<" vertices to remove" << std::endl;

    std::set_difference(vtxSurf1.begin(), vtxSurf1.end(), vtxRemove.begin(), vtxRemove.end(), std::inserter(vtxFinal, vtxFinal.begin()));

	/* Write .vtx file without common nodes */
	std::fstream fOutVtx;
	std::string outpath = path2vtx + "/" + vtxName + "_rm" + ".vtx";
	write2file(outpath , fOutVtx );
	fOutVtx << vtxFinal.size() << std::endl;
	fOutVtx << "intra" << std::endl;

    for( int i = 0; i < vtxFinal.size(); i++ )
    {
    	fOutVtx << vtxFinal[i] << std::endl;
    }

    fOutVtx.close();

}

void rmSeptum(const std::string& path_elem, const std::string& elem_file, const std::string& path_septum, const std::string& septum_file, const std::string& out_path, const std::string& out_file){

	/* Read in the files */
	elem whole_heart = ReadElem(path_elem,elem_file);
	surf septum_surf = ReadSurf(path_septum,septum_file);
	bool flag;
	int j;

	for(int i = 0; i < septum_surf.v1.size(); i++){
		flag = false;
		j = 0;
		while(flag == false){
			if(whole_heart.tag[j] == 2 || whole_heart.tag[j] == 26){
				if(whole_heart.v1[j] == septum_surf.v1[i] || whole_heart.v2[j] == septum_surf.v1[i]  || whole_heart.v3[j] == septum_surf.v1[i] || whole_heart.v4[j] == septum_surf.v1[i]
				|| whole_heart.v1[j] == septum_surf.v2[i] || whole_heart.v2[j] == septum_surf.v2[i]  || whole_heart.v3[j] == septum_surf.v2[i] || whole_heart.v4[j] == septum_surf.v2[i]
				|| whole_heart.v1[j] == septum_surf.v3[i] || whole_heart.v3[j] == septum_surf.v3[i]  || whole_heart.v3[j] == septum_surf.v3[i] || whole_heart.v4[j] == septum_surf.v3[i]){
					whole_heart.tag[j] = 1;
					flag = true;
				}
			}
			else 
				j++;
		}
	}

	/* Write .elem */
	std::fstream fOut;
	std::string outPath = out_path + "/" + out_file + ".elem";
	write2file(outPath, fOut);
	fOut << whole_heart.v1.size() << std::endl;
	
	for(int i = 0; i < whole_heart.v1[i]; i++)
		fOut << "Tt " << whole_heart.v1[i] << " " << whole_heart.v2[i] << " " << whole_heart.v3[i] << " " << whole_heart.v4[i] << whole_heart.tag[i] << std::endl;
}