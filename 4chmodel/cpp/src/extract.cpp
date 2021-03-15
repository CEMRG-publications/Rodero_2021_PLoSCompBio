/**
 * @brief Script with a function to extract a submesh.
 * 
 * @file extract.cpp
 * @author Marina Strocchi
 * @date 2018
 */

#include "../headers/extract.h"
/**
 * @brief Extracts a submesh from a given mesh.
 * 
 * The function extracts a submesh from a given three-dimensional mesh based on element tags. The submesh is 
 * saved in CARP format to a given path. Additionally, a <CODE>.vtx</CODE> file containing the map of the nodes
 * from the submesh to the original mesh and a <CODE>.elem.vtx</CODE> file containing the map of the elements
 * from the submesh to the original mesh are saved to make mapping easier.
 * 
 * @param mshPath(<CODE>string</CODE>) path to the mesh to extract tags from.
 * @param mshName (<CODE>string</CODE>) name of the mesh to extract tags from.
 * @param tags (<CODE>std::vector<int></CODE>) vector containing the tags to extract.
 * @param outPath (<CODE>string</CODE>) path to the output submesh.
 * @param outName (<CODE>string</CODE>) name of the output submesh.
 * @param flag (<CODE>int</CODE>) extract fibres as well (2 directions expected) with <VAR>flag</VAR>=1 or create
 * a default <CODE>.lon</CODE> file with <VAR>flag</VAR>=0/
 */


void extract(const std::string &mshPath,const std::string &mshName,const std::vector<int> &tags,const std::string &outPath,const std::string &outName,const int &flag)
{
    std::string outPath2;
	/* Read input mesh */
	points mshPts = ReadPts(mshPath,mshName);
    elem mshElem = ReadElem(mshPath,mshName);
    lon mshLon = ReadLon(mshPath,mshName);

    /* Initialize new output elements */
    std::vector<int> v1(mshElem.v1.size()),v2(mshElem.v1.size()),v3(mshElem.v1.size()),v4(mshElem.v1.size()),tag(mshElem.v1.size());
    std::vector<double> fx(mshElem.v1.size()),fy(mshElem.v1.size()),fz(mshElem.v1.size()),sx(mshElem.v1.size()),sy(mshElem.v1.size()),sz(mshElem.v1.size());
    
    elem submshElem;
    lon submshLon;

    std::vector<int> mapElem(mshElem.v1.size());

    int count = 0;

    /* Extract elements (and fibers) based on the element tags */
    for(int i = 0; i < mshElem.v1.size(); i++)
    {
    	if( std::find(tags.begin(), tags.end(), mshElem.tag[i]) != tags.end() )
    	{
    		v1[count] = mshElem.v1[i];
    		v2[count] = mshElem.v2[i];
    		v3[count] = mshElem.v3[i];
    		v4[count] = mshElem.v4[i];
    		tag[count] = mshElem.tag[i];

    		fx[count] = mshLon.x1[i];
			fy[count] = mshLon.y1[i];
			fz[count] = mshLon.z1[i];
			sx[count] = mshLon.x2[i];
			sy[count] = mshLon.y2[i];
			sz[count] = mshLon.z2[i];

    		mapElem[count] = i;

            count = count + 1;
    	}
    }

    v1.resize(count);
    v2.resize(count);
    v3.resize(count);
    v4.resize(count);
    tag.resize(count);

    fx.resize(count);
    fy.resize(count);
    fz.resize(count);
    sx.resize(count);
    sy.resize(count);
    sz.resize(count);

    mapElem.resize(count);

    submshElem.tag = tag;

    submshLon.x1 = fx;
    submshLon.y1 = fy;
    submshLon.z1 = fz;
    submshLon.x2 = sx;
    submshLon.y2 = sy;
    submshLon.z2 = sz;

    /* -----------------------------------------------------------------------------------------------

     Form a vector containing the indices of the nodes of the submesh */
    std::vector<int> vtx;

    /* Concatenate the vectors defining the four nodes of the elements */
	vtx.insert( vtx.end(), v1.begin(), v1.end() );
	vtx.insert( vtx.end(), v2.begin(), v2.end() );
	vtx.insert( vtx.end(), v3.begin(), v3.end() );
	vtx.insert( vtx.end(), v4.begin(), v4.end() );

    /* Sort and remove duplicates */
	std::sort( vtx.begin(), vtx.end() );
    std::vector<int>::iterator ip;
    ip = std::unique(vtx.begin(), vtx.end());
 
    vtx.resize(std::distance(vtx.begin(), ip));

    /* -------------------------------------------------------------------------------------------- */

    std::vector<int> v1I(v1.size()),v2I(v1.size()),v3I(v1.size()),v4I(v1.size());

    /* Re-index the elements */
    #pragma omp parallel for
    for( int i = 0; i < v1.size(); i++ )
    {	
    	std::vector<int>::iterator it;
   		it = find(vtx.begin(), vtx.end(), v1[i]);
   		v1I[i] = distance(vtx.begin(),it);

   		it = find(vtx.begin(), vtx.end(), v2[i]);
   		v2I[i] = distance(vtx.begin(),it);

   		it = find(vtx.begin(), vtx.end(), v3[i]);
   		v3I[i] = distance(vtx.begin(),it);

   		it = find(vtx.begin(), vtx.end(), v4[i]);
		v4I[i] = distance(vtx.begin(),it);
    }

    submshElem.v1 = v1I;
    submshElem.v2 = v2I;
    submshElem.v3 = v3I;
	submshElem.v4 = v4I;

    /* Check if the output folder exists. If not, create it */
    struct stat st;
    if(stat(outPath.c_str(),&st) != 0)
    {   
        system(("mkdir "+outPath).c_str());
    }

  	/* Write CARP output */
    std::fstream fOutPts;

    outPath2 = outPath + "/" + outName + ".pts";
    write2file(outPath2, fOutPts );
    fOutPts << vtx.size() << std::endl;

	for( int i = 0; i < vtx.size(); i++ )
	{
		fOutPts << mshPts.x[vtx[i]] << " " << mshPts.y[vtx[i]] << " " << mshPts.z[vtx[i]] << std::endl;    
	}	

	fOutPts.close();

    std::fstream fOutElem;
    outPath2 = outPath + "/" + outName + ".elem";
    write2file(outPath2, fOutElem );
	fOutElem << submshElem.v1.size() << std::endl;

    std::fstream fOutLon;
    outPath2 = outPath + "/" + outName + ".lon";
    write2file(outPath2, fOutLon );
	fOutLon << flag + 1 << std::endl;

	for( int i = 0; i < submshElem.v1.size(); i++ )
	{
		fOutElem << "Tt " << submshElem.v1[i] << " " << submshElem.v2[i] << " " << submshElem.v3[i] << " " << submshElem.v4[i] << " " << submshElem.tag[i] << std::endl;
		
        if(flag == 1)
        {
            fOutLon << submshLon.x1[i] << " " << submshLon.y1[i] << " " << submshLon.z1[i] << " " << submshLon.x2[i] << " " << submshLon.y2[i] << " " << submshLon.z2[i] << std::endl;
        }
        if(flag == 0)
        {
            fOutLon << 1 << " " << 0 << " " << 0 << std::endl;
        }
	}

	fOutElem.close();  
	fOutLon.close();  

	/* Write elem and pts maps */
    std::fstream fOutVtx;
    outPath2 = outPath + "/" + outName + ".vtx";
    write2file(outPath2, fOutVtx );
    fOutVtx << vtx.size() << std::endl;
    fOutVtx << "intra" << std::endl;

    for( int i = 0; i < vtx.size(); i++ )
    {
        fOutVtx << vtx[i] << std::endl;
    }

    fOutVtx.close();

    std::fstream fOutMapE;
    outPath2 = outPath + "/" + outName + ".elem.vtx";
    write2file(outPath2, fOutMapE );
    fOutMapE << mapElem.size() << std::endl;

    for( int i = 0; i < mapElem.size(); i++ )
    {
        fOutMapE << mapElem[i] << std::endl;
    }

    fOutMapE.close();
 }

void midseptum_point(const std::string &BiV_folder, const std::string &UVC_folder, const std::string &outPath){
	
	elem AHA_tags = ReadElem(BiV_folder, "BiV.retag", 0);
	UVC BiV_UVC = ReadUVC(BiV_folder + "/UVC/UVC", "COMBINED_COORDS_Z_RHO_PHI_V");
	double min_PHI = 10, max_PHI = -10, min_Z = 10, max_Z = -10, ideal_PHI, ideal_Z;
	int point_vtx;

	for(int i=0; i < BiV_UVC.PHI.size(); i++){

		if(AHA_tags.tag[i] == 9){
			if(BiV_UVC.PHI[i] < min_PHI)
				min_PHI = BiV_UVC.PHI[i];
			if(BiV_UVC.PHI[i] > max_PHI)
				max_PHI = BiV_UVC.PHI[i];
			if(BiV_UVC.Z[i] < min_Z)
				min_Z = BiV_UVC.Z[i];
			if(BiV_UVC.Z[i] > max_Z)
				max_Z = BiV_UVC.Z[i];
        }
	}

	ideal_Z = (min_Z + max_Z)/2.0;
	ideal_PHI = (min_PHI + max_PHI)/2.0;

	point_vtx = UVC2vtx(UVC_folder, "COMBINED_COORDS_Z_RHO_PHI_V", ideal_Z, 1, ideal_PHI, -1);

	std::fstream fOut(outPath + "/BiV.midseptum.vtx",std::fstream::out);

	fOut << "1\nintra\n" << point_vtx << std::endl;
	
	fOut.close();
    
    std::cout << "\nIdeal PHI: " << ideal_PHI << "\nIdeal Z: " << ideal_Z;
}