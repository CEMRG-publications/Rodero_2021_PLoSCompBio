/**
 * @brief Script to map submeshes to the original mesh and the other way around.
 * 
 * @file map.cpp
 * @author Marina Strocchi
 * @date 2018
 */

#include "../headers/map.h"

struct mapNodes
{   
    int nPoints;
    std::vector<int> n1;
    std::vector<int> n2;
};
/**
 * @brief The function maps a set of three-dimensional points of a surface mesh to a 'bigger' set of three-dimensional
    points of the original mesh.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file of the original mesh the surface was extracted from.
 * @param submshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the submesh file.
 * @param submshName (<CODE>string</CODE>) name of the surface file.
 * @return (<CODE>std::vector<int></CODE>) <CODE>vtx</CODE> vector containing the index of the nodes of the submesh into the original mesh.
 */

std::vector<int> mapPoints(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName, const std::string &outPath, const std::string &outName)
{

    std::string output_name;
    /* Read .pts files */
    points meshPts = ReadPts(mshPath,mshName);
    points submeshPts = ReadPts(submshPath,submshName);

    /* Store the .vtx to not read it back in later */
    std::vector<int> vtxOut(submeshPts.x.size());
    std::cout << "Mapping points of the submesh in the original mesh...\n";

    #pragma omp parallel for
    for( int i = 0; i < submeshPts.x.size(); i++ ) 
    {   
        int flag = 1;
        #pragma omp parallel for
        for( int j = 0; j < meshPts.x.size(); j++ )
        {   

            /* .pts file points in the submesh have a different precision */
            double diffX = std::abs(submeshPts.x[i] - meshPts.x[j]);
            double diffY = std::abs(submeshPts.y[i] - meshPts.y[j]);
            double diffZ = std::abs(submeshPts.z[i] - meshPts.z[j]);

            if(diffX < 0.1)
                if(diffY < 0.1)
                    if(diffZ < 0.1){
                vtxOut[i] = j;
                flag = 0;
            }
        }

        // if(flag)
        // {
        //     std::cout << "For point " << i << " there was not a correspondent point on the original mesh" << std::endl;
        // }

    }

    std::sort(vtxOut.begin(),vtxOut.end());

    std::cout << " Done" <<  std::endl;

    /* Open .vtx output file */ 
    std::fstream fOutVtx;
    output_name = outPath + "/" + outName + ".vtx";
    write2file(output_name, fOutVtx );
    fOutVtx << vtxOut.size() << std::endl;
    fOutVtx << "intra" << std::endl;

    for( int k = 0; k < vtxOut.size(); k++ )
    {
        fOutVtx << vtxOut[k] << std::endl;
    }

    fOutVtx.close();

    return vtxOut;

}
/**
 * @brief Maps a set of triangles of a surface mesh to the original mesh - both <CODE>.elem</CODE> and <CODE>.neubc</CODE> files
 * are mapped to the original mesh.
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file of the original mesh the surface was extracted from.
 * @param submshName (<CODE>string</CODE>) name of the submesh.
 * @param neubcName (<CODE>string</CODE>) name of the corresponding <CODE>.neubc</CODE> file.
 */

void mapSurface(const std::string &path,const std::string &mshName,const std::string &submshPath, const std::string &submshName,const std::string &neubcName, const std::string &outPath, const std::string &outName)
{   


    std::string output_name;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in); // Open the vtx file of the submesh

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(path,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read .elem file */
    surf elemSubMsh = ReadSurfElem(submshPath,submshName);

    std::cout << "Mapping surface triangles on the original mesh..." << std::endl;

    surf elemMsh;

    std::vector<int> v1(elemSubMsh.v1.size());
    std::vector<int> v2(elemSubMsh.v2.size());
    std::vector<int> v3(elemSubMsh.v3.size());

    /* Open .surf output file */ 
    std::fstream fOutSurf;

    output_name = outPath + "/" + outName + ".surf";
    write2file(output_name,fOutSurf );
    fOutSurf << elemSubMsh.v1.size() << std::endl;

    /* Map the surface triangles in the original mesh */
    for( int i = 0; i < elemSubMsh.v1.size(); i++)
    {   

        v1[i] = vtx[elemSubMsh.v1[i]];
        v2[i] = vtx[elemSubMsh.v2[i]];
        v3[i] = vtx[elemSubMsh.v3[i]];

        fOutSurf << "Tr " << v1[i] << " " << v2[i] << " " << v3[i] << std::endl; 

    }

    std::cout << " Done" << std::endl;

    fOutSurf.close();

    elemMsh.v1 = v1;
    elemMsh.v2 = v2;
    elemMsh.v3 = v3;

    /* --------------------------------------------------------------------------------------------*/
    std::ifstream inputneubc(submshPath + "/" + neubcName + ".neubc");	

	if (inputneubc.is_open()) 
	{	
    /* Read .neubc file */
    neubc neubcIn = ReadNeubc(submshPath,neubcName);

    neubc neubcOut;
    std::vector<int> v1Out(elemMsh.v1.size()); 
    std::vector<int> v2Out(elemMsh.v1.size()); 
    std::vector<int> v3Out(elemMsh.v1.size()); 
    std::vector<int> v4Out(elemMsh.v1.size()); 
    std::vector<int> tetOut(elemMsh.v1.size()); 

    /* Extract a new .neubc file for the submesh */
    #pragma omp parallel for 
    for( int i = 0; i < elemMsh.v1.size(); i++ )
    {
        #pragma omp parallel for
        for( int j = 0; j < neubcIn.v1.size(); j++ )
        {
            if(elemMsh.v1[i]==neubcIn.v1[j])
                if(elemMsh.v2[i]==neubcIn.v2[j])
                    if(elemMsh.v3[i]==neubcIn.v3[j]){   
                v1Out[i] = neubcIn.v1[j]; 
                v2Out[i] = neubcIn.v2[j];
                v3Out[i] = neubcIn.v3[j];
                v4Out[i] = neubcIn.v4[j];
                tetOut[i] = neubcIn.tet[j];
            }
        }       
    }

    neubcOut.v1 = v1Out;
    neubcOut.v2 = v2Out;
    neubcOut.v3 = v3Out;
    neubcOut.v4 = v4Out;
    neubcOut.tet = tetOut;

    std::cout << "Writing .neubc submesh file..." << std::endl;

    /* Open output .neubc file */
    std::fstream fOutNeubc;
    
    output_name = outPath + "/" + outName + ".neubc";
    write2file(output_name, fOutNeubc );
    fOutNeubc << elemMsh.v1.size() << " # " << "mesh elast " << neubcIn.nElem << " " << neubcIn.nNodes << " :" << std::endl;

    for( int i = 0; i < neubcOut.v1.size(); i++ )
    {
        fOutNeubc << neubcOut.v1[i] << " " << neubcOut.v2[i] << " " << neubcOut.v3[i] << " 1 " << neubcOut.v4[i] << " " << neubcOut.tet[i] << std::endl;
    }

    std::cout<< " Done" << std::endl;

    fOutNeubc.close();
    }

}/**
 * @brief Maps <CODE>.elem</CODE> file of a submesh into the <CODE>.elem</CODE> file of the original mesh. 
 * 
 * @param mshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file of the original mesh the surface was extracted from.
 * @param submshPath (<CODE>string</CODE>) path to the submesh.
 * @param submshName (<CODE>string</CODE>) name of the submesh.
 * @return (<CODE>std::vector<int></CODE>) index vector containing the index of the element of the submesh into the original mesh.

 */

std::vector<int> mapElem(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName)
{   

    std::fstream inputFile;
    inputFile.open(mshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(mshPath,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read .elem files */
    elem mshElem = ReadElem(mshPath,mshName);
    elem submshElem = ReadElem(submshPath,submshName);

    elem mapElem = submshElem;

    /* Transform the submesh elements to the original mesh making use of the vtx file */
    #pragma omp parallel for
    for( int i = 0; i < submshElem.v1.size(); i++ )
    {
        mapElem.v1[i] = vtx[submshElem.v1[i]];
        mapElem.v2[i] = vtx[submshElem.v2[i]];
        mapElem.v3[i] = vtx[submshElem.v3[i]];
        mapElem.v4[i] = vtx[submshElem.v4[i]];
    }

    /* Generate the mapping vector */
    std::vector<int> index(mapElem.v1.size());

    #pragma omp parallel for
    for( int i = 0; i < mapElem.v1.size(); i++ )
    {   
        for( int j = 0; j < mshElem.v1.size(); j++ )
        {
            if(mapElem.v1[i] == mshElem.v1[j])
                if(mapElem.v2[i] == mshElem.v2[j])
                    if(mapElem.v3[i] == mshElem.v3[j])
                        if(mapElem.v4[i] == mshElem.v4[j])
                            index[i] = j;
        }
    }

    return index;
}
/**
 * @brief Maps a biventricular mesh to the whole heart mesh.
 * 
 * The function reads in the base, the apex, the endocardial and epicardial meshes (mapped to the original mesh),
 * maps them to the local mesh and writes in output all the files needed for the <CODE>GlRuleFibres</CODE>
 * (<CODE>.vtx</CODE> files for apex, base, epi, endoRV and endoLV).
 * 
 * @param path (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the original mesh. 
 * @param submshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the submesh file.
 * @param submshName (<CODE>string</CODE>) name of the CARP submesh.
 */

void mapBiV(const std::string &path,const std::string &mshName,const std::string &submshPath,const std::string &submshName)
{
    std::string outPath;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(path,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* -------------------------------------------------------------------------- */
    /* Read the vtx files */

    std::vector<int> apex = ReadVtx(path,"apex");
    std::vector<int> LV_base = ReadVtx(path,"LV_base.surf_rm_rm");
    std::vector<int> RV_base = ReadVtx(path,"RV_base.surf_rm_rm");
    std::vector<int> LV_epi = ReadVtx(path,"LV_epi.surfmesh_rm_rm_rm");
    std::vector<int> RV_epi = ReadVtx(path,"RV_epi.surfmesh_rm_rm");
    std::vector<int> LV_endo = ReadVtx(path,"LV_endo.surfmesh_rm_rm_rm");
    std::vector<int> RV_endo = ReadVtx(path,"RV_endo.surfmesh_rm_rm");
    std::vector<int> septum = ReadVtx(path,"septum.surfmesh_rm_rm");

    std::vector<int>::iterator ip;

    /* -------------------------------------------------------------------------- */
    /* Map the apex */
    std::vector<int> apexMap(apex.size());
    for( int i = 0; i < apex.size(); i++ )
    {
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), apex[i]);
        apexMap[i] = distance(vtx.begin(),it);
    }
    std::sort(apexMap.begin(),apexMap.end());
    ip = std::unique(apexMap.begin(), apexMap.end());
    apexMap.resize(std::distance(apexMap.begin(), ip));

    std::fstream fOutApex;
    outPath = submshPath + "/apex.vtx";
    write2file(outPath, fOutApex );
    fOutApex << apex.size() << std::endl;
    fOutApex << "intra" << std::endl;
    for( int i = 0; i < apexMap.size(); i++ )
    {
        fOutApex << apexMap[i] << std::endl;
    }

    fOutApex.close();

    /* -------------------------------------------------------------------------- */
    /* Map the base of the LV and the RV and combine them in a unique file*/
    std::vector<int> baseMap(LV_base.size()+RV_base.size());

    #pragma omp parallel for
    for( int i = 0; i < LV_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), LV_base[i]);
        baseMap[i] = distance(vtx.begin(),it);
    }

    #pragma omp parallel for
    for( int i = 0; i < RV_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), RV_base[i]);
        baseMap[i+LV_base.size()] = distance(vtx.begin(),it);
    } 
    std::sort(baseMap.begin(),baseMap.end());
    ip = std::unique(baseMap.begin(), baseMap.end());
    baseMap.resize(std::distance(baseMap.begin(), ip));

    std::fstream fOutBase;
    outPath = submshPath + "/base.vtx";
    write2file(outPath, fOutBase );
    fOutBase << baseMap.size() << std::endl;
    fOutBase << "intra" << std::endl;
    for( int i = 0; i < baseMap.size(); i++ )
    {   
        fOutBase << baseMap[i] << std::endl;
    } 

    fOutBase.close();

    /* -------------------------------------------------------------------------- */
    /* Maps the LV endocardial surface and writes it in the endo.vtx file (endoLV+endoRV),
       in the endo_LV.vtx file alone and in the epi_RV.vtx (everything but RV endo) */
    std::vector<int> LV_endoMap(LV_endo.size());
    std::vector<int> endoMap(LV_endo.size()+RV_endo.size()+septum.size());
    std::vector<int> RV_epiMap(LV_endo.size()+septum.size()+LV_epi.size()+RV_epi.size());

    #pragma omp parallel for
    for( int i = 0; i < LV_endo.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), LV_endo[i]);
        LV_endoMap[i] = distance(vtx.begin(),it);
        endoMap[i] = distance(vtx.begin(),it);
        RV_epiMap[i] = distance(vtx.begin(),it);
    }   

    /* -------------------------------------------------------------------------- */
    /* Maps the RV endocardial surface and writes it in the endo.vtx file (endoLV+endoRV),
       in the endo_RV.vtx file alone and in the epi_LV.vtx (everything but LV endo) */
    std::vector<int> RV_endoMap(RV_endo.size()+septum.size());
    std::vector<int> LV_epiMap(RV_endo.size()+septum.size()+LV_epi.size()+RV_epi.size());

    #pragma omp parallel for
    for( int i = 0; i < RV_endo.size(); i++ )
    {
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), RV_endo[i]);
        RV_endoMap[i] = distance(vtx.begin(),it);
        endoMap[i+LV_endo.size()] = distance(vtx.begin(),it);
        LV_epiMap[i] = distance(vtx.begin(),it);
    }   

    #pragma omp parallel for
    for( int i = 0; i < septum.size(); i++ )
    {
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), septum[i]);
        RV_endoMap[i+RV_endo.size()] = distance(vtx.begin(),it);
        endoMap[i+LV_endo.size()+RV_endo.size()] = distance(vtx.begin(),it);
        LV_epiMap[i+RV_endo.size()] = distance(vtx.begin(),it);
    }  

    /* -------------------------------------------------------------------------- */
    /* Maps the LV epicardial surface and writes it in the epi.vtx file (epiLV + epiRV),
       in the epi_RV.vtx file alone and in the epi_LV.vtx */
    std::vector<int> epiMap(LV_epi.size()+RV_epi.size());
    
    #pragma omp parallel for
    for( int i = 0; i < LV_epi.size(); i++ )
    { 
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), LV_epi[i]);
        epiMap[i] = distance(vtx.begin(),it);
        RV_epiMap[i+LV_endo.size()] = distance(vtx.begin(),it);
        LV_epiMap[i+RV_endo.size()+septum.size()] = distance(vtx.begin(),it);
    }       

    #pragma omp parallel for
    for( int i = 0; i < RV_epi.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), RV_epi[i]);
        epiMap[i+LV_epi.size()] = distance(vtx.begin(),it);
        RV_epiMap[i+LV_endo.size()+LV_epi.size()] = distance(vtx.begin(),it);
        LV_epiMap[i+RV_endo.size()+septum.size()+LV_epi.size()] = distance(vtx.begin(),it);
    }

    /* -------------------------------------------------------------------------- */

    std::sort(endoMap.begin(),endoMap.end());
    ip = std::unique(endoMap.begin(), endoMap.end());
    endoMap.resize(std::distance(endoMap.begin(), ip));

    std::fstream fOutEndo;
    outPath = submshPath + "/endo.vtx";
    write2file(outPath, fOutEndo );
    fOutEndo << endoMap.size() << std::endl;
    fOutEndo << "intra" << std::endl;
    for(int i = 0; i < endoMap.size(); i++)
    {
        fOutEndo << endoMap[i] << std::endl;
    }
    fOutEndo.close();

    std::sort(LV_endoMap.begin(),LV_endoMap.end());
    ip = std::unique(LV_endoMap.begin(), LV_endoMap.end());
    LV_endoMap.resize(std::distance(LV_endoMap.begin(), ip));

    std::fstream fOutLVEndo;
    outPath = submshPath + "/endo_LV.vtx";
    write2file(outPath, fOutLVEndo );
    fOutLVEndo << LV_endoMap.size() << std::endl;
    fOutLVEndo << "intra" << std::endl;
    for(int i = 0; i < LV_endoMap.size(); i++)
    {
        fOutLVEndo << LV_endoMap[i] << std::endl;
    }
    fOutLVEndo.close();

    std::sort(RV_endoMap.begin(),RV_endoMap.end());
    ip = std::unique(RV_endoMap.begin(), RV_endoMap.end());
    RV_endoMap.resize(std::distance(RV_endoMap.begin(), ip));

    std::fstream fOutRVEndo;
    outPath = submshPath + "/endo_RV.vtx";
    write2file(outPath, fOutRVEndo );
    fOutRVEndo << RV_endoMap.size() << std::endl;
    fOutRVEndo << "intra" << std::endl;
    for(int i = 0; i < RV_endoMap.size(); i++)
    {
        fOutRVEndo << RV_endoMap[i] << std::endl;
    }
    fOutRVEndo.close();

    std::sort(epiMap.begin(),epiMap.end());
    ip = std::unique(epiMap.begin(), epiMap.end());
    epiMap.resize(std::distance(epiMap.begin(), ip));    

    std::fstream fOutEpi;
    outPath = submshPath + "/epi.vtx"; 
    write2file(outPath, fOutEpi );
    fOutEpi << epiMap.size() << std::endl;
    fOutEpi << "intra" << std::endl;
    for(int i = 0; i <epiMap.size(); i++)
    {
        fOutEpi << epiMap[i] << std::endl;
    }
    fOutEpi.close();

    std::sort(LV_epiMap.begin(),LV_epiMap.end());
    ip = std::unique(LV_epiMap.begin(), LV_epiMap.end());
    LV_epiMap.resize(std::distance(LV_epiMap.begin(), ip));  

    std::fstream fOutLVEpi;
    outPath = submshPath + "/epi_LV.vtx";
    write2file(outPath, fOutLVEpi );
    fOutLVEpi << LV_epiMap.size() << std::endl;
    fOutLVEpi << "intra" << std::endl;
    for(int i = 0; i <LV_epiMap.size(); i++)
    {
        fOutLVEpi << LV_epiMap[i] << std::endl;
    }
    fOutLVEpi.close();

    std::sort(RV_epiMap.begin(),RV_epiMap.end());
    ip = std::unique(RV_epiMap.begin(), RV_epiMap.end());
    RV_epiMap.resize(std::distance(RV_epiMap.begin(), ip));  

    std::fstream fOutRVEpi;
    outPath = submshPath + "/epi_RV.vtx";
    write2file(outPath, fOutRVEpi );
    fOutRVEpi << RV_epiMap.size() << std::endl;
    fOutRVEpi << "intra" << std::endl;
    for(int i = 0; i <RV_epiMap.size(); i++)
    {
        fOutRVEpi << RV_epiMap[i] << std::endl;
    }
    fOutRVEpi.close();

}

void mapExtendedBase(const std::string &path, const std::string &mshName, const std::string &submshPath, const std::string &submshName){
        
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(path,mshName,submshPath,submshName,submshPath,submshName);
    }

    std::vector<int> LV_base = ReadVtx(path,"LV_base");
    std::vector<int> RV_base = ReadVtx(path,"RV_base");
    std::vector<int> LV_Ao_base = ReadVtx(path,"LV_Ao_base");
    std::vector<int> RV_PA_base = ReadVtx(path,"RV_PA_base");


    /* Map the base of the LV and the RV and combine them in a unique file*/
    std::vector<int> baseMap(LV_base.size()+RV_base.size() + LV_Ao_base.size() + RV_PA_base.size());

    #pragma omp parallel for
    for( int i = 0; i < LV_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), LV_base[i]);
        baseMap[i] = distance(vtx.begin(),it);
    }

    #pragma omp parallel for
    for( int i = 0; i < RV_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), RV_base[i]);
        baseMap[i+LV_base.size()] = distance(vtx.begin(),it);
    } 

    #pragma omp parallel for
    for( int i = 0; i < LV_Ao_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), LV_Ao_base[i]);
        baseMap[i+LV_base.size() + RV_base.size()] = distance(vtx.begin(),it);
    } 

    #pragma omp parallel for
    for( int i = 0; i < RV_PA_base.size(); i++ )
    {   
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), RV_PA_base[i]);
        baseMap[i+LV_base.size() + RV_base.size() + LV_Ao_base.size()] = distance(vtx.begin(),it);
    } 

    std::sort(baseMap.begin(),baseMap.end());
    std::vector<int>::iterator ip = std::unique(baseMap.begin(), baseMap.end());
    baseMap.resize(std::distance(baseMap.begin(), ip));

    std::fstream fOutBase;
    std::string outPath = submshPath + "/extended_base.vtx";
    write2file(outPath, fOutBase );
    fOutBase << baseMap.size() << std::endl;
    fOutBase << "intra" << std::endl;
    for( int i = 0; i < baseMap.size(); i++ )
    {   
        fOutBase << baseMap[i] << std::endl;
    } 

    fOutBase.close();
}

/**
 * @brief Map the fibres from a whole heart mesh to a submesh.
 * 
 * The function reads in the submesh files and maps the <CODE>.lon</CODE> fibres into the original file.
 * The function writes in output a new fibre <CODE>.lon</CODE> file for the original mesh containing as new
 * fibres the fibres of the submesh.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the CARP mesh files.
 * @param mshName (<CODE>string</CODE>) name of the CARP mesh. 
 * @param submshPath (<CODE>string</CODE>) path to the CARP submesh files.
 * @param submshName (<CODE>string</CODE>) name of the CARP submesh.
 * @param tags (<CODE>std::vector<int></CODE>) tags of the zones to map the fibres.
 * @param fileOut (<CODE>string</CODE>) name of the output file.
 * 
 * @attention The elements of the submesh were extracted from the original mesh in order.
 */
void mapFibres(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName,const std::vector<int> &tags,const std::string &fileOut )
{
    std::string outPath;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(mshPath,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read .lon file */
    elem submeshElem = ReadElem(submshPath,submshName);
    lon submshLon = ReadLon(submshPath,submshName);

    /* Read heart.elem for element tag */
    elem mshElem = ReadElem(mshPath,mshName);
    lon mshLon = ReadLon(mshPath,mshName);
    lon lonOut = mshLon;

    /* Based on the assumption that the elements in the original mesh and in the biventricular mesh come in the same order */
    int j = 0;

    for( int i = 0; i < mshElem.tag.size(); i++ )
    {
        if( std::find(tags.begin(), tags.end(), mshElem.tag[i]) != tags.end() )
        {   
            lonOut.x1[i] = submshLon.x1[j];
            lonOut.y1[i] = submshLon.y1[j];
            lonOut.z1[i] = submshLon.z1[j];
            lonOut.x2[i] = submshLon.x2[j];
            lonOut.y2[i] = submshLon.y2[j];
            lonOut.z2[i] = submshLon.z2[j];
            j = j+1;
        }else 
        {
            lonOut.x1[i] = 1.0;
            lonOut.y1[i] = 0.0;
            lonOut.z1[i] = 0.0;
            lonOut.x2[i] = 0.0;
            lonOut.y2[i] = 1.0;
            lonOut.z2[i] = 0.0;
        }
    }

    /* Write output */ 
    std::fstream fOutFibres;
    outPath = mshPath + "/" + fileOut + ".lon";
    write2file(outPath, fOutFibres );
    fOutFibres << 2 << std::endl;

    for( int i = 0; i < lonOut.x1.size(); i++ )
    {
        fOutFibres << lonOut.x1[i] << " " << lonOut.y1[i] << " " << lonOut.z1[i] << " " << lonOut.x2[i] << " " << lonOut.y2[i] << " " << lonOut.z2[i] << std::endl;
    }

    fOutFibres.close();

}

/**
 * @brief takes in input a mesh, a submesh and a .neubc file mapped in the original mesh and maps it to the submesh.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the original mesh. 
 * @param submshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the submesh file.
 * @param submshName (<CODE>string</CODE>) name of the CARP submesh.
 * @param filePath  (<CODE>string</CODE>) path to the <CODE>.neubc</CODE> file to map.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.neubc</CODE> file to map.
 */
void mapNeubcBack(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName,const std::string &filePath,const std::string &fileName )
{   
    std::string outPath;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(mshPath,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read .neubc that has to be mapped */
    neubc neubcIn = ReadNeubc(filePath,fileName);

    /* Initialize output */
    neubc neubcOut = neubcIn;

    /* Map the .neubc back to the submesh */
    #pragma omp parallel for 
    for( int i = 0; i < neubcIn.v1.size(); i++ )
    {
        std::vector<int>::iterator it;
        it = find(vtx.begin(), vtx.end(), neubcIn.v1[i]);
        neubcOut.v1[i] = distance(vtx.begin(),it);

        it = find(vtx.begin(), vtx.end(), neubcIn.v2[i]);
        neubcOut.v2[i] = distance(vtx.begin(),it);

        it = find(vtx.begin(), vtx.end(), neubcIn.v3[i]);
        neubcOut.v3[i] = distance(vtx.begin(),it);

        it = find(vtx.begin(), vtx.end(), neubcIn.v4[i]);
        neubcOut.v4[i] = distance(vtx.begin(),it);
    }

    /* Find the index of the element the triangles belong to in the submesh */ 
    std::fstream inputElemFile;
    inputElemFile.open(submshPath + "/" + submshName + "elem.vtx",std::ios_base::out | std::ios_base::in);

    int nElem;

    if (inputElemFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        std::vector<int> vtxElem = ReadVtx(submshPath,submshName+".elem");
        nElem = vtxElem.size();

        #pragma omp parallel for 
        for( int i = 0; i < neubcIn.v1.size(); i++ )
        {
            std::vector<int>::iterator it;
            it = find(vtxElem.begin(), vtxElem.end(), neubcIn.tet[i]);
            neubcOut.tet[i] = distance(vtxElem.begin(),it);
        }
    } else
    {
        elem submshElem = ReadElem(submshPath,submshName);
        nElem = submshElem.v1.size();

        #pragma omp parallel for 
        for( int i = 0; i < neubcOut.v1.size(); i++ )
        {   
            std::vector<int> temp_ne(4);    

            temp_ne[0] = neubcOut.v1[i];
            temp_ne[1] = neubcOut.v2[i];
            temp_ne[2] = neubcOut.v3[i];
            temp_ne[3] = neubcOut.v4[i];    

            std::sort(temp_ne.begin(),temp_ne.end());   

            #pragma omp parallel for 
            for( int j = 0; j<submshElem.v1.size(); j++ )
            {
                std::vector<int> temp_el(4);    

                temp_el[0] = submshElem.v1[j];
                temp_el[1] = submshElem.v2[j];
                temp_el[2] = submshElem.v3[j];
                temp_el[3] = submshElem.v4[j];  

                std::sort(temp_el.begin(),temp_el.end());   

                if(temp_el[0]==temp_ne[0])
                    if(temp_el[1]==temp_ne[1])
                        if(temp_el[2]==temp_ne[2])
                            if(temp_el[3]==temp_ne[3])
                                neubcOut.tet[i] = j;
            }
        }
    }

    /* Open output .neubc file */
    std::fstream fOutNeubc;
    outPath = submshPath + "/" + fileName + "_mapped_back_in_"+submshName +".neubc";
    write2file(outPath, fOutNeubc );
    fOutNeubc << neubcOut.v1.size() << " # " << "mesh elast " << nElem << " " << vtx.size() << " :" << std::endl;

    for( int i = 0; i < neubcOut.v1.size(); i++ )
    {   
        fOutNeubc << neubcOut.v1[i] << " " << neubcOut.v2[i] << " " << neubcOut.v3[i] << " " << 1 << " " << neubcOut.v4[i] << " " << neubcOut.tet[i] << std::endl;
    }

    fOutNeubc.close();

}


/**
 * @brief Takes in input a mesh, a submesh and a <CODE>.surf</CODE> file mapped in the original mesh and maps it to the submesh.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the original mesh.
 * @param submshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the submesh file.
 * @param submshName (<CODE>string</CODE>) name of the CARP submesh.
 * @param filePath (<CODE>string</CODE>) path to the <CODE>.surf</CODE> file to map.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.surf</CODE> file to map.
 */
void mapSurfBack(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName,const std::string &filePath,const std::string &fileName )
{   
    std::string outPath;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(mshPath,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read .surf file */
    surf surfIn = ReadSurf(filePath,fileName);

    /* Initialize the output */
    surf surfOut = surfIn;

    /* Map the .surf file back to the submesh */
    #pragma omp parallel for 
    for( int i = 0; i < surfIn.v1.size(); i++ )
    {
        int v1 = surfIn.v1[i];
        int v2 = surfIn.v2[i];
        int v3 = surfIn.v3[i];

        for( int j = 0; j < vtx.size(); j++ )
        {
            if( vtx[j] == v1 )
            {
                surfOut.v1[i] = j;
            } else if ( vtx[j] == v2 )
            {
                surfOut.v2[i] = j;
            } else if ( vtx[j] == v3 )
            {
                surfOut.v3[i] = j;
            } else
            {
                continue;
            }
        }
    }

    /* Open output .surf file */
    std::fstream fOutSurf;
    outPath = submshPath + "/" + fileName + "_mapped_back_in_"+submshName +".surf";
    write2file(outPath, fOutSurf );
    fOutSurf << surfOut.v1.size() << std::endl;

    for( int i = 0; i < surfOut.v1.size(); i++ )
    {   
        fOutSurf << "Tr " << surfOut.v1[i] << " " << surfOut.v2[i] << " " << surfOut.v3[i] << std::endl;
    }

    fOutSurf.close();

}

/**
 * @brief Function to map the UVCs from a submesh to the mesh they came from.
 * A value of -10 is assigned to the points outside of the ventricles.
 * It writes 4 files, one per UVC, with the suffix "_mapped".
 * 
 * @param vtxPath Path to the vtx file of the submesh.
 * @param vtxName Vtx file name.
 * @param UVCPath Path to the UVC file.
 * @param UVCName Name of the joint UVC file (with the four UVCs).
 * @param meshPath Path of the mesh we map the UVC to.
 * @param meshName Name of the mesh we map the UVC to.
 * @param outPath Path where the output will be written.
 */

void mapUVCBack(const std::string &vtxPath, const std::string &vtxName,
                const std::string &UVCPath, const std::string &UVCName,
                const std::string &meshPath, const std::string &meshName,
                const std::string &outPath){
    
    std::vector<int> BiV_vtx = ReadVtx(vtxPath,vtxName);
    UVC BiV_UVC = ReadUVC(UVCPath,UVCName);
    points mesh_pts = ReadPts(meshPath, meshName);

    std::vector<double> mapped_Z(mesh_pts.x.size(),-10); 
	std::vector<double> mapped_RHO(mesh_pts.x.size(),-10); 
	std::vector<double> mapped_PHI(mesh_pts.x.size(),-10); 
	std::vector<double> mapped_V(mesh_pts.x.size(),-10); 
    long int j = 0;
    for(long int i = 0; i < mesh_pts.x.size(); i++){
        if(j < BiV_vtx.size()){
            if(BiV_vtx[j] == i){
                mapped_Z[i] = BiV_UVC.Z[j]; 
                mapped_PHI[i] = BiV_UVC.PHI[j]; 
                mapped_RHO[i] = BiV_UVC.RHO[j]; 
                mapped_V[i] = BiV_UVC.V[j];

                j++;
            }
        }        
    }
    std::cout << "\nWriting output...\n";
    std::fstream fOut;

    std::string new_outPath = outPath + "/COORDS_Z_mapped.dat";
    write2file(new_outPath, fOut);
    for( int i = 0; i < mapped_Z.size(); i++ )   
        fOut << mapped_Z[i] << std::endl;
    fOut.close();

    new_outPath = outPath + "/COORDS_PHI_mapped.dat";
    write2file(new_outPath, fOut);
    for( int i = 0; i < mapped_PHI.size(); i++ )   
        fOut << mapped_PHI[i] << std::endl;
    fOut.close();
    
    new_outPath = outPath + "/COORDS_RHO_mapped.dat";
    write2file(new_outPath, fOut);
    for( int i = 0; i < mapped_RHO.size(); i++ )   
        fOut << mapped_RHO[i] << std::endl;
    fOut.close();

    new_outPath = outPath + "/COORDS_V_mapped.dat";
    write2file(new_outPath, fOut);
    for( int i = 0; i < mapped_V.size(); i++ )   
        fOut << mapped_V[i] << std::endl;
    fOut.close();



}

/**
 * @brief Takes in input a mesh, a submesh and a <CODE>.vtx</CODE> file mapped in the original mesh and maps it to the submesh.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the original mesh. 
 * @param submshPath (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the submesh file.
 * @param submshName (<CODE>string</CODE>) name of the CARP submesh.
 * @param filePath (<CODE>string</CODE>) path to the <CODE>.vtx</CODE> file to map.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.vtx</CODE> file to map.
 */
void mapVtxBack(const std::string &mshPath,const std::string &mshName,const std::string &submshPath,const std::string &submshName,const std::string &filePath,const std::string &fileName )
{   
    std::string outPath;
    std::fstream inputFile;
    inputFile.open(submshPath + "/" + submshName + ".vtx",std::ios_base::out | std::ios_base::in);

    std::vector<int> vtx;  

    if (inputFile.is_open()) 
    {   
        /* Read the .vtx file if it exists */
        vtx = ReadVtx(submshPath,submshName);
    } else
    {
        /* Map points to the original mesh */ 
        vtx = mapPoints(mshPath,mshName,submshPath,submshName,submshPath,submshName);
    }

    /* Read the .vtx file to map */
    std::vector<int> vtxIn = ReadVtx(filePath,fileName);

    /* Initialize output */
    std::vector<int> vtxOut = vtxIn;

    /* Map the .vtx file back to the submesh */
    #pragma omp parallel for 
    for( int i = 0; i < vtxIn.size(); i++ )
    {
        for( int j = 0; j < vtx.size(); j++ )
        {
            if( vtx[j] == vtxIn[i] )
            {
                vtxOut[i] = j;
            }else
            {
                continue;
            }
        }
    }

    /* Open output .vtx file */
    std::fstream fOutVtx;
    outPath = submshPath + "/" + fileName + "_mapped_back_in_"+submshName +".vtx";
    write2file(outPath, fOutVtx );
    fOutVtx << vtxOut.size() << std::endl;
    fOutVtx << "intra" << std::endl;

    for( int i = 0; i < vtxOut.size(); i++ )
    {   
        fOutVtx << vtxOut[i] << std::endl;
    }

    fOutVtx.close();

}

/**
 * @brief Locates duplicated nodes with indices and points files.
 * 
 * Takes in input a <CODE>.pts</CODE> file with duplicated nodes and the original <CODE>.pts</CODE> file and generates a map
 * containing the index of the duplicated node and the corresponding index among the non-duplicated nodes.
 * 
 * @param path2NewPts (<CODE>string</CODE>) path to the <CODE>.pts</CODE> containing the duplicated nodes.
 * @param newPtsName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> file containing the duplicated nodes. 
 * @param path2OldPts (<CODE>string</CODE>) path to the <CODE>.pts</CODE> without the duplicated nodes.
 * @param oldPtsName (<CODE>string</CODE>) name of the <CODE>.pts</CODE> without the duplicated nodes.
 * @return (<CODE>mapNodes</CODE> Structure with the indices of the repeated nodes (<CODE>mapNodes.n1</CODE>)
 * and the corresponding index among the non-duplicated nodes (<CODE>mapNodes.n2</CODE>).
 * 
 * @attention It assumes that the duplicated nodes are at the end of the <CODE>.pts</CODE> file. 
 */
mapNodes mapDuplicate(const std::string &path2NewPts,const std::string &newPtsName,const std::string &path2OldPts,const std::string &oldPtsName )
{
    /* Read the .pts files */
    points newPts = ReadPts(path2NewPts,newPtsName);
    points oldPts = ReadPts(path2OldPts,oldPtsName);

    int nPoints    = oldPts.x.size();
    int nDoublePts = newPts.x.size() - nPoints;

    mapNodes mapN;

    /* Map the duplicated nodes to the original nodes */
    for(int i = 0; i<nDoublePts; i++)
    {
        for(int j = 0; j<nPoints; j++)
        {
            if(newPts.x[i+nPoints] == oldPts.x[j])
                if(newPts.y[i+nPoints] == oldPts.y[j])
                    if(newPts.z[i+nPoints] == oldPts.z[j]){
                mapN.n1.push_back(i+nPoints);
                mapN.n2.push_back(j);
                break;
            }
        }
    }

    mapN.nPoints = nPoints;

    return mapN;

}   

/**
 * @brief Creates an activation sequence file without repeated nodes.
 * 
 * This function takes in input an activation sequence generated on a mesh with dubplicated nodes and a map of the duplicated
 * nodes back to the original mesh and generates a new activation sequence by erasing the duplicated node and assigning
 * the remaining ones with the largest activation time. 
 * 
 * @param mapN (<CODE>mapNodes</CODE>) the map of the duplicated nodes back to the original mesh.
 * @param path2ActSeq (<CODE>string</CODE>) path to the activation sequence that has to be modified.
 * @param actSeqName (<CODE>string</CODE>) name of the activation sequence that has to be modified.
 * @param actSeqOutName (<CODE>string</CODE>) name of the output activation sequence. 
 */
void newActSeq(const mapNodes &mapN,const std::string &path2ActSeq,const std::string &actSeqName,const std::string &actSeqOutName )
{   
    std::string outPath;
    /* Read activation sequence */
    std::ifstream actSeqFile(path2ActSeq + "/" + actSeqName + ".dat");   

    if (!actSeqFile.is_open()) 
    {   
        std::cout << "Error opening file " << path2ActSeq << "/" << actSeqName << ".dat" << std::endl;
        exit(0);
    }     

    std::string line;   
    double temp;
    std::vector<double> actSeq;

    while(std::getline(actSeqFile, line))
    {
        std::istringstream lineStream(line);        
        lineStream >> temp;
        actSeq.push_back(temp);
    }

    /* Assign the node with the new activation time */
    std::vector<double> actSeqOut = actSeq;
    actSeqOut.resize(mapN.nPoints);

    for(int i = 0; i < mapN.n1.size(); i++)
    {   
        if(actSeq[mapN.n1[i]] > actSeq[mapN.n2[i]])
        {
            actSeqOut[mapN.n2[i]] = actSeq[mapN.n1[i]];
        }
    }

    /* Write output */
    std::fstream fOutSeq;
    outPath = path2ActSeq + "/" + actSeqOutName + ".dat";
    write2file(outPath, fOutSeq );

    for(int i = 0; i < mapN.nPoints; i++)
    {   
        fOutSeq << actSeqOut[i] << std::endl;
    }

    fOutSeq.close();

}


/*  ------------------------------------------------------------------------------------------------------------

   input -- path2msh (string) path to the .pts of the mesh file
            mshName (string) name of the original mesh 
            path2vtx (string) path to the .vtx file defining the points of the surface
            vtxName (string) name of the .vtx file defining the points of the surface

    This function takes in input a mesh and a .vtx file defining points of a surface and gives in output the 
    .neubc file of the surface

/* ------------------------------------------------------------------------------------------------------------ */
/**
 * @brief Converts the <CODE>.vtx</CODE> of a surface into a <CODE>.neubc</CODE> file.
 * 
 * @param path2msh (<CODE>string</CODE>) path to the <CODE>.pts</CODE> of the mesh file.
 * @param mshName (<CODE>string</CODE>) name of the original mesh.
 * @param path2vtx (<CODE>string</CODE>)  path to the <CODE>.vtx</CODE> file defining the points of the surface.
 * @param vtxName (<CODE>string</CODE>)  name of the <CODE>.vtx</CODE> file defining the points of the surface.
 */
void vtx2neubc(const std::string &path2msh,const std::string &mshName,const std::string &path2vtx,const std::string &vtxName)
{  

    /* Read input mesh */
    points mshPoints = ReadPts(path2msh,mshName);
    elem mshElem = ReadElem(path2msh,mshName);

    /* Read input vtx */
    std::vector<int> vtxSurf = ReadVtx(path2vtx,vtxName);
    std::sort(vtxSurf.begin(),vtxSurf.end());

    /* Initialize neubc elements */
    std::vector<int> v1(mshElem.v1.size()), v2(mshElem.v1.size()), v3(mshElem.v1.size()), v4(mshElem.v1.size()), tetInd(mshElem.v1.size());
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

            v1[count] = tr[0]; v2[count] = tr[1]; v3[count] = tr[2]; v4[count] = diff[0]; tetInd[count] = i;

            count = count + 1;

        }

    }

    v1.resize(count); v2.resize(count); v3.resize(count); v4.resize(count); tetInd.resize(count);

    neubc neubcOut;
    neubcOut.v1 = v1; neubcOut.v2 = v2; neubcOut.v3 = v3; neubcOut.v4 = v4; neubcOut.tet = tetInd;
    std::vector<double> v2v1(3), v3v1(3), v4v1(3), n(3);
    double dotProd;

    /* Adjust normal direction (outward) */
    for( int i = 0; i < v1.size(); i++ )
    {   

        v2v1[0] = mshPoints.x[v2[i]]-mshPoints.x[v1[i]]; v2v1[1] = mshPoints.y[v2[i]]-mshPoints.y[v1[i]]; v2v1[2] = mshPoints.z[v2[i]]-mshPoints.z[v1[i]];
        v3v1[0] = mshPoints.x[v3[i]]-mshPoints.x[v1[i]]; v3v1[1] = mshPoints.y[v3[i]]-mshPoints.y[v1[i]]; v3v1[2] = mshPoints.z[v3[i]]-mshPoints.z[v1[i]];
        v4v1[0] = mshPoints.x[v4[i]]-mshPoints.x[v1[i]]; v4v1[1] = mshPoints.y[v4[i]]-mshPoints.y[v1[i]]; v4v1[2] = mshPoints.z[v4[i]]-mshPoints.z[v1[i]];

        n[0] = v2v1[1]*v3v1[2]-v2v1[2]*v3v1[1]; n[1] = v2v1[2]*v3v1[0]-v2v1[0]*v3v1[2]; n[2] = v2v1[0]*v3v1[1]-v2v1[1]*v3v1[0];
        dotProd = n[0]*v4v1[0] + n[1]*v4v1[1] + n[2]*v4v1[2];

        if( dotProd > 0 )
        {
            neubcOut.v2[i] = v3[i]; neubcOut.v3[i] = v2[i];
        }

    }

    /* Open output .neubc file */
    std::fstream fOutNeubc;
    std::string outName = path2vtx + "/" + vtxName + ".neubc";
    write2file(outName, fOutNeubc);
    fOutNeubc << count << " # " << "mesh elast " << mshElem.v1.size() << " " << mshPoints.x.size() << " :" << std::endl;

    for( int i = 0; i < v1.size(); i++ )
    {   
        fOutNeubc << neubcOut.v1[i] << " " << neubcOut.v2[i] << " " << neubcOut.v3[i] << " " << 1 << " " << neubcOut.v4[i] << " " << neubcOut.tet[i] << std::endl;
    }

    fOutNeubc.close();

}

void map2bivptsdatfile(const std::string& path2datfile, const std::string& datfileName, const std::string& path2vtxfile, const std::string& vtxfileName,const std::string& path2output, const std::string& outputfileName){

    std::vector<double> datfile = ReadActTime(path2datfile,datfileName);
    std::vector<int> vtx = ReadVtx(path2vtxfile,vtxfileName);
    std::vector<double> mapped_dat(datfile.size());

    for(int i = 0; i < vtx.size(); i++)
        mapped_dat[i] = datfile[vtx[i]];

    std::fstream fOut;
    std::string outName = path2output + "/" + outputfileName + ".dat";
    write2file(outName,fOut);

    for(int i = 0; i < mapped_dat.size(); i++)
        fOut << mapped_dat[i] << "\n";

    fOut.close();
}

void map_error_varifold(const std::string &case_num){
    elem heart;
    std::vector<std::vector<double>> mat;
    std::fstream fOut;
    std::string outName;

    // mat =  ReadErrorTable("/data","DeterministicAtlas__EstimatedParameters__Residuals");
    mat =  ReadErrorTable("/data","normalised_residuals");


    heart = ReadElem("/data/aligned_healthy","h_case" + case_num + "algn");
    std::vector<double> elemdata(heart.nElem);

    #pragma omp parallel
    for(int j=0; j < heart.nElem; j++){
        elemdata[j] = mat[std::stoi(case_num)][heart.tag[j]] ;
    }
    
    outName = "/data/aligned_healthy/h_case" + case_num + "algn_varifolderror.dat";
    write2file(outName,fOut);

    for(int j=0; j<elemdata.size(); j++)
        fOut << elemdata[j] << "\n";

    fOut.close();

}