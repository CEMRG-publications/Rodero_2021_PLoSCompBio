/**
 * @brief Main fuction to set up the EP simulation.
 * 
 * @file EP.cpp
 * @author Marina Strocchi
 * @author Cristobal Rodero
 * @date 2018
 */

#include "../headers/EP.h"

/**
 * @brief Creates a spherical electrode on a mesh.
 * 
 * 	The function takes in input the CARP mesh points and the index of the center of the electrode and saves a
 * <CODE>.vtx</CODE> file containing the nodes defining a spherical electrode of radius 5 mm.
 * 
 * @param path (<CODE>string</CODE>) path to the CARP mesh.
 * @param mshName (<CODE>string</CODE>) name of the CARP mesh.
 * @param IDnode (<CODE>int</CODE>) index of the node defining the center of the electrode.
 * @param outFile (<CODE>string</CODE>) name of the output <CODE>.vtx</CODE> file.
 */

void electrodeVtx(const std::string &path,const std::string &mshName,const int &IDnode,const std::string &outFile)
{

	std::string outPath;

	/* Read .pts file */
	points pts = ReadPts( path,mshName );

	double xC = pts.x[IDnode];
	double yC = pts.y[IDnode];
	double zC = pts.z[IDnode];

	double r = 2000.0; /* radius of the electrode in micrometers */ 

	/* store indices of the nodes in a vector */
	std::vector<int> vtx;

	for( int i = 0; i < pts.x.size(); i++ )
	{	
		double temp = pow( pts.x[i] - xC,2 ) + pow( pts.y[i] - yC,2 ) + pow( pts.z[i] - zC,2 );

		if( temp <= pow(r,2) ) 
		{
			vtx.push_back(i);
		}
	}

	/* Write output file */ 
    std::fstream fOutVtx;
	outPath = path + "/" + outFile + ".vtx";
	write2file(outPath, fOutVtx );
    fOutVtx << vtx.size() << std::endl;
    fOutVtx << "intra" << std::endl;
    
    std::fstream fOutPts;
	outPath =  path + "/" + outFile + ".pts";
    write2file(outPath, fOutPts );
    fOutPts << vtx.size() << std::endl;

    for( int i = 0; i < vtx.size(); i++ )
	{	
		fOutVtx << vtx[i] << std::endl;
		fOutPts << pts.x[vtx[i]] << " " << pts.y[vtx[i]] << " "<< pts.z[vtx[i]] << std::endl;
	}

	fOutVtx.close();
	fOutPts.close();

}
/**
 * @brief Creates a circular electrode in the endocardium.
 * 
 * The function takes in input the CARP mesh points, the list of nodes of a surface and the index of the center 
 * of the electrode and saves a <CODE>.vtx</CODE> file containing the nodes defining a surface circular electrode of radius 
 * <VAR>r</VAR> mm on the surface given in input.
 * 
 * @param path (<CODE>string</CODE>) path to the CARP mesh.
 * @param mshName (<CODE>string</CODE>) name of the CARP mesh.
 * @param surfaceVtx (<CODE>int</CODE>) name of the <CODE>.vtx</CODE> file defining the surface the electrode needs to be defined on.
 * @param IDnode (<CODE>int</CODE>) Index of the centre of the electrode. 
 * @param outFile (<CODE>string</CODE>) name of the output <CODE>.vtx</CODE> file.
 * @param r (<CODE>double</CODE>) Radius of the electrode.
 */

void electrodeEndoVtx(const std::string &path, const std::string &mshName, const std::string &surfaceVtx, const int &IDnode , const std::string &outFile, const double &r)
{

	std::string outPath;

	/* Read .pts file */
	points pts = ReadPts( path,mshName );
	std::vector<int> vtx = ReadVtx( path,surfaceVtx );

	double xC = pts.x[IDnode];
	double yC = pts.y[IDnode];
	double zC = pts.z[IDnode];

	/* store indices of the nodes in a vector */
	std::vector<int> vtxOut;

	for( int i = 0; i < vtx.size(); i++ )
	{	
		double temp = pow( pts.x[vtx[i]] - xC,2 ) + pow( pts.y[vtx[i]] - yC,2 ) + pow( pts.z[vtx[i]] - zC,2 );

		if( temp <= pow(r,2) ) 
		{
			vtxOut.push_back(vtx[i]);
		}
	}

	/* Write output file */ 
	std::fstream fOutVtx;
	outPath = path + "/" + outFile + ".vtx";
	write2file(outPath, fOutVtx );
    fOutVtx << vtxOut.size() << std::endl;
    fOutVtx << "intra" << std::endl;
    
    std::fstream fOutPts;
	outPath = path + "/" + outFile + ".pts";
    write2file(outPath, fOutPts );
    fOutPts << vtxOut.size() << std::endl;

    for( int i = 0; i < vtxOut.size(); i++ )
	{	
		fOutVtx << vtxOut[i] << std::endl;
		fOutPts << pts.x[vtxOut[i]] << " " << pts.y[vtxOut[i]] << " "<< pts.z[vtxOut[i]] << std::endl;
	}

	fOutVtx.close();
	fOutPts.close();

}

/**
 * @brief Takes in input a CARP mesh and the activation times on nodes and converts it into a <CODE>.vtk</CODE> file.
 * 
 * @param mshPath (<CODE>string</CODE>) path to the CARP mesh.
 * @param mshName (<CODE>string</CODE>) name of the CARP mesh.
 * @param pathFile (<CODE>string</CODE>) path to the file for the activation times.
 * @param fileName (<CODE>string</CODE>) name of the <CODE>.dat</CODE> file containing activation times.
 */
void vtkActivationTime(const std::string &mshPath,const std::string &mshName,const std::string &pathFile,const std::string &fileName )
{	
	std::string outPath;

	/* Read .dat file */ 
	std::vector<double> actTime = ReadActTime( pathFile, fileName );

	/* Read mesh file */ 
	points points = ReadPts( mshPath, mshName );
	elem elem = ReadElem( mshPath, mshName );

	/* Open .vtk file output */
	std::fstream fOutVtk;
	outPath = pathFile + "/" + mshName + ".vtk";
	write2file(outPath , fOutVtk );

	/* Header */
    fOutVtk << "# vtk DataFile Version 3.0" << std::endl;
    fOutVtk << "vtk output" << std::endl;
	fOutVtk << "ASCII" << std::endl;
	fOutVtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fOutVtk << " " << std::endl;

	fOutVtk << "POINTS " << points.x.size() << " float" << std::endl;

	for( int i = 0; i < points.x.size(); i++ )
	{
		fOutVtk << points.x[i] << " " << points.y[i] << " " << points.z[i] << std::endl;
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
	fOutVtk << "POINT_DATA " << points.x.size() << std::endl;
	fOutVtk << "SCALARS activation_time float 1" << std::endl;
	fOutVtk << "LOOKUP_TABLE default" << std::endl;
	for( int i = 0; i < actTime.size(); i++ )
	{
		if( actTime[i] > 600.0)
		{
			fOutVtk << 0.0 << std::endl;
		}else
		{
			fOutVtk << actTime[i] << std::endl;
		}
	}  	

	fOutVtk << " " << std::endl;
	fOutVtk << "CELL_DATA " << elem.v1.size() << std::endl;
	fOutVtk << "SCALARS tag float 1" << std::endl;
	fOutVtk << "LOOKUP_TABLE default" << std::endl;
	for( int i = 0; i < elem.tag.size(); i++ )
	{
		fOutVtk << elem.tag[i] << std::endl;
	}	

	fOutVtk.close();

}


/**
 * @brief Adds an additional tag for the fast endocardial conduction layer given in input the <CODE>.vtx</CODE> of the endocardial surface.
 * 
 * @param path2msh (<CODE>string</CODE>) path to the CARP mesh.
 * @param elemName (<CODE>string</CODE>) name of the CARP mesh.
 * @param path2vtx (<CODE>string</CODE>) path to the <CODE>.vtx</CODE> file defining the endocardial surface.
 * @param vtxName (<CODE>string</CODE>) name of the <CODE>.vtx</CODE> file defining the endocardial surface.
 * @param tag_endo (<CODE>int</CODE>) 1 if making the FEC layer of the LV, 2 if making the RV one.
 * @param path2phie (<CODE>string</CODE>) path to the Laplace solution file (epi to endo).
 * @param phieName (<CODE>string</CODE>) name of the Laplace solution file.
 * @param percVol (<CODE>double</CODE>) percentage of the myocardial volume to be assigned to the FEC layer.
 * @param fileOut (<CODE>string</CODE>) name of the <CODE>.elem</CODE> output file containing the new label for the FEC.
 */
void endoLayer(const std::string &path2msh,const std::string &elemName,const std::string &path2vtx,const std::string &vtxName,const int &tag_endo,const std::string &path2phie,const std::string &phieName,const double &percVol,const std::string &fileOut )
{

	std::string outPath;
	/* Read input file */
	
	std::vector<int> vtx = ReadVtx( path2vtx, vtxName );
	elem elemMsh = ReadElem( path2msh, elemName );
	std::vector<double> laplace = ReadActTime(path2phie,phieName);

	std::vector<int>::iterator maxTag = std::max_element( elemMsh.tag.begin(), elemMsh.tag.end() );

	/* Initialize output */
	elem elemOut = elemMsh;
	std::vector<int>inner_nodes (laplace.size(),-1);

	// Check which nodes are the inner nodes.
	std::cout<<"Finding inner nodes...\n";


	#pragma omp parallel for
	for(long int node = 0; node < laplace.size(); node++){
		if(laplace[node] <= percVol)
			inner_nodes[node] = vtx[node];
	}

	inner_nodes.erase(std::remove(inner_nodes.begin(),inner_nodes.end(),-1),inner_nodes.end()); // Remove all the -1 so we have a smaller vector.

	// Check which are the inner elements.
	std::cout<<"Finding inner elements...\n";
	
	#pragma omp parallel for
	for(long int el = 0; el < elemMsh.v1.size(); el++){
		if(elemMsh.tag[el] == tag_endo){ // We are in the correct tag
			if(std::find(inner_nodes.begin(), inner_nodes.end(), elemMsh.v1[el]) != inner_nodes.end()){
				if(std::find(inner_nodes.begin(), inner_nodes.end(), elemMsh.v2[el]) != inner_nodes.end()){
					if(std::find(inner_nodes.begin(), inner_nodes.end(), elemMsh.v3[el]) != inner_nodes.end()){
						if(std::find(inner_nodes.begin(), inner_nodes.end(), elemMsh.v4[el]) != inner_nodes.end())
							elemOut.tag[el] = maxTag[0] + 1;
					}
				}
			}
		}
	}

	std::cout<<"Writing output...\n";
	/* Write output */
	std::fstream fOutElem;
	outPath = path2msh + "/" + fileOut + ".elem";
	write2file(outPath, fOutElem );
	fOutElem << elemOut.v1.size() << std::endl;

	for(long int i = 0; i < elemOut.v1.size(); i++ )
	{
		fOutElem << "Tt " << elemOut.v1[i] << " " << elemOut.v2[i] << " " << elemOut.v3[i] << " " << elemOut.v4[i] << " " << elemOut.tag[i] << std::endl;
	}

	fOutElem.close();

}
/**
 * @brief Create a virtual thread of equispaced electrodes (one point per electrode).
 * 
 * @param path2mesh (<CODE>string</CODE>) Path to the full mesh.
 * @param mshName (<CODE>string</CODE>) Name of the <CODE>.elem</CODE> file of the mesh.
 * @param surfaceVtx (<CODE>string</CODE>) Name of the <CODE>.vtx</CODE> file of the surface where to put the thread.
 * @param elec_centre (<CODE>int</CODE>) Index of the first electrode.
 * @param elec_end (<CODE>int</CODE>) Index indicating the direction where the electrode thread should follow.
 * @param elec_distance (<CODE>double</CODE>) Inter-electrode distance.
 * @param num_electrodes (<CODE>int</CODE>) Number of electrodes in the thread.
 */
void createElectrodesThread(const std::string &path2mesh, const std::string &mshName, const std::string &surfaceVtx, const int &elec_centre, const int &elec_end, const double &elec_distance, const int &num_electrodes,const std::string& outPath, const std::string & outName){
	
	points pts = ReadPts(path2mesh, mshName);
	std::vector<int> vtx = ReadVtx(path2mesh,surfaceVtx);
	std::vector<double> distance_vec(vtx.size());
	std::vector<int> elec_centres;
	elec_centres;

	std::vector<double> dir_vec, ideal_pos;

	std::fstream fOutVtx;
	std::fstream fOutPts;

	std::string aux = outPath + "/" + outName + ".vtx";
	write2file(aux, fOutVtx );
	aux = outPath + "/" + outName + ".vtk";
	write2file(aux, fOutPts );

	fOutVtx << num_electrodes << std::endl << "intra" << std::endl;
	fOutPts << "# vtk DataFile Version 3.0 \nvtk output \nASCII \nDATASET UNSTRUCTURED_GRID \n\nPOINTS " << num_electrodes << " float\n";


	elec_centres.push_back(elec_centre);

	for(int i=0; i<num_electrodes; i++){

		if(i<1)
			dir_vec = pts2vec(pts,elec_end) - pts2vec(pts,elec_centres[0]);
		else
			dir_vec = pts2vec(pts,elec_centres[i]) - pts2vec(pts,elec_centres[i-1]);
			
		ideal_pos = pts2vec(pts,elec_centres[i]) + (elec_distance/norm(dir_vec))*dir_vec;

		
		for(int j = 0; j<vtx.size(); j++)
			distance_vec[j] = norm(abs(pts2vec(pts,vtx[j])-ideal_pos));

		elec_centres.push_back(vtx[min_pos(distance_vec)]);

		// std::cout << elec_centres[i] << std::endl;
	
		fOutVtx << std::to_string(elec_centres[i]) << std::endl;
		fOutPts << pts.x[elec_centres[i]] << " " << pts.y[elec_centres[i]] << " " << pts.z[elec_centres[i]] << std::endl;
	}

		fOutVtx.close();	
		fOutPts.close();	

}
/**
 * @brief Creates the bipole and tripole activation times files.
 * 
 * Reading the activation times given by CARP of the individual activation of the electrodes, it gathers them in the same folder
 * and creates the tripole ones.
 * 
 * @param path (<CODE>string</CODE>) Path to the current case.
 * @param commonName (<CODE>string</CODE>) Partial path between the case path and the electrode activation time folder.
 * @param num_electrodes (<CODE>int</CODE>) Number of electrodes of the thread.
 */

void createQuadpoles(const std::string &path, const std::string &commonName, const int &num_electrodes){

	std::string outPath;

	system(("mkdir -p " + path + "/results/AT").c_str());

	std::vector<std::string> comb_vec;

	for(int el_1 = 1; el_1 <= num_electrodes-1; el_1++){
		for(int el_2 = el_1+1; el_2 <= num_electrodes; el_2++){
					comb_vec.push_back(std::to_string(el_1) + std::to_string(el_2));
		}
	}

	std::vector<std::vector<double>> actTimes(comb_vec.size());

	for(int file_num = 1; file_num <= num_electrodes; file_num++){
		system(("cp " + path + "/" + commonName + std::to_string(file_num) + "/vm_act_seq.dat " + path + "/results/AT/" + std::to_string(file_num) + ".dat").c_str());
		actTimes[file_num - 1] = ReadActTime(path + "/results/AT", std::to_string(file_num));
	}

	std::fstream fOutDat;


	std::vector<double> v1, v2, v_min;

	for(int comb = 0; comb < comb_vec.size(); comb++){

		outPath = path + "/results/AT/" + comb_vec.at(comb) + ".dat";
		write2file(outPath,fOutDat);
		v1 = actTimes.at((int)((comb_vec[comb]).at(0))-'0'-1);
		v2 = actTimes.at((int)((comb_vec[comb]).at(1))-'0'-1);

		v_min = min_positive(v1,v2);

		for (int i = 0; i < v_min.size(); i++)
			fOutDat << v_min.at(i) << std::endl;
		
		
		fOutDat.close();
	}
	
}

void mergeTwoElectrodes(const std::vector<std::string> &paths, const std::vector<std::string> &electrodes, const std::string &outpath, const std::string &outname){

	system(("mkdir -p " + outpath).c_str());

	std::vector<std::vector<double>> actTimes(2);

	for(int file_num = 0; file_num < paths.size(); file_num++){
		actTimes[file_num] = ReadActTime(paths[file_num], electrodes[file_num]);
	}

	std::fstream fOutDat;

	std::vector<double> v1, v2, v_min;
	std::string outstring = outpath + "/" + outname + ".dat";
		write2file(outstring, fOutDat);

		v1 = actTimes.at(0);
		v2 = actTimes.at(1);

		v_min = min_positive(v1,v2);

		for (int i = 0; i < v_min.size(); i++)
			fOutDat << v_min.at(i) << std::endl;
		
		
		fOutDat.close();
	}

/**
 * @brief Removes the fast endocardial conduction tags.
 * 
 * @param path2msh (<CODE>string</CODE>) Path to the mesh.
 * @param elemName (<CODE>string</CODE>) Name of the <CODE>.elem<CODE> file of the mesh.
 * @param fileOut (<CODE>string</CODE>) Name of the new mesh without the tags.
 * @param tag_FEC_LV (<CODE>int</CODE>) Tag of the LV FEC.
 * @param tag_FEC_RV (<CODE>int</CODE>) Tag of the RV FEC.
 * 
 * @attention The new LV FEC tag will be 1 (LV) and 2 for the RV FEC (RV).
 */
void removeFEC(const std::string &path2msh,const std::string &elemName,const std::string &fileOut,const int &tag_FEC_LV,const int &tag_FEC_RV){
	
	std::string outPath;
	elem elemMsh = ReadElem( path2msh, elemName );
	std::vector<int> septum_vtx = ReadVtx(path2msh,"BiV.rvsept.surf");

	elem elemOut = elemMsh;



	
	#pragma omp parallel for
	for(long int el = 0; el < elemMsh.v1.size(); el++){
		if(elemMsh.tag[el] == tag_FEC_RV) // We are in the correct tag
			if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v1[el]) != septum_vtx.end())
					elemOut.tag[el] = 1;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v2[el]) != septum_vtx.end())
					elemOut.tag[el] = 1;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v3[el]) != septum_vtx.end())
					elemOut.tag[el] = 1;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v4[el]) != septum_vtx.end())
					elemOut.tag[el] = 1;
		}

	#pragma omp parallel for
	for(long int i = 0; i < elemMsh.v1.size(); i++ )
	{
		if(elemMsh.tag[i] == tag_FEC_RV)
		{	
			if( (std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v1[i]) != septum_vtx.end())||(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v2[i]) != septum_vtx.end()) ||
				(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v3[i]) != septum_vtx.end())||(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v4[i]) != septum_vtx.end()) ){
				elemOut.tag[i] = 1;
				}
		}
	}

	#pragma omp parallel for
	for (long int el = 0; el < elemOut.v1.size(); el++)
	{
		if(elemOut.tag[el] == tag_FEC_LV)
			elemOut.tag[el] = 1;
		
		if (elemOut.tag[el] == tag_FEC_RV)
			elemOut.tag[el] = 2;
			
	}

	std::fstream fOutElem;
	outPath = path2msh + "/" + fileOut + ".elem";
	write2file(outPath, fOutElem );
	fOutElem << elemOut.v1.size() << std::endl;

	for(long int i = 0; i < elemOut.v1.size(); i++ )
		fOutElem << "Tt " << elemOut.v1[i] << " " << elemOut.v2[i] << " " << elemOut.v3[i] << " " << elemOut.v4[i] << " " << elemOut.tag[i] << std::endl;

	fOutElem.close();

}
/**
 * @brief Creates a fast endocardial conduction layer as an endocardial shell.
 * 
 * Creates a fast endocardial conduction layer in the endocardia with a specified maximum thickness and height w.r.t. to the myocardial
 * dimensions. It puts one element at least.
 * 
 * @param path2msh (<CODE>string</CODE>) Path to the mesh.
 * @param elemName (<CODE>string</CODE>) Name of the <CODE>.elem</CODE> file.
 * @param path2epi2endo (<CODE>string</CODE>) Path to the <CODE>phie.dat</CODE> file corresponding to the Laplace of solucion of the epicardium.
 * @param epi2endoName (<CODE>string</CODE>) Name of the <CODE>phie.dat</CODE> file corresponding to the Laplace of solucion of the epicardium.
 * @param path2apba (<CODE>string</CODE>) Path to the <CODE>phie.dat</CODE> file corresponding to the Laplace of solucion of the apex to base.
 * @param apbaName (<CODE>string</CODE>) Name of the <CODE>phie.dat</CODE> file corresponding to the Laplace of solucion of the apex to base.
 * @param max_thick (<CODE>double</CODE>) Percentage of the myocardia thickness to assign to the FEC layer.
 * @param max_height (<CODE>double</CODE>) Percentage of the myocardia height to assign to the FEC layer (apex to base). 
 */
void createFEC(const std::string &path2msh,const std::string &elemName, const std::string &path2epi2endo,const std::string &epi2endoName, const std::string &path2apba, const std::string &apbaName, const double &max_thick, const double &max_height)
{
	// std::string outPath;
	/* Read input file */
	
	// std::vector<int> LV_vtx = ReadVtx( path2msh, "LV_endo.surfmesh" );
	// std::vector<int> RV_vtx = ReadVtx( path2msh, "RV_endo.surfmesh" );
	// std::vector<int> septum_vtx = ReadVtx(path2msh,"septum.surfmesh");

	std::vector<int> LV_vtx = ReadVtx( path2msh, "BiV.lvendo.surf" );
	std::vector<int> RV_vtx = ReadVtx( path2msh, "BiV.rvendo.surf" );
	std::vector<int> septum_vtx = ReadVtx(path2msh,"BiV.rvsept.surf");

	elem elemMsh = ReadElem( path2msh, elemName );

	scale_UVC_file(path2apba, apbaName);

	std::vector<double> epi2endo = ReadActTime(path2epi2endo,epi2endoName);
	std::vector<double> apba = ReadActTime(path2apba, apbaName + "_scaled");


	std::vector<int>::iterator maxTag = std::max_element( elemMsh.tag.begin(), elemMsh.tag.end() );
	
	/********************
	 * LV FEC
	 * *******************/

	/* Initialize output */
	elem elemOut = elemMsh;
	std::vector<int>inner_nodes_LV (epi2endo.size(),-1);


	// Check which nodes are the inner nodes.
	std::cout<<"Finding inner nodes of the LV...\n";

	/*
	std::cout << "\n Epi2endo size : " <<epi2endo.size() << "\n apba size: " << apba.size() << "\n inner nodes LV: " << inner_nodes_LV.size() << "\n lv vtx: " << LV_vtx.size();
	*/


	#pragma omp parallel for
	// for(long int node = 0; node < epi2endo.size(); node++){
	for(long int node = 0; node < LV_vtx.size(); node++){
		if((epi2endo[node] <= max_thick) && (apba[node] <= max_height))
			inner_nodes_LV[node] = LV_vtx[node];
	}

	inner_nodes_LV.erase(std::remove(inner_nodes_LV.begin(),inner_nodes_LV.end(),-1),inner_nodes_LV.end()); // Remove all the -1 so we have a smaller vector.

	// Check which are the inner elements.
	std::cout<<"Finding inner elements of the LV...\n";
	
	#pragma omp parallel for
	for(long int el = 0; el < elemMsh.v1.size(); el++){
		if(elemMsh.tag[el] == 1){ // We are in the correct tag
			if(std::find(inner_nodes_LV.begin(), inner_nodes_LV.end(), elemMsh.v1[el]) != inner_nodes_LV.end()){
				if(std::find(inner_nodes_LV.begin(), inner_nodes_LV.end(), elemMsh.v2[el]) != inner_nodes_LV.end()){
					if(std::find(inner_nodes_LV.begin(), inner_nodes_LV.end(), elemMsh.v3[el]) != inner_nodes_LV.end()){
						if(std::find(inner_nodes_LV.begin(), inner_nodes_LV.end(), elemMsh.v4[el]) != inner_nodes_LV.end())
							elemOut.tag[el] = maxTag[0] + 1;
					}
				}
			}
		}
	}

		/* Add one element in the most inner part */

	std::cout<<"Adding one layer in the LV...\n";

	#pragma omp parallel for
	for(long int i = 0; i < elemMsh.v1.size(); i++ )
	{
		if(elemMsh.tag[i] == 1)
		{	
			if( (std::find(LV_vtx.begin(), LV_vtx.end(), elemMsh.v1[i]) != LV_vtx.end())||(std::find(LV_vtx.begin(), LV_vtx.end(), elemMsh.v2[i]) != LV_vtx.end()) ||
				(std::find(LV_vtx.begin(), LV_vtx.end(), elemMsh.v3[i]) != LV_vtx.end())||(std::find(LV_vtx.begin(), LV_vtx.end(), elemMsh.v4[i]) != LV_vtx.end()) )
			{	
				if((apba[elemMsh.v1[i]] <= max_height) || (apba[elemMsh.v2[i]] <= max_height) || (apba[elemMsh.v3[i]] <= max_height) || (apba[elemMsh.v4[i]] <= max_height))	
					elemOut.tag[i] = maxTag[0] + 1;
			}
		}
	}

	/********************
	 * RV FEC
	 * *******************/

	std::vector<int>inner_nodes_RV (epi2endo.size(),-1);


	// Check which nodes are the inner nodes.
	std::cout<<"Finding inner nodes of the RV...\n";


	#pragma omp parallel for
	// for(long int node = 0; node < epi2endo.size(); node++){
	for(long int node = 0; node < RV_vtx.size(); node++){
		if((epi2endo[node] <= max_thick) && (apba[node] <= max_height))
			inner_nodes_RV[node] = RV_vtx[node];
	}

	inner_nodes_RV.erase(std::remove(inner_nodes_RV.begin(),inner_nodes_RV.end(),-1),inner_nodes_RV.end()); // Remove all the -1 so we have a smaller vector.

	// Check which are the inner elements.
	std::cout<<"Finding inner elements of the RV...\n";
	
	#pragma omp parallel for
	for(long int el = 0; el < elemMsh.v1.size(); el++){
		if(elemMsh.tag[el] == 2){ // We are in the correct tag
			if(std::find(inner_nodes_RV.begin(), inner_nodes_RV.end(), elemMsh.v1[el]) != inner_nodes_RV.end()){
				if(std::find(inner_nodes_RV.begin(), inner_nodes_RV.end(), elemMsh.v2[el]) != inner_nodes_RV.end()){
					if(std::find(inner_nodes_RV.begin(), inner_nodes_RV.end(), elemMsh.v3[el]) != inner_nodes_RV.end()){
						if(std::find(inner_nodes_RV.begin(), inner_nodes_RV.end(), elemMsh.v4[el]) != inner_nodes_RV.end())
							elemOut.tag[el] = maxTag[0] + 2;
					}
				}
			}
		}
	}

	// We add the septum.
	std::cout<<"Adding the septum...\n";


	
	#pragma omp parallel for
	for(long int el = 0; el < elemMsh.v1.size(); el++){
		if(elemMsh.tag[el] == 1) // We are in the correct tag
			if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v1[el]) != septum_vtx.end())
				if((apba[elemMsh.v1[el]] <= max_height) && (apba[elemMsh.v2[el]] <= max_height) && (apba[elemMsh.v3[el]] <= max_height) && (apba[elemMsh.v4[el]] <= max_height))
					elemOut.tag[el] = maxTag[0] + 2;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v2[el]) != septum_vtx.end())
				if((apba[elemMsh.v1[el]] <= max_height) && (apba[elemMsh.v2[el]] <= max_height) && (apba[elemMsh.v3[el]] <= max_height) && (apba[elemMsh.v4[el]] <= max_height))
					elemOut.tag[el] = maxTag[0] + 2;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v3[el]) != septum_vtx.end())
				if((apba[elemMsh.v1[el]] <= max_height) && (apba[elemMsh.v2[el]] <= max_height) && (apba[elemMsh.v3[el]] <= max_height) && (apba[elemMsh.v4[el]] <= max_height))
					elemOut.tag[el] = maxTag[0] + 2;
			else if(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v4[el]) != septum_vtx.end())
				if((apba[elemMsh.v1[el]] <= max_height) && (apba[elemMsh.v2[el]] <= max_height) && (apba[elemMsh.v3[el]] <= max_height) && (apba[elemMsh.v4[el]] <= max_height))
					elemOut.tag[el] = maxTag[0] + 2;
		}
			
			
		/* Add one element in the most inner part */

	std::cout<<"Adding one layer in the RV...\n";

	#pragma omp parallel for
	for(long int i = 0; i < elemMsh.v1.size(); i++ )
	{
		if(elemMsh.tag[i] == 2)
		{	
			if( (std::find(RV_vtx.begin(), RV_vtx.end(), elemMsh.v1[i]) != RV_vtx.end())||(std::find(RV_vtx.begin(), RV_vtx.end(), elemMsh.v2[i]) != RV_vtx.end()) ||
				(std::find(RV_vtx.begin(), RV_vtx.end(), elemMsh.v3[i]) != RV_vtx.end())||(std::find(RV_vtx.begin(), RV_vtx.end(), elemMsh.v4[i]) != RV_vtx.end()) ){
				if((apba[elemMsh.v1[i]] <= max_height) && (apba[elemMsh.v2[i]] <= max_height) && (apba[elemMsh.v3[i]] <= max_height) && (apba[elemMsh.v4[i]] <= max_height))
				elemOut.tag[i] = maxTag[0] + 2;
				}
		}
	}

	std::cout<<"Adding one layer in the septum...\n";

	#pragma omp parallel for
	for(long int i = 0; i < elemMsh.v1.size(); i++ )
	{
		if(elemMsh.tag[i] == 1)
		{	
			if( (std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v1[i]) != septum_vtx.end())||(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v2[i]) != septum_vtx.end()) ||
				(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v3[i]) != septum_vtx.end())||(std::find(septum_vtx.begin(), septum_vtx.end(), elemMsh.v4[i]) != septum_vtx.end()) ){
				if((apba[elemMsh.v1[i]] <= max_height) && (apba[elemMsh.v2[i]] <= max_height) && (apba[elemMsh.v3[i]] <= max_height) && (apba[elemMsh.v4[i]] <= max_height))
				elemOut.tag[i] = maxTag[0] + 2;
				}
		}
	}



	std::cout<<"Writing output...\n";
	/* Write output */
	std::fstream fOutElem;
	std::string outPath = path2msh + "/" + elemName + "_FEC_w" + std::to_string((int)round(100*max_thick)) + "_h" + std::to_string((int)round(100*max_height)) + ".elem";

	fOutElem.open(outPath,std::ios::out);
	fOutElem << elemOut.v1.size() << std::endl;

	for(long int i = 0; i < elemOut.v1.size(); i++ )
		fOutElem << "Tt " << elemOut.v1[i] << " " << elemOut.v2[i] << " " << elemOut.v3[i] << " " << elemOut.v4[i] << " " << elemOut.tag[i] << std::endl;

	fOutElem.close();

}


std::vector<double> scale_UVC(UVC BiV_UVC, std::string coordinate){


	struct wrt{
        const std::vector<double> & value_vector;

        wrt(const std::vector<double> & val_vec):
            value_vector(val_vec) {}

        bool operator()(int i1, int i2)
        {
            return value_vector[i1] < value_vector[i2];
        }
    };
	
	long int n_pts = BiV_UVC.Z.size();
	std::vector<double> indices(n_pts);
	std::vector<double> act_seq(n_pts);
	std::vector<double> coordinate_sorted(n_pts);

	#pragma omp parallel
	for(long int i = 0; i < n_pts; i++)
		indices[i] = i;
	
	std::vector<double> indices_2 = indices;
    std::vector<double>::iterator maxcoor, mincoor;

	if(coordinate == "Z"){
		act_seq = BiV_UVC.Z;
		maxcoor = std::max_element(BiV_UVC.Z.begin(), BiV_UVC.Z.end());
    	mincoor = std::min_element(BiV_UVC.Z.begin(), BiV_UVC.Z.end());
	}
	if(coordinate == "PHI"){
		act_seq = BiV_UVC.PHI;
		maxcoor = std::max_element(BiV_UVC.PHI.begin(), BiV_UVC.PHI.end());
    	mincoor = std::min_element(BiV_UVC.PHI.begin(), BiV_UVC.PHI.end());
	}


    std::sort(indices.begin(), indices.end(), wrt(act_seq));
	std::sort(indices_2.begin(), indices_2.end(), wrt(indices));
	
	for(int i=0; i < indices_2.size()-1; i++)
		coordinate_sorted[i] = mincoor[0] + (maxcoor[0]-mincoor[0])*indices_2[i]/double(n_pts-1);

	return coordinate_sorted;
}

void scale_UVC_file(const std::string &path, const std::string &filename){


	struct wrt{
        const std::vector<double> & value_vector;

        wrt(const std::vector<double> & val_vec):
            value_vector(val_vec) {}

        bool operator()(int i1, int i2)
        {
            return value_vector[i1] < value_vector[i2];
        }
    };

	std::vector<double> act_seq = ReadActTime(path,filename);
	
	long int n_pts = act_seq.size();

	std::vector<double> indices(n_pts);

	#pragma omp parallel
	for(long int i = 0; i < n_pts; i++)
		indices[i] = i;
	
	std::vector<double> indices_2 = indices;	
    std::sort(indices.begin(), indices.end(), wrt(act_seq));
	std::sort(indices_2.begin(), indices_2.end(), wrt(indices));

	std::fstream fOut(path + "/" + filename + "_scaled.dat", std::fstream::out);
	for(long int i = 0; i < act_seq.size(); i++)
		fOut << indices_2[i]/double(n_pts-1) << std::endl;

	fOut.close();
	
}

int single_electrode_UVC(const UVC &BiV_UVC, const std::vector<double> &z_scaled, const double &apba, const double &rho_goal, const double &rotational, const double &V_goal, const double &apba_tol, const double &rotational_tol){


    int i;
	std::vector<double> phi_scaled = scale_UVC(BiV_UVC,"PHI");
	// std::vector<double> z_scaled = scale_UVC(BiV_UVC,"Z");
	bool in_septum = true, nothing_worked = true;
	
	if(V_goal > 0){
    	for(i=0; i < z_scaled.size(); i++){
			if(BiV_UVC.V[i] == 1){
				if(BiV_UVC.RHO[i] == rho_goal){
					if(std::abs(z_scaled[i] - apba) < apba_tol){
						if(std::abs(phi_scaled[i] - rotational) < rotational_tol){
							nothing_worked = false;
							in_septum = false;
							break;
						}
                    }
                }
            }
        }
    }
	
	if(V_goal == -1 || in_septum == true){
    	for(i=0; i <z_scaled.size(); i++){
			if(BiV_UVC.V[i] == -1){
				if(BiV_UVC.RHO[i] == rho_goal){
					if(std::abs(z_scaled[i] - apba) < apba_tol){
						if(std::abs(phi_scaled[i] - rotational) < rotational_tol){
							nothing_worked = false;
							break;
						}
                    }
                }
            }
        }
	}
	
	if(nothing_worked){	
		std::cout << "Increasing tolerance...\n";
		single_electrode_UVC(BiV_UVC, z_scaled, apba, rho_goal, rotational, V_goal, 1.1*apba_tol, 1.1*rotational_tol);
	}
	return i;
}

void lead2electrodes(const std::string &path2leads, const std::string &leadName, const std::string &outPath, const std::string &heart){

	std::vector<int> lead;
	std::fstream fOutVtx;
	std::string aux;

	lead = ReadVtx(path2leads,leadName);

	for(int i = 0; i < lead.size(); i++){
	aux = outPath + "/" + leadName + "_" + std::to_string(i + 1) + ".vtx";
	write2file(aux, fOutVtx );
	fOutVtx << "1\nintra\n" << lead[i];
	fOutVtx.close();
	}


}

int findRVapex(const UVC &BiV_UVC){
	
	std::vector<double> Z_vec;
	std::vector<int> vtx;
	int aux_pos, final_vtx;

	for(int i = 0; i < BiV_UVC.Z.size(); i++){
		if((BiV_UVC.V[i] == 1) && (BiV_UVC.RHO[i] == 0)){
			Z_vec.push_back(BiV_UVC.Z[i]);
			vtx.push_back(i);
		}
	}

	aux_pos = min_pos(Z_vec);
	final_vtx = vtx[aux_pos];

	return final_vtx;
}

int findLVapex(const UVC &BiV_UVC){

	std::vector<double> Z_vec;
	std::vector<int> vtx;
	int aux_pos, final_vtx;

	for(int i = 0; i < BiV_UVC.Z.size(); i++){
		if((BiV_UVC.V[i] == -1) && (BiV_UVC.RHO[i] == 1)){
			Z_vec.push_back(BiV_UVC.Z[i]);
			vtx.push_back(i);
		}
	}

	aux_pos = min_pos(Z_vec);
	final_vtx = vtx[aux_pos];

	return final_vtx;
}
