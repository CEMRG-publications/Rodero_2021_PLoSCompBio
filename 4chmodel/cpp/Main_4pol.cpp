/*
g++ -std=c++11 ./src/ReadFiles.cpp ./src/vectorops.cpp ./src/EP.cpp Main_4pol.cpp ./src/convert.cpp -o ./bin/main_4pol_compiled.o
*/

#include "./headers/EP.h"
#include "./headers/convert.h"


int main(int argc,char* argv[]){
	std::string current_case = argv[1];
	std::string HHF = argv[2];
	std::string action = argv[3];

	std::string path2UVC, path2mesh, outPath, cmd, outPath_bash;

	if(HHF == "H"){
		path2UVC = "/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case" + current_case + "/meshing/1000um/BiV/UVC/UVC";
		path2mesh = "/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case" + current_case + "/meshing/1000um/BiV";
		outPath = "/media/crg17/Seagate Backup Plus Drive/CT_Segmentations/h_case" + current_case + "/simulations/multipole";
		outPath_bash = "/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/h_case" + current_case + "/simulations/multipole";
	}
	else if(HHF == "HF"){
		path2mesh = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV";
		path2UVC = path2mesh + "/UVC";
		outPath =  "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/simulations/multipole";
		outPath_bash = outPath;
	}

	std::string UVCname = "COMBINED_COORDS_Z_RHO_PHI_V";
	std::string mshName = "BiV";
	std::string surfaceVtx = "BiV.epi.surf";
	std::string aux_path = outPath + "/BiV.RV_endo.apex.vtx";

	std::vector<int> aux;

	double Z = 0.8;
	double RHO = 1;
	double PHI;
	double elec_distance = 7500;

	int V = -1;
	int first_electrode, second_electrode;
	int num_electrodes = 8;
	int RVapex;

	UVC BiV_UVC;

	std::fstream fOutVtx;

	cmd = "mkdir -p " + outPath_bash;
	system(cmd.c_str());

	aux = ReadVtx(path2mesh,"BiV.apex");


	if(action == "leads"){

 	/*      POSTERIOR      */
		PHI = 2.14;

		first_electrode = UVC2vtx(path2UVC,UVCname,Z,RHO,PHI,V);
		second_electrode = UVC2vtx(path2UVC,UVCname,Z-0.1,RHO,PHI,V);

		// second_electrode = aux[0];
		
		createElectrodesThread(path2mesh,mshName,surfaceVtx,first_electrode,second_electrode,elec_distance,num_electrodes,outPath,"PO");
	
	/*      POSTERIOLATERAL      */

		PHI = 2.64;

		first_electrode = UVC2vtx(path2UVC,UVCname,Z,RHO,PHI,V);
		second_electrode = UVC2vtx(path2UVC,UVCname,Z-0.1,RHO,PHI,V);
		// second_electrode = aux[0];
		
		createElectrodesThread(path2mesh,mshName,surfaceVtx,first_electrode,second_electrode,elec_distance,num_electrodes,outPath,"PL");
	
	/*      LATERAL      */

		
	 	PHI = -3.14;

	 	first_electrode = UVC2vtx(path2UVC,UVCname,Z,RHO,PHI,V);
	 	second_electrode = UVC2vtx(path2UVC,UVCname,Z-0.1,RHO,PHI,V);
		
	 	// std::cout << first_electrode << "\n" << second_electrode;
		
	 	createElectrodesThread(path2mesh,mshName,surfaceVtx,first_electrode,second_electrode,elec_distance,num_electrodes,outPath,"LA");
		
	/*      ANTEROLATERAL      */

		PHI = -2.64;

		first_electrode = UVC2vtx(path2UVC,UVCname,Z,RHO,PHI,V);
		second_electrode = UVC2vtx(path2UVC,UVCname,Z-0.1,RHO,PHI,V);
		
		createElectrodesThread(path2mesh,mshName,surfaceVtx,first_electrode,second_electrode,elec_distance,num_electrodes,outPath,"AL");
		
	/*      ANTERIOR      */

		PHI = -2.14;

		first_electrode = UVC2vtx(path2UVC,UVCname,Z,RHO,PHI,V);
		second_electrode = UVC2vtx(path2UVC,UVCname,Z-0.1,RHO,PHI,V);
		
		createElectrodesThread(path2mesh,mshName,surfaceVtx,first_electrode,second_electrode,elec_distance,num_electrodes,outPath,"AN");
	}
	else if(action == "electrodes"){
		lead2electrodes(outPath, "AN" , outPath, HHF + current_case);
		lead2electrodes(outPath, "AL" , outPath, HHF + current_case);
		lead2electrodes(outPath, "LA" , outPath, HHF + current_case);
		lead2electrodes(outPath, "PL" , outPath, HHF + current_case);
		lead2electrodes(outPath, "PO" , outPath, HHF + current_case);

		// BiV_UVC = ReadUVC(path2UVC,UVCname);
		// RVapex = findRVapex(BiV_UVC);
		// write2file(aux_path, fOutVtx );
		// fOutVtx << 1 << std::endl << "intra\n" << RVapex << std::endl;
		// fOutVtx.close();

	}

	return 0;
}

