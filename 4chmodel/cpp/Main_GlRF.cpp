#include "./headers/Main_GlRF.h"

// <--- From bivExtraction.sh
// --> To run_UVC_eikonal.sh

/* ------------------------------------------------------------------------------------------------------------------------- 
   To compile the code and execute it, cd to the folder containing the main script:

/home/crg17/Desktop/scripts/4chmodel/cpp/bin/Main_GlRF.o 00

cd /home/crg17/Desktop/scripts/4chmodel/cpp; g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./src/vectorops.cpp ./src/map.cpp ./src/common.cpp ./src/visualize.cpp ./src/EP.cpp ./src/GlFibreCorrection.cpp ./Main_GlRF.cpp -o ./bin/Main_GlRF.o
---------------------------------------------------------------------------------------------------------------------------- */

int main(int argc,char* argv[])
{

	system("clear");

	/* Set the path to files */

    std::string cmd, pathSubmsh, current_case = argv[1], f_endo = argv[2], f_epi = argv[3];
	std::string mshName = "h_case" + current_case; 
    std::string path = "/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/" + mshName + "/meshing/1000um";


    // // <----- From fibres.sh


    pathSubmsh = "/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case" + current_case + "/meshing/1000um/BiV";

    mshName = "h_case" + current_case;
    std::vector<int> tags{1,2};

    cmd = "cp " + pathSubmsh + "/fibres/fibres_bayer_" + f_epi + "_" + f_endo + ".lon " + pathSubmsh + "/BiV.lon";
    system(cmd.c_str());

    GlFibreCorrection(pathSubmsh,"BiV");

    cmd = "cp " + pathSubmsh + "/BiV.lon " + pathSubmsh + "/BiV_precorrected.lon";
    system(cmd.c_str());

    cmd = "cp " + pathSubmsh + "/BiV_corrected.lon " + pathSubmsh + "/BiV.lon";
    system(cmd.c_str());

    orthogonalise_lon(pathSubmsh,"BiV");

    cmd = "cp " + pathSubmsh + "/BiV.lon " + pathSubmsh + "/BiV_preortho.lon";
    system(cmd.c_str());

    cmd = "cp " + pathSubmsh + "/BiV_orto.lon " + pathSubmsh + "/BiV.lon";
    system(cmd.c_str());

    mapFibres(path, mshName, pathSubmsh, "BiV", tags, mshName);

    cmd = "/home/common/bin/meshtool convert -imsh=" + path + "/" + mshName + " -omsh=" + path + "/" + mshName + ".vtk -ifmt=carp_txt";
    system(cmd.c_str());

    return 0;
}








