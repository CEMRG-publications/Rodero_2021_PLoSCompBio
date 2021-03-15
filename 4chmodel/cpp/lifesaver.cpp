#include "./headers/lifesaver.h"
#include "./headers/GlFibreCorrection.h"

/*
cd /home/crg17/Desktop/scripts/4chmodel/cpp; g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./src/vectorops.cpp ./src/convert.cpp ./src/map.cpp ./src/common.cpp ./src/visualize.cpp ./src/extract.cpp ./src/heart_metrics.cpp ./src/EP.cpp ./src/GlFibreCorrection.cpp ./lifesaver.cpp -o ./bin/lifesiver.o
*/
bool compare(int a, int b, double* data){
    return data[a]<data[b];
    }
int main(int argc,char* argv[]){

    std::string heart = argv[1];
    std::string meshName = "h_case" + heart;
    std::string meshPath = "/media/crg17/Seagate Expansion Drive/" + meshName + "/meshing/1000um";
    std::string vtxPath = meshPath + "/BiV";
    std::string vtxName = "BiV";
    std::string UVCPath = vtxPath + "/UVC/UVC";
    std::string UVCName = "COMBINED_COORDS_Z_RHO_PHI_V";
    std::string outPath = meshPath + "/UVC";

    mapUVCBack(vtxPath, vtxName, UVCPath, UVCName, meshPath, meshName, outPath);

    return 0;
}
