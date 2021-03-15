#include "./headers/Main_merge_4PM.h"

/* ------------------------------------------------------------------------------------------------------------------------- 
   To execute :

/home/crg17/Desktop/scripts/4chmodel/cpp/bin/Main_merge_4PM.o 00


To compile:
g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./Main_merge_4PM.cpp -o ./bin/Main_merge_4PM.o
----------------------------------------------------------------------------------------------------------------------- */

int main(int argc,char* argv[])
{

	system("clear");

	/* Set the path to files */

    std::string current_case = argv[1];
    std::string action = argv[2];

	std::string mshName = "h_case" + current_case; 
    std::string path = "/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/" + mshName + "/meshing/1000um";

    // /* --------------------------------------------------------------------*/
    // // Create bases
    if(action == "ava"){
        mergeVtx(path,{"Ao_RV_base_ava.surf","AV_base_ava.surf","PV_base_ava.surf","Ao_base_ava.surf","LV_MV_base_ava.surf","RV_TV_base_ava.surf"},path,"PM_base_ava.surf");
        mergeSurf(path,{"Ao_RV_base_ava","AV_base_ava","PV_base_ava","Ao_base_ava","LV_MV_base_ava","RV_TV_base_ava"},path,"PM_base_ava");
        mergeNeubc(path,{"Ao_RV_base_ava","AV_base_ava","PV_base_ava","Ao_base_ava","LV_MV_base_ava","RV_TV_base_ava"},path,"PM_base_ava.surf");
    }
    else if (action == "PM"){
    mergeVtx(path,{"Ao_base.surf","AV_base.surf","LV_base.surf","RV_base.surf","Ao_RV_base.surf"},path,"PM_base.surf");
    mergeSurf(path,{"Ao_base","AV_base","LV_base","RV_base","Ao_RV_base"},path,"PM_base");
    mergeNeubc(path,{"Ao_base","AV_base","LV_base","RV_base","Ao_RV_base"},path,"PM_base.surf");
    }
    

    return 0;

}









