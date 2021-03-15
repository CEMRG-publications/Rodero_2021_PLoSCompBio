#include "./headers/Main_volumes.h"



//cd /home/crg17/Desktop/scripts/4chmodel/cpp; g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./src/visualize.cpp ./src/map.cpp ./src/common.cpp ./src/heart_metrics.cpp ./src/vectorops.cpp Main_volumes.cpp -o ./bin/Main_volumes_compiled; ./bin/Main_volumes_compiled; rm ./bin/Main_volumes_compiled

/* <------ From Main_surfs.cpp */

int main(){

    std::string command, current_case = "-1";
    std::ofstream fOut;
    double totVol, ao_diam;
    std::vector<double> vol_surf(2);
    int int_case;

    std::string path = "/data/harddrive/CT_cases/h_case" + current_case + "/meshing/1000um";
    std::string path2save = "/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/forall";
    std::string table_filename = "chambers_volumes_synth";

    fOut.open(path2save + "/" + table_filename + ".txt");
	fOut << "LV_mL LV_cm2 RV_mL RV_cm2 LA_mL LA_cm2 RA_mL RA_cm2 \n";
    fOut.close();
    
    for(int_case = 21; int_case <= 39; int_case++){

        if(int_case < 10)
            current_case = "0" + std::to_string(int_case);
        else
            current_case = std::to_string(int_case);

        path = "/data/harddrive/CT_cases/h_case" + current_case + "/meshing/1000um/cavities";

        vol_surf = surfaceVolume(path, "LV_endo_closed", "h_case" + current_case);

        fOut.open(path2save+"/"+table_filename+".txt",std::ios_base::app);
        fOut << vol_surf.at(0) << " " << vol_surf.at(1) << " ";
        fOut.close();

        vol_surf = surfaceVolume(path, "RV_endo_closed", "h_case" + current_case);

        fOut.open(path2save+"/"+table_filename+".txt",std::ios_base::app);
        fOut << vol_surf.at(0) << " " << vol_surf.at(1) << " ";
        fOut.close();

        vol_surf = surfaceVolume(path, "LA_endo_closed", "h_case" + current_case);

        fOut.open(path2save+"/"+table_filename+".txt",std::ios_base::app);
        fOut << vol_surf.at(0) << " " << vol_surf.at(1) << " ";
        fOut.close();

        vol_surf = surfaceVolume(path, "RA_endo_closed", "h_case" + current_case);

        fOut.open(path2save+"/"+table_filename+".txt",std::ios_base::app);
        fOut << vol_surf.at(0) << " " << vol_surf.at(1) << std::endl;
        fOut.close();

        // vol_surf = surfaceVolume(path,"AV_Ao.surfmesh","h_case" + current_case);
        // // vol_surf = surfaceVolume(path,"AV_LV.surfmesh","h_case" + current_case);
        // ao_diam = 2.0*(sqrt(vol_surf.at(1)/(4*atan(1))))*10; // In mm

        // fOut.open(path2save+"/"+table_filename+".txt",std::ios_base::app);
        // fOut << ao_diam << "\n";
        // fOut.close();
    }


    return 0;
}
