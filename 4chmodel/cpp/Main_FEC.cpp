#include "./headers/Main_FEC.h"

/*
<--- From run_UVC.sh

./bin/Main_FEC.o $heart $HHF $thickness $height $firsttime

cd /home/crg17/Desktop/scripts/4chmodel/cpp; g++ -std=c++11 -fopenmp ./src/ReadFiles.cpp ./src/vectorops.cpp ./src/convert.cpp ./src/EP.cpp ./Main_FEC.cpp -o ./bin/Main_FEC.o
*/

int main(int argc,char* argv[]){


    std::string current_case = argv[1];
    std::string thickness_str = argv[3];
    std::string HHF = argv[2];
    std::string height_str = argv[4];
    std::string firsttime_str = argv[5];

    double thickness = 0.01*(std::stod(thickness_str));
    double height    = 0.01*(std::stod(height_str));
    bool firsttime = str2bool(firsttime_str);

    std::string cmd;
    

    std::string path = "/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/" + HHF + "_case" + current_case + "/meshing/1000um/BiV" ;
   // std::string path = "/home/crg17/Desktop/Seg3DProjects/HF_CT/" + current_case + "HF/Mesh/800um/BiV" ;
    //std::string path = "/data/harddrive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV";
    std::string path2whole = path + "/..";
    std::string path2apba = path + "/UVC/UVC";
    //std::string path2apba = path + "/UVC";

    if(firsttime){
        cmd = "mkdir -p " + path + "/FEC";
        system(cmd.c_str());

        cmd = "cp " + path + "/BiV.elem " + path + "/FEC/BiV_FEC_w0_h0.elem";
        system(cmd.c_str());

        cmd = "mkdir -p " + path2whole + "/FEC";
        system(cmd.c_str());

        cmd = "cp " + path2whole + "/h_case" + current_case + ".elem " + path2whole + "/FEC/h_case" + current_case + "_FEC_w0_h0.elem";
        system(cmd.c_str());
    }

    cmd = "cp " + path + "/FEC/BiV_FEC_w0_h0.elem " + path + "/BiV.elem";
    system(cmd.c_str());

    // WATCH OOOOOOOOUT
    //cmd = "cp " + path + "/FEC/BiV_FEC_w" + thickness_str + "_h" + height_str + ".elem " + path + "/BiV.elem";
    //system(cmd.c_str());

    cmd = "cp " + path2whole + "/FEC/h_case" + current_case + "_FEC_w0_h0.elem " + path2whole + "/h_case" + current_case + ".elem";
    system(cmd.c_str());

    std::cout << "\n\n\n Creating FEC in BiV...";

    createFEC(path,"BiV",path2apba,"COORDS_RHO",path2apba, "COORDS_Z", thickness,height);

    cmd = "mv " + path + "/BiV_FEC_w" + thickness_str + "_h" + height_str + ".elem " + path + "/FEC/.";
    system(cmd.c_str());

    std::cout << "\n\n\n Changing tags in BiV...";

    change_tags(path + "/FEC","BiV_FEC_w" + thickness_str + "_h" + height_str, {3,4}, {25,26}); // BiV_w5_h70_retagged.elem

    cmd = "cp " + path + "/FEC/BiV_FEC_w" + thickness_str + "_h" + height_str  + "_retagged.elem " + path + "/BiV.elem";
    system(cmd.c_str());

   // std::cout << "\n\n\n Mapping...";

    //cmd = "meshtool insert submesh -submsh=" + path + "/BiV -msh=" + path2whole + "/h_case" + current_case + " -ifmt=carp_txt -ofmt=carp_txt -outmsh=" + path2whole + "/FEC/h_case" + current_case + "_FEC_w" + thickness_str + "_h" + height_str + ".elem";
    //system(cmd.c_str());


    return 0;
}
