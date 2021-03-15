#include "./headers/Main_create_lead_comb.h"

/* ------------------------------------------------------------------------------------------------------------------------- 
Example of use: /home/crg17/Desktop/scripts/4chmodel/cpp/bin/create_all_ATs.o 01 HF AN /home/crg17/Desktop/Seg3DProjects/multipole_AT/AN
g++ -fopenmp -std=c++11 /home/crg17/Desktop/scripts/4chmodel/cpp/src/ReadFiles.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/src/EP.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/src/vectorops.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/Main_create_lead_comb.cpp -o /home/crg17/Desktop/scripts/4chmodel/cpp/bin/create_all_ATs.o
---------------------------------------------------------------------------------------------------------------------------- */
std::string leads2name(std::vector<int> lead_pos){

    std::vector<int> vector_code = {0,0,0,0,0,0,0,0};
    std::string codename = "_";

    for(int i = 0; i < lead_pos.size(); i++)
        vector_code[lead_pos[i]] = 1;

    for(int i = 0; i < vector_code.size(); i++)
        codename += std::to_string(vector_code[i]);

    return codename;
}

void combine2Ats(const std::string &path2AT_1, const std::string &AT_1Name, const std::string &path2AT_2, const std::string &AT_2Name,const std::string &path2AT_res, const std::string &AT_resName){

    std::vector<double> AT_1, AT_2, AT_res;
    std::string fOut_name, cmd;
	std::fstream fOutVtx;

    AT_1 = ReadActTime(path2AT_1, AT_1Name);
    AT_2 = ReadActTime(path2AT_2, AT_2Name);
    
    AT_res = min_vec(AT_1,AT_2);

    cmd = "mkdir -p " + path2AT_res;
    system(cmd.c_str());
    
    fOut_name = path2AT_res + "/" + AT_resName + ".dat";

        std::cout << "\n Writing in " << fOut_name << std::endl;

    write2file(fOut_name, fOutVtx);

    for(int k = 0; k < AT_res.size(); k++)
        fOutVtx << AT_res[k] << std::endl;
    
    fOutVtx.close();
}

int main(int argc,char* argv[])
{   

    std::string HHF = argv[1]; // h
    std::string heart = argv[2]; // 24
    std::string mod_param = argv[3]; // kFEC
    std::string value_param = argv[4]; // 4.5
    std::string lead_pos = argv[5]; // AL
    std::string eikonal_path = argv[6]; // /media/crg17/Seagate Backup Plus Drive/CT_cases
    std::string ATs_folder = argv[7]; // /data/SA_multipole

    // if(HHF == "H" || HHF == "h")
    //     eikonal_path = "/media/crg17/Seagate Backup Plus Drive/CT_Segmentations/h_case" + heart + "/simulations/multipole/eikonal";
    // else if (HHF == "HF" || HHF == "hf")
    //     eikonal_path = "/data/harddrive/CT_cases/HF_case" + heart + "/simulations/multipole/eikonal";
    



    // std::vector<std::string> PHI_lead = {"AN","AL","LA","PL","PO"};
    // std::vector<std::string> PHI_lead = {"LA"};
    std::string RV_apex = "BiV.RV_endo.apex";
    int num_electrodes = 8;

    std::string path2AT_1, AT_1Name, path2AT_2, AT_2Name, path2AT_res, AT_resName, root_path, path2ATs;

    /* Monopoles = Free wall electrodes + RV apex electrode */

    // for(int i = 0; i < PHI_lead.size(); i++){
    //     ATs_folder = "/home/crg17/Desktop/Seg3DProjects/multipole_AT/" + PHI_lead[i];
    //     for(int j = 0; j < num_electrodes; j++){
    //         path2AT_1 = eikonal_path + "/" + PHI_lead[i] + "_" + std::to_string(j+1);
    //         AT_1Name = "vm_act_seq";
    //         path2AT_2 = eikonal_path + "/" + RV_apex;
    //         AT_2Name = AT_1Name;
    //         path2AT_res = ATs_folder;
    //         AT_resName = HHF + heart + "_" + PHI_lead[i] + leads2name({j});

    //         combine2Ats(path2AT_1,AT_1Name,path2AT_2,AT_2Name,path2AT_res,AT_resName);
    //     }
    // }

        // for(int j = 0; j < num_electrodes; j++){
        //     root_path = eikonal_path + "/" + HHF + "_case" + heart + "/simulations/multipole/eikonal_" + mod_param + "_" + value_param;

        //     path2AT_1 = root_path + "/" + lead_pos + "_" + std::to_string(j+1);
        //     AT_1Name = "vm_act_seq";

        //     path2AT_2 = root_path + "/" + RV_apex;
        //     AT_2Name = AT_1Name;

        //     path2AT_res = ATs_folder + "/" + mod_param + "_" + value_param + "/" +  HHF + "/" + lead_pos;
        //     AT_resName = HHF + heart + "_" + lead_pos + leads2name({j});

        //     combine2Ats(path2AT_1,AT_1Name,path2AT_2,AT_2Name,path2AT_res,AT_resName);
        // }

    /* Dipoles = Monopoles + Monopoles */

    // for(int lead_pos = 0; lead_pos < PHI_lead.size(); lead_pos++){
    //     for(int elec_1 = 0; elec_1 < num_electrodes-1; elec_1++){
    //         for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes; elec_2++){
    //             ATs_folder = "/home/crg17/Desktop/Seg3DProjects/multipole_AT/" + PHI_lead[lead_pos];
    //             path2AT_1 = ATs_folder;
    //             AT_1Name = HHF + heart + "_" + PHI_lead[lead_pos] + leads2name({elec_1});
    //             path2AT_2 = ATs_folder;
    //             AT_2Name = HHF + heart + "_" + PHI_lead[lead_pos] + leads2name({elec_2});
    //             path2AT_res = ATs_folder;
    //             AT_resName = HHF + heart + "_" + PHI_lead[lead_pos] + leads2name({elec_1,elec_2});

    //             combine2Ats(path2AT_1,AT_1Name,path2AT_2,AT_2Name,path2AT_res,AT_resName);
    //         }
    //     }
    // }

        // for(int elec_1 = 0; elec_1 < num_electrodes-1; elec_1++){
        //     for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes; elec_2++){

        //         path2ATs = ATs_folder + "/" + mod_param + "_" + value_param + "/" +  HHF + "/" + lead_pos;

        //         AT_1Name = HHF + heart + "_" + lead_pos + leads2name({elec_1});

        //         AT_2Name = HHF + heart + "_" + lead_pos + leads2name({elec_2});

        //         AT_resName = HHF + heart + "_" + lead_pos + leads2name({elec_1,elec_2});

        //         combine2Ats(path2ATs,AT_1Name,path2ATs,AT_2Name,path2ATs,AT_resName);
        //     }
        // }

    // /* Tripoles */

    // for(int lead_pos = 0; lead_pos < PHI_lead.size(); lead_pos++){

        // for(int elec_1 = 0; elec_1 < num_electrodes - 2; elec_1++){
        //     for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes - 1; elec_2++){
        //         for(int elec_3 = elec_2 + 1; elec_3 < num_electrodes; elec_3++){

        //             path2AT_1 = ATs_folder;
        //             AT_1Name = HHF + heart + "_" + lead_pos + leads2name({elec_1});
        //             path2AT_2 = ATs_folder;
        //             AT_2Name = HHF + heart + "_" + lead_pos + leads2name({elec_2,elec_3});
        //             path2AT_res = ATs_folder;
        //             AT_resName = HHF + heart + "_" + lead_pos + leads2name({elec_1,elec_2,elec_3});

        //             combine2Ats(path2AT_1,AT_1Name,path2AT_2,AT_2Name,path2AT_res,AT_resName);
        //         }
        //     }
        // }
    // }

        // for(int elec_1 = 0; elec_1 < num_electrodes - 2; elec_1++){
        //     for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes - 1; elec_2++){
        //         for(int elec_3 = elec_2 + 1; elec_3 < num_electrodes; elec_3++){

        //             path2ATs = ATs_folder + "/" + mod_param + "_" + value_param + "/" +  HHF + "/" + lead_pos;

        //             AT_1Name = HHF + heart + "_" + lead_pos + leads2name({elec_1});

        //             AT_2Name = HHF + heart + "_" + lead_pos + leads2name({elec_2,elec_3});

        //             AT_resName = HHF + heart + "_" + lead_pos + leads2name({elec_1,elec_2,elec_3});

        //             combine2Ats(path2ATs,AT_1Name,path2ATs,AT_2Name,path2ATs,AT_resName);
        //         }
        //     }
        // }

    // /* Quadripoles */

    // for(int lead_pos = 0; lead_pos < PHI_lead.size(); lead_pos++){
        // for(int elec_1 = 0; elec_1 < num_electrodes - 3; elec_1++){
        //     for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes - 2; elec_2++){
        //         for(int elec_3 = elec_2 + 1; elec_3 < num_electrodes - 1; elec_3++){
        //             for(int elec_4 = elec_3 + 1; elec_4 < num_electrodes; elec_4++){

        //                 path2AT_1 = ATs_folder;
        //                 AT_1Name = HHF + heart + "_" + lead_pos + leads2name({elec_1});
        //                 path2AT_2 = ATs_folder;
        //                 AT_2Name = HHF + heart + "_" + lead_pos + leads2name({elec_2,elec_3,elec_4});
        //                 path2AT_res = ATs_folder;
        //                 AT_resName = HHF + heart + "_" + lead_pos + leads2name({elec_1,elec_2,elec_3,elec_4});

        //                 combine2Ats(path2AT_1,AT_1Name,path2AT_2,AT_2Name,path2AT_res,AT_resName);
        //             }
        //         }
        //     }
        // }
    // }

        for(int elec_1 = 0; elec_1 < num_electrodes - 3; elec_1++){
            for(int elec_2 = elec_1 + 1; elec_2 < num_electrodes - 2; elec_2++){
                for(int elec_3 = elec_2 + 1; elec_3 < num_electrodes - 1; elec_3++){
                    for(int elec_4 = elec_3 + 1; elec_4 < num_electrodes; elec_4++){

                        path2ATs = ATs_folder + "/" + mod_param + "_" + value_param + "/" +  HHF + "/" + lead_pos;

                        AT_1Name = HHF + heart + "_" + lead_pos + leads2name({elec_1});
                        AT_2Name = HHF + heart + "_" + lead_pos + leads2name({elec_2,elec_3,elec_4});
                        AT_resName = HHF + heart + "_" + lead_pos + leads2name({elec_1,elec_2,elec_3,elec_4});

                        combine2Ats(path2ATs,AT_1Name,path2ATs,AT_2Name,path2ATs,AT_resName);
                    }
                }
            }
        }

    return 0;
}
