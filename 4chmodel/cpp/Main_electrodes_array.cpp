#include "./headers/Main_electrodes_array.h"

/* ------------------------------------------------------------------------------------------------------------------------- 
   clear; g++ -fopenmp -std=c++11 /home/crg17/Desktop/scripts/4chmodel/cpp/src/ReadFiles.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/src/EP.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/src/vectorops.cpp /home/crg17/Desktop/scripts/4chmodel/cpp/Main_electrodes_array.cpp -o /home/crg17/Desktop/scripts/4chmodel/cpp/bin/elec_array_compiled.o; /home/crg17/Desktop/scripts/4chmodel/cpp/bin/elec_array_compiled.o /data/01/Mesh/800um/BiV/UVC/UVC COMBINED_COORDS_Z_RHO_PHI_V /data/01/Mesh/800um/BiV/biv/laplace_full_base phie_scaled /data/01/Mesh/800um/BiV basal01 30 337611; rm /home/crg17/Desktop/scripts/4chmodel/cpp/bin/elec_array_compiled.o
---------------------------------------------------------------------------------------------------------------------------- */

int main(int argc,char* argv[])
{   

    std::string path2UVC = argv[1];
    std::string UVCfile = argv[2];
    std::string path2z = argv[3];
    std::string zfile = argv[4];
    std::string outPath = argv[5];
    std::string outName = argv[6];
    std::string n_el_str = argv[7];
    std::string el0_str = argv[8];

    int n_el = std::stoi(n_el_str);
    int el0 = std::stoi(el0_str);
    std::vector<int> vector_array(n_el);
    vector_array[0] = el0;
    
    UVC BiV_UVC = ReadUVC(path2UVC,UVCfile);
    std::vector<double> phi_normalised = scale_UVC(BiV_UVC,"PHI");
    // std::vector<double> z_normalised = scale_UVC(BiV_UVC,"Z");
    std::vector<double> z_normalised = ReadActTime(path2z,zfile);

   	std::vector<double>::iterator maxphi = std::max_element(phi_normalised.begin(), phi_normalised.end());
    std::vector<double>::iterator minphi = std::min_element(phi_normalised.begin(), phi_normalised.end());

    double el_dist = std::abs(maxphi[0]-minphi[0])/n_el;
    double z_fixed = z_normalised[el0];
    double z_tolerance = 1e-2;
    double phi_goal = phi_normalised[el0];
    double phi_tolerance = 1e-2;
    
    std::cout << "Apico-basal accuracy of " << z_tolerance << " and rotational tolerance of " << phi_tolerance << "." << std::endl;	
	
    for(int i = 1; i < n_el; i++){
        phi_goal += el_dist;

        if(phi_goal > maxphi[0])
            phi_goal = minphi[0] - maxphi[0] + phi_goal;
        
        vector_array[i] = single_electrode_UVC(BiV_UVC, z_normalised, z_fixed, 1, phi_goal, 2, z_tolerance, phi_tolerance);
    }

    std::fstream fOut;
    std::string output_name = outPath + "/" + outName + ".vtx";
    write2file(output_name, fOut );
    fOut << n_el << std::endl << "intra" << std::endl;

    for(int i = 0; i < n_el; i++)
        fOut << vector_array[i] << std::endl;

    fOut.close();
    
    return 0;

}
