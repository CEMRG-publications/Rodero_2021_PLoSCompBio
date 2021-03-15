/**
 * @brief Functions to compute some metrics such a surfaces and volumes.
 * 
 * @file heart_metrics.cpp
 * @author Cristobal Rodero
 * @date 2018
 */
#include "../headers/heart_metrics.h"

/**
 * @brief Computes the surface and the volume enclosed by a closed mesh using divergence theorem.
 * 
 * @attention If the mesh is not closed, reasonable results are not guaranteed.
 * 
 * @param path (<CODE>string</CODE>) Path of the mesh.
 * @param surf_filename (<CODE>string</CODE>) <CODE>.surf</CODE> file of the mesh.
 * @param pts_filename (<CODE>string</CODE>) <CODE>.pts</CODE> file of the mesh.
 * @return (<CODE>std::vector<double></CODE>) Vector whose first component is the volume enclosed by the mesh and the second
 * one is the area of the surface. 
 */
std::vector<double> surfaceVolume(const std::string &path,const std::string &surf_filename,const std::string &pts_filename){

    std::vector<double> d13(3), d12(3), cr(3), res(2);
    double totalVolume, crNorm, zMean, nz, totalarea;

    points p = ReadPts(path,pts_filename);
    surf t = ReadSurf(path,surf_filename);

    std::vector<double> volume(t.v1.size()), area(t.v1.size());
    
    for(int tr = 0; tr < t.v1.size(); tr++){
        
        d13.at(0) = (p.x).at(((t.v1).at(tr))) - (p.x).at(((t.v3).at(tr)));
        d13.at(1) = (p.y).at(((t.v1).at(tr))) - (p.y).at(((t.v3).at(tr)));
        d13.at(2) = (p.z).at(((t.v1).at(tr))) - (p.z).at(((t.v3).at(tr)));
        
        d12.at(0) = (p.x).at(((t.v1).at(tr))) - (p.x).at(((t.v2).at(tr)));
        d12.at(1) = (p.y).at(((t.v1).at(tr))) - (p.y).at(((t.v2).at(tr)));
        d12.at(2) = (p.z).at(((t.v1).at(tr))) - (p.z).at(((t.v2).at(tr)));


        cr = cross(d13,d12);
        crNorm = norm(cr);
        
        area.at(tr) = 0.5*crNorm;
        zMean = ((p.z).at(((t.v1).at(tr))) + (p.z).at(((t.v2).at(tr))) + (p.z).at(((t.v3).at(tr)))) / 3.0;
        nz = cr.at(2)/crNorm;

        volume.at(tr) = area.at(tr)*zMean*nz;
    }

    totalVolume = sum(volume)*1e-12; // In mL
    totalarea = sum(area)*1e-8; // In cm2

    res.at(0) = fabs(totalVolume);
    res.at(1) = totalarea;

    return res;
}
/**
 * @brief Computes the distance between a vector and a surface of points.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @param p (<CODE>points</CODE>) <CODE>points</CODE> structure describing the surface of points.
 * @return (<CODE>double</CODE>) Distance from the vector to the <CODE>points</CODE> surface. 
 */
double dist2pts(const std::vector<double> &v, const points &p){

    std::vector<double> dist_vec(p.x.size());

    #pragma parallel for
    for(int i=0; i < p.x.size(); i++)
        dist_vec.at(i) = norm(v-extract_point(p,i));
    
    return(min_positive(dist_vec));
}
/**
 * @brief Computes the volume of a mesh as the sum of the volume of its elements.
 * 
 * @param path2elem (<CODE>string</CODE>) Path to the <CODE>.elem</CODE> file of the mesh.
 * @param elem_name (<CODE>string</CODE>) Name of the <CODE>.elem</CODE> file of the mesh.
 * @param path2pts (<CODE>string</CODE>) Path to the <CODE>.pts</CODE> file of the mesh.
 * @param pts_name (<CODE>string</CODE>) Name to the <CODE>.pts</CODE> file of the mesh.
 * @param wanna_file (<CODE>bool</CODE>) Flag to save a <CODE>.dat</CODE> file with the volume of each element (1 if must be saved,
 * 0 otherwise).
 * @param outpath (<CODE>string</CODE>) Path of the elements volumes file, if <VAR>wanna_file</VAR>=1.
 * @param outname (<CODE>string</CODE>) Name of the elements volumes file, if <VAR>wanna_file</VAR>=1.
 * @return (<CODE>double</CODE>) Volume of the mesh.
 * 
 */

double volume_mesh(const std::string &path2elem, const std::string &elem_name, const std::string &path2pts, const std::string &pts_name,const bool &wanna_file, const std::string &outpath, const std::string &outname){

    std::string outPath;
    elem el = ReadElem(path2elem, elem_name);
    points pts = ReadPts(path2pts, pts_name);

    std::vector<double> vol_vec(el.v1.size());
    std::vector<std::vector<double>> a(vol_vec.size(),std::vector<double>(3)), b(vol_vec.size(),std::vector<double>(3)), c(vol_vec.size(),std::vector<double>(3)), d(vol_vec.size(),std::vector<double>(3));

    std::cout << "Computing the volume of each element...\n";
    
    #pragma omp parallel
    for(long int i = 0; i < vol_vec.size(); i++){

        // std::cout << " i = " << i << ", el = " << el.v1[i] << " " << el.v2[i] << " " << el.v3[i] << " " << el.v4[i] << ", pts_x = " << pts.x[el.v1[i]] <<"\n";
        
        // std::cout << i << std::endl;

        a[i] = {pts.x[el.v1[i]], pts.y[el.v1[i]], pts.z[el.v1[i]]};
        b[i] = {pts.x[el.v2[i]], pts.y[el.v2[i]], pts.z[el.v2[i]]};
        c[i] = {pts.x[el.v3[i]], pts.y[el.v3[i]], pts.z[el.v3[i]]};
        d[i] = {pts.x[el.v4[i]], pts.y[el.v4[i]], pts.z[el.v4[i]]};

        // vol_vec[i] = (abs(dot(a[i]-d[i],cross(b[i]-d[i],c[i]-d[i])))/6.0) * 1e-12;
        vol_vec[i] = (abs(dot(a[i]-d[i],cross(b[i]-d[i],c[i]-d[i])))/6.0) * 1e-09;


        // vol_vec[i] = abs(dot({pts.x[el.v1[i]], pts.y[el.v1[i]], pts.z[el.v1[i]]}-{pts.x[el.v4[i]], pts.y[el.v4[i]], pts.z[el.v4[i]]},cross({pts.x[el.v2[i]], pts.y[el.v2[i]], pts.z[el.v2[i]]}-{pts.x[el.v4[i]], pts.y[el.v4[i]], pts.z[el.v4[i]]},{pts.x[el.v3[i]], pts.y[el.v3[i]], pts.z[el.v3[i]]}-{pts.x[el.v4[i]], pts.y[el.v4[i]], pts.z[el.v4[i]]})))/6.0;
    }

    // std::cout << "\n VOLUME OF 1ST ELEMENT " << vol_vec[0] << std::endl; 

    if(wanna_file){

        std::cout << "Saving the volumes file...\n";

        std::fstream fOutVol;
        outPath = outpath + "/" + outname + ".dat";
        write2file(outPath, fOutVol);
        
        for(long int i = 0; i < vol_vec.size(); i++)
            fOutVol << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << vol_vec[i] << std::endl;
        
        fOutVol.close();
    }

    return(sum(vol_vec));
}
/**
 * @brief Approximates the height of a piece of the aorta.
 * 
 * It reads the mesh of the closest part of the aorta to the LV and the farthest part and compute the height between them.
 * These meshses must be disjoint and its recommended to be as thin as possible.
 * 
 * @param path (<CODE>string</CODE>) Path to the meshes.
 * @param pts_aorta_base (<CODE>string</CODE>) Name of the <CODE>.pts</CODE> file of the closest part of the aorta to the LV.
 * @param pts_aorta_top (<CODE>string</CODE>) Name of the <CODE>.pts</CODE> file of the farthest part of the aorta from the LV.
 * @return (<CODE>double</CODE>) Approximate height of the aorta piece. 
 */
double aortaHeight(const std::string &path, const std::string &pts_aorta_base, const std::string &pts_aorta_top){

    std::vector<double> dist_vec;

    points p1 = ReadPts(path,pts_aorta_base);
    points p2 = ReadPts(path,pts_aorta_top);

    #pragma parallel for
    for(int i = 0; i < p1.x.size(); i++)
        dist_vec.push_back(dist2pts(extract_point(p1,i),p2));

    return ((sum(dist_vec)/dist_vec.size())*1e-3);
}

/**
 * @brief Compute the total activation time of a biventricular mesh (<VAR>TAT</VAR>), the total activation time of the left ventricle
 * (<VAR>TAT_LV</VAR>) and the activation time between the 10 and the 90% of the volume is activated (<VAR>AT 10-90</VAR>).
 * 
 * @param path2elem (<CODE>string</CODE>) Path to the <CODE>.elem</CODE> file of the biventricular mesh.
 * @param elem_file (<CODE>string</CODE>) Name of the <CODE>.elem</CODE> file of the biventricular mesh.
 * @param path2dat (<CODE>string</CODE>) Path to the <CODE>.dat</CODE> file of the biventricular mesh activation times.
 * @param dat_file (<CODE>string</CODE>) Name of the <CODE>.dat</CODE> file of the biventricular mesh activation times.
 * @return (<CODE>std::vector<double></CODE>) Vector whose first component is <VAR>TAT</VAR>, the second one is <VAR>AT 10-90</VAR>
 * and the third one is <VAR>TAT_LV</VAR>.
 */

std::vector<double> AT_metrics(const std::string &path2elem, const std::string &elem_file, const std::string &path2dat, const std::string &dat_file, const std::string &path2biv, const std::string &biv_file){

    struct wrt{
        const std::vector<double> & value_vector;

        wrt(const std::vector<double> & val_vec):
            value_vector(val_vec) {}

        bool operator()(int i1, int i2)
        {
            return value_vector[i1] < value_vector[i2];
        }
    };

    struct by_AT {
        bool operator()(AT_vol a, AT_vol b) const noexcept { 
        return a.el_ActTime < b.el_ActTime;
    }
    };



    elem el = ReadElem(path2elem,elem_file);
    elem el_biv = ReadElem(path2biv,biv_file);
    std::vector<double> AT_4ch = ReadActTime(path2dat,dat_file);
    std::vector<int> vtx = ReadVtx(path2biv,biv_file);
    std::vector<double> AT(vtx.size());
    std::vector<double> metrics_vec(3);
    std::string command;
    double tot_vol, AT_10, AT_90;
    bool is_AT = 0;
    std::vector<double> temp, el_AT(el_biv.nElem);
    std::vector<AT_vol> elem_data(el_biv.nElem);

    #pragma omp parallel
    for(int i = 0; i < AT.size(); i++)
        AT[i] = AT_4ch[vtx[i]];

    metrics_vec[0] = max_value_vec(AT);

    std::vector<double> vol_vec = ReadActTime(path2biv,"BiV_mesh_volume");
    tot_vol = sum(vol_vec);

    std::cout << "Sorting elements by AT...";
    #pragma omp parallel
    for(long int i = 0; i < el_biv.nElem; i++){
        el_AT[i] = max_value_vec({AT[el_biv.v1[i]],AT[el_biv.v2[i]],AT[el_biv.v3[i]],AT[el_biv.v4[i]]});
    }


    std::vector<double> AT_sorted = el_AT;
    std::sort(vol_vec.begin(),vol_vec.end(),wrt(el_AT));

    std::sort(AT_sorted.begin(),AT_sorted.end());
    std::cout << "Accumulating volumes...\n";


    for(long int i = 1; i < vol_vec.size(); i++)
        vol_vec[i] += vol_vec[i-1];
    std::cout << "Computing AT_10...\n";

    long int i = 0;
    while(is_AT == 0){
        if(vol_vec[i]> 0.1*tot_vol){
            AT_10 = AT_sorted[i];
            is_AT = 1;
        }
        i++;
    }


    std::cout << "Computing AT_90...\n";
    is_AT = 0;
    i = vol_vec.size()-1;

    while(is_AT == 0){
        if(vol_vec[i] < 0.9*tot_vol){
            AT_90 = AT_sorted[i];
            is_AT = 1;
        }
        i--;
    }
    


    metrics_vec[1] = AT_90 - AT_10;

    #pragma omp parallel
    for(long int j = 0; j < el_AT.size(); j++){
        if(el_biv.tag[j] != 1 && el_biv.tag[j] != 25)
            el_AT[j] = 0;
    }

    metrics_vec[2] =  max_value_vec(el_AT) - min_positive(el_AT);

    return metrics_vec;
}

// OUTPUT: metrics_vec[0] TAT metrics_vec[1] TATLV metrics_vec[2] AT1090


std::vector<double> AT_metrics_BiV(const std::string &path2biv, const std::string &bivName, const std::string &path2volfile, const std::string &volfileName, const std::string &path2AT, const std::string &ATName){

    struct wrt{
        const std::vector<double> & value_vector;

        wrt(const std::vector<double> & val_vec):
            value_vector(val_vec) {}

        bool operator()(int i1, int i2)
        {
            return value_vector[i1] < value_vector[i2];
        }
    };

    struct by_AT {
        bool operator()(AT_vol a, AT_vol b) const noexcept { 
        return a.el_ActTime < b.el_ActTime;
    }
    };

    elem el_biv = ReadElem(path2biv,bivName);
    std::cout << "Reading act time from " << path2AT << " with the name " << ATName;
    std::vector<double> AT = ReadActTime(path2AT,ATName,1);
    std::vector<double> metrics_vec(3);
    double tot_vol, AT_10, AT_90;
    bool is_AT = 0;
    std::vector<double> el_AT(el_biv.nElem);
    std::vector<AT_vol> elem_data(el_biv.nElem);

    metrics_vec[0] = max_value_vec(AT);

    std::vector<double> vol_vec2;
    vol_vec2 = ReadActTime(path2volfile,volfileName,1);

    tot_vol = sum(vol_vec2);

    std::cout << "Sorting elements by AT...";
    #pragma omp parallel
    for(long int i = 0; i < el_biv.nElem; i++){
        el_AT[i] = max_value_vec({AT[el_biv.v1[i]],AT[el_biv.v2[i]],AT[el_biv.v3[i]],AT[el_biv.v4[i]]});
    }


    std::vector<double> AT_sorted = el_AT; // Initialisation
    
    std::sort(vol_vec2.begin(),vol_vec2.end(),wrt(el_AT)); // First the volume of the elements that activate the earliest

    std::sort(AT_sorted.begin(),AT_sorted.end()); // Sorted by volume

    // We accumulate the volumes
    for(long int i = 1; i < vol_vec2.size(); i++)
        vol_vec2[i] += vol_vec2[i-1];
        
    std::cout << "Computing AT_10...\n";

    long int i = 0;
    while(is_AT == 0){
        if(vol_vec2[i]> 0.1*tot_vol){
            AT_10 = AT_sorted[i];
            is_AT = 1;
        }
        i++;
    }


    std::cout << "Computing AT_90...\n";
    is_AT = 0;
    i = vol_vec2.size()-1;

    while(is_AT == 0){
        if(vol_vec2[i] < 0.9*tot_vol){
            AT_90 = AT_sorted[i];
            is_AT = 1;
        }
        i--;
    }
    


    metrics_vec[2] = AT_90 - AT_10;

    #pragma omp parallel
    for(long int j = 0; j < el_AT.size(); j++){
        if(el_biv.tag[j] != 1 && el_biv.tag[j] != 25 && el_biv.tag[j] != 3 && el_biv.tag[j] != 27 && el_biv.tag[j] != 29)
            el_AT[j] = 0;
    }

    metrics_vec[1] =  max_value_vec(el_AT) - min_positive(el_AT);

    return metrics_vec;
}

void BuildHClusterTable(const std::string &path2inout){

    std::vector<double> metrics(3);
    std::string current_case, path2biv, path2biv_vol, aux;
    std::ofstream outfile_QRS, outfile_TATLV, outfile_AT1090;
    EP_metrics_table tab;
    std::vector<std::string> rownames, colnames;
    std::cout << "path2simulation const is " << path2inout;
    std::string path2simulation = path2inout;
    std::string path2table = path2inout;

    std::vector<std::string> names_ATs = read_directory(path2simulation);


    for(int i = 0; i < names_ATs.size(); i++){

        std::cout << "\n\nProcessing file " << i+1 << "/" << names_ATs.size() << "\n\n";

        if( (names_ATs[i])[1] == 'F'){
            current_case = (names_ATs[i]).substr(2,2);
	    path2biv = "/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV";
            path2biv_vol = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV";

            tab.comb_elec.push_back((names_ATs[i]).substr(8,8));
            tab.HF.push_back("1");
	    //tab.lead_position.push_back((names_ATs[i]).substr(5,2));
        }
        else{
            current_case = (names_ATs[i]).substr(1,2);
	    path2biv = "/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/h_case" + current_case + "/meshing/1000um/BiV";
            path2biv_vol = "/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case" + current_case + "/meshing/1000um/BiV";
            
            tab.comb_elec.push_back((names_ATs[i]).substr(7,8));
            tab.HF.push_back("0");
            //tab.lead_position.push_back((names_ATs[i]).substr(4,2));
        }

        metrics = AT_metrics_BiV(path2biv, "BiV", path2biv_vol, "BiV_mesh_volume", path2simulation, (names_ATs[i]).substr(0,(names_ATs[i]).size() - 4)  );
        
        tab.TAT.push_back(metrics[0]);
        tab.TAT_LV.push_back(metrics[1]);
        tab.AT_10_90.push_back(metrics[2]);
        tab.heart.push_back(current_case);

    }

    std::cout << "\nWriting dataframe(s)...\n";

    rownames = tab.comb_elec;

    std::sort(rownames.begin(),rownames.end());
    auto it = std::unique(std::begin(rownames), std::end(rownames));
    rownames.erase(it,rownames.end());

    colnames = tab.heart;

    std::sort(colnames.begin(),colnames.end());
    it = std::unique(std::begin(colnames), std::end(colnames));
    colnames.erase(it,colnames.end());


    aux = path2table + "/multipole_QRS.dat";
    outfile_QRS.open(aux, std::ofstream::out);
    aux = path2table + "/multipole_TATLV.dat";
    outfile_TATLV.open(aux, std::ofstream::out);
    aux = path2table + "/multipole_AT1090.dat";
    outfile_AT1090.open(aux, std::ofstream::out);

    for(int i = 0; i < colnames.size(); i++){
        outfile_QRS << " " << colnames[i]; // First row the name of the hearts
        outfile_TATLV << " " << colnames[i]; 
        outfile_AT1090 << " " << colnames[i];
    }

    for(int i = 0; i < rownames.size(); i++){
        outfile_QRS << "\n" << rownames[i];
        outfile_TATLV << "\n" << rownames[i];
        outfile_AT1090 << "\n" << rownames[i];
        for(int j = 0; j < colnames.size(); j++){
            for(int k = 0; k < tab.TAT.size(); k++){
                if(tab.heart[k] == colnames[j] && tab.comb_elec[k] == rownames[i]){
                    outfile_QRS << " " << tab.TAT[k];
                    outfile_TATLV << " " << tab.TAT_LV[k];
                    outfile_AT1090 << " " << tab.AT_10_90[k];
                }
            }
        }
    }

    outfile_QRS.close();
    outfile_TATLV.close();
    outfile_AT1090.close();
}

void AT_metrics_midseptum(const std::string &path2table, const std::string &path2simulations, std::string &which_cases){

    std::vector<double> AT_vec;
    std::string path2biv, path2AT, aux;
    std::ofstream outTable;
    std::vector<std::string> hearts;
    
    if(which_cases == "h"){
        hearts = {"01","02","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20"};

        aux = path2table + "/AT_metric_midseptum_h.dat";
        outTable.open(aux, std::ofstream::out);

        outTable << " QRS TAT_LV AT_10_90";

        for(int i = 0; i < hearts.size(); i++){
            path2biv = path2simulations + "/h_case" + hearts[i] + "/meshing/1000um/BiV";
            path2AT = path2simulations + "/h_case" + hearts[i] + "/simulations/multipole/eikonal_midseptum/BiV.midseptum";
            AT_vec = AT_metrics_BiV(path2biv,"BiV",path2biv,"BiV_mesh_volume",path2AT,"vm_act_seq");

            outTable << "\n" << hearts[i] << " " << AT_vec[0] << " " << AT_vec[1] << " " << AT_vec[2];
        }

        outTable.close();
    }
    if(which_cases == "HF"){
        hearts = {"01","02","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"};
        aux = path2table + "/AT_metric_midseptum_HF.dat";

        outTable.open(aux, std::ofstream::out);

        outTable << " QRS TAT_LV AT_10_90";

        for(int i = 0; i < hearts.size(); i++){
            path2biv = path2simulations + "/HF_case" + hearts[i] + "/meshing/1000um/BiV";
            path2AT = path2simulations + "/HF_case" + hearts[i] + "/simulations/multipole/eikonal_midseptum/BiV.midseptum";
            AT_vec = AT_metrics_BiV(path2biv,"BiV",path2biv,"BiV_mesh_volume",path2AT,"vm_act_seq");

            outTable << "\n" << hearts[i] << " " << AT_vec[0] << " " << AT_vec[1] << " " << AT_vec[2];
        }

        outTable.close();
    }

}

std::vector<double> AT_metrics_BiV_optimised(const elem &el_biv,std::vector<double> vol_vec, const std::string &path2AT, const std::string &ATName){

    bool outflag=0;

    struct wrt{
        const std::vector<double> & value_vector;

        wrt(const std::vector<double> & val_vec):
            value_vector(val_vec) {}

        bool operator()(int i1, int i2)
        {
            return value_vector[i1] < value_vector[i2];
        }
    };

    struct by_AT {
        bool operator()(AT_vol a, AT_vol b) const noexcept { 
        return a.el_ActTime < b.el_ActTime;
    }
    };

    std::vector<double> AT = ReadActTime(path2AT,ATName,outflag);

    std::vector<double> metrics_vec(3);
    double tot_vol, AT_10, AT_90;
    bool is_AT = 0;
    std::vector<double> el_AT(el_biv.nElem);
    std::vector<AT_vol> elem_data(el_biv.nElem);

    metrics_vec[0] = max_value_vec(AT);

    tot_vol = sum(vol_vec);

    // std::cout << "Sorting elements by AT...";
    #pragma omp parallel
    for(long int i = 0; i < el_biv.nElem; i++){
        el_AT[i] = max_value_vec({AT[el_biv.v1[i]],AT[el_biv.v2[i]],AT[el_biv.v3[i]],AT[el_biv.v4[i]]});
    }
    

    std::sort(vol_vec.begin(),vol_vec.end(),wrt(el_AT)); // First the volume of the elements that activate the earliest

    // We accumulate the volumes
    for(long int i = 1; i < vol_vec.size(); i++)
        vol_vec[i] += vol_vec[i-1];

    std::vector<double> AT_sorted = el_AT; // Initialisation

    std::sort(AT_sorted.begin(),AT_sorted.end()); // Sorted by volume

        
    // std::cout << "Computing AT_10...\n";

    long int i = 0;
    while(is_AT == 0){
        if(vol_vec[i]> 0.1*tot_vol){
            AT_10 = AT_sorted[i];
            is_AT = 1;
        }
        i++;
    }


    // std::cout << "Computing AT_90...\n";
    is_AT = 0;
    i = vol_vec.size()-1;

    while(is_AT == 0){
        if(vol_vec[i] < 0.9*tot_vol){
            AT_90 = AT_sorted[i];
            is_AT = 1;
        }
        i--;
    }
    


    metrics_vec[2] = AT_90 - AT_10;

    #pragma omp parallel
    for(long int j = 0; j < el_AT.size(); j++){
        if(el_biv.tag[j] != 1 && el_biv.tag[j] != 25 && el_biv.tag[j] != 3 && el_biv.tag[j] != 27 && el_biv.tag[j] != 29)
            el_AT[j] = 0;
    }

    metrics_vec[1] =  max_value_vec(el_AT) - min_positive(el_AT);

    return metrics_vec;
}

// Ideally each folder will be of only one case.

void BuildHClusterTable_optimised(const std::string &path2inout, const std::string &path2biv){


    std::ofstream outfile_QRS, outfile_TATLV, outfile_AT1090;
    std::vector<std::string> rownames, colnames;

    std::string aux;
    std::string path2table = path2inout;
    std::string path2simulation = path2inout;
    std::vector<std::string> names_ATs = read_directory(path2inout);

    EP_metrics_table tab;
    std::vector<std::vector<double>> metrics(names_ATs.size(),std::vector<double>(3));
    std::vector<double> TAT_vec(names_ATs.size());
    std::vector<double> TAT_LV_vec = TAT_vec;
    std::vector<double> AT_10_90_vec = TAT_vec;
    std::vector<std::string> current_case(names_ATs.size());
    std::vector<std::string> comb_elec_vec = current_case;
    std::vector<std::string> lead_position_vec = current_case;
    bool outflag = false;



    // Read the volume file
    std::vector<double> vol_vec1 = ReadActTime(path2biv,"BiV_mesh_volume",outflag);

    // Read the elem file
    elem el_biv = ReadElem(path2biv,"BiV",outflag);

    //#pragma omp parallel
    for(int i = 0; i < names_ATs.size(); i++){

        if(i%10 == 0){
            std::cout << "\nProcessing file " << i << "\\" << names_ATs.size() << "...\n";
        }

        if( (names_ATs[i])[1] == 'F'){
            current_case[i] = (names_ATs[i]).substr(2,2);
            comb_elec_vec[i] = (names_ATs[i]).substr(8,8);
            lead_position_vec[i] = (names_ATs[i]).substr(5,2);


        }
        else if( ((names_ATs[i])[0] == 'H') || (names_ATs[i])[0] == 'h'){
            current_case[i] = (names_ATs[i]).substr(1,2);
            comb_elec_vec[i] = (names_ATs[i]).substr(4,8);
            comb_elec_vec[i] = (names_ATs[i]).substr(7,8);

        }
        else
        continue;



        metrics[i] = AT_metrics_BiV_optimised(el_biv, vol_vec1, path2simulation, (names_ATs[i]).substr(0,(names_ATs[i]).size() - 4)  );
        
        TAT_vec[i] = (metrics[i])[0];
        TAT_LV_vec[i] = (metrics[i])[1];
        AT_10_90_vec[i] = (metrics[i])[2];

    }

    tab.AT_10_90 = AT_10_90_vec;
    tab.comb_elec = comb_elec_vec;
    tab.heart = current_case;
    tab.TAT = TAT_vec;
    tab.TAT_LV = TAT_LV_vec;
    tab.lead_position = lead_position_vec;

    // std::cout << "\nWriting dataframe(s)...\n";

    rownames = tab.comb_elec;

    std::sort(rownames.begin(),rownames.end());
    auto it = std::unique(std::begin(rownames), std::end(rownames));
    rownames.erase(it,rownames.end());

    colnames = tab.lead_position;

    std::sort(colnames.begin(),colnames.end());
    it = std::unique(std::begin(colnames), std::end(colnames));
    colnames.erase(it,colnames.end());


    aux = path2table + "/multipole_QRS.dat";
    outfile_QRS.open(aux, std::ofstream::out);
    aux = path2table + "/multipole_TATLV.dat";
    outfile_TATLV.open(aux, std::ofstream::out);
    aux = path2table + "/multipole_AT1090.dat";
    outfile_AT1090.open(aux, std::ofstream::out);

    for(int i = 0; i < colnames.size(); i++){
        outfile_QRS << " " << colnames[i]; // First row the name of the lead positions
        outfile_TATLV << " " << colnames[i]; 
        outfile_AT1090 << " " << colnames[i];
    }

    for(int i = 0; i < rownames.size(); i++){
        outfile_QRS << "\n" << rownames[i];
        outfile_TATLV << "\n" << rownames[i];
        outfile_AT1090 << "\n" << rownames[i];
        for(int j = 0; j < colnames.size(); j++){
            for(int k = 0; k < tab.TAT.size(); k++){
                if(tab.lead_position[k] == colnames[j] && tab.comb_elec[k] == rownames[i]){
                    outfile_QRS << " " << tab.TAT[k];
                    outfile_TATLV << " " << tab.TAT_LV[k];
                    outfile_AT1090 << " " << tab.AT_10_90[k];
                }
            }
        }
    }

    outfile_QRS.close();
    outfile_TATLV.close();
    outfile_AT1090.close();
}


void BuildHClusterTable_RV(const std::string &path2hearts, const std::string &suffix, const std::string &path2out, const std::string &condition){
    
    std::vector<double> AT_vec;
    std::string path2biv, path2AT, aux, heart_str;
    std::ofstream outTable;
    int num_hearts = 0;

    if(condition == "h")
        num_hearts = 20;
    else if(condition == "HF")
        num_hearts = 24;

    aux = path2out + "/" + suffix + "/" + condition + "/multipole_RVapex.dat";
    outTable.open(aux, std::ofstream::out);

    outTable << " QRS TAT_LV AT_10_90";

    for(int i = 0; i < num_hearts; i++){
        if(i < 9)
            heart_str = "0" + std::to_string(i+1);
        else
            heart_str = std::to_string(i+1);

        path2biv = path2hearts + "/" + condition + "_case" + heart_str + "/meshing/1000um/BiV";
        path2AT = path2hearts + "/" + condition + "_case" + heart_str +  "/simulations/multipole/eikonal_" + suffix + "/BiV.RV_endo.apex";

        AT_vec = AT_metrics_BiV(path2biv,"BiV",path2biv,"BiV_mesh_volume",path2AT,"vm_act_seq");

        outTable << "\n" << heart_str << " " << AT_vec[0] << " " << AT_vec[1] << " " << AT_vec[2];
    }

    outTable.close();
}

void get_RVapex_UVC(){
    
    int RVapex;
    std::string heart, path2apex, path2UVC;
    UVC RV_apex_UVC;
    std::ofstream outfile_UVC;
    std::vector<int> RVapex_vec;

    outfile_UVC.open("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/RVapex_UVC.dat");

    for(int i = 1; i <= 20; i++){
        if(i < 10)
            heart = "0" + std::to_string(i);
        else
            heart = std::to_string(i);
        
        path2apex = "/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case" + heart + "/simulations/multipole";
        RVapex_vec = ReadVtx(path2apex,"BiV.RV_endo.apex");
        RVapex = RVapex_vec[0];

        path2UVC = "/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case" + heart + "/meshing/1000um/BiV/UVC/UVC";
        RV_apex_UVC = ReadUVC(path2UVC, "COMBINED_COORDS_Z_RHO_PHI_V");

        outfile_UVC << heart + "HC " << RV_apex_UVC.Z[RVapex] << " " << RV_apex_UVC.RHO[RVapex] << " " << RV_apex_UVC.PHI[RVapex]<< " " << RV_apex_UVC.V[RVapex] << std::endl; 

    }

    for(int i = 1; i <= 24; i++){
        if(i < 10)
            heart = "0" + std::to_string(i);
        else
            heart = std::to_string(i);
        
        path2apex = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + heart + "/simulations/multipole";
        RVapex_vec = ReadVtx(path2apex,"BiV.RV_endo.apex");
        RVapex = RVapex_vec[0];

        path2UVC = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + heart + "/meshing/1000um/BiV/UVC";
        RV_apex_UVC = ReadUVC(path2UVC, "COMBINED_COORDS_Z_RHO_PHI_V");

        outfile_UVC << heart + "HF " << RV_apex_UVC.Z[RVapex]<< " "  << RV_apex_UVC.RHO[RVapex]<< " "  << RV_apex_UVC.PHI[RVapex]<< " "  << RV_apex_UVC.V[RVapex] << std::endl; 

    }

    outfile_UVC.close();
}

void get_midseptum_AHA(const std::string &path2mesh, const std::string &path2UVC, const std::string & UVC_name, const std::string & AHA_name, const std::string & path2output, const std::string & output_name){
    
    UVC BiV_UVC;
    std::vector<double> AHA_map, AHA9_PHI, AHA9_Z;
    points LV_pts, BiV_pts;
    double avg_Z, max_PHI;
    int midseptum_vtx;
    std::string aux;
    std::ofstream output;
    
    
    BiV_UVC = ReadUVC(path2UVC, UVC_name);


    // We create a sequential vector of indices.
    std::vector<int> AHA9_idx(BiV_UVC.Z.size());
    std::iota (std::begin(AHA9_idx), std::end(AHA9_idx), 0);

    AHA9_PHI = BiV_UVC.PHI;
    AHA9_Z = BiV_UVC.Z;

    AHA_map = ReadActTime(path2mesh,AHA_name);
    
    #pragma omp parallel
    for(int i = 0; i < AHA_map.size(); i++){
        if((AHA_map[i] != 9) || (BiV_UVC.RHO[i] != 1)){
            AHA9_PHI[i] = -10;
            AHA9_Z[i] = -10;
            AHA9_idx[i] = -10;
        }
    }

    // We get only the points in the epicardium within the tag 9

    AHA9_PHI.erase(std::remove(AHA9_PHI.begin(), AHA9_PHI.end(), -10), AHA9_PHI.end());
    AHA9_Z.erase(std::remove(AHA9_Z.begin(), AHA9_Z.end(), -10), AHA9_Z.end());
    AHA9_idx.erase(std::remove(AHA9_idx.begin(), AHA9_idx.end(), -10), AHA9_idx.end());

    avg_Z = average(AHA9_Z);

    #pragma omp parallel
    for(int i = 0; i < AHA9_Z.size(); i++){
        if((AHA9_Z[i] > (avg_Z + 0.01)) || (AHA9_Z[i] < (avg_Z - 0.01))){
            AHA9_PHI[i] = -10;
            AHA9_Z[i] = -10;
            AHA9_idx[i] = -10;
        }
    }

    midseptum_vtx = AHA9_idx[0];
    max_PHI = AHA9_PHI[0];

    for(int i = 1; i < AHA9_Z.size(); i++){
        if(AHA9_PHI[i] > max_PHI){
            max_PHI = AHA9_PHI[i];
            midseptum_vtx = AHA9_idx[i];
        }
    }

    aux = path2output + "/" + output_name + ".vtx";

    output.open(aux);

    output << "1\nintra\n" << midseptum_vtx;

    output.close();
}
