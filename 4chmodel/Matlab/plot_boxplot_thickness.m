function stats = stats_thickness(threshold)
    
    stats = [];
    
    for i=1:24
        if(i < 10)
            heart = "0" + num2str(i);
        end

        BiV_thickness = dlmread("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + heart + "/meshing/1000um/BiV/BiV_thickness.dat" ,' ',0,0);


        stats(end+1) = invprctile(BiV_thickness(BiV_thickness > 0),threshold);
    end
    
end