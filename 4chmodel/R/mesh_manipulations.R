Visualise_AT1090 <- function(heart){

  warning("The elem file must be without the header and the activation file must be element-wise (use meshtool)")
  
  aux <- read.table("~/Desktop/transfer/vm_act_seq_11HC.dat", quote="\"", comment.char="")
  AT_12 <- as.vector(aux$V1)
  rm(aux)
  aux <- read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case11/meshing/1000um/BiV/BiV_mesh_volume.dat", quote="\"", comment.char="")
  BiV_mesh_volume <- as.vector(aux$V1)
  
  # We sort the volumes according to the ones activated first
  
  idx2sort <- order(AT_12)
  vol_sorted <- BiV_mesh_volume[idx2sort]
  
  # We get the percentage of the total volume that it's each element activated
  
  tot_vol = sum(vol_sorted)
  
  accumul_vol = 100*cumsum(vol_sorted)/tot_vol
  
  # We recover the original sorting
  
  idx2unsort <- match(AT_12,sort(AT_12))
  
  final_vec <- accumul_vol[idx2unsort]
  
  write.table(final_vec,"/home/crg17/Desktop/transfer/h11_percentages_AT.dat",row.names = FALSE,col.names = FALSE)
}


Map_thickness_transmurally <- function(case_number){
  
  require("svMisc")
  require("dplyr")
  
  print(paste0("Case number ", case_number))
  # flush.console()
  
  BiV_thickness <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_thickness.dat"), quote="\"", comment.char="")
  BiV_thickness <- as.vector(BiV_thickness$V1)
  UVC <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/UVC/COMBINED_COORDS_Z_RHO_PHI_V.dat"), quote="\"", comment.char="")
  colnames(UVC) <- c("Z","RHO","PHI","V")
  
  # We find the indices of the surface, so with possitive thickness
  mesh_idx <- c(1:length(BiV_thickness))
  surface_idx <- which(BiV_thickness != -1)
  
  p <- progress_estimated(length(mesh_idx[!(mesh_idx%in%surface_idx)]))
  for (i in mesh_idx[!(mesh_idx%in%surface_idx)]) {
    if(UVC[i,"V"] == -1){
      p$tick()
      p$print()
      # if(i %% 10 == 0){
      #   print(paste0("Iteration ",toString(i)))
      # }
      diffPHI <- abs(UVC[surface_idx,"PHI"] - UVC[i,"PHI"])
      diffZ <- abs(UVC[surface_idx,"Z"] - UVC[i,"Z"])
      
      mean_vec <- apply(cbind(diffPHI,diffZ),1,'mean')
      
      idx_in_surface <- which.min(mean_vec)[1]
      idx_in_mesh <- surface_idx[idx_in_surface]
      
      BiV_thickness[i] <- BiV_thickness[idx_in_mesh]
    }
  }
  
  print("Writing results...")
  
  write.table(BiV_thickness,file = paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_thickness_transmural.dat"),quote = FALSE,row.names = FALSE,col.names = FALSE)
}

Select_scar_thickness <- function(case_number){
  print(paste0("Case number ", case_number, "\n"))
  
  BiV_thickness <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_thickness_transmural.dat"), quote="\"", comment.char="")
  BiV_thickness <- as.vector(BiV_thickness$V1)
  
  BiV_AHA <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_AHA.dat"), quote="\"", comment.char="")
  BiV_AHA <- as.vector(BiV_AHA$V1)
  
  BiV_septum <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV.rvsept.surf.vtx"), quote="\"", comment.char="")
  ll <- dim(BiV_septum)[1]
  septum_vtx <- as.double(BiV_septum$V1[3:ll])
  
  thin_vtx <- which(BiV_thickness < 5000)
  
  noapex_vtx <- which(BiV_AHA < 17)
  nobase_vtx <- which(BiV_AHA > 0)
  
  chosen_vtx_septum <- Reduce(intersect, list(thin_vtx,noapex_vtx,nobase_vtx))
  chosen_vtx <- chosen_vtx_septum[!(chosen_vtx_septum %in% septum_vtx)]
  
  is_scar <- rep(0,length(BiV_AHA))
  is_scar[chosen_vtx] <- 1
  
  write.table(is_scar,file = paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_scar_5mm_no_septum.dat"),quote = FALSE,row.names = FALSE,col.names = FALSE)
}

Change_scar_tag <- function(case_number){
  print(paste0("Case number ", case_number))
  
  BiV_noheader <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_FEC70_noheader.elem"), quote="\"", comment.char="")
  colnames(BiV_noheader) <- c("el","p1","p2","p3","p4","tag")
  
  BiV_scar <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_scar_5mm_elem.dat"), quote="\"", comment.char="")
  
  BiV_noheader$tag[which(BiV_scar > 0)] <- 100
  
  write.table(BiV_noheader,file = paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",case_number,"/meshing/1000um/BiV/BiV_FEC70_scar5mm_noheader.elem"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  
}

Extract_endo_PN <- function(case_number, flag_debuggin = FALSE){
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  
  Load_Install_Packages("dplyr")
  
  
  path2mesh <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",
                      case_number,"/meshing/1000um/")
  
  mesh_elem <- Read_table(file = paste0(path2mesh,"h_case",case_number,".elem"),
                          skip = 1)
  colnames(mesh_elem) <- c("Tt","p1",)
  mesh_pts <- Read_table(file = paste0(path2mesh,"h_case",case_number,".pts"),
                          skip = 1)
  mesh_lon <- Read_table(file = paste0(path2mesh,"h_case",case_number,".lon"),
                          skip = 1)
  
  
}

Merge_fibre_files <- function(path2files, LVfibfile, RVfibfile, tagfilename,
                              output_name){
  
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  Load_Install_Packages("dplyr")
  
  LVfibres <- paste0(path2files, "/", LVfibfile, ".lon") %>%
              Read_table(file = ., skip = 1)
  
  RVfibres <- paste0(path2files, "/", RVfibfile, ".lon") %>%
              Read_table(file = ., skip = 1)
  
  tagfile <- paste0(path2files, "/", tagfilename, ".dat") %>%
             Read_table(file = ., skip = 1)
  tagfile[which(tagfile  == 25),] <- 1
  tagfile[which(tagfile  == 26),] <- 2
  
  outfibres <- LVfibres
  
  outfibres[which(tagfile  == 2),] <- RVfibres[which(tagfile  == 2),]
  
  formatted_file <- apply(outfibres, 1, function(x) gsub(",", "", toString(x)))
  formatted_file <- c(2,formatted_file)
  
  Write_table(x = formatted_file,
              file = paste0(path2files, "/", output_name, ".lon"),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}