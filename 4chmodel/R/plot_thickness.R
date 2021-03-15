plot_thickness <- function(){
  
  hearts <- c(paste0("0",1:2))

  
  for(heart in hearts){
    
    BiV_thickness <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",heart,"/meshing/1000um/BiV/BiV_thickness.dat"), quote="\"", comment.char="")
    BiV_thickness[BiV_thickness < 0] = NA
    
    if(heart == hearts[1])
      df_boxplot <- as.data.frame(BiV_thickness)
    
   cbind(df_boxplot,BiV_thickness)
    
  }
  
  return(df_boxplot)
}