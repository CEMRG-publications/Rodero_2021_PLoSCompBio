# list.of.packages <- c("nnet","corrplot","berryFunctions","plot.matrix","RColorBrewer","ggplot2","treemap","latex2exp")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# flags<-lapply(list.of.packages, require, character.only = TRUE)


merge_restored_cavities <- function(heart_case, flag_read = FALSE){
  
  if(flag_read)
    print(paste0("You need the file /media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.LV_original.csv"))
  
  cav.LV_original <- read.csv(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.LV_original.csv"), stringsAsFactors = FALSE)
  
  if(flag_red)
    print(paste0("You need the file /media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk350_restored/cav.LV.csv"))
  cav.LV_restored <- read.csv(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk350_restored/cav.LV.csv"), comment.char = "#")
  
  cav.LV_final <- cav.LV_original
  
  row_to_delete_from <- as.numeric(cav.LV_restored$Time[2])+51
  cav.LV_final <- cav.LV_final[-c((row_to_delete_from):nrow(cav.LV_final)),]
  cav.LV_final <- rbind(cav.LV_final,cav.LV_restored)
  
  write.csv(cav.LV_final, file = paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.LV.csv"),row.names = FALSE )
  
  print(paste0("You need the file /media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.RV_original.csv"))
  cav.RV_original <- read.csv(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.RV_original.csv"), stringsAsFactors = FALSE)
  print(paste0("You need the file /media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk350_restored/cav.RV.csv"))
  cav.RV_restored <- read.csv(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk350_restored/cav.RV.csv"), comment.char = "#")
  
  cav.RV_final <- cav.RV_original
  
  row_to_delete_from <- as.numeric(cav.RV_restored$Time[2])+51
  cav.RV_final <- cav.RV_final[-c((row_to_delete_from):nrow(cav.RV_final)),]
  cav.RV_final <- rbind(cav.RV_final,cav.RV_restored)
  
  write.csv(cav.RV_final, file = paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart_case,"/simulations/cohort/",heart_case,"HC_wk3_tpeak_120/cav.RV.csv"),row.names = FALSE )
  
}
#' @description Function to create (either printing or returning) a table with
#' the correlations between modes and phenotypes.
#' @param cohort "CT", "EXTREME3", "EXTREME2", "EXTREME1", "PROJECTED5",
#'  "PROJECTED9", "PROJECTED14", "PROJECTED18"
#' @param with_EP TRUE, FALSE
#' @param action "savefen", "savecorr", "saveplot", "plot", "return"
#' @param ventricle "LV", "RV"
#' @param remove_phenotypes Vector containing the names of the phenotypes you
#' don't want to take into account in the correlation.
#' @param read_fen_from_files TRUE, FALSE. If in the laptop, use TRUE
#' @param flag_read If TRUE prints the directory and name of all the files read
#' and written.
plot_correlation <- function(cohort, with_EP = TRUE, action = "plot",
                             ventricle = "LV", remove_phenotypes = c(), 
                             read_fen_from_file = FALSE, flag_read = FALSE){
  
  set.seed(1)
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  Load_Install_Packages(c("berryFunctions","corrplot"))
  #Creating the data
  if(flag_read)
    print(paste0("You need file /data/modes.csv"))
  modes <- as.data.frame(read.csv("/data/modes.csv", header=TRUE, sep = "\t"))
  modes <- modes[,2:ncol(modes)]
  # colnames(modes) <- paste0("Mode",c(1:ncol(modes)))
  
  if(cohort == "CT"){
    hearts <- c("01","02","03","04","05","06","07","08","09",10:19)
    failure <- c(9,10)
    hearts <- hearts[-failure]
    harddrive <- "/Seagate Backup Plus Drive/CT_cases"
  }
  else if(cohort == "EXTREME3"){
    hearts <- c(21:39)
    failure <- c(5,6,13,18) 
    hearts <- hearts[-failure]
    harddrive <- "/Seagate Backup Plus Drive/CT_cases"
  }
  else if(cohort == "EXTREME2"){
    hearts <- c(43,44,45,40,41,46:56,42,57)
    failure <- c(4)
    hearts <- hearts[-failure]
    harddrive <- "/Seagate Backup Plus Drive/CT_cases"
  }
  else if(cohort == "EXTREME1"){
    hearts <- c(134,135)
    failure <- c()
    harddrive <- "/Seagate Backup Plus Drive/CT_cases"
  }
  else if(cohort == "PROJECTED9"){
    hearts <- c(58:76)
    failure <- c()
    harddrive <- "/Seagate Expansion Drive"
  }
  else if(cohort == "PROJECTED5"){
    hearts <- c(77:84,86:95)
    failure <- 9
    harddrive <- "/Seagate Expansion Drive"
  }
  else if(cohort == "PROJECTED14"){
    hearts <- c(98:101,103:114)
    failure <- c(1,2,7)
    harddrive <- "/Seagate Expansion Drive"
  }
  else if(cohort == "PROJECTED18"){
    hearts <- c(116:120,122:133)
    failure <- c(1,7)  
    harddrive <- "/Seagate Expansion Drive"
  }
  
  
  if(!read_fen_from_file){
    # We initialise the variables:
    P_vent <- matrix(nrow = 801, ncol = length(hearts))
    dP_vent <- matrix(nrow = 800, ncol = length(hearts))
    V_vent <- P_vent
    Flux_vent <- P_vent
    DFlux_vent <- P_vent
    State_vent <- P_vent
    Pout_vent <- P_vent
    Qout_vent <- P_vent
    DPout_vent <- P_vent
    DPin_vent <- P_vent
    Pin_vent <- P_vent
    Qin_vent <- P_vent
    
    
    # We fill the matrices in:
    
    for(heart_num in c(1:length(hearts))){
      heart <- hearts[heart_num]
      
      if(flag_read)
        print(paste0("You need the file","/media/crg17",harddrive,"/h_case",heart,"/simulations/cohort/",heart,"HC_wk3_tpeak_120/cav.",ventricle,".csv"))
      cav.vent <- read.csv(paste0("/media/crg17",harddrive,"/h_case",heart,"/simulations/cohort/",heart,"HC_wk3_tpeak_120/cav.",ventricle,".csv"),stringsAsFactors = FALSE)
      new_phase <- Correct_phase_name(cav.vent = cav.vent[2:nrow(cav.vent),])
      aux <- cav.vent$Pressure[52:nrow(cav.vent)]
      
      P_vent[,heart_num] <- as.double(aux)
      dP_vent[,heart_num] <- P_vent[-1,heart_num] - P_vent[-nrow(P_vent),heart_num]
      V_vent[,heart_num] <- as.double(cav.vent$Volume[52:nrow(cav.vent)])
      Flux_vent[,heart_num] <- as.double(cav.vent$Flux[52:nrow(cav.vent)])
      DFlux_vent[,heart_num] <- as.double(cav.vent$D_Flux[52:nrow(cav.vent)])
      # State_vent <- mapply(trimws,cav.vent$State[52:nrow(cav.vent)])
      State_vent[,heart_num] <- new_phase[51:length(new_phase)]
      Pout_vent[,heart_num] <- as.double(cav.vent$P_out[52:nrow(cav.vent)])
      Qout_vent[,heart_num] <- as.double(cav.vent$Q_out[52:nrow(cav.vent)])
      DPout_vent[,heart_num] <- as.double(cav.vent$D_P_out[52:nrow(cav.vent)])
      DPin_vent[,heart_num] <- as.double(cav.vent$D_P_in[52:nrow(cav.vent)])
      Pin_vent[,heart_num] <- as.double(cav.vent$P_in[52:nrow(cav.vent)])
      Qin_vent[,heart_num] <- as.double(cav.vent$Q_in[52:nrow(cav.vent)])
      
    }
    
    # We convert it to dataframe columns:
    
    fenotypes_vent <- data.frame("EDP" = rep(NA,length(hearts)),
                                 "EDV" = rep(NA,length(hearts)),
                                 "Myo_vol" = rep(NA,length(hearts)),
                                 "ESV" = rep(NA,length(hearts)),
                                 "SV" = rep(NA,length(hearts)),
                                 "EF" = rep(NA,length(hearts)),
                                 "V1" = rep(NA,length(hearts)), 
                                 "EF1" = rep(NA,length(hearts)), 
                                 "ESP" = rep(NA,length(hearts)),
                                 "dPdtmax" = rep(NA,length(hearts)),
                                 "dPdtmin" = rep(NA,length(hearts)),
                                 "PeakP" = rep(NA,length(hearts)),
                                 "tpeak" = rep(NA,length(hearts)), 
                                 "ET" = rep(NA,length(hearts)), 
                                 "ICT" = rep(NA,length(hearts)),
                                 "IRT" = rep(NA,length(hearts))
    )
    
    fenotypes_vent$PeakP <- apply(P_vent,2,max)
    fenotypes_vent$tpeak <- apply(P_vent,2,'which.max')
    fenotypes_vent$dPdtmax <- apply(dP_vent,2,max)
    fenotypes_vent$dPdtmin <- apply(dP_vent,2,min)
    
    # fenotypes_vent$ICT <- apply(State_vent,2,function(x) (which(x == "IVC")[1]-1))
    fenotypes_vent$ICT <- apply(State_vent,2,function(x) (length(which(x == "IVC"))))
    fenotypes_vent$ET <- apply(State_vent,2,function(x) (length(which(x == "ejec"))))
    fenotypes_vent$IRT <- apply(State_vent,2,function(x) (length(which(x == "IVR"))))
    fenotypes_vent$tsys <- apply(State_vent,2,function(x) (length(which(x == "IVC")) + length(which(x == "ejec"))))
    
    fenotypes_vent$EDV <- V_vent[1,] 
    fenotypes_vent$EDP <- P_vent[1,]
    fenotypes_vent$ESV <- apply(V_vent,2,min)
    fenotypes_vent$ESP <- apply(-dP_vent,2,'which.max')
  aux <- apply(State_vent,2,function(x2) (which(x2 == "IVR")[1]-1))
  fenotypes_vent$ESP <- P_vent[aux]
  
  fenotypes_vent$SV <- fenotypes_vent$EDV - fenotypes_vent$ESV
  fenotypes_vent$EF <- 100 * fenotypes_vent$SV / fenotypes_vent$EDV
  
  for(heart_num in c(1:length(hearts))){
    fenotypes_vent$V1[heart_num] <- V_vent[fenotypes_vent$tpeak[heart_num],heart_num]
    if(flag_read){
      print(paste0("You need the file ","/media/crg17",harddrive,"/h_case",hearts[heart_num],"/meshing/1000um/BiV/",ventricle,"_mesh_volume.dat"))
    }
    vent_volume <- read.delim(paste0("/media/crg17",harddrive,"/h_case",hearts[heart_num],"/meshing/1000um/BiV/",ventricle,"_mesh_volume.dat"),stringsAsFactors = FALSE,header = FALSE)
    fenotypes_vent$Myo_vol[heart_num] <- sum(as.numeric(vent_volume$V1))/1000
  }
  
  fenotypes_vent$EF1 <- 100 * (fenotypes_vent$EDV - fenotypes_vent$V1) / fenotypes_vent$EDV
  
  
  
  if(with_EP){
  # We add the electrophysiology:
    if(flag_read)
      print("You need the file /media/crg17/Seagate Backup Plus Drive/CT_cases/forall/ATs_values_healthy_synthetic.txt")
  ATs_values_temp <- read.delim("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/ATs_values_healthy_synthetic.txt")
  rownames(ATs_values_temp)[1:9] <- c("01","02","03","04","05","06","07","08","09")
  
  ATs_values_healthy_synthetic <- ATs_values_temp[hearts,]
    fenotypes_vent$QRS <- as.numeric(ATs_values_healthy_synthetic$TAT)
    fenotypes_vent$AT1090 <- as.numeric(ATs_values_healthy_synthetic$AT.10.90)
    
    if(ventricle == "LV"){
      fenotypes_vent$AT <- as.numeric(ATs_values_healthy_synthetic$TAT.LV)
    }

  }
  
  
  if(action != "savefen"){
  corr_fen_vent <- fenotypes_vent[,2:ncol(fenotypes_vent)]
  }
  else{
    corr_fen_vent <- fenotypes_vent[,1:ncol(fenotypes_vent)]
  }
  }

  else{
    if(flag_read)
      print(paste0("You need the file /data/fenotypes_",ventricle,"_",cohort,".csv"))
    corr_fen_vent <- read.csv2(paste0("/data/fenotypes_",ventricle,"_",cohort,".csv"))
    if(action != "savefen"){
      corr_fen_vent <- corr_fen_vent[,-1]
    }
  }
  
  # We discard the variables we are not interested in:
  if(action != "savefen"){
  if(cohort == "CT"){
    corr_modes <- modes[1:19,]
  }
  else if(cohort == "PROJECTED9"){
    corr_modes <- modes[1:19,1:9]
  }
  else if(cohort == "PROJECTED5"){
    corr_modes <- modes[1:19,1:5]
  }
  else if(cohort == "PROJECTED14"){
    corr_modes <- modes[1:19,1:14]
  }
  else if(cohort == "PROJECTED18"){
    corr_modes <- modes[1:19,1:18]
  }
  else{
    corr_modes <- modes[hearts,]
  }
  # We create the (linear) correlation matrices:
    if(!read_fen_from_file){
      # We discard the cases that diverged
      if(length(failure) > 0 && cohort != "EXTREME3" && cohort != "EXTREME2"){
        
        corr_modes <- corr_modes[-failure,]
      }
    corr_vent <- sapply(1:ncol(corr_fen_vent), function(i,j) cor(corr_fen_vent[,i], corr_modes[,j]))
    }
    else{
      corr_vent <- sapply(as.vector(1:ncol(corr_fen_vent)), function(i,j) {
        if(length(which(is.na(corr_fen_vent[,i]))) == 0){
          cor(corr_fen_vent[,i],corr_modes[,j])
        }
        else{
          cor(corr_fen_vent[-which(is.na(corr_fen_vent[,i])),i],corr_modes[-which(is.na(corr_fen_vent[,i])),j])
        }
        }
        )
    }
    colnames(corr_vent) <- colnames(corr_fen_vent)
    rownames(corr_vent) <- paste0("Mode ",c(1:nrow(corr_vent)))
    
  if(cohort == "EXTREME2" || cohort == "EXTREME3" || cohort == "PROJECTED9")
    corr_vent = corr_vent[1:9,]
  if(cohort == "PROJECTED5")
    corr_vent = corr_vent[1:5,]
  }
  
  if(length(remove_phenotypes) > 0){
    corr_vent <- corr_vent[,-which(colnames(corr_vent) %in% remove_phenotypes)]
    }
  
  if(action != "savefen" && action != "savecorr"){
  if(action == "saveplot"){
    if(cohort == "CT"){
      png(filename=paste0("/home/crg17/Pictures/",cohort,"_",ventricle,".png"), width = 1770, height = 1600,
          units = "px", pointsize = 12, bg = "white", res = 120)
    }
    if(cohort == "PROJECTED18"){
      png(filename=paste0("/home/crg17/Pictures/",cohort,"_",ventricle,".png"), width = 1700, height = 1700,
          units = "px", pointsize = 12, bg = "white", res = 120)
    }
    if(cohort == "PROJECTED14"){
      png(filename=paste0("/home/crg17/Pictures/",cohort,"_",ventricle,".png"), width = 1700, height = 1400,
          units = "px", pointsize = 12, bg = "white", res = 120)
    }
    if(cohort == "EXTREME2" || cohort == "EXTREME3" || cohort == "PROJECTED9"){
      png(filename=paste0("/home/crg17/Pictures/",cohort,"_",ventricle,".png"), width = 1600, height = 1000,
          units = "px", pointsize = 12, bg = "white", res = 120)
    }
    if(cohort == "PROJECTED5"){
      png(filename=paste0("/home/crg17/Pictures/",cohort,"_",ventricle,".png"), width = 1600, height = 700,
          units = "px", pointsize = 12, bg = "white", res = 120)
    }
  }
  if(action == "plot" || action == "saveplot"){
  # We plot the matrices (for better results zoom it in in a new window):
    mar_left = 4
    mar_bottom = 0
    mar_right = 4
    mar_top = 8
    
    delta_left = 0
    delta_bottom = 0
    delta_right = 2
    delta_top = 0
  if(ventricle == "LV")
    corrplot(corr_vent,"square",mar=c(mar_bottom,mar_left,mar_top,mar_right),tl.cex=3,cl.cex=3,tl.col="red",tl.srt=60,cl.pos = "n")
  else if(ventricle == "RV")
    corrplot(corr_vent[,1:(ncol(corr_vent)-2)],"square",mar=c(mar_bottom + delta_bottom,mar_left + delta_left,mar_top + delta_top,mar_right + delta_right),tl.cex=3,cl.cex=3,tl.col="red",tl.srt=60)
    
  if(cohort == "CT")
    mtext(paste0("Correlation matrix for the output of\n the ", cohort," cohort for the ", ventricle), at=7, cex=4, line = -3)
  if(cohort == "EXTREME2" || cohort == "EXTREME3" || cohort == "PROJECTED9"|| cohort == "PROJECTED5")
    mtext(paste0("Correlation matrix for the ", ventricle, " of the ", cohort," cohort"), at=9, cex=3)
  if(cohort == "PROJECTED18")
    mtext(paste0("Correlation matrix for the ", ventricle, " of the ", cohort," cohort"), at=9, cex=3)
    if(cohort == "PROJECTED14")
      mtext(paste0("Correlation matrix for the ", ventricle, " of the ", cohort," cohort"), at=9.5, cex=3)
  }
  if(action == "saveplot"){
      dev.off()
  } 
  
  if(action == "return")
     return(corr_vent)
  }
  else{
    if(length(failure) > 0 && !read_fen_from_file){
      corr_fen_vent <- insertRows(corr_fen_vent,failure)
    }
    if(action == "savefen"){
      # if(cohort == "CT"){
          write.table(corr_fen_vent,paste0("/data/files4matlab/fenotypes_",ventricle,"_",cohort,".csv"),row.names = FALSE,sep = ",", dec = ".", col.names = TRUE)
      # }
      # else{
        # write.table(corr_fen_vent,paste0("/data/fenotypes_",ventricle,"_",cohort,".csv"),row.names = FALSE,sep = ";", dec = ",", col.names = TRUE)
        # }
    }
    if(action == "savecorr")
    write.table(corr_vent,paste0("/data/correlations_",ventricle,"_",cohort,".csv"),row.names = FALSE,sep = ",", dec = ".", col.names = TRUE)
    }

}

plot_difference <- function(cohort1,cohort2,with_EP = TRUE, ventricle, flag_read = FALSE){
  
  corr_plot_original <- plot_correlation(cohort = cohort1, with_EP = with_EP, action = "return",  ventricle = ventricle)
  corr_plot_projected <- plot_correlation(cohort = cohort2, with_EP = with_EP, action = "return", ventricle = ventricle)
  
  png(filename=paste0("/home/crg17/Pictures/similarity_original_projected_",ventricle,".png"), width = 1600, height = 1600,
      units = "px", pointsize = 12, bg = "white", res = 120)
  corrplot(is.corr = FALSE, 1-abs(corr_plot_original-corr_plot_projected)/2,"pie",mar=c(0,0,0.25,0),tl.cex=1.5,cl.cex=1.5,tl.col="red",tl.srt=60,cl.lim=c(0,1),col=colorRampPalette(c("white","white","black"))(200))
  mtext(paste0("Similarity of the ",ventricle," correlation matrices"), at=10, cex=3.5)
  dev.off()

}

plot_difference_original_18 <- function(flag_read = FALSE){
  
  if(flag_read){
    print("You need the file /data/fenotypes_LV_CT.csv")
    print("You need the file /data/fenotypes_RV_CT.csv")
    print("You need the file /data/fenotypes_LV_PROJECTED18.csv")
    print("You need the file /data/fenotypes_RV_PROJECTED18.csv")
  }
  fenotypes_LV_CT <- read.csv2("/data/fenotypes_LV_CT.csv")
  fenotypes_RV_CT <- read.csv2("/data/fenotypes_RV_CT.csv")
  fenotypes_LV_PROJECTED18 <- read.csv2("/data/fenotypes_LV_PROJECTED18.csv")
  fenotypes_RV_PROJECTED18 <- read.csv2("/data/fenotypes_RV_PROJECTED18.csv")
  
  colnames(fenotypes_LV_CT) <- paste0(colnames(fenotypes_LV_CT)," LV")
  colnames(fenotypes_RV_CT) <- paste0(colnames(fenotypes_RV_CT)," RV")
  colnames(fenotypes_LV_PROJECTED18) <- paste0(colnames(fenotypes_LV_PROJECTED18)," LV")
  colnames(fenotypes_RV_PROJECTED18) <- paste0(colnames(fenotypes_RV_PROJECTED18)," RV")
  
  fenotypes_CT <- cbind(fenotypes_LV_CT[,2:ncol(fenotypes_LV_CT)],fenotypes_RV_CT[,2:(ncol(fenotypes_RV_CT-2))])
  fenotypes_PROJECTED18 <- cbind(fenotypes_LV_PROJECTED18[,2:ncol(fenotypes_LV_PROJECTED18)],fenotypes_RV_PROJECTED18[,2:(ncol(fenotypes_RV_PROJECTED18-2))])
  
  fenotypes_diff <- data.frame()
  
  for (i in c(1:nrow(fenotypes_CT))) {
    for (j in c(1:ncol(fenotypes_CT))) {
      fenotypes_diff[i,j] <- 1-abs((fenotypes_CT[i,j]-fenotypes_PROJECTED18[i,j])/max(fenotypes_PROJECTED18[i,j],fenotypes_CT[i,j]))
    }
  }
  
  colnames(fenotypes_diff) <- colnames(fenotypes_PROJECTED18)
  
  fenotypes_diff <- as.matrix(fenotypes_diff)
  fenotypes_diff[is.na(fenotypes_diff)] <- 0
  
  png(filename=paste0("/home/crg17/Pictures/similarity_original_projected.png"), width = 3000, height = 2000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  corrplot(is.corr = FALSE, fenotypes_diff,"pie",mar=c(0,0,0.25,0),tl.cex=3,cl.cex=1.5,tl.col="red",tl.srt=60,cl.lim=c(0,1))
  mtext(paste0("Similarity of the original and PCA fenotypes"), at=18, cex=5)
  dev.off()
  
}

#' @description Function to create the grid-like/tables for different 
#' measurements.
#' @param ventricle "LV","RV"
#' @param table "x2", "return_x2", "range", "range_norm", "median_residuals",
#'  "cv_residuals", "sensitivity", "mean_residuals", "all"
#' @param flag_read If TRUE prints all the files read and written.
#' @param flag_return TRUE if you want to get the table rather than printing it.
#' Not available for all tables.
#' @param remove_phenotypes Vector containing the names of the phenotypes you
#' want to exclude from the table (it removes these columns).
plot_tables <- function(ventricle, table, flag_read = FALSE, 
                        flag_return = FALSE,remove_phenotypes=c()){
  
  if(table == "all"){
    plot_correlation(cohort = "CT", action = "saveplot", ventricle = "LV", read_fen_from_file = TRUE)
    plot_correlation(cohort = "CT", action = "saveplot", ventricle = "RV", read_fen_from_file = TRUE)
    plot_tables(ventricle = "LV", table = "x2")
    plot_tables(ventricle = "RV", table = "x2")
    plot_tables(ventricle = "LV", table = "mean_residuals")
    plot_tables(ventricle = "RV", table = "mean_residuals")
    plot_tables(ventricle = "LV", table = "range_norm")
    plot_tables(ventricle = "RV", table = "range_norm")
    plot_tables(ventricle = "LV", table = "sensitivity")
    plot_tables(ventricle = "RV", table = "sensitivity")
    }

  if(flag_read)
    print("You need the file /data/modes.csv")
  # modes_whole <- read.csv("/data/modes.csv", header = FALSE)
  modes_whole <- read.csv("/data/modes.csv",header = TRUE,sep = "\t")
  modes <- mapply(sd,modes_whole[1:19,2:10])
  
  if(flag_read)
    print(paste0("You need the file /data/fenotypes_",ventricle,"_CT.csv"))
  
  feno_CT <- read.csv(paste0("/data/fenotypes_",ventricle,"_CT.csv"), header = TRUE, sep = ";", dec = ",")
  feno_CT <- feno_CT[,2:ncol(feno_CT)]
  if(flag_read)
    print(paste0("You need the file /data/fenotypes_",ventricle,"_EXTREME2.csv"))
  feno_extreme2 <- read.csv(paste0("/data/fenotypes_",ventricle,"_EXTREME2.csv"), header = TRUE, sep = ";", dec = ",")
  feno_extreme2 <- feno_extreme2[,2:ncol(feno_extreme2)]
  if(flag_read)
    print(paste0("You need the file /data/fenotypes_",ventricle,"_EXTREME3.csv"))
  feno_extreme3 <- read.csv(paste0("/data/fenotypes_",ventricle,"_EXTREME3.csv"), header = TRUE, sep = ";", dec = ",")
  feno_avg <- feno_extreme3[1,2:ncol(feno_extreme3)]
  feno_extreme3 <- feno_extreme3[2:nrow(feno_extreme3),2:ncol(feno_extreme3)]
  if(flag_read)
    print(paste0("You need the file /data/fenotypes_",ventricle,"_EXTREME1.csv"))
  feno_extreme1 <- read.csv(paste0("/data/fenotypes_",ventricle,"_EXTREME1.csv"), header = TRUE, sep = ";", dec = ",")
  feno_extreme1 <- feno_extreme1[,2:ncol(feno_extreme1)]
  
  quad_terms <- NA*feno_extreme2[1:9,]
  diverged <- 0*quad_terms
  range_extreme2 <- quad_terms
  range_extreme <- quad_terms
  range_normalised <- quad_terms
  median_residuals <- quad_terms
  mean_residuals <- quad_terms
  mean_residuals_no_norm <- quad_terms
  cv_residuals <- quad_terms
  
  sensitivity <- NA*feno_CT[1:18,]
  predictions_feno <- c()
  residuals_no_norm <- c()
  
  range_CT <- lapply(feno_CT, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  range_CT <- as.data.frame(range_CT)
  avg_CT <- lapply(feno_CT, mean, na.rm = TRUE)
  
  for(mode in c(1:9)){

    for(fenotype in c(1:ncol(feno_extreme2))){

      x <- c(3*modes[mode],-3*modes[mode],0,2*modes[mode],-2*modes[mode])
      y <- c(feno_extreme3[c((mode*2)-1,mode*2),fenotype],feno_avg[1,fenotype],feno_extreme2[c((mode*2)-1,mode*2),fenotype])
      
      if(mode == 2){
        x <- c(x,modes[mode],-modes[mode])
        y <- c(y,feno_extreme1[1,fenotype],feno_extreme1[2,fenotype])
      }
      
      if(sum(is.na(y)) > 0){
        #Remove the missing values
        x <- x[-which(is.na(y))]
        y <- y[-which(is.na(y))]
      }
      
      if(table == "x2" || table == "return_x2"){
      #Non linear least squares
      nonlin_mod=nls(y~a*x^2+b*x+c,start=list(a=0,b=0,c=0),control = list(warnOnly=TRUE), na.action = na.exclude) 
      
      quad_terms[mode,fenotype] <-  environment(nonlin_mod[["m"]][["resid"]])[["env"]][["a"]]
      if(nonlin_mod$convInfo$stopCode != 0)
        diverged[mode,fenotype] <- 1
      }
      
      if(table == "range" || table == "range_norm"){
      range_extreme[mode,fenotype] <- max(y) - min(y)
      # range_normalised[mode,fenotype] <- range_extreme2[mode,fenotype]/range_CT[fenotype]
      range_normalised[mode,fenotype] <- range_extreme[mode,fenotype]/as.numeric(avg_CT[fenotype])
      }
      
      
      x <- modes_whole[1:19,mode]
      y <- feno_CT[,fenotype]
      
      if(table == "sensitivity"){
      lm_temp <- lm(y ~ x)
      slope <- lm_temp$coefficients["x"]
      # sensitivity[mode,fenotype] <- slope/(as.numeric(avg_CT[fenotype])) # Slope over the average of the sample
      sensitivity[mode,fenotype] <- slope
      }
      # x <- c(0,2*modes[mode],-2*modes[mode])
      # y <- c(feno_avg[1,fenotype],feno_extreme2[c((mode*2)-1,mode*2),fenotype])
      
      if(table %in% c("mean_residuals","median_residuals","cv_residuals")){
        x <- c(3*modes[mode],-3*modes[mode],0,2*modes[mode],-2*modes[mode])
        y <- c(feno_extreme3[c((mode*2)-1,mode*2),fenotype],feno_avg[1,fenotype],feno_extreme2[c((mode*2)-1,mode*2),fenotype])
        
        if(mode == 2){
          x <- c(x,modes[mode],-modes[mode])
          y <- c(y,feno_extreme1[1,fenotype],feno_extreme1[2,fenotype])
        }
        
        if(sum(is.na(y)) > 0){
          #Remove the missing values
          x <- x[-which(is.na(y))]
          y <- y[-which(is.na(y))]
        }
        
        lm_temp <- lm(y ~ x)
        new_y <- predict.lm(object = lm_temp, newdata = data.frame(x = modes_whole[1:19,mode]))
        predictions_feno <- cbind(predictions_feno,new_y)
        residuals_no_norm <- cbind(residuals_no_norm,abs(new_y-feno_CT[,fenotype]))
        mean_residuals_no_norm[mode,fenotype] <- mean(abs(new_y - feno_CT[,fenotype]),na.rm = TRUE)
        mean_residuals[mode,fenotype] <- mean(abs(new_y - feno_CT[,fenotype]), na.rm = TRUE)/as.numeric(avg_CT[fenotype])
        median_residuals[mode,fenotype] <- median(abs(new_y - feno_CT[,fenotype]), na.rm = TRUE)
        cv_residuals[mode,fenotype] <- 100*sd(abs(new_y - feno_CT[,fenotype]), na.rm = TRUE)/mean(abs(new_y - feno_CT[,fenotype]), na.rm = TRUE)
      }
    }
  }
  
  if(table == "sensitivity"){
    for(mode in c(10:18)){
      
      for(fenotype in c(1:ncol(feno_extreme2))){
        x <- modes_whole[1:19,mode]
        y <- feno_CT[,fenotype]
        
        lm_temp <- lm(y ~ x)
        slope <- lm_temp$coefficients["x"]

        sensitivity[mode,fenotype] <- slope
      }
    }
  }
  
  if(table == "x2"){
    
    rownames(quad_terms) <- paste0("Mode",rownames(quad_terms))
    
    png(paste0('/home/crg17/Pictures/coeffx2_',ventricle,'.png'),width=2300,height=1000,res=100)
    
    cex_value = 3
    corrplot(as.matrix(quad_terms),
             "square",
             mar=c(0,0,6,0),
             tl.cex=cex_value,
             cl.cex=cex_value,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    mtext(paste0("Quadratic term of the fit for the synthetic cohorts in the ",ventricle),
          side=3,
          cex=cex_value+1)
    
    dev.off()
  }
  if(table == "return_x2"){
    return(list(quad_terms,diverged))
  }
  if(table == "range"){
    rownames(range_extreme) <- paste0("Mode",rownames(range_extreme))
    
    png(paste0('/home/crg17/Pictures/range_',ventricle,'.png'),width=2300,height=1000,res=100)
    
    # par(mar=c(11.1, 10.1, 4.1, 4.1))
    cex_value = 3
    corrplot(as.matrix(range_extreme),
             "square",
             mar=c(0,0,8,0),
             tl.cex=cex_value,
             cl.cex=cex_value,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    mtext(paste0("Range of the synthetic cohort for the ",ventricle),
          side=3,
          cex=cex_value+1)
    
    dev.off()
  }
  if(table == "range_norm"){
    rownames(range_normalised) <- paste0("Mode ",rownames(range_normalised))
    
    if(length(remove_phenotypes) > 0){
      range_normalised <- range_normalised[,-which(colnames(range_normalised) %in% remove_phenotypes)]
    }
    
    
    if(!flag_return){
    png(paste0('/home/crg17/Pictures/range_norm_',ventricle,'.png'),width=2300,height=1000,res=100)
    
    # par(mar=c(11.1, 10.1, 4.1, 4.1))
    cex_value = 3
    corrplot(as.matrix(range_normalised),
             "square",
             mar=c(0,0,8,0),
             tl.cex=cex_value,
             cl.cex=cex_value,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    mtext(paste0("Range of the synthetic cohort normalised\n by the avg. of the CT cohort for the ",ventricle),
          side=3,
          cex=cex_value+1,
          line=-3)
    dev.off()
    }
    else{
      return(range_normalised)
    }
  }
  if(table == "mean_residuals"){
    rownames(mean_residuals) <- paste0("Mode",rownames(mean_residuals))
    
    png(paste0('/home/crg17/Pictures/mean_residuals_',ventricle,'.png'),width=2300,height=1000,res=100)
    
    cex_value = 3
    corrplot(as.matrix(mean_residuals),
             "square",
             mar=c(0,0,8,0),
             tl.cex=cex_value,
             cl.cex=cex_value,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    mtext(paste0("Normalised mean of the absolute residuals between the CT\n cohort's output and the fitting of the synthetic cohort for the ",ventricle),
          side=3,
          cex=cex_value+1,
          line=-3)
    dev.off()
    return(list(mean_residuals,feno_CT,predictions_feno,mean_residuals_no_norm,residuals_no_norm))
  }
  if(table == "median_residuals"){
    png(paste0('/home/crg17/Pictures/median_residuals_',ventricle,'.png'),width=2300,height=1000,res=100)
    par(mar=c(11.1, 10.1, 4.1, 4.1))
    cex_value = 2
    plot(as.matrix(median_residuals),
         main = "Median of the residuals of the EM output with the extreme fitting",
         fmt.cell='%.2f',
         col=rev(brewer.pal(name = 'Greens', n = 9)),
         key = list(cex.axis=cex_value),
         fmt.key = '%.0f',
         xlab = '',
         ylab = '',
         axis.col = NULL,
         axis.row = NULL,
         cex = cex_value,
         cex.main = cex_value)
    
    lablist<-paste(colnames(quad_terms),"     ")
    text(seq(1, ncol(quad_terms), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 45, pos = 1, xpd = TRUE, cex = cex_value)
    lablist <- paste0("0",1:9)
    text(0.3, 9:1, labels = lablist, srt = 0, pos = 2, xpd = TRUE, cex = cex_value)
    text(11, -1.3, labels = paste0("Functional parameters of the ", ventricle), srt = 0, pos = 2, xpd = TRUE, cex = cex_value)
    text(-0.6, 5.6, labels = "Modes", srt = 90, pos = 2, xpd = TRUE, cex = cex_value)
    
    dev.off()
  }
  if(table == "cv_residuals"){
    png(paste0('/home/crg17/Pictures/cv_residuals_',ventricle,'.png'),width=2300,height=1000,res=100)
    par(mar=c(11.1, 10.1, 4.1, 4.1))
    cex_value = 2
    plot(as.matrix(cv_residuals),
         main = "Coefficient of variation of the residuals of the EM output with the extreme fitting",
         fmt.cell='%.2f',
         col=rev(brewer.pal(name = 'Greens', n = 9)),
         key = list(cex.axis=cex_value),
         fmt.key = '%.0f',
         xlab = '',
         ylab = '',
         axis.col = NULL,
         axis.row = NULL,
         cex = cex_value,
         cex.main = cex_value)
    
    lablist<-paste(colnames(quad_terms),"     ")
    text(seq(1, ncol(quad_terms), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 45, pos = 1, xpd = TRUE, cex = cex_value)
    lablist <- paste0("0",1:9)
    text(0.3, 9:1, labels = lablist, srt = 0, pos = 2, xpd = TRUE, cex = cex_value)
    text(11, -1.3, labels = paste0("Functional parameters of the ", ventricle), srt = 0, pos = 2, xpd = TRUE, cex = cex_value)
    text(-0.6, 5.6, labels = "Modes", srt = 90, pos = 2, xpd = TRUE, cex = cex_value)
    
    dev.off()
  }
  if(table == "sensitivity"){
    
    rownames(sensitivity) <- paste0("Mode",rownames(sensitivity))
    png(filename=paste0('/home/crg17/Pictures/sensitivity_',ventricle,'.png'),
        width = 1500,
        height = 1350,
        units = "px",
        pointsize = 12,
        bg = "white",
        res = 120)
    
    cex_value = 1.8
    
    corrplot(as.matrix(sensitivity),
             "square",
             mar=c(0,0,8,0),
             tl.cex=cex_value,
             cl.cex=cex_value,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    
    mtext(paste0("Slope of the linear regression of\n the CT cohort's output for the ",ventricle),
          side=3,
          cex=cex_value+1.5,
          line=-2)
    
    dev.off()
  }
}

preprocess_fenotypes_GP <- function(action = "save", remove_phenotypes = c(),flag_read = FALSE){
  
  cohorts <- c("CT","EXTREME3","EXTREME2","EXTREME1")
  modes <- read.table(file = "/data/modes.csv",header = TRUE)

  GP_output <- data.frame()
  x_output <- data.frame()
  
  for(coh in cohorts){
    # We read the LV dataframe
    if(flag_read)
      print(paste0("You need the file /data/fenotypes_LV_",coh,".csv"))
    fen_df <- read.csv2(paste0("/data/fenotypes_LV_",coh,".csv"))
    # Remove the EDP
    fen_df <- fen_df[,-1]
    # We read the RV dataframe
    if(flag_read)
      print(paste0("You need the file /data/fenotypes_RV_",coh,".csv"))
    fen_df_RV <- read.csv2(paste0("/data/fenotypes_RV_",coh,".csv"))
    # We remove the EDP, QRS and AT1090
    fen_df_RV <- fen_df_RV[,-c(1,18,19)]
    
    colnames(fen_df) <- paste0(colnames(fen_df),"_LV")
    colnames(fen_df_RV) <- paste0(colnames(fen_df_RV),"_RV")
    
    # We paste them together horizontally
    fen_df <- cbind(fen_df,fen_df_RV)
    
    #We name the rows
    rownames(fen_df) <- paste0(coh,"_",rownames(fen_df))
    
    #We delete the ones with NAs
    rows2delete <- apply(fen_df,1,function(x) sum(is.na(x)) > 0)
    rows2delete <- which(rows2delete == TRUE)
    
    #Subset of the modes
    x_temp <- modes[modes$Cohort == coh,2:ncol(modes)]
    
    if(length(rows2delete) > 0){
      fen_df <- fen_df[-rows2delete,]
      x_temp <- x_temp[-rows2delete,]
    }
    
    if(length(remove_phenotypes) > 0)
      fen_df <- fen_df[,-which(names(fen_df) %in% remove_phenotypes)]
  
    GP_output <- rbind(GP_output,fen_df)
    x_output <- rbind(x_output,x_temp)
    
  }
  if(action == "return")
    return(GP_output)
  if(action == "save"){
    write.table(GP_output,"/data/PCA_GP/y.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
    write.table(colnames(GP_output),"/data/PCA_GP/y_labels.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(rownames(GP_output),"/data/PCA_GP/heart_labels.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(x_output,"/data/PCA_GP/x.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
    write.table(colnames(x_output),"/data/PCA_GP/modes_labels.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
}

difference_boxplot <- function(flag_read = FALSE){
  
  if(flag_read)
    print("/data/fenotypes_LV_CT.csv")
  fenotypes_LV_CT <- read.csv2("/data/fenotypes_LV_CT.csv")
  
  df_CT_LV <- as.data.frame(fenotypes_LV_CT[,2])
  colnames(df_CT_LV) <- "Value"
  df_CT_LV$Fenotype <- paste0(colnames(fenotypes_LV_CT)[2],"_LV")
  
  for(j in c(3:ncol(fenotypes_LV_CT))){
    temp_df <- as.data.frame(fenotypes_LV_CT[,j])
    colnames(temp_df) <- "Value"
    temp_df$Fenotype <- paste0(colnames(fenotypes_LV_CT)[j],"_LV")
    
    df_CT_LV <- rbind(df_CT_LV,temp_df)
  }
  df_CT_LV$Cohort <- "CT"
  
  if(flag_read){
    print("You need the file /data/fenotypes_LV_PROJECTED18.csv")
  }
  fenotypes_LV_PROJECTED18 <- read.csv2("/data/fenotypes_LV_PROJECTED18.csv")
  
  df_18_LV <- as.data.frame(fenotypes_LV_PROJECTED18[,2])
  colnames(df_18_LV) <- "Value"
  df_18_LV$Fenotype <- paste0(colnames(fenotypes_LV_PROJECTED18)[2],"_LV")
  
  for(j in c(3:ncol(fenotypes_LV_PROJECTED18))){
    temp_df <- as.data.frame(fenotypes_LV_PROJECTED18[,j])
    colnames(temp_df) <- "Value"
    temp_df$Fenotype <- paste0(colnames(fenotypes_LV_PROJECTED18)[j],"_LV")
    
    df_18_LV <- rbind(df_18_LV,temp_df)
  }
  df_18_LV$Cohort <- "PROJECTED18"
  
  if(flag_read)
    print("You need the file /data/fenotypes_RV_CT.csv")
  fenotypes_RV_CT <- read.csv2("/data/fenotypes_RV_CT.csv")
  
  df_CT_RV <- as.data.frame(fenotypes_RV_CT[,2])
  colnames(df_CT_RV) <- "Value"
  df_CT_RV$Fenotype <- paste0(colnames(fenotypes_LV_CT)[2],"_RV")
  
  for(j in c(3:(ncol(fenotypes_RV_CT)-2))){
    temp_df <- as.data.frame(fenotypes_RV_CT[,j])
    colnames(temp_df) <- "Value"
    temp_df$Fenotype <- paste0(colnames(fenotypes_RV_CT)[j],"_RV")
    
    df_CT_RV <- rbind(df_CT_RV,temp_df)
  }
  df_CT_RV$Cohort <- "CT"
  
  if(flag_read)
    print("You need the file /data/fenotypes_RV_PROJECTED18.csv")
  fenotypes_RV_PROJECTED18 <- read.csv2("/data/fenotypes_RV_PROJECTED18.csv")
  
  df_18_RV <- as.data.frame(fenotypes_RV_PROJECTED18[,2])
  colnames(df_18_RV) <- "Value"
  df_18_RV$Fenotype <- paste0(colnames(fenotypes_RV_PROJECTED18)[2],"_RV")
  
  for(j in c(3:(ncol(fenotypes_RV_PROJECTED18)-2))){
    temp_df <- as.data.frame(fenotypes_RV_PROJECTED18[,j])
    colnames(temp_df) <- "Value"
    temp_df$Fenotype <- paste0(colnames(fenotypes_RV_PROJECTED18)[j],"_RV")
    
    df_18_RV <- rbind(df_18_RV,temp_df)
  }
  df_18_RV$Cohort <- "PROJECTED18"
  
  df_boxplot <- rbind(df_CT_LV,df_CT_RV,df_18_LV,df_18_RV)
  
  p <- ggplot(data = df_boxplot, aes(x=Fenotype, y=Value)) + geom_boxplot(aes(fill=Cohort))
  p <- p + facet_wrap( ~ Fenotype, scales="free")
  
  show(p)
  
  return(df_boxplot)
}

Show_EF <- function(path_to_cavfile, flag_read = FALSE){
    
  if(flag_read)
    print(paste0("You need the file ",path_to_cavfile,"/cav.LV.csv"))
  
    cav.vent <- read.csv(paste0(path_to_cavfile,"/cav.LV.csv"),stringsAsFactors = FALSE)
    volumes <- as.double(cav.vent$Volume[52:nrow(cav.vent)])
    EDV <- volumes[1]
    ESV <- min(volumes)
    EF <- 100*(EDV - ESV)/EDV
    
    print("LV EF:")
    print(EF)
    
    if(flag_read)
      print(paste0("You need the file ",path_to_cavfile,"/cav.RV.csv"))
    
    cav.vent <- read.csv(paste0(path_to_cavfile,"/cav.RV.csv"),stringsAsFactors = FALSE)
    volumes <- as.double(cav.vent$Volume[52:nrow(cav.vent)])
    EDV <- volumes[1]
    ESV <- min(volumes)
    EF <- 100*(EDV - ESV)/EDV
    
    print("RV EF:")
    print(EF)
    
}

Check_negative_R2 <- function(GP_folder,ksplit=5){
  
  r2scores <- read.table(paste0("/data/PCA_GP/",GP_folder,"/r2scores.txt"), quote="\"", comment.char="")
  r2scores <- r2scores$V1
  
  y_labels <- read.table(paste0("/data/PCA_GP/",GP_folder,"/y_labels.txt"), quote="\"", comment.char="")
  y_labels <- y_labels$V1
  
  r2_split <- split(r2scores, ceiling(seq_along(r2scores)/ksplit))
  r2_matrix<-as.matrix(as.data.frame(r2_split))
  
  indices <- which(apply(r2_matrix, 2, mean)<0)
  return(list(y_labels[indices],paste0(length(indices==TRUE),"/",length(y_labels))))
}

#' @description Script to plot the results of the GSA as rectangles where
#' the size of each rectange is proportional to the percentage of variance
#' explained by each mode.
#' 
#' @param GP_folder Name of the folder of the GSA.
#' @param r2score If TRUE, it reads the values of the r2 score, otherwise it
#' uses the MSE.
#' @param flag_debugging If TRUE, it plots whatever is reading and writing.

Plot_GSA_tree <- function(GP_folder,r2score=TRUE,flag_debugging = FALSE){

  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  Load_Install_Packages(c("dplyr","treemap"))
  
  y_labels <- paste0("/data/PCA_GP/",GP_folder,"/y_labels.txt") %>%
                Read_table(., quote="\"", comment.char="", as.is = TRUE,
                           flag_debugging = flag_debugging)
  y_labels <- y_labels$V1
  
  if(r2score){
    r2scores <- paste0("/data/PCA_GP/",GP_folder,"/r2scores.txt") %>%
    Read_table(., quote="\"", comment.char="",
               flag_debugging = flag_debugging)
  }
  else{
    r2scores <- paste0("/data/PCA_GP/",GP_folder,"/mse_scores.txt") %>%
      Read_table(., quote="\"", comment.char="", flag_debugging = flag_debugging)
  }
  r2scores <- r2scores$V1
  
  r2_split <- split(r2scores, ceiling(seq_along(r2scores)/round(length(r2scores)/length(y_labels))))
  r2_matrix<-as.matrix(as.data.frame(r2_split))
  
  for(i in 1:length(y_labels)){
  
  fo_effects <- paste0("/data/PCA_GP/",GP_folder,"/",i-1,"/Si.txt") %>%
    Read_table(., quote="\"", comment.char="", flag_debugging = flag_debugging)

    total_effects <- paste0("/data/PCA_GP/",GP_folder,"/",i-1,"/STi.txt") %>%
      Read_table(., quote="\"", comment.char="", flag_debugging = flag_debugging)
    
  all_effects_df <- as.data.frame(matrix(nrow=19,ncol=2)) 
  colnames(all_effects_df) <- c("Mode","Mean")
  
  all_effects_df[1:18,1] <- paste0("Mode ",1:18)
  all_effects_df[19,1] <- "Higher order effects"
  
  all_effects_df[1:18,2] <- apply(fo_effects,2,function(x) max(0,mean(x)))
  # all_effects_df[19,2] <- max(0,1-sum(all_effects_df[1:18,2]))
  all_effects_df[19,2] <- max(0,sum(apply(total_effects,2,mean)-all_effects_df[1:18,2]))
  # Plot
  
  # original_palette <- heat.colors(19,rev = TRUE)
  original_palette_fun <- colorRampPalette(c("#ffd040", "#ed2424"))
  original_palette <- original_palette_fun(19)
  actual_order <- c(19,1,10:18,2:9)
  color_palette <- original_palette[actual_order]
    
    
  reorder <- c(2,4,5,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
  if(r2score)
    r2_value <- round(max(r2_matrix[,i]),2)
  else
    r2_value <- round(min(r2_matrix[,i]),2)
  if(mean(r2_matrix[,i]) < 0)
    r2_value <- round(mean(r2_matrix[,i]),2)
  
  if(r2score)
    title_plot <- paste0("Sensitivity of ",substr(y_labels[i],1,nchar(y_labels[i])-3)," in the ",substr(y_labels[i],nchar(y_labels[i])-1,nchar(y_labels[i]))," (r2=",r2_value,")")
  else
    title_plot <- paste0("Sensitivity of ",substr(y_labels[i],1,nchar(y_labels[i])-3)," in the ",substr(y_labels[i],nchar(y_labels[i])-1,nchar(y_labels[i]))," (MSE=",r2_value,")")
  
  png(paste0('/home/crg17/Pictures/GSA_',GP_folder,'_',y_labels[i],'.png'),width=1600,height=1000,res=100)
  
  treemap(all_effects_df,
          
          # data
          index=c("Mode","Mode"),
          vSize="Mean",
          type="index",
          
          # Main
          title=title_plot,
          fontsize.title = 60,
          aspRatio = 1,
          palette=color_palette,
          # palette=rev(heat.colors(19)),
          
          # Borders:
          border.col=c("black"),             
          border.lwds=1.5,                         
          
          # Labels
          fontsize.labels=0.7,
          fontcolor.labels="black",
          fontface.labels=1,            
          bg.labels="transparent", # Transparent background labels
          align.labels=c("centre","centre"),                                  
          overlap.labels=0.5,
          inflate.labels=T                        # If true, labels are bigger when rectangle is bigger.
  )
  
  
  dev.off()
  }
}

#' @description Script to plot the results of the GSA as a corrplot.
#' 
#' @param GP_folder Name of the folder of the GSA.
#' @param ventricle "LV" or "RV".
#' @param MSE_flag If FALSE, it reads the values of the r2 score, otherwise it
#' uses the MSE.
#' @param show_MSE If TRUE, it creates a row with the values of the MSE/r2score.
#' @param flag_return If TRUE, it also returns the matrix.

Plot_GSA_table <- function(GP_folder,ventricle,MSE_flag=TRUE,show_MSE=TRUE,flag_return=FALSE){
  
  y_labels_original <- read.table(paste0("/data/PCA_GP/",GP_folder,"/y_labels.txt"),
                                  quote="\"", comment.char="", as.is = TRUE)
  y_labels_original <- y_labels_original$V1
  
  if(MSE_flag)
    MSE_values <- read.table(paste0("/data/PCA_GP/",GP_folder,"/mse_scores.txt"), quote = "\"", comment.char = "")
  else
    MSE_values <- read.table(paste0("/data/PCA_GP/",GP_folder,"/r2scores.txt"), quote = "\"", comment.char = "")
  
  MSE_values <- MSE_values$V1
  
  MSE_split <- split(MSE_values, ceiling(seq_along(MSE_values)/round(length(MSE_values)/length(y_labels_original))))
  MSE_matrix<-as.matrix(as.data.frame(MSE_split))
  
  if(MSE_flag)
    final_MSE <- as.vector(apply(MSE_matrix, 2, min))
  else
    final_MSE <- as.vector(apply(MSE_matrix, 2, max))
  
  idx <- which(as.vector(sapply(y_labels_original,function(x) substr(x,nchar(x)-1,nchar(x)))) == ventricle)
  y_labels <- as.vector(sapply(y_labels_original[idx],function(x) substr(x,1,nchar(x)-3)))
  final_MSE <- final_MSE[idx]
  
  all_effects_df <- as.data.frame(matrix(nrow=20,ncol=length(y_labels)))
  if(MSE_flag)
    rownames(all_effects_df)[1] <- "MSE"
  else
    rownames(all_effects_df)[1] <- "r2"
  rownames(all_effects_df)[2:19] <- paste0("Mode ",1:18)
  rownames(all_effects_df)[20] <- "Higher order effects"
  colnames(all_effects_df) <- y_labels
  
  for(i in 1:length(idx)){
    
    fo_effects <- read.table(paste0("/data/PCA_GP/",GP_folder,"/",idx[i]-1,"/Si.txt"), quote="\"", comment.char="")
    total_effects <- read.table(paste0("/data/PCA_GP/",GP_folder,"/",idx[i]-1,"/STi.txt"),quote="\"",comment.char = "")
    
    all_effects_df[1,i] <- 0
    all_effects_df[2:19,i] <- apply(fo_effects,2,function(x) max(0,mean(x)))
    all_effects_df[20,i] <- max(0,sum(apply(total_effects,2,mean)-all_effects_df[2:19,i]))
  }
  
  MSE_matrix <- as.matrix(0*all_effects_df)
  MSE_matrix[1,] <- final_MSE
  
  if(!flag_return){
  if(show_MSE){
  
  png(paste0('/home/crg17/Pictures/GSA_',ventricle,'_',GP_folder,'.png'),
      width=1800,height=1500,res=100)
  
  # cex_value = 3
  # corrplot(as.matrix(all_effects_df),
  #          "square",
  #          p.mat = MSE_matrix,
  #          insig = "p-value",
  #          sig.level = min(MSE_matrix[1,])/2,
  #          mar=c(0,0,8,0),
  #          tl.cex=cex_value,
  #          cl.cex=cex_value,
  #          number.cex = 10,
  #          cl.lim = c(0,1),
  #          cl.ratio = 0.3,
  #          tl.col="red",
  #          tl.srt=60,
  #          is.corr = FALSE)
  # mtext(paste0("Variance explained by each one of the modes \nusing Gaussian processes in the ",ventricle),
  #       side=3,
  #       cex=cex_value+1,
  #       line=-3,
  #       at=7)
  
  cex_value <- 3
  mag.factor <- 1.3
  cex.before <- par("cex")
  par(cex = 1.4)
  corrplot(as.matrix(all_effects_df),
           "square",
           p.mat = MSE_matrix,
           insig = "p-value",
           sig.level = min(MSE_matrix[1,])/2,
           mar=c(0,0,8,0),
           tl.cex=par("cex") * mag.factor,
           cl.cex=par("cex") * mag.factor,
           number.cex = 10,
           cl.lim = c(0,1),
           cl.ratio = 0.2,
           tl.col="red",
           tl.srt=60,
           is.corr = FALSE)
  mtext(paste0("Variance explained by each one of the modes \nusing Gaussian processes in the ",ventricle),
        side=3,
        cex=cex_value+1,
        line=-3,
        at=7)
  par(cex = cex.before)
  dev.off()
  }
  
  else{
    png(paste0('/home/crg17/Pictures/GSA_',ventricle,'_',GP_folder,'.png'),
        width=1800,height=1500,res=100)
    
    cex_value <- 3
    mag.factor <- 1.3
    cex.before <- par("cex")
    par(cex = 1.4)
    corrplot(as.matrix(all_effects_df[-1,]),
             "square",
             mar=c(0,0,8,0),
             tl.cex=par("cex") * mag.factor,
             cl.cex=par("cex") * mag.factor,
             number.cex = 10,
             cl.lim = c(0,1),
             cl.ratio = 0.2,
             tl.col="red",
             tl.srt=60,
             is.corr = FALSE)
    mtext(paste0("Variance explained by each one of the modes \nusing Gaussian processes in the ",ventricle),
          side=3,
          cex=cex_value+1,
          line=-3,
          at=7)
    par(cex = cex.before)
    dev.off()
  }
  }
  if(flag_return){
    return(all_effects_df)
  }
    
}

#' @description Function to create the grid-like images of the manuscript with
#' the same definition and image parameters.
#' @param output: "GSA", "range_norm", "correlation", "LSA" or "all"
#' @param venricle: "LV" or "RV"

Print_PNG <- function(output = "GSA",ventricle = "LV"){
  
  Load_Install_Packages("corrplot")
  
  if(output == "all"){
    Print_PNG(output = "range_norm", ventricle = "LV")
    Print_PNG(output = "range_norm", ventricle = "RV")
    Print_PNG(output = "correlation", ventricle = "LV")
    Print_PNG(output = "correlation", ventricle = "RV")
    Print_PNG(output = "GSA", ventricle = "LV")
    Print_PNG(output = "GSA", ventricle = "RV")
  }
  else{
  set_width = 4000
  set_height = 4000
  set_res = 300
  
  if(output == "GSA"){
  file_name <- paste0('/home/crg17/Pictures/GSA_',ventricle,'.png')
  title_name <- paste0("Variance explained by each one of the modes \nusing Gaussian processes in the ",ventricle)
  
  matrix <- Plot_GSA_table(GP_folder="tov7_final",ventricle=ventricle,
                           MSE_flag=TRUE,show_MSE=FALSE,flag_return=TRUE)
  matrix <- as.matrix(matrix[2:(nrow(matrix)),])
  
  is_corr_flag <- FALSE
  colour_lim <- c(0,1)
  if(ventricle == "LV"){
    low_title <- -6
    left_title <- 7
  }
  else{
    low_title <- -5
    left_title <- 6
  }
  }
  else if(output == "correlation"){
  file_name <- paste0('/home/crg17/Pictures/CT_',ventricle,'.png')
  title_name <- paste0("Correlation matrix for the output of\n the CT cohort for the ", ventricle)
  
  if(ventricle == "LV"){
    matrix <- plot_correlation(cohort = "CT", action = "return",
                             ventricle = ventricle, read_fen_from_file = TRUE,
                             remove_phenotypes = c("EDV","Myo_vol","ICT","IRT","ESP","EF"))
  }
  else if(ventricle == "RV"){
    matrix <- plot_correlation(cohort = "CT", action = "return",
                               ventricle = ventricle, read_fen_from_file = TRUE,
                               remove_phenotypes = c("EDV","Myo_vol","ICT","IRT"))
    matrix <-  matrix[,1:(ncol(matrix)-2)]
  }
  is_corr_flag <- FALSE
  colour_lim <- c(-1,1)
  low_title <- -3
  left_title <- 7
  }
  else if(output == "range_norm"){
    file_name <- paste0('/home/crg17/Pictures/range_norm_',ventricle,'.png')
    title_name <- paste0("Range of the synthetic cohort normalised\n by the avg. of the CT cohort for the ",ventricle)

    if(ventricle == "LV"){
      matrix <- plot_tables(ventricle = ventricle, table = "range_norm",flag_return = T,
                            remove_phenotypes = c("EDV","Myo_vol","ICT","IRT","ESP","EF"))
      matrix <- as.matrix(matrix)
      low_title <- -13
      left_title <- 7
    }
    if(ventricle == "RV"){
      matrix <- plot_tables(ventricle = ventricle,table = "range_norm",flag_return = T,
                            remove_phenotypes = c("EDV","Myo_vol","ICT","IRT"))
      matrix <- as.matrix(matrix[,1:(ncol(matrix)-2)])
      low_title <- -12
      left_title <- 6
    }
    is_corr_flag <- F
    colour_lim <- c(-0.24,1.18)

    if(ventricle == "LV"){
      low_title <- -14
      left_title <- 7
    }
  }
  else if(output == "LSA"){
    file_name <- paste0('/home/crg17/Pictures/LSA_',ventricle,'.png')
    title_name <- paste0("Sensitivity coefficients for the ",ventricle)
    
    if(ventricle == "LV"){
      remove_phenotypes <- c("EDV","Myo_vol","ICT","IRT","ESP","EF")
      low_title <- -8
      left_title <- 5
    }
    else if(ventricle == "RV"){
      remove_phenotypes <- c("EDV","Myo_vol","ICT","IRT")
      low_title <- -9
      left_title <- 6
    }
    matrix <- Plot_LSA(ventricle = ventricle,
                       remove_phenotypes = remove_phenotypes)
    matrix <- as.matrix(matrix)
    
    is_corr_flag <- F
  }
  
  colour_lim <- c(min(-1,min(matrix)),max(1,max(matrix)))
  
  png(file_name,width=set_width,height=set_height,res=set_res)
  
  label_size <- 2.5
  title_size <- 3.5

 set_margin = c(0,0,10,10) # bottom, left, top, right

   corrplot(matrix,"square",mar=set_margin,
           tl.cex=label_size, # Labels
           tl.col="black", # Labels colors
           tl.srt=60, # Label angle
           cl.cex=label_size, # Colorbar
           cl.lim = colour_lim, # Colorbar limits
           cl.ratio = 0.2, # Colorbar width
           is.corr = is_corr_flag)
  mtext(title_name, # Title text
        side=3, # Margin of the text
        cex=title_size, # Size 
        line=low_title, # How low
        at=left_title) # How to the left

    dev.off()
  }
 
}

#' @description Function to create a corrplot of the LSA.
#' 
#' @param ventricle "LV" or "RV".
#' @param flag_debugging If TRUE, it prints whatever is reading and writing.
#' @param remove_phenotypes Vector with the phenotypes not wanted included in
#' the matrix for returning or plotting.
#' 
#' @return The matrix of the LSA.

Plot_LSA <- function(ventricle, flag_debugging = FALSE, remove_phenotypes = c()){
  
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  
  Read_phenotypes <- function(path2LSA, ventricle, flag_debugging = parent.frame()$flag_debugging){
    
      cav.vent <- Read_csv(paste0(path2LSA,"/cav.",ventricle,".csv"),
                           stringsAsFactors = FALSE,
                           flag_debugging = parent.frame()$flag_debugging)
      new_phase <- Correct_phase_name(cav.vent = cav.vent[2:nrow(cav.vent),])
      aux <- cav.vent$Pressure[52:nrow(cav.vent)]
      
      P_vent <- as.double(aux)
      dP_vent <- P_vent[-1] - P_vent[-length(P_vent)]
      V_vent <- as.double(cav.vent$Volume[52:nrow(cav.vent)])
      Flux_vent <- as.double(cav.vent$Flux[52:nrow(cav.vent)])
      DFlux_vent <- as.double(cav.vent$D_Flux[52:nrow(cav.vent)])
      # State_vent <- mapply(trimws,cav.vent$State[52:nrow(cav.vent)])
       State_vent <- new_phase[51:length(new_phase)]
      Pout_vent <- as.double(cav.vent$P_out[52:nrow(cav.vent)])
      Qout_vent <- as.double(cav.vent$Q_out[52:nrow(cav.vent)])
      DPout_vent <- as.double(cav.vent$D_P_out[52:nrow(cav.vent)])
        DPin_vent <- as.double(cav.vent$D_P_in[52:nrow(cav.vent)])
      Pin_vent <- as.double(cav.vent$P_in[52:nrow(cav.vent)])
      Qin_vent <- as.double(cav.vent$Q_in[52:nrow(cav.vent)])
      
    
    # We convert it to dataframe columns:
    
    phenotypes_vent <- data.frame("EDV" = NA,
                                  "ESV" = NA,
                                  "SV" = NA,
                                  "EF" = NA,
                                  "V1" = NA,
                                  "EF1" = NA,
                                  "ESP" = NA,
                                  "dPdtmax" = NA,
                                  "dPdtmin" = NA,
                                  "PeakP" = NA,
                                  "tpeak" = NA,
                                  "ET" = NA,
                                  "ICT" = NA,
                                  "IRT" = NA)
    
    phenotypes_vent$EDV <- V_vent[1] 
    phenotypes_vent$PeakP <- max(P_vent)
    phenotypes_vent$tpeak <- which.max(P_vent)
    phenotypes_vent$dPdtmax <- max(dP_vent)
      phenotypes_vent$dPdtmin <- min(dP_vent)
    
    # phenotypes_vent$ICT <- which(State_vent == "IVC")[1]-1
    phenotypes_vent$ICT <- length(which(State_vent == "IVC"))
    phenotypes_vent$ET <- length(which(State_vent == "ejec"))
    phenotypes_vent$IRT <- length(which(State_vent == "IVR"))
    phenotypes_vent$tsys <- length(which(State_vent == "IVC")) +
                            length(which(State_vent == "ejec"))
    
    phenotypes_vent$ESV <- min(V_vent)
    phenotypes_vent$ESP <- which.max(-dP_vent)
    aux <- which(State_vent == "IVR")[1]-1
    phenotypes_vent$ESP <- P_vent[aux]
    
    phenotypes_vent$SV <- phenotypes_vent$EDV - phenotypes_vent$ESV
    phenotypes_vent$EF <- 100 * phenotypes_vent$SV / phenotypes_vent$EDV
    
    phenotypes_vent$V1 <- V_vent[phenotypes_vent$tpeak]
    phenotypes_vent$EF1 <- 100 * (phenotypes_vent$EDV - phenotypes_vent$V1) / phenotypes_vent$EDV
    
    return(phenotypes_vent)
      
  }
  
  LSA_path <- paste0("/media/crg17/Seagate Expansion Drive/h_case21/",
                     "simulations/cohort/LSA")
  
  phenotypes <- Read_phenotypes(paste0(LSA_path,"/../21HC_wk3_tpeak_120"), "LV")
  avg_parameters <- data.frame("CV" = 0.8, "CVFEC" = 5.6, 
                               "AV_resistance" = 0.03,
                               "PV_resistance" = 0.015,
                               "lvEDP_kPa" = 1.6, "rvEDP_kPa" = 0.8, 
                               "aortic_press" = 77, "PA_press_mmHg" = 17.4,
                               "BC_stiffness" = 0.001, "t_peak" = 120,
                               "t_dur" = 550, "scaling_Guccione" = 1.7, 
                               "c_neohook" = 7.45)
  
  
  for(param in colnames(avg_parameters)){
    folder_name <- paste0(LSA_path,"/LSA_wk3_",param,"_",
                          format(0.9*avg_parameters[1,param], scientific = FALSE))
    phenotypes[nrow(phenotypes)+1,] <- Read_phenotypes(folder_name,
                                                       ventricle = ventricle)
    
    folder_name <- paste0(LSA_path,"/LSA_wk3_",param,"_",
                          format(1.1*avg_parameters[1,param], scientific = FALSE))
    phenotypes[nrow(phenotypes)+1,] <- Read_phenotypes(folder_name,
                                                       ventricle = ventricle)
  }
  
  LSA_matrix <- matrix(nrow = length(avg_parameters), ncol = ncol(phenotypes))
  
  
  yprev <- phenotypes[seq(2, nrow(phenotypes), by = 2),]
  ynext <- phenotypes[seq(3, nrow(phenotypes), by = 2),]
  y0 <- t(phenotypes[1,])
  x0 <- t(avg_parameters[1,])
  dx <- 0.1*x0
  # Derivative with central difference scheme
  LSA_matrix <- (ynext - yprev)/(2*dx)
  #Normalise with x0/y0
  LSA_matrix <- x0*LSA_matrix
  LSA_matrix <- t(apply(as.matrix(LSA_matrix), 1, function(xx) xx/as.matrix(y0)))
  
  rownames(LSA_matrix) <- c("CV", "CV FEC","AV resistance", "PV resistance",
                            "LV pressure", "RV pressure",
                            "Aortic pressure", "PA pressure",
                            "BC stiffness", "Peak isometric tension",
                            "Duration of transient", "a Guccione", "c Neohookean")
  colnames(LSA_matrix) <- colnames(phenotypes)
  rownames(yprev) <- colnames(avg_parameters)
  rownames(ynext) <- colnames(avg_parameters)
  
  
  
  if(length(remove_phenotypes) > 0){
    LSA_matrix <- LSA_matrix[,-which(colnames(LSA_matrix) %in% remove_phenotypes)]
  }
  
  return(LSA_matrix)
  
  }

#' @description Function to return a matrix of the values of the LSA with the 
#' parameters changed qualitatively.
#' 
#' @param ventricle "LV" or "RV".
#' @param flag_debugging If TRUE, it prints whatever is reading and writing.
#' @param remove_phenotypes Vector with the phenotypes not wanted included in
#' the matrix for returning or plotting.
#' 
#' @return The matrix of the LSA.

get_matrix_LSA2 <- function(ventricle, flag_debugging = FALSE,
                            remove_phenotypes = c()){
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  
  Read_phenotypes <- function(path2LSA, ventricle,
                              flag_debugging = parent.frame()$flag_debugging){
    
    cav.vent <- Read_csv(paste0(path2LSA,"/cav.",ventricle,".csv"),
                         stringsAsFactors = FALSE,
                         flag_debugging = parent.frame()$flag_debugging)
    new_phase <- Correct_phase_name(cav.vent = cav.vent[2:nrow(cav.vent),])
    aux <- cav.vent$Pressure[52:nrow(cav.vent)]
    
    P_vent <- as.double(aux)
    dP_vent <- P_vent[-1] - P_vent[-length(P_vent)]
    V_vent <- as.double(cav.vent$Volume[52:nrow(cav.vent)])
    Flux_vent <- as.double(cav.vent$Flux[52:nrow(cav.vent)])
    DFlux_vent <- as.double(cav.vent$D_Flux[52:nrow(cav.vent)])
    # State_vent <- mapply(trimws,cav.vent$State[52:nrow(cav.vent)])
    State_vent <- new_phase[51:length(new_phase)]
    Pout_vent <- as.double(cav.vent$P_out[52:nrow(cav.vent)])
    Qout_vent <- as.double(cav.vent$Q_out[52:nrow(cav.vent)])
    DPout_vent <- as.double(cav.vent$D_P_out[52:nrow(cav.vent)])
    DPin_vent <- as.double(cav.vent$D_P_in[52:nrow(cav.vent)])
    Pin_vent <- as.double(cav.vent$P_in[52:nrow(cav.vent)])
    Qin_vent <- as.double(cav.vent$Q_in[52:nrow(cav.vent)])
    
    
    # We convert it to dataframe columns:
    
    phenotypes_vent <- data.frame("EDV" = NA,
                                  "ESV" = NA,
                                  "SV" = NA,
                                  "EF" = NA,
                                  "V1" = NA,
                                  "EF1" = NA,
                                  "ESP" = NA,
                                  "dPdtmax" = NA,
                                  "dPdtmin" = NA,
                                  "PeakP" = NA,
                                  "tpeak" = NA,
                                  "ET" = NA,
                                  "ICT" = NA,
                                  "IRT" = NA)
    
    phenotypes_vent$EDV <- V_vent[1] 
    phenotypes_vent$PeakP <- max(P_vent)
    phenotypes_vent$tpeak <- which.max(P_vent)
    phenotypes_vent$dPdtmax <- max(dP_vent)
    phenotypes_vent$dPdtmin <- min(dP_vent)
    
    # phenotypes_vent$ICT <- which(State_vent == "IVC")[1]-1
    phenotypes_vent$ICT <- length(which(State_vent == "IVC"))
    phenotypes_vent$ET <- length(which(State_vent == "ejec"))
    phenotypes_vent$IRT <- length(which(State_vent == "IVR"))
    phenotypes_vent$tsys <- length(which(State_vent == "IVC")) +
      length(which(State_vent == "ejec"))
    
    phenotypes_vent$ESV <- min(V_vent)
    phenotypes_vent$ESP <- which.max(-dP_vent)
    aux <- which(State_vent == "IVR")[1]-1
    phenotypes_vent$ESP <- P_vent[aux]
    
    phenotypes_vent$SV <- phenotypes_vent$EDV - phenotypes_vent$ESV
    phenotypes_vent$EF <- 100 * phenotypes_vent$SV / phenotypes_vent$EDV
    
    phenotypes_vent$V1 <- V_vent[phenotypes_vent$tpeak]
    phenotypes_vent$EF1 <- 100 * (phenotypes_vent$EDV - phenotypes_vent$V1) / phenotypes_vent$EDV
    
    return(phenotypes_vent)
    
  }
  
  
  path2LSA <- paste0("/media/crg17/Seagate Expansion Drive/h_case21/",
                     "simulations/cohort/LSA")
  total_path <- paste0(path2LSA, "/LSA_wk3_CV_0.72")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- phen_result
  
  total_path <- paste0(path2LSA, "/LSA_wk3_CV_0.88")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_CVFEC_5.04")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_CVFEC_6.16")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_-55_75")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_-65_85")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_L0RN")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_L0RW")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_LNR0")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_fibres_LWR0")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/21HC_biv_wk3")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  total_path <- paste0(path2LSA, "/LSA_wk3_purkinje_on")
  phen_result <- Read_phenotypes(total_path, ventricle = ventricle)
  whole_df <- rbind(whole_df, phen_result)
  
  rownames(whole_df) <- c("Lower CV", "Higher CV","Lower CVFEC", "Higher CVFEC",
                          "LNRN", "LWRW", "L0RN", "L0RW", "LNR0", "LWR0", "BiV",
                          "purkinje")
  
  if(ventricle == "LV"){
  path2volumefile <-  paste0("/media/crg17/Seagate Expansion Drive/h_case21/",
                             "meshing/1000um/BiV")
  BiVvol <- "BiV_mesh_volume.dat"
  ATfilename <- "vm_act_seq_biv_elem.dat"
  tagfile <- "tags.dat"
  
  total_path <- paste0(path2LSA, "/LSA_EP_CV_0.72")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["Lower CV","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)

  whole_df["Lower CV","AT1090"] <- EP_list[[2]]
  whole_df["Lower CV","AT"] <- EP_list[[1]]
  
  
  total_path <- paste0(path2LSA, "/LSA_EP_CV_0.88")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["Higher CV","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["Higher CV","AT1090"] <- EP_list[[2]]
  whole_df["Higher CV","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_CVFEC_5.04")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["Lower CVFEC","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["Lower CVFEC","AT1090"] <- EP_list[[2]]
  whole_df["Lower CVFEC","AT"] <- EP_list[[1]]
  
  
  total_path <- paste0(path2LSA, "/LSA_EP_CVFEC_6.16")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["Higher CVFEC","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["Higher CVFEC","AT1090"] <- EP_list[[2]]
  whole_df["Higher CVFEC","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_-55_75")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["LNRN","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["LNRN","AT1090"] <- EP_list[[2]]
  whole_df["LNRN","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_-65_85")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["LWRW","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["LWRW","AT1090"] <- EP_list[[2]]
  whole_df["LWRW","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_L0RN")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["L0RN","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["L0RN","AT1090"] <- EP_list[[2]]
  whole_df["L0RN","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_L0RW")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["L0RW","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["L0RW","AT1090"] <- EP_list[[2]]
  whole_df["L0RW","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_LNR0")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["LNR0","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["LNR0","AT1090"] <- EP_list[[2]]
  whole_df["LNR0","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/LSA_EP_fibres_LWR0")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["LWR0","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["LWR0","AT1090"] <- EP_list[[2]]
  whole_df["LWR0","AT"] <- EP_list[[1]]
  
  total_path <- paste0(path2LSA, "/../../../meshing/1000um/purkinje")
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = "both",
                                              flag_debugging = flag_debugging)
  whole_df["purkinje","QRS"] <- EP_list[[1]]
  
  EP_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                              volumefilename = BiVvol,
                                              path2ATfile = total_path,
                                              ATfilename = ATfilename,
                                              path2tagfile = path2volumefile,
                                              tagfilename = tagfile,
                                              ventricle = ventricle,
                                              flag_debugging = flag_debugging)
  
  whole_df["purkinje","AT1090"] <- EP_list[[2]]
  whole_df["purkinje","AT"] <- EP_list[[1]]
  }
  
  total_path <- paste0(path2LSA,"/../21HC_wk3_tpeak_120")
  avg_pheno <- Read_phenotypes(total_path, ventricle = ventricle)
  
  if(ventricle == "LV"){
    avg_pheno$QRS <- 76.0493
    avg_pheno$AT1090 <- 35.9436
    avg_pheno$AT <- 74.0493
  }
  
  percentage_change <-t(apply(whole_df,1,function(x) 100*(x - t(avg_pheno))/t(avg_pheno)))
  colnames(percentage_change) <- colnames(whole_df)
  
  if(length(remove_phenotypes) > 0){
    percentage_change <- percentage_change[,-which(colnames(percentage_change) %in% remove_phenotypes)]
  }
  return(percentage_change)
}