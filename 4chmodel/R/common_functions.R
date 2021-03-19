#' @description Loads the needed packages and install them if they are not downloaded.
#' @param list.of.packages Vector with the strings of the name of the packages
Load_Install_Packages <- function(list.of.packages){
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
    flags<-lapply(list.of.packages, require, character.only = TRUE)
}


#' @description It prints a message of the file reading or writing.
#' @param string Directory and name of the file.
#' @param io "r" if reading the file, "w" if writing it.
#' @param flag_debugging If FALSE it does not run anything.
Debug_message <- function(string,io,flag_debugging=TRUE){
  if(flag_debugging){
    Load_Install_Packages(c("dplyr"))
    
    if(io == "r"){
      paste0("You need the file ",string) %>%
        print(.)
    }
    else if(io == "w"){
      paste0("You are getting the file ",string) %>%
        print(.)
    }
  }
}

#' @description Same as read.csv function with the Debug_message function added.
#' 
#' @param ... Same as the ones in read.csv
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.
Read_csv <- function(file, header = TRUE, sep = ",", quote = "\"",dec = ".",
                     fill = TRUE, comment.char = "", flag_debugging = FALSE,
                     stringsAsFactors = default.stringsAsFactors()){
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  
  Debug_message(string = file,io = "r",flag_debugging = parent.frame()$flag_debugging)
  
  file_read <- read.csv(file = file, header = header, sep = sep, quote = quote,
                        dec = dec, fill = fill, comment.char = comment.char,
                        stringsAsFactors = stringsAsFactors)
  return(file_read)
}


#' @description The phase naming of CARP didn't seem to be correct (easy to 
#' check by plotting the volume and pressure curves), so this functions is to
#' correct that.
#' 
#' @param cav.vent Data frame from the output of the CARP simulations (cav.LV
#' or cav.RV) to correct the phase.
#' 
#' @return The data frame with the phase naming corrected. 
Correct_phase_name <- function(cav.vent){
    volume <- as.double(cav.vent$Volume)
    phase <- NA*volume
    
    phase[1:50] <- "load"
    phase[51:length(phase)] <- "IVC"
    
    position <- 52
    
    while(abs(volume[position + 10] - volume[position])/10 < 0.1){
      position <- position + 1
    }
    
    phase[position:length(phase)] <- "ejec"
    
    while(abs(volume[position + 10] - volume[position])/10 > 0.1){
      position <- position + 1
    }
    
    phase[position:length(phase)] <- "IVR"
    
    while(abs(volume[position + 10] - volume[position])/10 <= 0.1){
      position <- position + 1
    }
    
    phase[position:length(phase)] <- "fill"
    

    return(phase)
    
}