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
  source("./common_functions.R")
  
  Debug_message(string = file,io = "r",flag_debugging = parent.frame()$flag_debugging)
  
  file_read <- read.csv(file = file, header = header, sep = sep, quote = quote,
                        dec = dec, fill = fill, comment.char = comment.char,
                        stringsAsFactors = stringsAsFactors)
  return(file_read)
}

#' @description Same as write.csv function with the Debug_message function
#' added.
#' 
#' @param ... Same as the ones in write.csv
#' @param flag_debugging If TRUE, prints whatever is writtig
#' 
Write_csv <- function(x, file = "", append = FALSE, quote = TRUE, sep = " ",
                      eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                      col.names = TRUE, qmethod = c("escape", "double"),
                      fileEncoding = "", flag_debugging = FALSE){
  source("./common_functions.R")
  
  Debug_message(string = file,io = "w",
                flag_debugging = parent.frame()$flag_debugging)
  
  write.csv(x = x, file = file, append = append, quote = quote, sep = sep,
            eol = eol, na = na, dec = dec, row.names = row.names,
            col.names = col.names, qmethod = qmethod,
            fileEncoding = fileEncoding)
}

#' @description Same as write.table function with the Debug_message function
#' added.
#' 
#' @param ... Same as the ones in write.table
#' @param flag_debugging If TRUE, prints whatever is writtig
#' 
Write_table <- function(x, file = "", append = FALSE, quote = TRUE, sep = " ",
                        eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "", flag_debugging = FALSE){
  source("./common_functions.R")
  
  Debug_message(string = file,io = "w",
                flag_debugging = parent.frame()$flag_debugging)
  
  write.table(x = x, file = file, append = append, quote = quote, sep = sep,
            eol = eol, na = na, dec = dec, row.names = row.names,
            col.names = col.names, qmethod = qmethod,
            fileEncoding = fileEncoding)
}

#' @description Same as read.delim function with the Debug_message function
#' added.
#' 
#' @param ... Same as the ones in read.delim
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.
Read_delim <- function(file, header = TRUE, sep = "\t", quote = "\"",
                       dec = ".", fill = TRUE, comment.char = "",
                       flag_debugging = FALSE){
  source("./common_functions.R")
  
  Debug_message(string = file,io = "r",
                flag_debugging = parent.frame()$flag_debugging)
  
  file_read <- read.delim(file = file, header = header, sep = sep,
                          quote = quote, dec = dec, fill = fill,
                          comment.char = comment.char)
  return(file_read)
}

#' @description Same as read.csv2 function with the Debug_message function
#' added.
#' 
#' @param ... Same as the ones in read.csv2
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.
Read_csv2 <- function(file, header = TRUE, sep = ";", quote = "\"",
                      dec = ",", fill = TRUE, comment.char = "",
                      flag_debugging = FALSE){
  source("./common_functions.R")
  
  Debug_message(string = file,io = "r",
                flag_debugging = parent.frame()$flag_debugging)
  
  file_read <- read.csv2(file = file, header = header, sep = sep,
                          quote = quote, dec = dec, fill = fill,
                          comment.char = comment.char)
  return(file_read)
}

#' @description Same as read.table function with the Debug_message function
#' added.
#' 
#' @param ... Same as the ones in read.csv2
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.
Read_table <- function(file, header = FALSE, sep = "", quote = "\"'", dec = ".",
                       numerals = c("allow.loss", "warn.loss", "no.loss"),
                       row.names, col.names, as.is = !stringsAsFactors,
                       na.strings = "NA", colClasses = NA, nrows = -1,
                       skip = 0, check.names = TRUE, fill = !blank.lines.skip,
                       strip.white = FALSE, blank.lines.skip = TRUE,
                       comment.char = "#",
                       allowEscapes = FALSE, flush = FALSE,
                       stringsAsFactors = default.stringsAsFactors(),
                       fileEncoding = "", encoding = "unknown", text,
                       skipNul = FALSE, flag_debugging = FALSE){
  source("./common_functions.R")
  
  Debug_message(string = file,io = "r",
                flag_debugging = parent.frame()$flag_debugging)
  
  file_read <- read.table(file = file, header = header, sep = sep,
                          quote = quote, dec = dec, numerals = numerals,
                          row.names = row.names, col.names = col.names,
                          as.is = as.is, na.strings = na.strings,
                          colClasses = colClasses, nrows = nrows, skip = skip,
                          check.names = check.names, fill = fill,
                          strip.white = strip.white,
                          blank.lines.skip = blank.lines.skip,
                          comment.char = comment.char,
                          allowEscapes = allowEscapes, flush = flush,
                          stringsAsFactors = stringsAsFactors,
                          fileEncoding = fileEncoding, encoding = encoding,
                          text = text, skipNul = skipNul)
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