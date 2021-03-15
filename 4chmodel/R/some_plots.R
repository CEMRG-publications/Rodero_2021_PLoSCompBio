#' @description Function to plot the individual and explained variance by the 
#' PCA modes. The input data is hardcoded.
#' @return A png image. 

Plot_variance_modes <- function(){
  
  Individual_explained_variance = c(32.6772271,19.0222101,13.5231487,7.24375708, 
                                   5.17672572,3.89260612,3.6237595,2.35576423, 
                                   2.27320087,1.81989127,1.63818581,1.50912978, 
                                   1.2496146,0.986255296,0.946770655,0.790044095, 
                                   0.668580993,0.603128055)
  
  Cumulative_explained_variance = c(32.677227,51.699437,65.222586,72.466343, 
                                   77.643069,81.535675,85.159434,87.515199, 
                                   89.788399,91.608291,93.246477,94.755606, 
                                   96.005221,96.991476,97.938247,98.728291, 
                                   99.396872,100)

 png(filename=paste0("/home/crg17/Pictures/explained_variance.png"),
      width = 1920, height = 1003,
     units = "px", pointsize = 12, bg = "white", res = 120)
  twoord.plot(1:18,Individual_explained_variance,
              1:18,Cumulative_explained_variance,
              type = c("bar","b"),
              xlab = "Principal component",
              xtickpos = 1:18,
              xticklab = 1:18,
              main = list("PCA explained variance",cex=4),
              ylab = "Individual (%)",
              lcol = "#808080",
              lylim = c(0,1.1*max(Individual_explained_variance)),
              lytickpos = seq(0,35,by=5),
              lwd = 4,
              rylab = "Cumulative (%)",
              rcol = "#b32430",
              rylim = c(30,105),
              rpch = 4,
              rytickpos = seq(30,100,by=10),
              cex = 2, # Size of the crosses
              axislab.cex = 2.1)
  dev.off()
  
}

#' @description Script to plot as doughnut charts the results of the GSA. 
#' Aesthetically they look like core-less pie charts, but it's basically stacked
#' rectangles in polar coordinates.
#' @param gp_folder Name of the folder of the GP results. 
#' @return A png image for each phenotype + a plot for the legend
Plot_doughnut_chart <- function(gp_folder = "tov7_final", flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  Load_Install_Packages(c("dplyr","ggplot2","magrittr"))
  
  y_labels <- Read_table(paste0("/data/PCA_GP/",gp_folder,"/y_labels.txt"),
                         quote="\"", comment.char="",
                         flag_debugging = flag_debugging, as.is = TRUE)
  y_labels <- y_labels$V1
  
  
  for(i in 1:length(y_labels)){
  # i=18

    fo_effects <- Read_table(paste0("/data/PCA_GP/",gp_folder,"/",i-1,
                            "/Si.txt"), quote="\"", comment.char="", 
                            flag_debugging = flag_debugging)
    total_effects <- Read_table(paste0("/data/PCA_GP/",gp_folder,"/",i-1,
                            "/STi.txt"),quote="\"",comment.char = "",
                            flag_debugging = flag_debugging)
    
    all_effects_df <- as.data.frame(matrix(nrow=11,ncol=3)) 
    colnames(all_effects_df) <- c("Mode","Mean","Percentage")
    
    all_effects_df[1:9,1] <- paste0("Mode ",1:9)
    all_effects_df[10,1] <- "Modes 10-18"
    all_effects_df[11,1] <- "Multifactorial effects"
    
    all_modes_effects <- apply(fo_effects,2,function(x) max(0,mean(x)))
    
    all_effects_df[1:9,2] <- all_modes_effects[1:9]
    all_effects_df[10,2] <- sum(all_modes_effects[10:18])
    all_effects_df[11,2] <- (apply(total_effects,2,mean)-all_modes_effects) %>%
      sum(.) %>% max(0,.)
    
    all_effects_df[,3] <- all_effects_df[,2]/sum(all_effects_df[,2])
    
    # In previous iterations we neededto plot only a subset (chunks)
    
    indices_sorted <- order(all_effects_df$Percentage, decreasing = TRUE)
    # doughnut_dataframe <- all_effects_df[indices_sorted[1:(chunks-1)],
    #                                      c("Mode","Percentage")]
    doughnut_dataframe <- all_effects_df[indices_sorted,c("Mode","Percentage")]
    # doughnut_dataframe[chunks,1] <- "Others"
    # doughnut_dataframe[chunks,2] <- all_effects_df[indices_sorted[chunks:19],
    #                                                "Percentage"] %>% sum(.)
    # For colour issues, the names need to be in alphabetical order so we add
    # a 0 to all the single cipher
    single_cipher_idx <- as.numeric(rownames(doughnut_dataframe)) < 10
        # single_cipher_idx <- c(single_cipher_idx,FALSE) # So we don't change the last name
    doughnut_dataframe$Mode[single_cipher_idx] <-  rownames(doughnut_dataframe)[single_cipher_idx] %>%
                                                  paste0("Mode 0",.)

    
    # We fixed the colours for each mode 
    # all_colours <- c("#67001f",
    #                  "#8e0c25",
    #                  "#b0212f",
    #                  "#c7433e",
    #                  "#d96754",
    #                  "#e98b6f",
    #                  "#f6ae8d",
    #                  "#facab1",
    #                  "#fde1d2",
    #                  "#fff5f0",
    #                  "#f6f6f6",
    #                  "#e6e6e6",
    #                  "#d4d4d4",
    #                  "#c0c0c0",
    #                  "#a7a7a7",
    #                  "#8c8c8c",
    #                  "#6e6e6e",
    #                  "#515151",
    #                  "#353535",
    #                  "#1a1a1a")
    
    all_colours <- c("#67001f",
                     "#8e0c25",
                     "#b0212f",
                     "#c7433e",
                     "#d96754",
                     "#e98b6f",
                     "#f6ae8d",
                     "#facab1",
                     "#fde1d2",
                     "#c0c0c0",
                     "#1a1a1a")
    
  
    # final_indices <- rownames(doughnut_dataframe)[1:(nrow(doughnut_dataframe)-1)] %>%
    #                  as.numeric(.) %>% as.vector(.)
    final_indices <- rownames(doughnut_dataframe) %>% as.numeric(.) %>%
                      as.vector(.)

     final_colours <- all_colours[final_indices] # We select the colors
    final_colours <- final_colours[order(final_indices)]# We sort them alphabetically
    final_colours <- c(final_colours,all_colours[length(all_colours)])

  # Make the plot

    
ggplot(doughnut_dataframe,
       aes(x = 1, y = Percentage, fill = Mode)) + # x is size of the circle
  geom_bar(stat = "identity", color = "white", size = 3) + # Contour
  coord_polar(theta = "y", start = 0)+ # Changes from rectangles to circle
  scale_fill_manual(values = final_colours) + # Colours
  theme_void()+ # Remove background
  xlim(-2, 2) + # Lower inferior limit is thinner. Upper has to be bigger than
                # x in aes
  labs(title = paste0(substr(y_labels[i],
                             1,
                             nchar(y_labels[i])-3),
                      " (",
                      substr(y_labels[i],
                             nchar(y_labels[i])-1,
                             nchar(y_labels[i])),
                      ")")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 100, face = "bold", hjust = 0.5))

  ggsave(
    filename = paste0("/home/crg17/Pictures/doughnut_",y_labels[i],".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 30,
    height = 30,
    units = "cm",
    dpi = 300
  )

   }
  legend_df <- data.frame(Mode = c(paste0("0",c(1:9)),"10-18","Multifactorial"),
                          Percentage = rep(1/11,11))
  
  ggplot(legend_df,
         aes(x = 1, y = Percentage, fill = Mode)) + # x is size of the circle
    geom_bar(stat = "identity", color = "white", size = 3) + # Contour
    coord_polar(theta = "y", start = 0)+ # Changes from rectangles to circle
    scale_fill_manual(values = all_colours) + # Colours
    theme_void()+ # Remove background
    xlim(-2, 2) + # Lower inferior limit is thinner. Upper has to be bigger than
    # x in aes
    labs(title = "Relative position", fill = "Shape Mode") +
    theme(plot.title = element_text(size = 100, face = "bold", hjust = 0.5))+
    theme(legend.text = element_text(size = 85), 
          legend.key.width = unit(2.6,"cm"),
          legend.title = element_text(size = 85))
  
  ggsave(
    filename = paste0("/home/crg17/Pictures/legend.png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 60,
    height = 60,
    units = "cm",
    dpi = 300
  )
  
  }