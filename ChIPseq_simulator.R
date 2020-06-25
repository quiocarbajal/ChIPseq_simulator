library(tidyverse)
library(gridExtra)
library(ggpmisc)


simulate_chip <- function(prot_mean=500, prot_sd=30, prot_length=40, no_cut_l=0, no_cut_r=0,
                          cut_l=0, cut_r=0, repetitions=1000, length=1000, read_length=75,
                          dna_top_size=250, dna_bottom_size=75, plot = TRUE, intensity=1) { 
  # prot_mean       : the coordinate in which the protein center can be found (mean of normal distribution)
  # prot_length     : length of DNA covered by the protein, this is a "no-cut" zone.
  # no_cut_l        : coordinates, relative to the left side of the protein, in which no cuts are generated  (input as "c(-10,0)", beeing 0 the left side of the protein)
  # no_cut_r        : the same as above but coodinates are positive, with 0 being the rigt side of the protein
  # cut_l / cut_r   : coordinates (relative to the left-most or right-most side of the protein) of a fixed cut to the left or the right of the protein respectively
  # repetitions     : number of reads that are going to be included in the final table (discarded reads are...discarded, hence, they are not included in the final table)
  # length          : length of the coordinate system
  # read_length     : length of sequencing, determined by the sequencer (75, 150, etc)
  # dna_top_size    : max DNA length. Any random piece of DNA above this threshold will be discarded
  # dna_bottom_size : minimal DNA length. Any random piece of DNA below this threshold will be discarded. If you are purifying your adapter-ligated chipped DNA with ampure beads using 1X buffer, the max size you keep is about 200-150 bp, minus the adapter (75 bp), you get DNA pieces of at least 75-125 bp long 
  # plot            : Do you want to get the plots automatically? You might want to change it to "FALSE" for composite peaks
  # intensity       : Intensity of the peaks, for doing composite peaks
  
  i <- 1
  #reads_table <-  as_tibble(matrix(ncol = length, nrow = repetitions)) #This was before it gave Fw and Rv strand information
  reads_table_F <- matrix(nrow = repetitions, ncol = length)
  reads_table_R <- matrix(nrow = repetitions, ncol = length)
  reads_table_F[] <- 0
  reads_table_R[] <- 0
  coordinates <- 1:length
  colnames(reads_table_F) <- coordinates
  colnames(reads_table_R) <- coordinates
  
  while (i < repetitions) {
    #determine prot coordinates
    prot_center <- as.integer(rnorm(1, mean=prot_mean, sd=prot_sd))
    prot_left <-  prot_center - (prot_length/2)
    prot_right <- prot_left + prot_length
    
    # determine left and right cuts
    
    left_cut <- sample(1:(prot_left-1),1)
    if (cut_l > 0) {
        stop('cut_l should be a negative integer, representing the position of the cut relative to the leftmost protein boundary')
      }
    if (prot_left + cut_l > left_cut) {
        left_cut <- prot_left + cut_l
    }
    
    
    right_cut <- sample((prot_right+1):length,1)
    if (cut_r < 0) {
      stop('cut_r should be a positive integer, representing the position of the cut relative to the leftmost protein boundary')
      }
    if (prot_right + cut_r < right_cut) {
        right_cut <- prot_right + cut_r
      }
    
    #Determine if cuts are in no cut zone
    
    cut_restriction_left <- TRUE

    if (length(no_cut_l)>1) {
      if (no_cut_l[1] > no_cut_l[2]) {
        stop('no_cut_l should have the following format: \'c(-x,y)\', being x smaller than y,
             and both are the coordinate of the no cut zone relative to the protein position')
        }
      no_cut_l_distal <- prot_left + no_cut_l[1]
      no_cut_l_proximal <- prot_left + no_cut_l[2]
      if ((no_cut_l_distal < left_cut) & (left_cut < no_cut_l_proximal)) {
      cut_restriction_left <- FALSE
      }
    }
    
    cut_restriction_right <- TRUE  
    if (length(no_cut_r)>1) {
      if (no_cut_r[1] > no_cut_r[2]) {
          stop('no_cut_r should have the following format: \'c(x,y)\', being x smaller than y, 
               and both are the coordinate of the no cut zone relative to the protein position')
      }
      no_cut_r_distal <- prot_right + no_cut_r[1]
      no_cut_r_proximal <- prot_right + no_cut_r[2]
      
      if (no_cut_r_distal < right_cut & right_cut < no_cut_r_proximal) {
      cut_restriction_right <- FALSE
      }
    }
    
    # Keep read if DNA size < top_DNA_size and if left and right cut are not within no-cut zones, if so, 
    # fill read length with intensity values
    if ( (right_cut - left_cut) < dna_top_size & 
         (right_cut - left_cut) > dna_bottom_size &
         left_cut > 1 &
         right_cut < length &
         cut_restriction_left & 
         cut_restriction_right) {
      reads_table_F[i,left_cut:(left_cut+read_length)] <- intensity
      reads_table_R[i,(right_cut-read_length):right_cut] <- intensity
      i = i+1
    }
    # If left cut is of boundaries...
    if ( (right_cut - left_cut) < dna_top_size & 
         (right_cut - left_cut) > dna_bottom_size &
         left_cut < 1 &
         left_cut + read_length > 1 &
         right_cut < length &
         cut_restriction_left & 
         cut_restriction_right) {
      
      reads_table_F[i,1:(left_cut+read_length)] <- intensity
      reads_table_R[i,(right_cut-read_length):right_cut] <- intensity
      i = i+1
    }
    # If right cut is of boundaries...
    if ( (right_cut - left_cut) < dna_top_size & 
         (right_cut - left_cut) > dna_bottom_size &
         left_cut > 1 &
         right_cut > length &
         right_cut - read_length < length &
         cut_restriction_left & 
         cut_restriction_right) {
      
      reads_table_F[i,left_cut:(left_cut+read_length)] <- intensity
      reads_table_R[i,(right_cut-read_length):length] <- intensity
      i = i+1
    }
  }
  
  #Calculate average value per coordinate
  i <- 1
  average_table_F <- matrix(nrow = 2, ncol = length)
  average_table_R <- matrix(nrow = 2, ncol = length)
  
  rownames(average_table_F) <- c("coordinates","Fw")
  rownames(average_table_R) <- c("coordinates","Rv")
  
  average_table_F[1,] <- coordinates
  average_table_R[1,] <- coordinates
  
  while (i<(length+1)) {
    average_table_F[2,i] <- mean(reads_table_F[,i])
    average_table_R[2,i] <- mean(reads_table_R[,i])
    i = i+1
  }
  
  # Transform matrixes to tibbles so it can be used in ggplot (I am not sure if this is actually necessary)
  average_table_F <- as_tibble(t(average_table_F))
  average_table_R <- as_tibble(t(average_table_R))
  
  # merge both tibbles, calculate sum of Fw and Rv, re-shape tibble so that it is useful for ggplot
  average_table <- left_join(average_table_F, average_table_R) %>% mutate("Fw_plus_Rv" = Fw + Rv) %>%
                  pivot_longer(col= -coordinates, names_to = "Strand", values_to = "Average_signal")
  #Change the order of factors so when graphed in ggplot Fw_plus_Rv comes first (and does not cover the other 2)
  average_table$Strand <- factor(average_table$Strand, levels =c("Fw_plus_Rv","Fw","Rv"))

  #Get parameters in a table
  parameters <- tibble("prot_mean"= prot_mean, "prot_sd"= prot_sd, "prot_length"= prot_length, 
                       "no_cut_l" = no_cut_l, "no_cut_r"= no_cut_r, "cut_l"= cut_l,"cut_r"=cut_r,
                       "repetitions"= repetitions, "length"= length, "read_length"= read_length, 
                       "dna_top_size"= dna_top_size, "dna_bottom_size" = dna_bottom_size)
  
  #Create a graph
  if (plot) {
    
  
  y_axis_max <- average_table %>% filter(Strand=="Fw_plus_Rv") %>% select("Average_signal") %>% max()
  
  
  table1 <- tibble(x=0,y= y_axis_max + (y_axis_max*0.2), tb=list(parameters))
  
  print(ggplot(data = average_table, aes(x = coordinates, y = Average_signal, color = Strand)) +
    geom_line() +
    geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  }

  # Return average and parameters tables
    return(list("average_table" =average_table, "parameters" = parameters))
  
    # Add the line below to return also the matrixes
    # "reads_matrix_F" = reads_table_F, "reads_matrix_R" = reads_table_R,
}

#################################################################
############### SIMULATE COMPOSITE PEAKS FUNCTION ###############
#################################################################

simulate_composite_peaks <- function(parameters){
  #parameters is a list containing all the parameters that can be modified in the simulate_peak function. Use following model (Cmnd+Shift+c to uncomment all the lines together):
    # a <- list(prot_mean = c(500,550), prot_sd = c(30,30), prot_length = c(35,35), 
    # no_cut_l = c(c(0,0),c(0,0)), no_cut_r = c(c(0,0),c(0,0)), cut_l = c(0,0), cut_r = c(0,1),
    # repetitions = c(2000,2000), length = c(2000,2000), read_length = c(75,75), 
    # dna_top_size = c(250,250), dna_bottom_size = c(85,85), plot = c(FALSE,FALSE),
    # number_peaks = 2, names=c("Peak1", "Peak2"))
  
  # Call "simulate_peak" function for every peak in the parameters list
  results <- list()
  for (i in (1:parameters$number_peaks[[1]])) {
    results[[ parameters$names[[i]] ]] <- simulate_chip(prot_mean = parameters$prot_mean[[i]],
                                               prot_sd = parameters$prot_sd[[i]],
                                               prot_length = parameters$prot_length[[i]],
                                               no_cut_l = parameters$no_cut_l[[i]],
                                               no_cut_r = parameters$no_cut_r[[i]],
                                               cut_l = parameters$cut_l[[i]],
                                               cut_r = parameters$cut_r[[i]],
                                               repetitions = parameters$repetitions[[i]],
                                               length = parameters$length[[i]],
                                               read_length = parameters$read_length[[i]],
                                               dna_top_size = parameters$dna_top_size[[i]],
                                               dna_bottom_size =parameters$dna_bottom_size[[i]],
                                               plot = parameters$plot[[i]],
                                               intensity = parameters$intensity[[i]])
    # Add the peak name info to each result table
    results[[c(i,1)]] <- results[[c(i,1)]] %>% mutate(peak_name = parameters$names[[i]])
  }
  # Create tibbles to combine peaks information
  composite_result <- tibble()
  composite_parameters <- tibble()
  
  for (n in (1:parameters$number_peaks[[1]]) ) {
    composite_result <- bind_rows(composite_result, results[[c(n,1)]])
    composite_parameters <- bind_rows(composite_parameters, results[[c(n,2)]])
    n <- n+1
  }
  
  # Get average signal for composite peak
  composite_peak <- composite_result %>% group_by(coordinates, Strand) %>% 
                    summarise("Average_signal" = sum(Average_signal), peak_name = "composite_peak") %>% ungroup() #Is ungroup necessary?
  composite_result <- bind_rows(composite_result, composite_peak)
  
  # Print plots
  y_axis_max <- composite_result %>% filter(Strand=="Fw_plus_Rv") %>% select("Average_signal") %>% max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb=list(composite_parameters))
  print(ggplot(data = composite_result, aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
    geom_line() +
    geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  
  y_axis_max <- composite_result %>% filter(Strand=="Fw_plus_Rv", peak_name == "composite_peak") %>% select("Average_signal") %>% max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb=list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name == "composite_peak"), 
               aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
          geom_line() +
          geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  
  y_axis_max <- composite_result %>% filter(Strand=="Fw_plus_Rv", peak_name != "composite_peak") %>% select("Average_signal") %>% max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb=list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name != "composite_peak"), 
               aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
          geom_line() +
          geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  
  #geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  
  # Return table
  return(composite_result)
}

############################################################################################################
### EXAMPLE OF SINGLE PEAK
############################################################################################################
result <- simulate_chip(repetitions=1000, prot_mean = 1000, prot_sd =40, prot_length = 35,
                        read_length = 75, dna_top_size=300, dna_bottom_size = 85, length=2000, intensity=1)

############################################################################################################
### EXAMPLE OF COMPOSITE PEAK
############################################################################################################
####### 2 peaks #############
a <- list(prot_mean = c(500,520), prot_sd = c(70,70), prot_length = c(50,50),
          no_cut_l = c(c(0,0),c(0,0)), no_cut_r = c(c(0,0),c(0,0)), cut_l = c(0,-10), cut_r = c(10,0),
          repetitions = c(2000,2000), length = c(1000,1000), read_length = c(75,75),
          dna_top_size = c(300,300), dna_bottom_size = c(85,85), plot = c(FALSE,FALSE),
          number_peaks = 2, names=c("Peak1", "Peak2"), intensity=c(1,2))


x <- simulate_composite_peaks(a)

#trim33 profile is reproduced when the distance between the peak is slightly smaller than the size of the protein

####### 3 peaks #############
a <- list(prot_mean = c(2000,3000,4000), prot_sd = c(200,200,200), prot_length = c(40,40,40), 
          no_cut_l = c(0,0,0), no_cut_r = c(0,0,0), cut_l = c(0,0,-1), cut_r = c(1,0,0),
          repetitions = c(2000,2000,2000), length = c(6000,6000,6000), read_length = c(75,75,75), 
          dna_top_size = c(600,600,600), dna_bottom_size = c(85,85,85), plot = c(FALSE,FALSE,FALSE),
          number_peaks = 3, names=c("Peak1", "Peak2","a"), intensity=c(4,1,4))
x <- simulate_composite_peaks(a)


####### 4 peaks #############
a <- list(prot_mean = c(1000,1500,1500,2000), prot_sd = c(100,200,200,100), prot_length = c(40,40,40,40), 
          no_cut_l = c(0,0,0,0), no_cut_r = c(0,0,0,0), cut_l = c(0,0,0,-1), cut_r = c(1,0,0,0),
          repetitions = c(2000,2000,2000,2000), length = c(6000,6000,6000,6000), read_length = c(75,75,75,75), 
          dna_top_size = c(350,300,300,300), dna_bottom_size = c(85,85,85,85), plot = c(FALSE,FALSE,FALSE,FALSE),
          number_peaks = 4, names=c("Peak1", "Peak2","a","b"), intensity=c(4,1,1,4))
x <- simulate_composite_peaks(a)
