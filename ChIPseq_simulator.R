library(tidyverse)
library(gridExtra)
library(ggpmisc)
# Added protein position matrix

simulate_chip <- function(prot_mean = 500, prot_sd = 30, prot_length = 15, no_cut_l = c(0,0), no_cut_r = c(0,0),
                          cut_l = 0, cut_r = 0, repetitions = 1000, length = 1000, read_length = 75,
                          number_cuts_per_kb = 4, dna_bottom_size = 85, plot = TRUE, intensity = 1, PE_full_length_read = FALSE) { 
  # prot_mean       : the coordinate in which the protein center can be found (mean of normal distribution)
  # prot_length     : length of DNA covered by the protein, this is a "no-cut" zone.
  # no_cut_l        : coordinates, relative to the left side of the protein, in which no cuts are generated  (input as "c(-10,0)", being 0 the left side of the protein)
  # no_cut_r        : the same as above but coordinates are positive, with 0 being the right side of the protein
  # cut_l / cut_r   : coordinates (relative to the left-most or right-most side of the protein) of a fixed cut to the left or the right of the protein respectively. 0 = No cut
  # repetitions     : number of reads that are going to be included in the final table (discarded reads are...discarded, hence, they are not included in the final table)
  # length          : length of the coordinate system
  # read_length     : length of sequencing, determined by the sequencer (75, 150, etc)
  # dna_mean_size   : DNA mean size (normal distribution)
  # dna_size_sd     : DNA size standard deviation (normal distribution)
  # dna_bottom_size : minimal DNA length. Any random piece of DNA below this threshold will be discarded. If you are purifying your adapter-ligated chipped DNA with ampure beads using 1X buffer, the max size you keep is about 200-150 bp, minus the adapter (75 bp), you get DNA pieces of at least 75-125 bp long 
  # plot            : Do you want to get the plots automatically? You might want to change it to "FALSE" for composite peaks
  # intensity       : Intensity of the peaks, for doing composite peaks
  # PE_full_length_read : If TRUE, then the matrix is filled with "intensity" values all over the length of the piece of DNA, not just at "read_length" on each side. With this, Forward, reverse and Fw+Rv strands profiles will be equal.
  
  # Set index to count succesful iterations
  i <- 1
  
  # Set matrixes to store results
  reads_table_F <- matrix(nrow = repetitions, ncol = length)
  reads_table_R <- matrix(nrow = repetitions, ncol = length)
  reads_table_protein <- matrix(nrow = repetitions, ncol = length)
  reads_dna_size <- matrix(nrow = repetitions, ncol = 1)
  reads_table_F[] <- 0
  reads_table_R[] <- 0
  reads_table_protein[] <- 0
  coordinates <- 1:length
  colnames(reads_table_F) <- coordinates
  colnames(reads_table_R) <- coordinates
  colnames(reads_table_protein) <- coordinates
  
  #########
  while (i < repetitions + 1) {
    # Prot. coordinates
    ## Check prot mean is within boundaries
    if (prot_mean < 0 | prot_mean > length) {
      stop('Prot_mean out of boundaries')
    }
    ## Set prot. coordinates
    prot_center <- as.integer(rnorm(1, mean = prot_mean, sd = prot_sd))
    prot_left <-  as.integer(prot_center - (prot_length/2))
    prot_right <- prot_left + prot_length
    if (prot_left < 1 | prot_right > length) {
      next
    }
    #########
    # determine left and right cuts
    ##cuttable_dna_size <- length - no_cut_l[1] + no_cut_l[2] - no_cut_r[2] + no_cut_r[1]
    number_of_cuts <- as.integer(length * number_cuts_per_kb/1000)
    cuts <- as.integer(runif(number_of_cuts, min = 1, max = length))
    
    ## Check if cut_l should be used or not
    if (cut_l > 0) {
      stop('cut_l should be a negative integer, representing the position of the cut relative to the leftmost protein boundary')
    }
    if (cut_l < 0)  {
      cuts <- c(cuts,cut_l + prot_left)
    }
    
    ## Check if cut_r should be used or not
    if (cut_r < 0) {
      stop('cut_r should be a positive integer, representing the position of the cut relative to the leftmost protein boundary')
    }
    if (cut_r > 0)  {
      cuts <- c(cuts,cut_r + prot_right)
    }
    
    # Eliminate cutswithin no cut zones
    cuts <-  cuts[cuts[] > prot_right | cuts[] < prot_left]
    if (no_cut_l[1] > no_cut_l[2]) {
      stop('no_cut_l should have the following format: \'c(-x,y)\', being x smaller than y,
             and both are the coordinate of the no cut zone relative to the protein position')
    }
    
    if (no_cut_l[1] < 0) {
      no_cut_l_distal <- prot_left + no_cut_l[1]
      no_cut_l_proximal <- prot_left + no_cut_l[2]
      cuts <-  cuts[cuts[] > no_cut_l_proximal | cuts[] < no_cut_l_distal]
    }
    
    if (no_cut_r[1] > no_cut_r[2]) {
      stop('no_cut_r should have the following format: \'c(x,y)\', being x smaller than y, 
               and both are the coordinate of the no cut zone relative to the protein position')
    }
    if (no_cut_r[2] > 0) {
      no_cut_r_distal <- prot_right + no_cut_r[2]
      no_cut_r_proximal <- prot_right + no_cut_r[1]
      cuts <-  cuts[cuts[] < no_cut_r_proximal | cuts[] > no_cut_r_distal]
    }
    if (length(cuts) == 0) {
      next
    }
    
    ## Choose smallest dna size
    left_cuts <- cuts[cuts[] < prot_left]
    right_cuts <- cuts[cuts[] > prot_right]
    if (length(left_cuts) == 0 | 
        length(right_cuts) == 0) {
      next
    }
    
    left_cut <- max(cuts[cuts[] < prot_left])
    right_cut <- min(cuts[cuts[] > prot_right])
  
    #Check it is bigger than bottom DNA size
    if (right_cut - left_cut < dna_bottom_size) {
      next
    }
    
    #########
    # Fill the matrix according to left/right cuts and intensity of the signal
    if (PE_full_length_read) {
      reads_table_F[i,left_cut:right_cut] <- intensity
      reads_table_R[i,left_cut:right_cut] <- intensity
      reads_table_protein[i,prot_left:prot_right] <- intensity
      reads_dna_size[[i]] <- right_cut - left_cut
      i = i + 1
    } else {
    reads_table_F[i,left_cut:(left_cut + read_length)] <- intensity
    reads_table_R[i,(right_cut - read_length):right_cut] <- intensity
    reads_table_protein[i,prot_left:prot_right] <- intensity
    reads_dna_size[i,1] <- right_cut - left_cut
    i = i + 1
    }
  }
  
  #Calculate average value per coordinate
  i <- 1
  average_table_F <- matrix(nrow = 2, ncol = length)
  average_table_R <- matrix(nrow = 2, ncol = length)
  average_table_protein <- matrix(nrow = 2, ncol = length)
  
  rownames(average_table_F) <- c("coordinates","Fw")
  rownames(average_table_R) <- c("coordinates","Rv")
  rownames(average_table_protein) <- c("coordinates","Protein")
  
  average_table_F[1,] <- coordinates
  average_table_R[1,] <- coordinates
  average_table_protein[1,] <- coordinates
  
  while (i < (length + 1)) {
    average_table_F[2,i] <- mean(reads_table_F[,i])
    average_table_R[2,i] <- mean(reads_table_R[,i])
    average_table_protein[2,i] <-  mean(reads_table_protein[,i])
    i = i + 1
  }
  
  
  # Transform matrixes to tibbles so it can be used in ggplot (I am not sure if this is actually necessary)
  average_table_F <- as_tibble(t(average_table_F))
  average_table_R <- as_tibble(t(average_table_R))
  average_table_protein <- as_tibble(t(average_table_protein))
  # Calculate DNA-size mean and sd
  dna_statistics <- tibble("mean_DNA_size" = mean(reads_dna_size[,1]), "DNA_size_sd" = sd(reads_dna_size[,1])) %>%
    mutate(mean_DNA_size = sprintf("%0.0f", mean_DNA_size), DNA_size_sd = sprintf("%0.0f", DNA_size_sd))
  
  # merge tibbles, calculate sum of Fw and Rv, 
  average_table <- left_join(average_table_F, average_table_R) %>% 
    mutate("Fw_plus_Rv" = Fw + Rv)
  # Re size protein average so that max value is the same as Fw_plus_Rv average
  max_fw_rv <- max(average_table$Fw_plus_Rv)
  max_prot <- max(average_table_protein$Protein)
  average_table_protein <- average_table_protein %>% 
                            mutate("Protein" = max_fw_rv * Protein / max_prot)
  
  # Merge tables and re-shape so that it can be plotted using gpglot2
  average_table <- average_table %>% 
    left_join(average_table_protein) %>% 
    pivot_longer(col = -coordinates, names_to = "Strand", values_to = "Average_signal")
  #Change the order of factors so when graphed in ggplot Fw_plus_Rv comes first (and does not cover the other 2)
  average_table$Strand <- factor(average_table$Strand, levels = c("Fw_plus_Rv","Fw","Rv","Protein"))
  
  #Get parameters in a table
  parameters <- tibble("prot_mean " = prot_mean,
                       "prot_sd " = prot_sd,
                       "prot_length" = prot_length, 
                       "no_cut_l" = paste(no_cut_l[1],no_cut_l[2], sep = ","), 
                       "no_cut_r" = paste(no_cut_r[1],no_cut_r[2], sep = ","), 
                       "cut_l" = cut_l, 
                       "cut_r" = cut_r,
                       "reps" = repetitions, 
                       "read_length" = read_length, 
                       "DNA_btm_size" = dna_bottom_size, 
                       "number_cuts_per_kb" = number_cuts_per_kb, 
                       "PE_FL_read" = PE_full_length_read,
                       "DNA_mean_size" = dna_statistics$mean_DNA_size,
                       "DNA_size_sd" = dna_statistics$DNA_size_sd)
  
  #Create a graph
  if (plot) {
    y_axis_max <- average_table %>% filter(Strand == "Fw_plus_Rv") %>% select("Average_signal") %>% max()
    
    table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max * 0.2), tb = list(parameters))
    # Agregate profile
    print(ggplot(data = average_table, aes(x = coordinates, y = Average_signal, color = Strand)) +
            geom_line() +
            geom_table(data = table1, aes(x, y, label = tb), size = 3.5))
    # DNA size plot
    table2 <- tibble(x = 1, y = 1 , tb = list(dna_statistics))
    reads_dna_size <- as_tibble_col(reads_dna_size[,1], column_name = "DNA_size" )
    print(ggplot(data = reads_dna_size) +
          geom_freqpoly(mapping = aes(x = DNA_size), binwidth = 3) +
          geom_table_npc(data = dna_statistics,
                         label = list(dna_statistics),
                         npcx = 1, npcy = 1,
                         size = 5))
  }
  # Return average and parameters tables
  return(list("average_table" = average_table, "parameters" = parameters, "dna_statistics" = dna_statistics, "reads_dna_size" = reads_dna_size))
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

  # Call "simulate_chip" function for every peak in the parameters list
  # Make list to store resuts
  results <- list()
  # For parameters that are only input once, repeat them as many times as peaks there are
  if (length(parameters$names) < parameters$number_peaks[[1]]) {
  for (i in length(parameters$names):(parameters$number_peaks[[1]] + length(parameters$names))) {
    parameters$names[[i]] <- str_c("peak_",i)  
    }
  }
  
  for (i in 1:17) { #17 is the number of parameters to input to the simulate_chip function
    if (length(parameters[[i]]) == 1) {
      while (length(parameters[[i]]) < parameters$number_peaks[[1]]) {
        parameters[[i]] <- c(parameters[[i]],parameters[[i]])
      }
      a <- parameters[[i]]
      parameters[[i]] <- a[1:parameters$number_peaks[[1]]]
     # parameters[[i]] <- c(parameters[[i]],parameters[[i]],parameters[[i]])
    }
    
  }
  
  # Call simulate_peak() for all the peaks
  for (w in (1:parameters$number_peaks[[1]])) {
    results[[ parameters$names[[w]] ]] <- simulate_chip(prot_mean = parameters$prot_mean[[w]],
                                               prot_sd = parameters$prot_sd[[w]],
                                               prot_length = parameters$prot_length[[w]],
                                               no_cut_l = parameters$no_cut_l[[w]],
                                               no_cut_r = parameters$no_cut_r[[w]],
                                               cut_l = parameters$cut_l[[w]],
                                               cut_r = parameters$cut_r[[w]],
                                               repetitions = parameters$repetitions[[w]],
                                               length = parameters$length[[w]],
                                               read_length = parameters$read_length[[w]],
                                               number_cuts_per_kb = parameters$number_cuts_per_kb[[w]],
                                               dna_bottom_size = parameters$dna_bottom_size[[w]],
                                               plot = parameters$plot[[w]],
                                               intensity = parameters$intensity[[w]],
                                               PE_full_length_read = parameters$PE_full_length_read[[w]])
    # Add the peak name info to each result table
    results[[c(w,1)]] <- results[[c(w,1)]] %>% mutate(peak_name = parameters$names[[w]])
  }
  # Create tibbles to combine peaks information
  composite_result <- tibble()
  composite_parameters <- tibble()
  
  for (n in (1:parameters$number_peaks[[1]]) ) {
    composite_result <- bind_rows(composite_result, results[[c(n,1)]])
    composite_parameters <- bind_rows(composite_parameters, results[[c(n,2)]])
    n <- n + 1
  }
  
  # Get average signal for composite peak...it is actually the sum of the averages of the two strands
  composite_peak <- composite_result %>% group_by(coordinates, Strand) %>% 
                    summarise("Average_signal" = sum(Average_signal), peak_name = "composite_peak") %>%
                    ungroup() #Is ungroup necessary?
  composite_result <- bind_rows(composite_result, composite_peak)
  composite_result$peak_name <- as.factor(composite_result$peak_name)
  # Print plots
  # All together (composite peaks plus individual peaks)
   y_axis_max <- composite_result %>% 
                filter(Strand == "Fw_plus_Rv") %>% select("Average_signal") %>% 
                max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb = list(composite_parameters))
  print(ggplot(data = composite_result, 
               aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
    geom_line() +
    geom_table(data = table1, aes(x, y, label = tb), size = 3.5))
  # Individual Peaks only (no composite peak)
  y_axis_max <- composite_result %>% filter(Strand == "Fw_plus_Rv", peak_name != "composite_peak") %>% 
                select("Average_signal") %>% max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb = list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name != "composite_peak"), 
               aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
          geom_line() +
          geom_table(data = table1, aes(x, y, label = tb), size = 3.5))
  # Composite peak only
  y_axis_max <- composite_result %>% 
    filter(Strand == "Fw_plus_Rv", peak_name == "composite_peak") %>% 
    select("Average_signal") %>% 
    max()
  table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max*0.2), tb = list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name == "composite_peak"), 
               aes(x = coordinates, y = Average_signal, color = Strand, linetype = peak_name)) +
          geom_line() +
          geom_table(data= table1, aes(x, y, label = tb), size = 3.5))
  
  return(composite_result)
}

#######################################
###### EXAMPLE OF SINGLE PEAK ########
######################################
# simulate_chip(repetitions = 10000,
#               prot_mean = 1500,
#               prot_sd = 200,
#               prot_length = 20,
#               read_length = 20,
#               dna_bottom_size = 75,
#               length = 3000,
#               number_cuts_per_kb = 10,
#               intensity = 7,
#               PE_full_length_read = FALSE,
#               cut_r = 0,
#               cut_l = -0,
#               no_cut_l = c(0,0))
################################################
###### EXAMPLE OF COMPOSITE PEAKS ##############
################################################
####### 2 peaks #############
# a <- list(prot_mean = c(990,1010),
#           prot_sd = c(60,60),
#           prot_length = c(30,30),
#           no_cut_l = list(c(0,0),c(0,0)),
#           no_cut_r = list(c(0,0),c(0,0)),
#           cut_l = c(0,-30),
#           cut_r = c(30,0),
#           repetitions = c(3000,3000),
#           length = c(2000,2000),
#           read_length = c(75,75),
#           number_cuts_per_kb = c(2,2),
#           dna_bottom_size = c(50,50),
#           plot = c(FALSE,FALSE),
#           number_peaks = c(2,2),
#           names = c("Peak1", "Peak2"),
#           intensity = c(1,1),
#           PE_full_length_read = c(FALSE,FALSE))
# x <- simulate_composite_peaks(a)

####### 3 peaks #############
# a <- list(prot_mean = c(2100,3000,3900),
#           prot_sd = c(250,100,250), 	          
#           prot_length = c(20,20,20),           
#           no_cut_l = list(c(0,0)),
#           no_cut_r = list(c(0,0)), 
#           cut_l = c(0,0,-1), 	          
#           cut_r = c(1,0,0),
#           repetitions = c(3000),
#           length = c(6000),
#           read_length = c(75),
#           number_cuts_per_kb = c(3.8,2.3,3.8),
#           dna_bottom_size = c(85),
#           plot = c(FALSE),
#           number_peaks = c(3),
#           names = c("Peak1", "Peak2","Peak3"),
#           intensity = c(7,2,7),
#           PE_full_length_read = c(FALSE))
# 
# x <- simulate_composite_peaks(a)

# combined <- bind_rows(filter(x, Strand == "Fw_plus_Rv" & peak_name == "composite_peak"), 
#                       filter(PE, Strand == "Fw_plus_Rv" & peak_name == "composite_peak"))
# 

# #
# # # ####### 4 peaks #############
# a <- list(prot_mean = c(2100,3000,3000,3900),
#           prot_sd = c(220,100,220,100),
#           prot_length = c(30,30,30),
#           no_cut_l = list(c(0,0),c(0,0), c(0,0)),
#           no_cut_r = list(c(0,0),c(0,0), c(0,0)),
#           cut_l = c(0,0,-1),
#           cut_r = c(1,0,0),
#           repetitions = c(2000,2000,2000),
#           length = c(6000,6000,6000,6000),
#           read_length = c(75,75,75),
#           number_cuts_per_kb = c(3.8,2.3,3.8),
#           dna_bottom_size = c(85,85,85),
#           plot = c(FALSE,FALSE,FALSE),
#           number_peaks = c(3,3,3),
#           names = c("Peak1", "Peak2","Peak3"),
#           intensity = c(7,2,7),
#           PE_full_length_read = c(FALSE,FALSE,FALSE))
# x <- simulate_composite_peaks(a)
