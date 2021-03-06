library(tidyverse)
library(gridExtra)
library(ggpmisc)
# Added protein position matrix

simulate_chip <- function(prot_mean = 500, prot_sd = 30, prot_length = 15, 
                          no_cut = c(0,0), cut_l = 0, cut_r = 0, repetitions = 1000,
                          length = 1000, read_length = 75, number_cuts_per_kb = 4, 
                          dna_bottom_size = 85, plot = TRUE, intensity = 1, 
                          fixed_pos_cut = 0) { 

# parameters exlplained ---------------------------------------------------
  # prot_mean       : the Coordinate in which the protein center can be found (mean of normal distribution)
  # prot_length     : length of DNA covered by the protein, this is a "no-cut" zone.
  # no_cut        : Coordinates, relative to the ends of the protein, in which no cuts are generated.  
  #                 eg for composite peak: "list(list(c(100,200),c(250,400)))". The first list is for each of the                       peaks and the second for each of the no_cut zones
  # cut_l / cut_r   : Coordinates (relative to the left-most or right-most side of the protein) of a fixed cut to the left or the right of the protein respectively. 0 = No cut
  # repetitions     : number of reads that are going to be included in the final table (discarded reads are...discarded, hence, they are not included in the final table)
  # length          : length of the Coordinate system
  # read_length     : length of sequencing, determined by the sequencer (75, 150, etc)
  # dna_mean_size   : DNA mean size (normal distribution)
  # dna_size_sd     : DNA size standard deviation (normal distribution)
  # dna_bottom_size : minimal DNA length. Any random piece of DNA below this threshold will be discarded. If you are purifying your adapter-ligated chipped DNA with ampure beads using 1X buffer, the max size you keep is about 200-150 bp, minus the adapter (75 bp), you get DNA pieces of at least 75-125 bp long 
  # plot            : Do you want to get the plots automatically? You might want to change it to "FALSE" for composite peaks
  # intensity       : Intensity of the peaks, for doing composite peaks

# functions ---------------------------------------------------------------
  delete_no_cuts <- function(i, cuts){
    if (i[1] > i[2]) {
      stop(str_c(i,' -->no_cut should have the following format: \'c(x,y)\', being x smaller than y,
                 and both are the Coordinate of the no cut zone relative to the protein position'))
    }
    if (i[1] >= 0 ) { # Then both not-cut values are to the right
      no_cut_left <- prot_right + i[1]
      no_cut_right <- prot_right + i[2]
    }
    if (i[2] <= 0 ) { # Then both not-cut values are to the left
      no_cut_left <- prot_left + i[1]
      no_cut_right <- prot_left + i[2]
    } else { # Then 1st is on the left and second is on the right
      no_cut_left <- prot_left + i[1]
      no_cut_right <- prot_right + i[2]
    }
    
    CUTS <-  cuts[cuts[] < no_cut_left | cuts[] > no_cut_right]
    return(CUTS)
    }

# Simulation --------------------------------------------------------------
  
  # Set index to count succesful iterations
  i <- 1
  
  # Set matrixes to store results -----
  reads_table_F <- matrix(nrow = repetitions, ncol = length)
  reads_table_R <- matrix(nrow = repetitions, ncol = length)
  reads_table_full_length <- matrix(nrow = repetitions, ncol = length)
  reads_table_protein <- matrix(nrow = repetitions, ncol = length)
  reads_dna_size <- matrix(nrow = repetitions, ncol = 1)
  reads_table_F[] <- 0
  reads_table_R[] <- 0
  reads_table_full_length[] <- 0
  reads_table_protein[] <- 0
  Coordinates <- 1:length
  colnames(reads_table_F) <- Coordinates
  colnames(reads_table_R) <- Coordinates
  colnames(reads_table_full_length) <- Coordinates
  colnames(reads_table_protein) <- Coordinates
  
  while (i < repetitions + 1) {
    # Prot. Coordinates ----
    ## Check prot mean is within boundaries
    if (prot_mean < 0 | prot_mean > length) {
      stop('Prot_mean out of boundaries')
    }
    ## Set prot. Coordinates
    prot_center <- as.integer(rnorm(1, mean = prot_mean, sd = prot_sd))
    prot_left <-  as.integer(prot_center - (prot_length/2))
    prot_right <- prot_left + prot_length
    if (prot_left < 1 | prot_right > length) {
      next
    }
    
    # determine left and right cuts ----

        number_of_cuts <- as.integer(length * number_cuts_per_kb/1000)
    if (fixed_pos_cut > 0) {
    cuts <- as.integer(c(runif(number_of_cuts, min = 1, max = length), as.integer(fixed_pos_cut)))  
    } else {
      cuts <- as.integer(runif(number_of_cuts, min = 1, max = length))
    }
    
    if (fixed_pos_cut > length) {
      stop('fixed_pos_cut is outside of boundaries (too high)')
    }
    ## Check if cut_l should be used or not ----
    if (cut_l > 0) {
      stop('cut_l should be a negative integer, representing the position of the cut relative to the leftmost protein boundary')
    }
    if (cut_l < 0)  {
      cuts <- c(cuts,cut_l + prot_left)
    }
    
    ## Check if cut_r should be used or not ----
    if (cut_r < 0) {
      stop('cut_r should be a positive integer, representing the position of the cut relative to the leftmost protein boundary')
    }
    if (cut_r > 0)  {
      cuts <- c(cuts,cut_r + prot_right)
    }
    
    # Delete cuts in no-cut zone ----
    if (length(no_cut) > 1) {
      for (s in seq_along(no_cut)) {
      cuts <- delete_no_cuts(no_cut[[s]], cuts = cuts)
      }
    }
    else{
      cuts <- delete_no_cuts(no_cut[[1]], cuts = cuts)
    }
    if (length(cuts) == 0) {
      next
    }
    
    ## Choose smallest dna size -----
    left_cuts <- cuts[cuts[] < prot_left]
    right_cuts <- cuts[cuts[] > prot_right]
    if (length(left_cuts) == 0 | 
        length(right_cuts) == 0) {
      next
    }
    
    left_cut <- max(cuts[cuts[] < prot_left])
    right_cut <- min(cuts[cuts[] > prot_right])
    
    #Check it is bigger than bottom DNA size -----
    if (right_cut - left_cut < dna_bottom_size) {
      next
    }
    
    #########
    # Fill the matrix according to left/right cuts and intensity of the signal
    reads_table_full_length[i,left_cut:right_cut] <- intensity
    reads_table_F[i,left_cut:(left_cut + read_length)] <- intensity
    reads_table_R[i,(right_cut - read_length):right_cut] <- intensity
    reads_table_protein[i,prot_left:prot_right] <- intensity
    reads_dna_size[i,1] <- right_cut - left_cut
    i = i + 1
    
    }
  
  #Calculate average value per Coordinate
  i <- 1
  average_table_F <- matrix(nrow = 2, ncol = length)
  average_table_R <- matrix(nrow = 2, ncol = length)
  average_table_full_length <- matrix(nrow = 2, ncol = length)
  average_table_protein <- matrix(nrow = 2, ncol = length)
  
  rownames(average_table_F) <- c("Coordinates","Fw")
  rownames(average_table_R) <- c("Coordinates","Rv")
  rownames(average_table_full_length) <- c("Coordinates","Full_length")
  rownames(average_table_protein) <- c("Coordinates","Protein")
  
  average_table_F[1,] <- Coordinates
  average_table_R[1,] <- Coordinates
  average_table_full_length[1,] <- Coordinates
  average_table_protein[1,] <- Coordinates
  
  while (i < (length + 1)) {
    average_table_F[2,i] <- mean(reads_table_F[,i])
    average_table_R[2,i] <- mean(reads_table_R[,i])
    average_table_full_length[2,i] <- mean(reads_table_full_length[,i])
    average_table_protein[2,i] <-  mean(reads_table_protein[,i])
    i = i + 1
  }
  
  
  # Transform matrixes to tibbles so it can be used in ggplot (I am not sure if this is actually necessary)
  average_table_F <- as_tibble(t(average_table_F))
  average_table_R <- as_tibble(t(average_table_R))
  average_table_full_length <- as_tibble(t(average_table_full_length))
  average_table_protein <- as_tibble(t(average_table_protein))
  # Calculate DNA-size mean and sd
  dna_statistics <- tibble("mean_DNA_size" = mean(reads_dna_size[,1]), "DNA_size_sd" = sd(reads_dna_size[,1])) %>%
    mutate(mean_DNA_size = sprintf("%0.0f", mean_DNA_size), DNA_size_sd = sprintf("%0.0f", DNA_size_sd))
  
  # merge tibbles, calculate sum of Fw and Rv, 
  average_table <- left_join(average_table_F, average_table_R) %>% left_join(average_table_full_length)
  # Re size protein average so that max value is the same as Fw_plus_Rv average
  max_fw_rv <- max(average_table$Full_length)
  max_prot <- max(average_table_protein$Protein)
  average_table_protein <- average_table_protein %>% 
    mutate("Protein" = max_fw_rv * Protein / max_prot)
  
  # Merge tables and re-shape so that it can be plotted using gpglot2
  average_table <- average_table %>% 
    left_join(average_table_protein) %>% 
    pivot_longer(col = -Coordinates, names_to = "Strand", values_to = "Average_signal")
  #Change the order of factors so when graphed in ggplot Fw_plus_Rv comes first (and does not cover the other 2)
  average_table$Strand <- factor(average_table$Strand, levels = c("Full_length","Fw","Rv","Protein"))
  
  #Get parameters in a table
  parameters <- tibble("prot_mean " = prot_mean,
                       "prot_sd " = prot_sd,
                       "prot_length" = prot_length, 
                       "no_cut" = paste(no_cut[1],no_cut[2], sep = ","), 
                       "cut_l" = cut_l, 
                       "cut_r" = cut_r,
                       "reps" = repetitions, 
                       "read_length" = read_length, 
                       "DNA_btm_size" = dna_bottom_size, 
                       "number_cuts_per_kb" = number_cuts_per_kb,
                       "DNA_mean_size" = dna_statistics$mean_DNA_size,
                       "DNA_size_sd" = dna_statistics$DNA_size_sd,
                       "intensity" = intensity,
                       "fixed_pos_cut" = fixed_pos_cut)
  
  #Create a graph
  if (plot) {
    y_axis_max <- average_table %>% filter(Strand == "Full_length") %>% select("Average_signal") %>% max()
    
    table1 <- tibble(x = 0, y = y_axis_max + (y_axis_max * 0.2), tb = list(parameters))
    # Agregate profile
    print(ggplot(data = average_table, aes(x = Coordinates, y = Average_signal, color = Strand)) +
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
  # no_cut = c(c(0,0),c(0,0)), cut_l = c(0,0), cut_r = c(0,1),
  # repetitions = c(2000,2000), length = c(2000,2000), read_length = c(75,75), 
  # dna_top_size = c(250,250), dna_bottom_size = c(85,85), plot = c(FALSE,FALSE),
  # number_peaks = 2, names=c("Peak1", "Peak2"))
  
  # Call "simulate_chip" function for every peak in the parameters list
  # Make list to store resuts
  results <- list()
  number_peaks <- length(parameters$prot_mean)
  # For parameters that are only input once, repeat them as many times as peaks there are
  # For peaks names
  if (length(parameters$names) < number_peaks) {
    for (i in length(parameters$names):(number_peaks + length(parameters$names))) {
      parameters$names[[i]] <- str_c("peak_",i)  
    }
  }
  #For the rest of parameters
  for (i in 1:length(parameters)) { # number of parameters to input to the simulate_chip function
    if (length(parameters[[i]]) == 1) {
      while (length(parameters[[i]]) < number_peaks) {
        parameters[[i]] <- c(parameters[[i]],parameters[[i]])
      }
      a <- parameters[[i]]
      parameters[[i]] <- a[1:number_peaks]
      # parameters[[i]] <- c(parameters[[i]],parameters[[i]],parameters[[i]])
    }
    
  }
  
  # Call simulate_peak() for all the peaks
  for (w in (1:number_peaks)) {
    results[[ parameters$names[[w]] ]] <- simulate_chip(prot_mean = parameters$prot_mean[[w]],
                                                        prot_sd = parameters$prot_sd[[w]],
                                                        prot_length = parameters$prot_length[[w]],
                                                        no_cut = parameters$no_cut[[w]],
                                                        repetitions = parameters$repetitions[[w]],
                                                        length = parameters$length[[w]],
                                                        read_length = parameters$read_length[[w]],
                                                        number_cuts_per_kb = parameters$number_cuts_per_kb[[w]],
                                                        dna_bottom_size = parameters$dna_bottom_size[[w]],
                                                        plot = parameters$plot[[w]],
                                                        fixed_pos_cut = parameters$fixed_pos_cut[[w]],
                                                        intensity = parameters$intensity[[w]])
    # Add the peak name info to each result table
    results[[c(w,1)]] <- results[[c(w,1)]] %>% mutate(peak_name = parameters$names[[w]])
  }
  # Create tibbles to combine peaks information
  composite_result <- tibble()
  composite_parameters <- tibble()
  
  for (n in (1:number_peaks) ) {
    composite_result <- bind_rows(composite_result, results[[c(n,1)]])
    composite_parameters <- bind_rows(composite_parameters, results[[c(n,2)]])
    n <- n + 1
  }
  
  # Get average signal for composite peak..
  composite_peak <- composite_result %>% group_by(Coordinates, Strand) %>% 
    summarise("Average_signal" = sum(Average_signal), peak_name = "composite_peak") %>%
    ungroup() #Is ungroup necessary?
  composite_result <- bind_rows(composite_result, composite_peak)
  composite_result$peak_name <- as.factor(composite_result$peak_name)
  max_Fw_Rv <- composite_result %>% filter(Strand == "Fw" | Strand == "Rv") %>% select(Average_signal) %>% max()
  max_Prot <- composite_result %>% filter(Strand == "Protein") %>% select(Average_signal) %>% max()
  max_FL <- composite_result %>% filter(Strand == "Full_length") %>% select(Average_signal) %>% max()
  norm <- list()
  for (strd in c("Fw", "Rv")) {
    norm[[as.character(strd)]] <- composite_result %>% 
      filter(Strand == strd) %>% 
      mutate("Average_signal" = Average_signal / max_Fw_Rv)
  }
  norm[["Protein"]] <- composite_result %>% 
      filter(Strand == "Protein") %>% 
      mutate("Average_signal" = 2 * Average_signal / max_Prot)
  norm[["Full_length"]] <- composite_result %>% 
    filter(Strand == "Full_length") %>% 
    mutate("Average_signal" = 2 * Average_signal / max_FL)
  composite_result <- tibble()
  for (i in seq_along(norm)) {
    composite_result <- bind_rows(composite_result, norm[[i]])
  }
  
  # Print plots
  # All together (composite peaks plus individual peaks)
  y_axis_max <- composite_result %>% 
    filter(Strand == "Full_length") %>% select("Average_signal") %>% 
    max()
  table1 <- tibble(x = 0, y = y_axis_max * (1 + .07 * number_peaks), tb = list(composite_parameters))
  print(ggplot(data = composite_result, 
               aes(x = Coordinates, y = Average_signal, color = Strand, linetype = peak_name, alpha = Strand)) +
          geom_line() +
          scale_alpha_manual(values = c(1,1,1,0.3)) +
          geom_table(data = table1, aes(x, y, label = tb), size = 3) +
          geom_vline(xintercept = parameters$prot_mean, alpha = 0.25, linetype = "dashed"))
  # Individual Peaks only (no composite peak)
  y_axis_max <- composite_result %>% filter(peak_name != "composite_peak") %>% 
    select("Average_signal") %>% max()
  table1 <- tibble(x = 0, y = y_axis_max * (1 + .07 * number_peaks), tb = list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name != "composite_peak"), 
               aes(x = Coordinates, y = Average_signal, color = Strand, linetype = peak_name,alpha = Strand)) +
          geom_line() +
          scale_alpha_manual(values = c(1,1,1,0.3)) +
          geom_table(data = table1, aes(x, y, label = tb), size = 3) +
          geom_vline(xintercept = parameters$prot_mean, alpha = 0.25, linetype = "dashed"))
  # Composite peak only
  y_axis_max <- composite_result %>% 
    filter(Strand == "Full_length", peak_name == "composite_peak") %>% 
    select("Average_signal") %>% 
    max()
  table1 <- tibble(x = 0, y = y_axis_max * (1 + .07 * number_peaks), tb = list(composite_parameters))
  print(ggplot(data = filter(composite_result, peak_name == "composite_peak"), 
               aes(x = Coordinates, y = Average_signal, color = Strand, linetype = peak_name, alpha = Strand)) +
          geom_line() +
          scale_alpha_manual(values = c(1,1,1,0.3)) +
          geom_table(data = table1, aes(x, y, label = tb), size = 3) +
          geom_vline(xintercept = parameters$prot_mean, alpha = 0.25, linetype = "dashed"))
  
  return(composite_result)
}