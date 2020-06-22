library(tidyverse)


simulate_chip <- function(prot_l=990, prot_r=1010, prot_length=40, no_cut_l=NULL, no_cut_r=NULL,
                          cut_l=NULL, cut_r=NULL, repetitions=10, length=2000, read_length=75,
                          dna_top_size=500) { 
  # prot_l and prot_r form the range of coordinates in which the protein center can be found
  # prot_length is the length of DNA that the protein covers, this is a "no-cut"zone
  # no_cut_l are the coordinates, relative to the left side of the protein, in which no cuts are generated  (input as "c(-10,0), beeing 0 the left of the protein)
  # no_cut_l are the coordinates, relative to the left side of the protein, in which no cuts are generated  (input as "c(-10,0), beeing 0 the left of the protein)
  # no_cut_r is the same as above but coodinates are positive, with 0 being the rigt side of the protein
  # cut_l and cut_r are the coordinates (relative to the left-most or right most side of the protein) of a fixed cut on the left or right of the protein
  # repetitions is the number of reads that are going to be included in the final table, default =1000
  # length of the coordinate system
  # read length is the length of the reads, determined by the sequencer (75, 150, etc) 
  
  i <- 1
  #reads_table <-  as.tibble(matrix(ncol = length, nrow = repetitions))
  reads_table_F <- matrix(nrow = repetitions, ncol = length)
  reads_table_R <- matrix(nrow = repetitions, ncol = length)
  reads_table_F[] <- 0
  reads_table_R[] <- 0
  coordinates <- 1:length
  colnames(reads_table_F) <- coordinates
  colnames(reads_table_R) <- coordinates
  
  while (i<repetitions) {
    #determine prot coordinates
    prot_left <- sample(prot_l:prot_r,1) #sample(a:b,N) makes N random integers between a and b 
    prot_right <- prot_left + prot_length
    
    #determine left and right cuts
    
    ####################################################
    ####################################################
    ####################################################
    left_cut <- sample(1:(prot_left-1),1)
     
    if (!is.null(cut_l)) {
      if (cut_l > 0) {
        stop('cut_l should be a negative integer, representing the position of the cut relative to the leftmost protein boundary')
      }
       if ((prot_left + cut_l) > left_cut) {
        left_cut <- prot_left + cut_l
       }
    }
    
    right_cut <- sample((prot_left+1):length,1)
    if (!is.null(cut_r)) {
      if (cut_r < 0) {
        stop('cut_r should be a positive integer, representing the position of the cut relative to the leftmost protein boundary')
      }
      if ((prot_right + cut_r) < right_cut) {
        right_cut <- prot_right + cut_r
      }
    }
    ####################################################
    ####################################################
    ####################################################
    #Determine if cuts are in no cut zone
    cut_restriction_left <- TRUE
    if (!is.null(no_cut_l)) {
      #Check that the first coordinate is smaller than the second
      if (no_cut_l[1] > no_cut_l[2]) {
        stop('no_cut_l should have the following format: \'c(-x,y)\', being x smaller than y, and both are the coordinate of the no cut zone relative to the protein position')
      }
      no_cut_l_distal <- prot_left + no_cut_l[1]
      no_cut_l_proximal <- prot_left + no_cut_l[2]
      if ((no_cut_l_distal < left_cut) & (left_cut < no_cut_l_proximal)) {
      cut_restriction_left <- FALSE
      }
    }
    
    cut_restriction_right <- TRUE  
    if (!is.null(no_cut_r)) {
      if (no_cut_r[1] > no_cut_r[2]) {
        stop('no_cut_r should have the following format: \'c(x,y)\', being x smaller than y, and both are the coordinate of the no cut zone relative to the protein position')
      }
      no_cut_r_distal <- prot_right + no_cut_r[1]
      no_cut_r_proximal <- prot_right + no_cut_r[2]
      if (no_cut_r_distal < right_cut & right_cut < no_cut_r_proximal) {
      cut_restriction_right <- FALSE
      }
    }
    
    # Keep read if DNA size <top_DNA_size and if left and right cut are not within no-cut zones
    if ( (right_cut - left_cut) < dna_top_size & cut_restriction_left & cut_restriction_right) {
      reads_table_F[i,left_cut:(left_cut+read_length)] <- 1
      reads_table_R[i,(right_cut-read_length):right_cut] <- 1
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
  average_table_F <- as.tibble(t(average_table_F))
  average_table_R <- as.tibble(t(average_table_R))
  
  # merge both tibbles, calculate sum of Fw and Rv, re-shape tibble so that it is useful for ggplot
  average_table <- left_join(average_table_F, average_table_R) %>% mutate("Fw_plus_Rv" = Fw + Rv) %>%
                  pivot_longer(col= -coordinates, names_to = "Strand", values_to = "Average_signal")
  #Change the order of factors so when graphed in ggplot Fw_plus_Rv comes first (and does not cover the other 2)
  average_table$Strand <- factor(average_table$Strand, levels =c("Fw_plus_Rv","Fw","Rv"))

  #Get parameters in a table
  parameters <- tibble("prot_l"= prot_l, "prot_r"= prot_r, "prot_length"= prot_length, 
                       "no_cut_l" = no_cut_l, "no_cut_r"= no_cut_r, "cut_l"= cut_l,"cut_r"=cut_r,
                       "repetitions"= repetitions, "length"= length, "read_length"= read_length, "dna_top_size"= dna_top_size)
  
  parameters[[parameters == "NULL"]] <- "NULL"
  
  
return(list("reads_matrix_F" = reads_table_F, "reads_matrix_R" = reads_table_R, "average_table" =average_table, "parameters" = parameters))
}


#rm(result)
result <- simulate_chip(repetitions=2000, prot_l = 950, prot_r =1050, prot_length = 45, 
                        read_length = 75, dna_top_size=300, cut_r = 10)


ggplot(data = result$average_table, aes(x = coordinates, y = Average_signal, color = Strand )) +
  geom_line()

# aggregate_profile <- function(table){
#   average_table <- matrix(nrow = 1, ncol = length)
#   while (i<(length+1)) {
#   average_table[1,i] <- mean(reads_table[,i])
#     }    
# }


