#######################################
###### EXAMPLE OF SINGLE PEAK ########
######################################
simulate_chip(repetitions = 1000,
              prot_mean = 1500,
              prot_sd = 200,
              prot_length = 20,
              read_length = 20,
              dna_bottom_size = 75,
              length = 3000,
              number_cuts_per_kb = 10,
              intensity = 7,
              cut_r = 0,
              cut_l = -0,
              no_cut_l = c(0,0))



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
a <- list(prot_mean = c(2654,2854,3146,3346),
          prot_sd = c(100),
          prot_length = c(30),
          no_cut_l = list(c(0,0)),
          no_cut_r = list(c(0,0)),
          cut_l = c(0),
          cut_r = c(0),
          repetitions = c(500),
          length = c(6000),
          read_length = c(50),
          number_cuts_per_kb = c(6),
          dna_bottom_size = c(85),
          plot = c(FALSE),
          number_peaks = c(4),
          names = c("Peak1"),
          intensity = c(4,7,7,4),
          PE_full_length_read = c(FALSE,FALSE,FALSE))
x <- simulate_composite_peaks(a)
