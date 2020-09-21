#######################################
###### EXAMPLE OF SINGLE PEAK ########
# ######################################
# simulate_chip(repetitions = 1000,
#               prot_mean = 1500,
#               prot_sd = 200,
#               prot_length = 20,
#               read_length = 20,
#               dna_bottom_size = 75,
#               length = 3000,
#               number_cuts_per_kb = 10,
#               intensity = 7,
#               cut_r = 0,
#               cut_l = -0,
#               no_cut_l = c(0,0))



################################################
###### EXAMPLE OF COMPOSITE PEAKS ##############
################################################
####### 2 peaks #############
a <- list(prot_mean = c(1000,2000),
          prot_sd = c(20),
          prot_length = c(10),
          no_cut_l = list(list(c(100,200),c(250,400))),
          no_cut_r = c(0,0),
          cut_l = c(0,0),
          cut_r = c(0,0),
          repetitions = c(10),
          length = c(3000),
          read_length = c(75),
          number_cuts_per_kb = c(4),
          dna_bottom_size = c(80),
          plot = c(FALSE),
          names = c("Peak1"),
          intensity = c(1),
          fixed_pos_cut = 0)
x <- simulate_composite_peaks(a)

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
# # # ####### 2 peaks #############
single_peak <- list(prot_mean = c(1000),
          prot_sd = c(50),
          prot_length = c(10),
          no_cut_l = list(c(0,0)),
          no_cut_r = list(c(0,0)),
          cut_l = c(-250),
          cut_r = c(0),
          repetitions = c(1000),
          length = c(2000),
          read_length = c(50),
          number_cuts_per_kb = c(6),
          dna_bottom_size = c(85),
          plot = c(FALSE),
          names = c("Peak1"),
          intensity = c(1),
          fixed_pos_cut = (0))
# x <- simulate_composite_peaks(single_peak)
# 
# ggplot(data = filter(x, peak_name == "composite_peak", Strand != "Protein", Strand != "Full_length")) + 
#   geom_line(mapping = aes(x = Coordinates, y = Average_signal, 
#                           color = Strand)) +
#   scale_color_manual(values = c("green", "blue"))


# #
# # # ####### pT28 #############
pT28 <- list(prot_mean = c(2100,2654,2854,3146,3346,3900),
          prot_sd = c(220,100,100,100,100,220),
          prot_length = c(20),
          no_cut_l = list(c(0,0)),
          no_cut_r = list(c(0,0)),
          cut_l = c(-0,0,0,0,0,-1),
          cut_r = c(1,0,0,0,0,0),
          repetitions = c(3000),
          length = c(6000),
          read_length = c(50),
          number_cuts_per_kb = c(2.4,4.5,4.5,4.5,4.5,2.4),
          dna_bottom_size = c(85),
          plot = c(FALSE),
          names = c("Peak1"),
          intensity = c(7,0.5,1,1,0.5,7),
          PE_full_length_read = c(FALSE,FALSE,FALSE))
# pT28 <- simulate_composite_peaks(pT28)
# pT28_b <- filter(pT28, Strand != "Protein", Strand != "Full_length") %>% 
#   mutate("Average_signal" = Average_signal / 1.65) %>% 
#   bind_rows(filter(pT28, Strand != "Protein", Strand == "Full_length"))

# ggplot(data = filter(pT28_b, peak_name == "composite_peak", Strand != "Protein")) + 
  # geom_line(mapping = aes(x = Coordinates, y = Average_signal, 
  #                         color = Strand)) 
  # scale_color_manual(values = c("green", "blue"))
  
# # # ####### T33 #############
T33 <- list(prot_mean = c(2654,2854,3146,3346),
             prot_sd = c(100,100,100,100),
             prot_length = c(20),
             no_cut_l = list(c(0,0)),
             no_cut_r = list(c(0,0)),
             cut_l = c(0,0,0,0),
             cut_r = c(0,0,0,0),
             repetitions = c(4000),
             length = c(6000),
             read_length = c(50),
             number_cuts_per_kb = c(5),
             dna_bottom_size = c(85),
             plot = c(FALSE),
             number_peaks = c(4),
             names = c("Peak1"),
             intensity = c(0.4,1,1,0.4),
             PE_full_length_read = c(FALSE),
             fixed_pos_cut = 0)

# simT33 <- simulate_composite_peaks(T33)

# ggplot(data = filter(T33, peak_name == "composite_peak", Strand != "Protein")) + 
  geom_line(mapping = aes(x = Coordinates, y = Average_signal, 
                          color = Strand)) 

# T28 zyg lib5 ------------------------------------------------------------
T28_zyg_rep3 <- list(prot_mean = c(2450,2654,2854,3146,3346,3550),
                     prot_sd = c(10),
                     prot_length = c(10),
                     no_cut_l = list(c(-73,0),c(-73,0),c(-73,0),c(-73,0),c(-73,0),c(-73,0)),
                     no_cut_r = list(c(0,73),c(0,73),c(0,73),c(0,73),c(0,73),c(0,73)),
                     cut_l = c(0),
                     cut_r = c(0),
                     repetitions = c(1000),
                     length = c(6000),
                     read_length = c(1),
                     number_cuts_per_kb = c(13),
                     dna_bottom_size = c(140),
                     plot = c(FALSE),
                     names = c("Peak1"),
                     intensity = c(2,6,8,8,6,2),
                     PE_full_length_read = c(FALSE),
                     fixed_pos_cut = 0)
# simT28_zyg_rep3 <- simulate_composite_peaks(T28_zyg_rep3)
# 
# simT28_zyg_rep3_smothed <- smooth(simT28_zyg_rep3, value = "Average_signal",
#        protein = "peak_name")
# ggplot(data = filter(simT28_zyg_rep3_smothed, peak_name == "composite_peak")) +
#   geom_line(mapping = aes(x = Coordinates, y = Average_signal,color = Strand))

# SMARCA5 -----------------------------------------------------------------
SMARCA5 <- list(prot_mean = c(1900,2411,2662,2861,3139,3338,3589,4100),
                     prot_sd = c(320,20,20,20,20,20,20,320),
                     prot_length = c(30),
                     no_cut_l = list(c(0,0)),
                     no_cut_r = list(c(0,0)),
                     cut_l = c(0,0,0,0,0,0,0,-1),
                     cut_r = c(1,0,0,0,0,0,0,0),
                     repetitions = c(4000),
                     length = c(6000),
                     read_length = c(1),
                     number_cuts_per_kb = c(2.3,10,10,10,10,10,10,2.3),
                     dna_bottom_size = c(80,140,140,140,140,140,140,80),
                     plot = c(FALSE),
                     names = c("Peak1"),
                     intensity = c(25,2.3,3,4.4,4.4,3,2.3,25),
                     PE_full_length_read = c(FALSE),
                     fixed_pos_cut = 0)
# simSMARCA5 <- simulate_composite_peaks(SMARCA5)
# r <- smooth(simSMARCA5, value = "Average_signal", protein = "peak_name")
# ggplot(data = filter(r, peak_name=="composite_peak")) +
#   geom_line(mapping = aes(x = Coordinates, y=Average_signal, color = Strand))
