#RSCIPT FOR PLATEREADER DATA ANALYSIS - FERROZINE ASSAY
# Will need to base script off of the way the Excel file is formated
# when it comes off of the plate reader

# 1. Write code for standard curve
# 2. Reference Excelsheet from Plate reader
# 3. Figure out how to assign data to specific dataset >> This may need adjusted each time
# 4. Calculation Fe3+ (Total detectable iron - Fe2+)
# 5. Convert to PPM based on Standard Curve
# 6. Average Data
# 7. Calculation standard deviation (SD)
# 8. Create bar graph for each time point PPM with SD
# 9. Create line graph for change over time PPM with SD
# 10. Calc t-test, 1-way ANOVA, 2-way ANOVA PPM
# 11. Calc Fe2+ / Total detectable
# 12. Calc Fe3+ / Total detectable
# 13. Create bar graph for each time point Fe2+/ Total detectable && Fe3+ / total detectable with SD
# 14. Create line graph for change over time for Fe2+ / Total detectable iron with SD
# 15. Create line graph for change over time for Fe3+ / Total detectable iron  with SD
# 16. Cal t-test, 1-way ANOVA, 2-way ANOVA for % Fe2+/Tot & Fe3+/Tot


# ?? When should I normalize my data? Should I normalize each day?
# Or should I only normalize day 0?


setwd("~/Documents/Duquesne/Trun Lab/Bioinformatics/R/Rscripts")
#same location as data file (excel file)


# library(ggpubr)
# install.packages("Hmisc")
# install.packages("broom")
# library(broom) 
# library(lattice)
# library(Hmisc)
# install.packages("ggpmisc")
# library(ggpmisc)
# 
# 
# library(dplyr)
# 
# library(growthcurver)
# library(purrr)
# library(matrixStats)
# library(reshape2)
# 
# library(ggpubr)
# library(broom) 
# library(lattice)
# library(Hmisc)
# library(ggpmisc)

library(readr) 
library(ggplot2)

EXP137 <- read_csv("AV-137-Ferrozine-RScript.csv", col_names = TRUE) #load in datafile. Skip 1st 17 rows. Column names listed
#Latin1 was added, as I was getting an error code later in the doc



#create standard curve to calculate y=mx+b so we can measure the OD readings at PPM

graph <- ggplot(EXP137, aes(x = EXP137$STD_Value, y = EXP137$Absorbancy_562nm)) +  
  geom_point()+

  labs(x = "Fe(II) conc. PPM", y = "Absorbancy Fe(II) (OD600)") +
  
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.x = element_line(colour = "grey85", linetype = "dashed"),
  ) +
  ggtitle("Fe(II) Conc. (PPM) vs. Absorbancy Fe(II)") +
  
  theme(plot.title = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = .5)) +
  theme(axis.text = element_text(face = "plain", size = 12)) + 
  theme(axis.title = element_text(size = 15)) + 
  geom_smooth(method = "lm")

print(graph)  


# y = mx+b

# PPM = (OD-b)/m

#Would it be easier to pull the final numbers from excel into a more organized spreadsheet
#then use that spreadsheet to make these graphs? 

##load the data into R
EXP137PPM <- read_csv("AV-137-Ferrozine-PPM-RScript.csv", col_names = TRUE) #load in datafile. Skip 1st 17 rows. Column names listed
#Latin1 was added, as I was getting an error code later in the doc

##take the mean for the replicates
Fe2_G33_Day0_Aug_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Aug_AV1-W_G#33`)
#beacuse I am missing values in this column, it will not give me an average
#will need to remove the "NA" and leave in blank
#this still doesn't resolve the issue when it's left blank
Fe2_G33_Day0_Oct_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Oct_AV1-W_G#33`)
Fe2_G33_Day0_1xSx_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_1xSx_AV1-W_G#33`)

Fe2_G31_Day0_Aug_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Aug_AV1-W_G#31`)
#this does give me the average for the PPM
Fe2_G31_Day0_Oct_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Oct_AV1-W_G#31`)
Fe2_G31_Day0_1xSx_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_1xSx_AV1-W_G#31`)

Fe2_Sterile_Day0_Aug_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Aug_Sterile`)
Fe2_Sterile_Day0_Oct_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Oct_Sterile`)
Fe2_Sterile_Day0_1xSx_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_1xSx_Sterile`)
#take the st dev for the replicates

Fe2_MMC_Day0_Aug_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Aug_MMC`)
Fe2_MMC_Day0_Oct_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_S.B1_Oct_MMC`)
Fe2_MMC_Day0_1xSx_Avg <- mean(EXP137PPM$`Fe2+_PPM_DAY0_1xSx_MMC`)

## Bar graph for single time point 
# so I need to take the averages for each condition and make those they y values
# and need to make the x values the type of condition 


# gc_R2A <- gc %>%
#   rowwise() %>%
#   mutate(
#     m.R2A = mean(c(A8,B8,C8,D8,E8,F8,G8,H8)),
#     sd.R2A = sd(c(A8,B8,C8,D8,E8,F8,G8,H8))
#   )

graphPPM <- ggplot(EXP137PPM, aes(y = Fe2_G31_Day0_1xSx_Avg)) +  
  geom_bar()

print(graphPPM)

#this will graph just 1 condition mean, we want to graph all conditions on the
# same graph as individual bars
  



labs(x = "Fe(II) conc. PPM", y = "Absorbancy Fe(II) (OD600)") +
  
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.x = element_line(colour = "grey85", linetype = "dashed"),
  ) +
  ggtitle("Fe(II) Conc. (PPM) vs. Absorbancy Fe(II)") +
  
  theme(plot.title = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = .5)) +
  theme(axis.text = element_text(face = "plain", size = 12)) + 
  theme(axis.title = element_text(size = 15)) + 
  geom_smooth(method = "lm")

print(graph)  


# Tex for subscripts and superscripts
# $u_{i}$ will give you a subscript
# $Fe^{2+} will superscript the 2+ together
# $\mathrm{Fe^{2+}}$ will remove the italices that is automatically formated for Tex 

