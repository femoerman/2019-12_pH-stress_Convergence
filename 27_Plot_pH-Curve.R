#1) Clear memory
rm(list=ls())

#2) Set working directory
#setwd("F:/Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#3) Load data
library(readxl)
library(tidyverse)
pHcurve_medium_BIO309 <- read_excel("2_data/pHcurve_medium_BIO309.xlsx")

#Rename daya
dd <- pHcurve_medium_BIO309
dd$total <- dd$`total amount HCl ((μl)`

#4) Plot data
fig1 <- ggplot(dd, aes(x = total, y  =pH)) + geom_point(size=2) + geom_line(size=1) + xlab("HCl concentration (μl/100mL)") + ylab("pH") + 
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.text.y=element_text(size=16), axis.text.x=element_text(size=16),
         plot.title = element_text(size=20, hjust = 0.5)) + theme_light()
ggsave(plot = fig1,filename = "FigS1.png", path = "4_results/Figures", width = 175, height = 125, units = "mm", device = "png")
