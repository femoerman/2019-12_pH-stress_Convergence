#1) clear memory
{
  rm(list=ls())
}

#2) Set working directory and load data Import RData files
{
  setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")
  #setwd("F:/Documenten/PhD/05_student_project_2017")
  load("2_data/DensitiesDuringEvolution.RData")
  dd <- evolve
  
  #Define if replicate was extinct or survived during experiment
  dd$survive <- ifelse(dd$file %in% c("sample_00015", "sample_00005", "sample_00001"), " -extinct", " -surviving")
  dd$colour <- paste(dd$pH, dd$survive)
  
  #Define strainnames
  dd$strainname <- ifelse(dd$strain=="B2086.2", "Strain 1", ifelse(dd$strain=="CU427.4", "Strain 2", ifelse(dd$strain=="CU428.2", "Strain 3", "Strain 4")))
  
}

#3) Load packages
{
  library(tidyverse)
  library(rethinking)
  library(glmer2stan)
  }

#4) Plot figures of evolution dynamics
{
  #4.1) Over all replicates per treatment
  
     library(scales)
   fig1 <- ggplot(filter(dd, strainname=="Strain 1"),aes(x=DaysSinceStart,y=indiv_per_ml, colour=colour))+ geom_line( aes(group=unique))+ geom_point()+
     ylab("Population density, N (indiv/mL)") + scale_y_continuous(trans=pseudo_log_trans(base=10, sigma=10), labels=scientific, breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6), limits = c(0, 1.5e6)) +
     theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                             axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                             axis.title=element_text(size=16), strip.text = element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))+
     theme(strip.background =element_rect(fill="white"))+
     theme(strip.text = element_text(colour = 'black')) +
     scale_color_manual(values=c("#2d33f9", "#000000", "#e41a1c"), breaks=c("high  -surviving", "low  -surviving", "low  -extinct"), name="Evolution treatment", labels=c("NpH", "LpH (surviving)", "LpH (extinct)"))+ xlim(0, 45)
  
   fig2 <- ggplot(filter(dd, strainname=="Strain 2"),aes(x=DaysSinceStart,y=indiv_per_ml, colour=colour))+ geom_line( aes(group=unique))+ geom_point() +
     scale_y_continuous(trans=pseudo_log_trans(base=10, sigma=10), labels=scientific, breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6), limits = c(0, 1.5e6)) +
     theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                             axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                             axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                             axis.title=element_text(size=16), strip.text = element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))+
     theme(strip.background =element_rect(fill="white"))+
     theme(strip.text = element_text(colour = 'black')) +
     scale_color_manual(values=c("#2d33f9", "#000000", "#e41a1c"), breaks=c("high  -surviving", "low  -surviving", "low  -extinct"), name="Evolution treatment", labels=c("NpH", "LpH (surviving)", "LpH (extinct)"))+ xlim(0, 45)
   
   fig3 <- ggplot(filter(dd, strainname=="Strain 3"),aes(x=DaysSinceStart,y=indiv_per_ml, colour=colour))+ geom_line( aes(group=unique))+ geom_point() +
     xlab("Time (d)") +  ylab("Population density, N (indiv/mL)") + scale_y_continuous(trans=pseudo_log_trans(base=10, sigma=10), labels=scientific, breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6), limits = c(0, 1.5e6)) +
     theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                             axis.title=element_text(size=16), strip.text = element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))+
     theme(strip.background =element_rect(fill="white"))+
     theme(strip.text = element_text(colour = 'black')) +
     scale_color_manual(values=c("#2d33f9", "#000000", "#e41a1c"), breaks=c("high  -surviving", "low  -surviving", "low  -extinct"), name="Evolution treatment", labels=c("NpH", "LpH (surviving)", "LpH (extinct)"))+ xlim(0, 45)
   

   fig4 <- ggplot(filter(dd, strainname=="Strain 4"),aes(x=DaysSinceStart,y=indiv_per_ml, colour=colour))+ geom_line( aes(group=unique))+ geom_point() +
     ylab("") + xlab("Time (d)")+ scale_y_continuous(trans=pseudo_log_trans(base=10, sigma=10), labels=scientific, breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6), limits = c(0, 1.5e6)) +
     theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                             axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                             axis.title=element_text(size=16), strip.text = element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))+
     theme(strip.background =element_rect(fill="white"))+
     theme(strip.text = element_text(colour = 'black')) +
     scale_color_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("high  -surviving", "low  -surviving"), name="Evolution treatment", labels=c("NpH", "LpH (surviving)"))+ xlim(0, 45)
   
}

#Combine figures
library(ggpubr)
ggarrange(fig1, fig2, fig3, fig4,  nrow=2, ncol=2, align = "hv", common.legend = T, legend =  "bottom",
                  labels =list("A", "B", "C", "D"), hjust = -3)

ggsave(filename = "Fig1.png", path = "4_results/Figures", width = 300, height = 200, units = "mm", device = "png")

