#1) clear memory
{
  rm(list=ls())
}

#2) import libraries to use
{
  library(ggfortify)
  library(tidyverse)
}

#3) set working directory
{
  setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")
}

#4) import data, add hours to first one and collect in one data.frame
{
   #4.1) import simplified and regular trait data
  {
    traits <- read.csv("2_data/04_strain_traits")
    simplTraits <- read.csv("2_data/04_strain_traits_simplified.csv")
    
    #load posterior data
    load("2_data/4_PCA_Analysis/posteriorBH_means.RData")
    
    #Prepare data for use
    dd <- filter(sumdata, timing=="ancestral") #filter ancesor data
    dd2 <- data.frame(row.names = as.character(unique(as.character(dd$ident2)))) #Create data frame with identifier as rowname
    dd2 <- mutate(dd2, strain = NA, r0low = NA, r0neut = NA, Klow = NA, Kneut = NA, dlow = NA, dneut = NA, alphalow = NA, alphaneut = NA)
    row.names(dd2) <- unique(as.character(dd$ident2))
    
    #Copy life history traits to new data frame
    for (i in rownames(dd2)){
      templow <- filter(dd, ident2==i & pHmedium==5)
      tempneut <- filter(dd, ident2==i & pHmedium==6.5)
      dd2[i, "strain"] <- templow$strain
      dd2[i, "r0low"] <- templow$log_r0mean
      dd2[i, "Klow"] <- templow$log_Kmean
      dd2[i, "dlow"] <- templow$log_dmean
      dd2[i, "alphalow"] <- templow$log_alfamean
      dd2[i, "r0neut"] <- tempneut$log_r0mean
      dd2[i, "Kneut"] <- tempneut$log_Kmean
      dd2[i, "dneut"] <- tempneut$log_dmean
      dd2[i, "alphaneut"] <- tempneut$log_alfamean
    }
  }
}
# set strain names
dd2$strain <- ifelse(dd2$strain =="B2086.2", "Genotype 1", ifelse(dd2$strain =="CU427.4", "Genotype 2", ifelse(dd2$strain=="CU428.2", "Genotype 3", "Genotype 4") ))

#5) Perform PCA
{
  #5.1) for traits including dispersal
  {
      df <- traits[2:18]
      autoplot(prcomp(df, scale=TRUE), data = traits, colour = 'strain', loadings=TRUE, loadings.label=TRUE)
  }
  
  #5.2) for simplified traits
  {
    df <- simplTraits[2:6]
    autoplot(prcomp(df, scale=TRUE), data = traits, colour = 'strain', loadings=TRUE, loadings.label=TRUE)
  }
  
}

xy.list <- as.list(as.data.frame(t(dd2)))
#Perform PCA for ancestor data
{
  df <- dd2[, c(2, 3, 5, 7)]
  autoplot(prcomp(df, scale=TRUE), data = dd2, colour = 'strain', loadings=TRUE, loadings.label=TRUE, size=2) + theme_light() +   
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
    axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5))
  
  df <- dd2[, c(2:9)]
  dd3 <- rename(dd2, Genotype=strain)
  figX <- autoplot(prcomp(df, scale=TRUE), data = dd3, colour = 'Genotype', loadings=TRUE, loadings.label=TRUE, size=2) + theme_light() +   
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5))
  ggsave(figX, filename = "4_results/Figures/FigS2.png", device="png", width = 255, height = 204, units = "mm")
  
  dd3 <- dd2 %>% select(-"r0low", -"Klow", -"alphalow", -"dlow") %>% gather("Variable", "Value", -strain)
  ggplot(dd3, aes(x = strain, y = exp(Value), colour = strain)) + geom_point() + facet_wrap(~Variable, scales = "free")
  
  dd3 <- dd2  %>% gather("Variable", "Value", -strain)
  dd3$var <- ifelse(dd3$Variable %in% c("r0low", "r0neut"), "r0", ifelse(dd3$Variable %in% c("Klow", "Kneut"), "K", ifelse(dd3$Variable %in% c("dlow", "dneut"), "d", "alpha")))
  dd3$pH <- ifelse(dd3$Variable %in% c("r0low", "alphalow", "Klow", "dlow"), "low", "neutral")
  
  
  dd3$Genotype <- as.factor(ifelse(dd3$strain=="B2086.2", "Genotype 1", ifelse(dd3$strain=="CU427.4", "Genotype 2", ifelse(dd3$strain == "CU428.2", "Genotype 3", "Genotype 4"))))
  
  Fig.resp <- ggplot(dd3, aes(x = Genotype, y = exp(Value), colour = Genotype)) + geom_point() + facet_wrap(~Variable, scales = "free") + facet_grid(var~pH, scales = "free_y") + 
    xlab("Genotype") + ylab("Life-history trait") + 
    theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
           plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) )  +
    theme(strip.text = element_text(colour = 'black', size=16))
  
    
  ggsave(filename = "Figresp.png", path = "4_results/Figures", width = 255, height = 204, units = "mm", device = "png")
}