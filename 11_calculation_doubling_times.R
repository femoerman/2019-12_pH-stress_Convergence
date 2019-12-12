# clear memory
rm(list=ls())

#Set working directory
#setwd("D://Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(nlme)
library(coefplot2)

#Load dataset
load("2_data/posteriorBH_means.RData")

#1) Calcultation of doubling times
{
  sumdata$doublingTimes <- log(2)/log(1 +exp(sumdata$log_r0mean))
}

#2) Plot doubling times for ancestral populations
{
  plotdata <- filter(sumdata, timing=="ancestral")
  plotdata$Genotype <- as.factor(ifelse(plotdata$strain=="B2086.2", "Genotype 1", ifelse(plotdata$strain=="CU427.4", "Genotype 2", ifelse(plotdata$strain == "CU428.2", "Genotype 3", "Genotype 4"))))
  plotdata$pH <- ifelse(plotdata$pHmedium==5, "4.5", "6.5")
  
  ggplot(plotdata, aes(y=doublingTimes, x=Genotype, colour=pH)) + geom_boxplot() + ylab("Doubling time (h)") + xlab("Genotype")+ 
    theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
           plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) ) + labs(colour = "pH of the assay medium")+ theme(legend.position="top")
  ggsave(filename = "Figresp2.png", path = "4_results/Figures", width = 255, height = 204, units = "mm", device = "png")
  
  ggplot(plotdata, aes(y = exp(log_r0mean)-exp(log_dmean), x = Genotype)) + geom_point()
}

#3) Plot doubling time for all populations
{
  plotdata <- filter(sumdata, timing!="nocg")
  ggplot(plotdata, aes(y=doublingTimes, x=strain, colour=as.factor(pHmedium))) + geom_boxplot(position="identity") + facet_wrap(~paste(timing, pHorigin))
}

ggplot(plotdata, aes(x=alfamean, y=doublingTimes, colour=as.factor(pHmedium), shape=strain))+ geom_point() + facet_wrap(~paste(timing, pHorigin))
ggplot(plotdata, aes(x=r0mean, y=doublingTimes, colour=as.factor(pHmedium), shape=strain))+ geom_point() + facet_wrap(~paste(timing, pHorigin))

#4) Visualize density over time for treatments of ancestors
{
  times <- seq(from=1, to=14*7, by=1)
  ode.model_BH = function(t,N,p){
    with(as.list(p),{
      dNdt = ((p[1] + p[2])/(1 + ((p[1]/(p[3]*p[2])) * N)) - p[2])*N
      return(list(dNdt))
    })
  }
  
  ancestors <- filter(sumdata, timing == "ancestral")
  datatimeseries <- data.frame(row.names=1:(32*49*24))
  dens=vector("double")
  datatimeseries$times <- rep(1:(49*24), 32)
  strain=vector("character")
  pHmedium <- vector("double")
  ident <- vector("character")
  for (i in 1:32){
    temp <- ancestors[i, ]
    dens <- c(dens, c( as.numeric(lsoda(y=temp$Kmean/2, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 72), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    for (i in 1:6){
    dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 72), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
    }
    
    strain <- c(strain, rep(temp$strain, (49*24)))
    pHmedium <- c(pHmedium, rep(temp$pHmedium, (49*24)))
    ident <- c(ident, rep(temp$ident, (49*24)))
  }
  datatimeseries$pHmedium <- pHmedium
  datatimeseries$strain <- strain
  datatimeseries$dens <- as.numeric(dens)
  datatimeseries$ident <- ident
  
  #Plot this data
  ggplot(datatimeseries, aes(x=times, y=dens,  colour=as.factor(pHmedium), shape=ident)) + geom_line() + facet_wrap(~strain)
  ggplot(filter(datatimeseries, times<(24*14)), aes(x=times, y=dens,  colour=as.factor(pHmedium), shape=ident)) + geom_line() + facet_wrap(~strain) + ggtitle("Ancestral")
}

#5) Visualize doubling times for evolved strains
{
  {
    times <- seq(from=1, to=14*7, by=1)
    ode.model_BH = function(t,N,p){
      with(as.list(p),{
        dNdt = ((p[1] + p[2])/(1 + ((p[1]/(p[3]*p[2])) * N)) - p[2])*N
        return(list(dNdt))
      })
    }
    
    evolved <- filter(sumdata, timing == "cg")
    evolved <- evolved[complete.cases(evolved), ]
    datatimeseries2 <- data.frame(row.names=1:(58*49*24))
    dens=vector("double")
    datatimeseries2$times <- rep(1:(49*24), 58)
    strain=vector("character")
    pHmedium <- vector("double")
    pHorigin <- vector("double")
    ident <- vector("character")
    for (i in 1:58){
      temp <- evolved[i, ]
      dens <- c(dens, c( as.numeric(lsoda(y=temp$Kmean/2, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
      dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
      dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 72), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
      for (i in 1:6){
        dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
        dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 48), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
        dens <- c(dens, c( as.numeric(lsoda(y=tail(dens, 1)/2*temp$Kmean, times=seq(1, 72), func=ode.model_BH, parms=c(temp$r0mean, temp$dmean, temp$Kmean))[,2]))/temp$Kmean)
      }
      
      strain <- c(strain, rep(temp$strain, (49*24)))
      pHmedium <- c(pHmedium, rep(temp$pHmedium, (49*24)))
      pHorigin <- c(pHorigin, rep(temp$pHorigin, (49*24)))
      ident <- c(ident, rep(temp$ident, (49*24)))
    }
    datatimeseries2$pHmedium <- pHmedium
    datatimeseries2$pHorigin <- pHorigin
    datatimeseries2$strain <- strain
    datatimeseries2$dens <- as.numeric(dens)
    datatimeseries2$ident <- ident
    
    #Plot this data
    ggplot(datatimeseries2, aes(x=times, y=dens,  colour=as.factor(pHmedium), shape=ident)) + geom_line() + facet_wrap(~strain)
    ggplot(filter(datatimeseries2, pHorigin==5, times<(24*14)), aes(x=times, y=dens,  colour=as.factor(pHmedium), shape=ident)) + geom_line() + facet_wrap(~strain) + ggtitle("low-evolved")
    ggplot(filter(datatimeseries2,pHorigin==6.5,  times<(24*14)), aes(x=times, y=dens,  colour=as.factor(pHmedium), shape=ident)) + geom_line() + facet_wrap(~strain) + ggtitle("high-evolved")
  }
}

