# clear memory
rm(list=ls())

#Set working directory
#setwd("D://Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(rethinking)
library(gridExtra)

#Load data and prepare for usage (make variables for pH numerical, filter out data from before common garden)
load("2_data/posteriorBH_means.RData")

sumdata$strainval <- sumdata$strain
sumdata <- mutate(sumdata, strain=ifelse(sumdata$strain=="CU427.4", 1, ifelse(sumdata$strain=="CU428.2", 2, ifelse(sumdata$strain=="SB3539", 3, 4))))
sumdata$evolved <- ifelse(sumdata$timing=="cg", 1, 0)
sumdata$pHorigin <- ifelse(sumdata$pHorigin==6.5, 6.5, 4.5)
sumdata$highevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==6.5, 1, 0)
sumdata$lowevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==4.5, 1, 0)
sumdata$pHmedium <- ifelse(sumdata$pHmedium==6.5, 1, 0)
data.all <- sumdata %>% filter(timing != "nocg")
colorsch <- unique(paste(data.all$timing, data.all$pHorigin))

#Order dataset
data.all <- arrange(data.all, pHmedium, timing, pHorigin)

#2) Look at CU427.4 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #Get data for CU427.4
  dataCU427.4 <- filter(data.all, strain=="CU427.4")
  
  #Compare models
  models <- compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                    WAIC=FALSE, func=DIC)
  models
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16)
  best_model <- mr0.10
  precis(best_model, depth = 2, prob = 0.95)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e3)
  
  #1.3) Calculate the posteriors with 95% probability interval
  dataCU427.4 <- filter(data.all, strainval=="CU427.4")
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  dataCU427.4 <- filter(dataCU427.4, timing != "nocg")
  
  #Create weighted posteriors based on DIC criterion
  post.pred <- ensemble(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                       func=DIC, data=data, n=16e3)
  
  mu.pred.mean <- apply( post.pred$link , 2 , mean )
  mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred1 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred1 <- mutate(posteriorpred1, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  posteriorpred1 <- filter(posteriorpred1, evolved!= 1)
  posteriorpred1$colorsch <- paste(posteriorpred1$timingval, posteriorpred1$pHoriginval)
  
  #1.5) Plot reaction norms
  r0CU427.4 <- ggplot(dataCU427.4, aes(y=exp(log_r0mean), x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
    geom_line(data=posteriorpred1, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=colorsch), size=2) + 
    geom_ribbon(data=posteriorpred1, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ylim(0, 0.5) + 
    theme_light() + 
    theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
}

#3) Look at B2086.2 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/B2086.2/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.B2086.2 <- filter(data.all, strainval=="B2086.2")
  
  #Compare models
  models <- compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                    WAIC=FALSE, func=DIC)
  models
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16)
  best_model <- mr0.07
  precis(best_model, depth = 2, prob = 0.95)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e3)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data.B2086.2 <- filter(data.all, strainval=="B2086.2")
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                     lowevolved=unique(data.all$lowevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
  data.B2086.2 <- filter(data.B2086.2, timing != "nocg")
  
  #Create weighted posteriors based on DIC criterion
  post.pred <- ensemble(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                        func=DIC, data=data, n=16e3)
  
  mu.pred.mean <- apply( post.pred$link , 2 , mean )
  mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred4 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred4 <- mutate(posteriorpred4, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  r0B2086.2 <- ggplot(data.B2086.2, aes(y=exp(log_r0mean),  x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
    geom_line(data=posteriorpred4, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred4, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ylab(expression("r"[0]*" (1/h)")) + ylim(0, 0.5)     + 
    theme_light() + 
    theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
          axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
}

#4) Look at SB3539 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/SB3539/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.SB3539 <- filter(data.all, strainval=="SB3539")
  
  #Compare models
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16)
  models <- compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- mr0.07
  precis(best_model, depth = 2, prob = 0.95)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e3)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                     lowevolved=unique(data.all$lowevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
  data.SB3539 <- filter(data.SB3539, timing != "nocg")
  
  #Create weighted posteriors based on DIC criterion
  post.pred <- ensemble(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                        func=DIC, data=data, n=16e3)
  
  mu.pred.mean <- apply( post.pred$link , 2 , mean )
  mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred2 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred2 <- mutate(posteriorpred2, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  r0SB3539 <- ggplot(data.SB3539, aes(y=exp(log_r0mean),  x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
    geom_line(data=posteriorpred2, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred2, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ylim(0, 0.5)    +  theme_light() + 
    theme(axis.text=element_text(size=12), legend.text=element_text(size=16), legend.title=element_text(size=20),strip.text.x=element_text(20),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
}

#5) Look at CU428.2 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU428.2/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.CU428.2 <- filter(data.all, strainval=="CU428.2")
  
  #Compare models
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16)
  models <- compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- mr0.07
  precis(best_model, depth = 2, prob = 0.95)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e3)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                     lowevolved=unique(data.all$lowevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
  data.CU428.2 <- filter(data.CU428.2, timing != "nocg")
  
  #Create weighted posteriors based on DIC criterion
  post.pred <- ensemble(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, mr0.11, mr0.12, mr0.13, mr0.14, mr0.15, mr0.16,
                        func=DIC, data=data, n=16e3)
  
  mu.pred.mean <- apply( post.pred$link , 2 , mean )
  mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred3 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred3 <- mutate(posteriorpred3, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  r0CU428.2 <- ggplot(data.CU428.2, aes(y=exp(log_r0mean),  x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) +  
    geom_line(data=posteriorpred3, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred3, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ylim(0, 0.5)    + 
   theme_light() + 
    theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
}


#Do now all for K
{
  #2) Look at CU427.4 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/CU427.4/", full.names=T))
    lapply(file_names[2:17],load,.GlobalEnv)
    
    #Get data for CU427.4
    dataCU427.4 <- filter(data.all, strainval=="CU427.4")
    
    #Compare models
    coeftab(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16)
    models <- compare(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- mK.02
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    dataCU427.4 <- filter(data.all, strainval=="CU427.4")
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    dataCU427.4 <- filter(dataCU427.4, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred8 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred8 <- mutate(posteriorpred8, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    posteriorpred8 <- filter(posteriorpred8, !(lowevolved==0 & evolved == 1))
    posteriorpred8 <- filter(posteriorpred8, !(highevolved==0 & evolved == 1))
    
    #1.5) Plot reaction norms
    KCU427.4 <- ggplot(dataCU427.4, aes(y=exp(log_Kmean)/1e5, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred8, inherit.aes = F, aes(x=pHmediumval, y=mean/1e5, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred8, inherit.aes = F, aes(ymax=PI97.5/1e5, ymin=PI2.5/1e5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
      ylim(2, 9) + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #3) Look at B2086.2 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/B2086.2/", full.names=T))
    lapply(file_names[2:17],load,.GlobalEnv)
    
    #get data
    data.B2086.2 <- filter(data.all, strainval=="B2086.2")
    
    #Compare models
    coeftab(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16)
    models <- compare(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- mK.09
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    precis(best_model, digits=4, depth = 2, prob = 0.95)
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data.B2086.2 <- filter(data.all, strainval=="B2086.2")
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data.B2086.2 <- filter(data.B2086.2, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred7 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred7 <- mutate(posteriorpred7, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    posteriorpred7 <- filter(posteriorpred7, !(lowevolved==0 & evolved == 1))
    posteriorpred7 <- filter(posteriorpred7, !(highevolved==0 & evolved == 1))
    
    #1.5) Plot reaction norms
    KB2086.2 <- ggplot(data.B2086.2, aes(y=exp(log_Kmean)/1e5, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred7, inherit.aes = F, aes(x=pHmediumval, y=mean/1e5, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred7, inherit.aes = F, aes(ymax=PI97.5/1e5, ymin=PI2.5/1e5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
      ylab("K (1e5 indiv/mL)") + ylim(2, 9)   + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #4) Look at SB3539 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/", full.names=T))
    lapply(file_names[2:17],load,.GlobalEnv)
    
    #get data
    data.SB3539 <- filter(data.all, strain=="SB3539")
    
    #Compare models
    coeftab(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16)
    models <- compare(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- mK.09
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data.SB3539 <- filter(data.all, strainval=="SB3539")
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data.SB3539 <- filter(data.SB3539, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred6 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred6 <- mutate(posteriorpred6, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    posteriorpred6 <- filter(posteriorpred6, !(lowevolved==0 & evolved == 1))
    posteriorpred6 <- filter(posteriorpred6, !(highevolved==0 & evolved == 1))
    
    #1.5) Plot reaction norms
    KSB3539 <- ggplot(data.SB3539, aes(y=exp(log_Kmean)/1e5, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred6, inherit.aes = F, aes(x=pHmediumval, y=mean/1e5, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred6, inherit.aes = F, aes(ymax=PI97.5/1e5, ymin=PI2.5/1e5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
      ylim(2, 9) + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #5) Look at CU428.2 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/CU428.2/", full.names=T))
    lapply(file_names[2:17],load,.GlobalEnv)
    
    #get data
    data.CU428.2 <- filter(data.all, strainval=="CU428.2")
    
    #Compare models
    coeftab(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16)
    models <- compare(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- mK.08
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data.CU428.2 <- filter(data.all, strainval=="CU428.2")
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data <- filter(data, evolved != 1)
    data <- mutate(data, evolved = ifelse(lowevolved==1 | highevolved==1, 1, 0))
    data.CU428.2 <- filter(data.CU428.2, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(mK.01, mK.02, mK.03, mK.04, mK.05, mK.06, mK.07, mK.08, mK.09, mK.10, mK.11, mK.12, mK.13, mK.14, mK.15, mK.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred5 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred5 <- mutate(posteriorpred5, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    
    
    
    
    #1.5) Plot reaction norms
    KCU428.2 <- ggplot(data.CU428.2, aes(y=exp(log_Kmean)/1e5, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred5, inherit.aes = F, aes(x=pHmediumval, y=mean/1e5, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred5, inherit.aes = F, aes(ymax=PI97.5/1e5, ymin=PI2.5/1e5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
      ylim(2, 9) + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
}

#And now for alpha
{
  #2) Look at CU427.4 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU427.4/", full.names=T))
    lapply(file_names[1:16],load,.GlobalEnv)
    
    #Get data for CU427.4
    dataCU427.4 <- filter(data.all, strain=="CU427.4")
    
    #Compare models
    coeftab(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16)
    models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- malfa.10
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    dataCU427.4 <- filter(data.all, strainval=="CU427.4")
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    dataCU427.4 <- filter(dataCU427.4, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred9 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred9 <- mutate(posteriorpred9, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    
    #1.5) Plot reaction norms
    aCU427.4 <- ggplot(dataCU427.4, aes(y=exp(log_alfamean)*1e7, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred9, inherit.aes = F, aes(x=pHmediumval, y=mean*1e7, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred9, inherit.aes = F, aes(ymax=PI97.5*1e7, ymin=PI2.5*1e7, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
       xlab("pH of assay medium") + ylim(0, 12.5) + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title=element_text(size=16), axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(size=20, hjust = 0.5)) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #3) Look at B2086.2 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/", full.names=T))
    lapply(file_names[2:17],load,.GlobalEnv)
    
    #get data
    data.B2086.2 <- filter(data.all, strainval=="B2086.2")
    
    #Compare models
    coeftab(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16)
    models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- malfa.09
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data.B2086.2 <- filter(data.B2086.2, timing != "nocg")
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred10 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred10 <- mutate(posteriorpred10, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    
    #1.5) Plot reaction norms
    aB2086.2 <- ggplot(data.B2086.2, aes(y=exp(log_alfamean)*1e7, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred10, inherit.aes = F, aes(x=pHmediumval, y=mean*1e7, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred10, inherit.aes = F, aes(ymax=PI97.5*1e7, ymin=PI2.5*1e7, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
      ylab("α (1e7 mL/(h indiv))") + xlab("pH of assay medium") + ylim(0, 12.5) + 
      theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5)) +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #4) Look at SB3539 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/SB3539/", full.names=T))
    lapply(file_names[1:16],load,.GlobalEnv)
    
    #get data
    data.SB3539 <- filter(data.all, strainval=="SB3539")
    
    #Compare models
    coeftab(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16)
    models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- malfa.11
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    precis(best_model, digits=4, depth = 2, prob = 0.95)
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                       lowevolved=unique(data.all$lowevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred11 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred11 <- mutate(posteriorpred11, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    
    #1.5) Plot reaction norms
    aSB3539 <- ggplot(data.SB3539, aes(y=exp(log_alfamean)*10e6, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred11, inherit.aes = F, aes(x=pHmediumval, y=mean*10e6, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred11, inherit.aes = F, aes(ymax=PI97.5*10e6, ymin=PI2.5*10e6, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
       xlab("pH of assay medium") + ylim(0, 12.5) + 
      guides(color = FALSE, size = FALSE)+ theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title=element_text(size=16), axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(size=20, hjust = 0.5), legend.position = "none") +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))
  }
  
  #5) Look at CU428.2 models
  {
    file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU428.2/", full.names=T))
    lapply(file_names[1:16],load,.GlobalEnv)
    
    #get data
    data.CU428.2 <- filter(data.all, strainval=="CU428.2")
    
    #Compare models
    coeftab(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16)
    models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                      WAIC=FALSE, func=DIC)
    models
    best_model <- malfa.08
    precis(best_model, depth = 2, prob = 0.95)
    
    #1.2) Extract samples
    post <- extract.samples(best_model, n=16e3)
    
    #1.3) Calculate the posteriors with 95% probability interval
    data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                       lowevolved=unique(data.all$lowevolved))
    data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
    data <- filter(data, !(lowevolved==1 & highevolved == 1))
    data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
    
    #Create weighted posteriors based on DIC criterion
    post.pred <- ensemble(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                          func=DIC, data=data, n=16e3)
    
    mu.pred.mean <- apply( post.pred$link , 2 , mean )
    mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
    
    #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
    posteriorpred12 <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
    posteriorpred12 <- mutate(posteriorpred12, pHmediumval=ifelse(pHmedium==1, 6.5, 4.5), pHoriginval=ifelse(lowevolved==1, 4.5, 6.5), 
                            timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
    
    #1.5) Plot reaction norms
    
    aCU428.2 <- ggplot(data.CU428.2, aes(y=exp(log_alfamean)*1e7, x=ifelse(pHmedium==1, 6.5, 4.5), colour=paste(timing, pHorigin), group = paste(timing, pHorigin))) + geom_jitter(size=3, position = position_dodge2(width = 0.4)) + 
      geom_line(data=posteriorpred12, inherit.aes = F, aes(x=pHmediumval, y=mean*1e7, colour=paste(timingval, pHoriginval)), size=2) + 
      geom_ribbon(data=posteriorpred12, inherit.aes = F, aes(ymax=PI97.5*1e7, ymin=PI2.5*1e7, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
       xlab("pH of assay medium") + ylim(0, 12.5) + 
      guides(color = FALSE, size = FALSE)+ theme_light() + 
      theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
            axis.title=element_text(size=16), axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(size=20, hjust = 0.5), legend.position = "none") +
      scale_color_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + 
      scale_fill_manual(values=c("#2db3f9", "#e41a1c", "#2d33f9"), breaks=colorsch, name= "Origin", labels=c("NpH", "LpH", "ANC")) + scale_x_continuous(breaks = c(4.5, 6.5))

    }
}

library(ggpubr)
#Combine all figures
figure1 <- ggarrange(r0B2086.2, r0CU427.4, r0CU428.2, r0SB3539, 
          KB2086.2, KCU427.4, KCU428.2, KSB3539,
          aB2086.2, aCU427.4, aCU428.2, aSB3539,
          nrow=3, ncol=4, common.legend = T, legend="bottom",
          labels =list("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), align = 'v', hjust = -3, vjust = 1)

figure1
ggsave(plot = figure1,filename = "Fig2.png", path = "4_results/Figures", width = 350, height = 250, units = "mm", device = "png")
