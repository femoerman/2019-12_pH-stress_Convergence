# clear memory
rm(list=ls())

#Set working directory
setwd("F:/Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(rethinking)

#prepare data
sumdata$strainval <- sumdata$strain
#MAke factors numerical (needed for Bayesian model fitting with rethinking package)
sumdata <- mutate(sumdata, strain=ifelse(sumdata$strain=="CU427.4", 1, ifelse(sumdata$strain=="CU428.2", 2, ifelse(sumdata$strain=="SB3539", 3, 4))))
sumdata$evolved <- ifelse(sumdata$timing=="cg", 1, 0)
sumdata$highevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==6.5, 1, 0)
sumdata$lowevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==5, 1, 0)
sumdata$pHmedium <- ifelse(sumdata$pHmedium==6.5, 1, 0)

#Filter out data to only keep ancestors and evolved strains after common garden
data.all <- sumdata %>% filter(timing != "nocg")

#1) Look at global models
{load("2_data/14,16,18_EvolutionModelsBayesianContinued/14,16,18_EvolutionModelsBayesianContinued/")
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/global/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #Compare models
  models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                    WAIC=FALSE, func=DIC)
 
  modeltorun <- malfa.07
  precis(modeltorun, depth=2)
  
  #1.2) Extract samples
  precis(modeltorun,  depth=2)
  post <- extract.samples(modeltorun, n=8e3)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(strain=unique(data.all$strain), pHorigin=unique(data.all$pHorigin), pHmedium=unique(data.all$pHmedium), timing=unique(data.all$timing))
  data <- filter(data, (timing==1 | pHorigin == 1))
  post.pred=link(modeltorun, data, n=8e3)
  
  mu.pred.mean <- apply( post.pred , 2 , mean )
  mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred <- mutate(posteriorpred, pHmediumval=ifelse(pHmedium==1, 6.5, 5), pHoriginval=ifelse(pHorigin==1, 6.5, 5), timingval=ifelse(timing==1, "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  ggplot(data.all, aes(y=exp(log_alfamean), x=pHmediumval, colour=paste(timingval, pHoriginval))) + geom_point()  + 
    geom_line(data=posteriorpred, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval))) + 
    facet_wrap(~strain) + geom_ribbon(data=posteriorpred, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ggtitle("alfa model for all strains") + ylab("mean alfa") + xlab("pH of medium") + 
    scale_fill_discrete(name="Origin", breaks=c("ancestral 6.5", "cg 5", "cg 6.5"), labels=c("ancestral", "low-evolved", "high-evolved")) + 
    guides(color = FALSE, size = FALSE)+ theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5))
  ggsave(plot=last_plot(), filename = "2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/global/plotalfa.png")
  
}

#2) Look at CU427.4 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU427.4/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #Get data for CU427.4
  dataCU427.4 <- filter(data.all, strain=="CU427.4")
  
  #Compare models
  models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- malfa.10
  precis(best_model, depth=2)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e4)
  
  #1.3) Calculate the posteriors with 95% probability interval
  dataCU427.4 <- filter(data.all, strainval=="CU427.4")
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  dataCU427.4 <- filter(dataCU427.4, timing != "nocg")
  
  post.pred=link(best_model, data, n=16e4)
  
  mu.pred.mean <- apply( post.pred , 2 , mean )
  mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred <- mutate(posteriorpred, pHmediumval=ifelse(pHmedium==1, 6.5, 5), pHoriginval=ifelse(lowevolved==1, 5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  ggplot(dataCU427.4, aes(y=exp(log_alfamean), x=ifelse(pHmedium==1, 6.5, 5), colour=paste(timing, pHorigin))) + geom_point(size=2)  + 
    geom_line(data=posteriorpred, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ggtitle("Strain CU427.4") + ylab("Intraspecific competitive ability alpha") + xlab("pH of test-medium") + 
    scale_fill_discrete(name="Origin", breaks=c("ancestral 6.5", "cg 5", "cg 6.5"), labels=c("ancestral", "low-evolved", "high-evolved")) + 
    guides(color = FALSE, size = FALSE)+ theme_light() + 
    theme(axis.text=element_text(size=20), legend.text=element_text(size=16), strip.text.x=element_text(20),
          axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))
  ggsave(plot=last_plot(), filename = "2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU427.4/plotalfa.png")
}

#3) Look at B2086.2 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.B2086.2 <- filter(data.all, strainval=="B2086.2")
  
  #Compare models
  models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- malfa.09
  precis(best_model, depth=2)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e4)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), lowevolved=unique(data.all$lowevolved), highevolved=unique(data.all$highevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data.B2086.2 <- filter(data.B2086.2, timing != "nocg")
  
  post.pred=link(best_model, data, n=16e4)
  
  mu.pred.mean <- apply( post.pred , 2 , mean )
  mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred <- mutate(posteriorpred, pHmediumval=ifelse(pHmedium==1, 6.5, 5), pHoriginval=ifelse(lowevolved==1, 5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  ggplot(data.B2086.2, aes(y=exp(log_alfamean), x=ifelse(pHmedium==1, 6.5, 5), colour=paste(timing, pHorigin))) + geom_point(size=2)  + 
    geom_line(data=posteriorpred, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ggtitle("Strain B2086.2") + ylab("Intraspecific competitive ability alpha") + xlab("pH of test-medium") + 
    scale_fill_discrete(name="Origin", breaks=c("ancestral 6.5", "cg 5", "cg 6.5"), labels=c("ancestral", "low-evolved", "high-evolved")) + 
    guides(color = FALSE, size = FALSE)+ theme_light() + 
    theme(axis.text=element_text(size=20), legend.text=element_text(size=16), strip.text.x=element_text(20),
          axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))
  
  ggsave(plot=last_plot(), filename = "2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/plotalfa.png")
}

#4) Look at SB3539 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/SB3539/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.SB3539 <- filter(data.all, strainval=="SB3539")
  
  #Compare models
  models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- malfa.11
  precis(best_model, depth=2)
  
  #1.2) Extract samples
  precis(best_model, digits=4, depth=2)
  post <- extract.samples(best_model, n=16e4)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                     lowevolved=unique(data.all$lowevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
  
  post.pred=link(best_model, data, n=16e4)
  mu.pred.mean <- apply( post.pred , 2 , mean )
  mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred <- mutate(posteriorpred, pHmediumval=ifelse(pHmedium==1, 6.5, 5), pHoriginval=ifelse(lowevolved==1, 5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  ggplot(data.SB3539, aes(y=exp(log_alfamean), x=ifelse(pHmedium==1, 6.5, 5), colour=paste(timing, pHorigin))) + geom_point(size=2)  + 
    geom_line(data=posteriorpred, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ggtitle("Strain SB3539") + ylab("Intraspecific competitive ability alpha") + xlab("pH of test-medium") + 
    scale_fill_discrete(name="Origin", breaks=c("ancestral 6.5", "cg 5", "cg 6.5"), labels=c("ancestral", "low-evolved", "high-evolved")) + 
    guides(color = FALSE, size = FALSE)+ theme_light() + 
    theme(axis.text=element_text(size=20), legend.text=element_text(size=16), strip.text.x=element_text(20),
          axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))
  
  ggsave(plot=last_plot(), filename = "2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/SB3539/plotalfa.png")
}

#5) Look at CU428.2 models
{
  file_names <- paste(list.files("2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU428.2/", full.names=T))
  lapply(file_names[1:16],load,.GlobalEnv)
  
  #get data
  data.CU428.2 <- filter(data.all, strainval=="CU428.2")
  
  #Compare models
  models <- compare(malfa.01, malfa.02, malfa.03, malfa.04, malfa.05, malfa.06, malfa.07, malfa.08, malfa.09, malfa.10, malfa.11, malfa.12, malfa.13, malfa.14, malfa.15, malfa.16,
                    WAIC=FALSE, func=DIC)
  models
  best_model <- malfa.08
  precis(best_model, depth=2)
  
  #1.2) Extract samples
  post <- extract.samples(best_model, n=16e4)
  
  #1.3) Calculate the posteriors with 95% probability interval
  data<- expand.grid(evolved=unique(data.all$evolved), pHmedium=unique(data.all$pHmedium), highevolved=unique(data.all$highevolved), 
                     lowevolved=unique(data.all$lowevolved))
  data <- filter(data, (evolved==1 & (lowevolved==0 & highevolved==0)) | (evolved == 0))
  data <- filter(data, !(lowevolved==1 & highevolved == 1))
  data$evolved <- ifelse(data$lowevolved == 1 | data$highevolved == 1, 1, 0)
  
  post.pred=link(best_model, data, n=16e4)
  mu.pred.mean <- apply( post.pred , 2 , mean )
  mu.pred.PI <- apply(post.pred, 2, PI, prob=0.95)
  
  #1.4) Backtransform r0 predictions from log scale, and create new variables with original pH and timing values
  posteriorpred <- cbind(data, data.frame(mean=exp(mu.pred.mean), PI97.5=exp(mu.pred.PI[2, ]), PI2.5=exp(mu.pred.PI[1, ])))
  posteriorpred <- mutate(posteriorpred, pHmediumval=ifelse(pHmedium==1, 6.5, 5), pHoriginval=ifelse(lowevolved==1, 5, 6.5), 
                          timingval=ifelse(evolved==1 | (lowevolved == 1 | highevolved == 1), "cg", "ancestral"))
  
  #1.5) Plot reaction norms
  ggplot(data.CU428.2, aes(y=exp(log_alfamean), x=ifelse(pHmedium==1, 6.5, 5), colour=paste(timing, pHorigin))) + geom_point(size=2)  + 
    geom_line(data=posteriorpred, inherit.aes = F, aes(x=pHmediumval, y=mean, colour=paste(timingval, pHoriginval)), size=2) + 
    geom_ribbon(data=posteriorpred, inherit.aes = F, aes(ymax=PI97.5, ymin=PI2.5, x=pHmediumval, fill=paste(timingval, pHoriginval)),alpha=0.3) + 
    ggtitle("Strain CU428.2") + ylab("Intraspecific competitive ability alpha") + xlab("pH of test-medium") + 
    scale_fill_discrete(name="Origin", breaks=c("ancestral 6.5", "cg 5", "cg 6.5"), labels=c("ancestral", "low-evolved", "high-evolved")) + 
    guides(color = FALSE, size = FALSE)+ theme_light() + 
    theme(axis.text=element_text(size=20), legend.text=element_text(size=16), strip.text.x=element_text(20),
          axis.title=element_text(size=20), plot.title = element_text(size=20, hjust = 0.5))
  
  ggsave(plot=last_plot(), filename = "2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/CU428.2/plotalfa.png")
}
