# clear memory
rm(list=ls())

#Set working directory
# setwd("F:/Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(rethinking)
library(ggfortify)
library(ggpubr)
library(car)

#Load data and prepare for usage
load("2_data/posteriorBH_means.RData")
sumdata$strainval <- sumdata$strain
sumdata$evolved <- ifelse(sumdata$timing=="cg", 1, 0)
sumdata$pHmedium <- ifelse(sumdata$pHmedium==6.5, 6.5, 4.5)
sumdata$pHorigin <- ifelse(sumdata$pHorigin==6.5, 6.5, 4.5)
sumdata$highevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==6.5, 1, 0)
sumdata$lowevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==4.5, 1, 0)
sumdata$pHmedium <- as.factor(sumdata$pHmedium)
data.all <- sumdata %>% filter(timing != "nocg")
data.all <- mutate(data.all, origin=ifelse(timing!="cg", "ancestor", ifelse(pHorigin==6.5, "high-evolved", "low-evolved")))

#Order the strains and turn them into a factor
data.all <- arrange(data.all, strain)
data.all$strain <- as.factor(data.all$strain)


#1) Prepare data for variance comparison
{
  #1.1) Calculate absolute difference between r0, K, alpha of observations and overall mean
  data.low <- filter(data.all, pHmedium==4.5)
  
  data.low$r0mean <- ifelse(data.low$origin =='ancestor', mean(filter(data.low, data.low$origin == "ancestor")$log_r0mean), ifelse(data.low$origin == "low-evolved",
                                                                                                                                   mean(filter(data.low, data.low$origin == "low-evolved")$log_r0mean),mean(filter(data.low, data.low$origin == "high-evolved")$log_r0mean)))
  data.low$Kmean <- ifelse(data.low$origin =='ancestor', mean(filter(data.low, data.low$origin == "ancestor")$log_Kmean), ifelse(data.low$origin == "low-evolved",
                                                                                                                                 mean(filter(data.low, data.low$origin == "low-evolved")$log_Kmean),mean(filter(data.low, data.low$origin == "high-evolved")$log_Kmean)))
  data.low$alphamean <- ifelse(data.low$origin =='ancestor', mean(filter(data.low, data.low$origin == "ancestor")$log_alfamean), ifelse(data.low$origin == "low-evolved",
                                                                                                                                        mean(filter(data.low, data.low$origin == "low-evolved")$log_alfamean),mean(filter(data.low, data.low$origin == "high-evolved")$log_alfamean)))
  data.low$log_r0diff <- abs(data.low$log_r0mean-data.low$r0mean)
  data.low$log_Kdiff <-  abs(data.low$log_Kmean-data.low$Kmean)
  data.low$log_alphadiff <-  abs(data.low$log_alfamean-data.low$alphamean)
  
  
  data.neutral <- filter(data.all, pHmedium==6.5)
  data.neutral$r0mean <- ifelse(data.neutral$origin =='ancestor', mean(filter(data.neutral, data.neutral$origin == "ancestor")$log_r0mean), ifelse(data.neutral$origin == "low-evolved",
                                                                                                                                                   mean(filter(data.neutral, data.neutral$origin == "low-evolved")$log_r0mean),mean(filter(data.neutral, data.neutral$origin == "high-evolved")$log_r0mean)))
  data.neutral$Kmean <- ifelse(data.neutral$origin =='ancestor', mean(filter(data.neutral, data.neutral$origin == "ancestor")$log_Kmean), ifelse(data.neutral$origin == "low-evolved",
                                                                                                                                                 mean(filter(data.neutral, data.neutral$origin == "low-evolved")$log_Kmean),mean(filter(data.neutral, data.neutral$origin == "high-evolved")$log_Kmean)))
  data.neutral$alphamean <- ifelse(data.neutral$origin =='ancestor', mean(filter(data.neutral, data.neutral$origin == "ancestor")$log_alfamean), ifelse(data.neutral$origin == "low-evolved",
                                                                                                                                                        mean(filter(data.neutral, data.neutral$origin == "low-evolved")$log_alfamean),mean(filter(data.neutral, data.neutral$origin == "high-evolved")$log_alfamean)))
  data.neutral$log_r0diff <- abs(data.neutral$log_r0mean-data.neutral$r0mean)
  data.neutral$log_Kdiff <-  abs(data.neutral$log_Kmean-data.neutral$Kmean)
  data.neutral$log_alphadiff <-  abs(data.neutral$log_alfamean-data.neutral$alphamean)
  colorsch <- unique(data.all$origin)
}

#2) Load models from convergence
{
  load("2_data/convergenceModelsr0Low.RData")
  load("2_data/convergenceModelsr0Neutral.RData")
  load("2_data/convergenceModelsKLow.RData")
  load("2_data/convergenceModelsKneutral.RData")
  load("2_data/convergenceModelsAlphaLow.RData")
  load("2_data/convergenceModelsAlphaNeutral.RData")
}

#3) Make figures showing raw data
{
  #1.1) r0
  f1.1<- ggplot(data.low, aes(y=exp(log_r0mean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab(expression("r"[0]*" (1/h)")) + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(0.03, 0.4)) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
  f1.2 <-ggplot(data.neutral, aes(y=exp(log_r0mean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab(expression("r"[0]*" (1/h)")) + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(0.03, 0.4)) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
  #1.2) K
  f2.1 <-ggplot(data.low, aes(y=exp(log_Kmean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab("K (indiv/mL)") + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(exp(12.5), exp(13.7))) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
  f2.2 <-ggplot(data.neutral, aes(y=exp(log_Kmean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab("K (indiv/mL)") + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(exp(12.5), exp(13.7))) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
  #1.3) alpha
  f3.1 <-ggplot(data.low, aes(y=exp(log_alfamean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab("α (mL/(h indiv))") + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(exp(-16.5), exp(-13.5))) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
  f3.2 <- ggplot(data.neutral, aes(y=exp(log_alfamean), x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
    ylab("α (mL/(h indiv))") + xlab("Origin of population")  + 
    coord_trans(y = "log", limy = c(exp(-16.5), exp(-13.5))) +
    theme_light() + 
    theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
          axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("NpH", "LpH", "ANC"))
  
  
}

#4) Draw predictions based on weighted models
{
  #4.1) At low pH
  {
    #4.1.1) r0
    {
      data<- expand.grid(evolved=unique(data.low$evolved), lowevolved=unique(data.low$lowevolved), highevolved=unique(data.low$highevolved), strain = unique(data.low$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsr0low$m1, modelsr0low$m2, modelsr0low$m3, modelsr0low$m4, modelsr0low$m5, 
                            modelsr0low$m6, modelsr0low$m7, modelsr0low$m8, modelsr0low$m9, modelsr0low$m10, data=data, n=16e3)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred1 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred1 <- mutate(posteriorpred1, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      

      fp1 <- ggplot(posteriorpred1, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) + ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.low, aes(x = paste(origin, strain), y = log_r0diff, group=strain)) + 
        theme_light() +  ylab(expression("Diff. to mean r"[0]*" (ln)")) + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
    
    #4.1.2) K
    {
      data<- expand.grid(evolved=unique(data.low$evolved), lowevolved=unique(data.low$lowevolved), highevolved=unique(data.low$highevolved), strain = unique(data.low$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsKlow$m1, modelsKlow$m2, modelsKlow$m3, modelsKlow$m4, modelsKlow$m5, 
                            modelsKlow$m6, modelsKlow$m7, modelsKlow$m8, modelsKlow$m9, modelsKlow$m10, data=data, n=16e3)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred2 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred2 <- mutate(posteriorpred2, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      
      
      fp2 <- ggplot(posteriorpred2, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) +ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.low, aes(x = paste(origin, strain), y = log_Kdiff, group=strain)) + 
        theme_light() +  ylab("Diff. to mean K (ln)") + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
    
    #4.1.3) alpha
    {
      data<- expand.grid(evolved=unique(data.low$evolved), lowevolved=unique(data.low$lowevolved), highevolved=unique(data.low$highevolved), strain = unique(data.low$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsalphalow$m1, modelsalphalow$m2, modelsalphalow$m3, modelsalphalow$m4, modelsalphalow$m5, 
                            modelsalphalow$m6, modelsalphalow$m7, modelsalphalow$m8, modelsalphalow$m9, modelsalphalow$m10, data=data, n=16e3)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred3 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred3 <- mutate(posteriorpred3, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      
      
      fp3 <- ggplot(posteriorpred3, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) + ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.low, aes(x = paste(origin, strain), y = log_alphadiff, group=strain)) + 
        theme_light() +  ylab("Diff. to mean α (ln)") + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
  }
  
  #4.2) At neutral pH
  {
    #4.1.1) r0
    {
      data<- expand.grid(evolved=unique(data.neutral$evolved), lowevolved=unique(data.neutral$lowevolved), highevolved=unique(data.neutral$highevolved), strain = unique(data.neutral$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsr0neut$m1, modelsr0neut$m2, modelsr0neut$m3, modelsr0neut$m4, modelsr0neut$m5, 
                            modelsr0neut$m6, modelsr0neut$m7, modelsr0neut$m8, modelsr0neut$m9, modelsr0neut$m10, data=data, n=16e3)
      compare(modelsr0neut$m1, modelsr0neut$m2, modelsr0neut$m3, modelsr0neut$m4, modelsr0neut$m5, 
                            modelsr0neut$m6, modelsr0neut$m7, modelsr0neut$m8, modelsr0neut$m9, modelsr0neut$m10)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred4 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred4 <- mutate(posteriorpred4, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      posteriorpred4$strain <- as.character(posteriorpred4$strain)
      
      fp4 <- ggplot(posteriorpred4, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) + ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.neutral, aes(x = paste(origin, strain), y = log_r0diff, group=strain)) + 
        theme_light() +  ylab(expression("Diff. to mean r"[0]*" (ln)")) + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
    
    #4.1.2) K
    {
      data<- expand.grid(evolved=unique(data.neutral$evolved), lowevolved=unique(data.neutral$lowevolved), highevolved=unique(data.neutral$highevolved), strain = unique(data.neutral$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsKneut$m1, modelsKneut$m2, modelsKneut$m3, modelsKneut$m4, modelsKneut$m5, 
                            modelsKneut$m6, modelsKneut$m7, modelsKneut$m8, modelsKneut$m9, modelsKneut$m10, data=data, n=16e3)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred5 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred5 <- mutate(posteriorpred5, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      
      
      fp5 <- ggplot(posteriorpred5, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) +ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.neutral, aes(x = paste(origin, strain), y = log_Kdiff, group=strain)) + 
        theme_light() +  ylab("Diff. to mean K (ln)") + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
    
    #4.1.2) alpha
    {
      data<- expand.grid(evolved=unique(data.neutral$evolved), lowevolved=unique(data.neutral$lowevolved), highevolved=unique(data.neutral$highevolved), strain = unique(data.neutral$strain))
      data <- filter(data, !(lowevolved==1 & highevolved==1))
      data <- filter(data, !(evolved==0 & highevolved == 1))
      data <- filter(data, !(evolved==0 & lowevolved == 1))
      data <- filter(data, !(evolved==1 & (lowevolved==0 & highevolved==0)))
      post.pred <- ensemble(modelsalphaneut$m1, modelsalphaneut$m2, modelsalphaneut$m3, modelsalphaneut$m4, modelsalphaneut$m5, 
                            modelsalphaneut$m6, modelsalphaneut$m7, modelsalphaneut$m8, modelsalphaneut$m9, modelsalphaneut$m10, data=data, n=16e3)
      
      mu.pred.mean <- apply( post.pred$link , 2 , mean )
      mu.pred.PI <- apply(post.pred$link, 2, PI, prob=0.95)
      
      posteriorpred6 <- cbind(data, data.frame(mean=mu.pred.mean, PI97.5=mu.pred.PI[2, ], PI2.5=mu.pred.PI[1, ]))
      posteriorpred6 <- mutate(posteriorpred6, origin=ifelse(evolved==0, "ancestor", ifelse(lowevolved == 1, "low-evolved", "high-evolved")))
      posteriorpred6$xpos <- paste(posteriorpred6$origin, posteriorpred6$strain)
      data.neutral$xpos <- paste(data.neutral$origin, data.neutral$strain)
      
      fp6 <- ggplot(posteriorpred6, aes(x=paste(origin, strain), y=mean, fill=origin, color=origin)) +
        guides(fill=F) +ylim(-0.1, 1.2) + 
        geom_boxplot(aes(ymin = PI2.5, ymax = PI97.5, lower = PI2.5, upper=PI97.5, middle=mean, color=NA), stat="identity", alpha = 0.3)+ 
        geom_boxplot(aes(ymin = mean, ymax = mean, lower = mean, upper=mean, middle=mean, fill=NA), stat="identity", alpha = 0.3)+
        geom_point(inherit.aes = F, data=data.neutral, aes(x = paste(origin, strain), y = log_alphadiff, group=strain)) + 
        theme_light() + ylab("Diff. to mean α (ln)") + 
        theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(size=20, hjust = 0.5) ) +
        scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC")) + 
        scale_fill_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=colorsch, name="origin", labels=c("NpH", "LpH", "ANC"))
      
      
    }
  }
}

#5) Combine figures
{
  flow <- ggarrange(f1.1, fp1, f2.2, fp2, f3.1, fp3,  nrow=3, ncol=2, align = "hv", common.legend = T, legend =  "bottom",
                    labels =list("A", "B", "E", "F", "I", "J"))
  fneut <- ggarrange(f1.2, fp4, f2.2, fp5, f3.2, fp6, nrow=3, ncol=2, align = "hv", common.legend = T, legend =  "bottom",
                     labels =list("C", "D", "G", "H", "K", "L"))
  
  ggarrange(flow, "blank",  fneut, ncol = 3, nrow = 1, align = "hv", common.legend = T, legend =  "bottom", widths = c(5, 1, 5))
  ggsave(filename = "4_results/Figures/FigS3.png", device = "png", height = 9, width = 12, units = "in")
  }
