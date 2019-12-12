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

#1) Prepare data for variance comparison
{
  #1.1) Calculate absolute difference between r0, K, alpha of observations and overall mean
  data.low <- filter(data.all, pHmedium==4.5)
  data.low$log_r0diff <- abs(data.low$log_r0mean-mean(data.low$log_r0mean))
  data.low$log_Kdiff <- abs(data.low$log_Kmean-mean(data.low$log_Kmean))
  data.low$log_alphadiff <- abs(data.low$log_alfamean-mean(data.low$log_alfamean))
  data.neutral <- filter(data.all, pHmedium==6.5)
  data.neutral$log_r0diff <- abs(data.neutral$log_r0mean-mean(data.neutral$log_r0mean))
  data.neutral$log_Kdiff <- abs(data.neutral$log_Kmean-mean(data.neutral$log_Kmean))
  data.neutral$log_alphadiff <- abs(data.neutral$log_alfamean-mean(data.neutral$log_alfamean))
  
  #Plot the distributions
  ggplot(data.low, aes(y=(log_r0diff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.low, aes(y=(log_Kdiff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.low, aes(y=(log_alphadiff), colour = origin, x = origin)) + geom_boxplot()
  
  ggplot(data.neutral, aes(y=(log_r0diff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.neutral, aes(y=(log_Kdiff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.neutral, aes(y=(log_alphadiff), colour = origin, x = origin)) + geom_boxplot()
}

#2) Do the models for r0
{
  #2.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-0.7, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )

  
  mr0.03 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )


  mr0.04 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestr0low <- mr0.03
  precis(bestr0low, depth=2, prob = 0.95)
  
  #2.2) At neutral pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-0.7, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dnorm(0.5, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-0.7, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestr0neutral <- mr0.01
  precis(bestr0neutral, depth=2, prob=0.95)
}

#3) Do the models for K
{
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-1.9, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestKlow <- mr0.02
  precis(bestKlow, depth=2, prob=0.95)
  
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-1.9, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-1.9, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )

  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestKneutral <- mr0.03
  precis(bestKneutral, depth=2, prob = 0.95)
}

#4) Do the models for alpha
{
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-1.4, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestalphalow <- mr0.03
  precis(bestalphalow, depth=2)
  
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a  ,
      a ~dlnorm(-1.4, 0.5),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bev ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bl ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bh ~ dnorm(0, 1),
      
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved  + bl*lowevolved ,
      a ~dlnorm(-1.4, 0.5),
      
      bh ~ dnorm(0, 0.5),
      
      bl ~ dnorm(0, 0.5),
      
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05)
  bestalphaneutral <- mr0.02
  precis(bestalphaneutral, depth=2)
}

#5) Save the selected models
{
  output_rKalphaconvWO <- list(r0low = bestr0low, r0neutral = bestr0neutral, Klow = bestKlow, Kneutral = bestKneutral, alphalow = bestalphalow, alphaneutral = bestalphaneutral)
  save(output_rKalphaconvWO, file="2_data/rKalphaConvergenceModelsWO.RData")
  load("2_data/rKalphaConvergenceModels.RData")
}
typeof(output_rKalphaconv[[1]])

#5) Look at distribution of r0, K and alpha
{
  #3.2) At low pH
  {
    #r0
    f1 <- ggplot(data.low, aes(y=log_r0diff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ggtitle("Low pH of medium") + ylab("Diff. to mean r0 (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #K
    f2 <- ggplot(data.low, aes(y=log_Kdiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3)+ylab("Diff. to mean K (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #alpha
    f3 <- ggplot(data.low, aes(y=log_alphadiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) +ylab("Diff. to mean alpha (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.low$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
  }
  
  #3.2) At neutral pH
  {
    #r0
    f4 <-ggplot(data.neutral, aes(y=log_r0diff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ggtitle("Neutral pH of medium") + ylab("Diff. to mean r0 (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #K
    f5 <- ggplot(data.neutral, aes(y=log_Kdiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean K (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #alpha
    f6 <- ggplot(data.neutral, aes(y=log_alphadiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean alpha (log)") + 
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#e41a1c", "#2db3f9",  "#2d33f9"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
  }
  
ggarrange(f1, f4, f2, f5, f3, f6, nrow=3, ncol=2, align = "v", common.legend = T, legend =  "bottom")
}