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
  
  #First get data at low pH of assay medium
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
  
  #Secondly, get data at neutral pH of assay medium
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
  
  
  #Plot the distributions
  ggplot(data.low, aes(y=exp(log_r0diff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.low, aes(y=exp(log_Kdiff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.low, aes(y=exp(log_alphadiff), colour = origin, x = origin)) + geom_boxplot()
  
  ggplot(data.neutral, aes(y=exp(log_r0diff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.neutral, aes(y=exp(log_Kdiff), colour = origin, x = origin)) + geom_boxplot()
  ggplot(data.neutral, aes(y=exp(log_alphadiff), colour = origin, x = origin)) + geom_boxplot()
  

}

#2) Do the models for r0
{
  #2.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )

  
  mr0.03 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )


  mr0.04 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #without random variables
  mr0.06 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a,
      a ~dlnorm(-1.2, 0.3),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-1.2, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsr0low <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsr0low, file = "2_data/convergenceModelsr0Low.RData")
  bestr0low <- mr0.08
  precis(bestr0low, depth=2, prob = 0.95)
  
  #2.2) At neutral pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #Without random variables
  mr0.06 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-1.6, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-1.6, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_r0diff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #Compare models; and save them
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsr0neut <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsr0neut, file = "2_data/21, 23_r-K-alpha_Convergence/convergenceModelsr0Neutral.RData")
  bestr0neutral <- mr0.02
  precis(bestr0neutral, depth=2, prob=0.95)
}

#3) Do the models for K
{
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-2.3, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-2.3, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-2.3, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-2.3, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-2.3, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #without random variables
  mr0.06 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-2.3, 0.3),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-2.3, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-2.3, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-2.3, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-2.3, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsKlow <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsKlow, file = "2_data/21, 23_r-K-alpha_Convergence/convergenceModelsKLow.RData")
  bestKlow <- mr0.07
  precis(bestKlow, depth=2, prob=0.95)
  
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )

  #Without random variables
  mr0.06 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      a ~dlnorm(-1.6, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-1.6, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-1.6, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_Kdiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-1.6, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsKneut <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsKneut, file = "2_data/21, 23_r-K-alpha_Convergence/convergenceModelsKneutral.RData")
  bestKneutral <- mr0.08
  precis(bestKneutral, depth=2, prob = 0.95)
}

#4) Do the models for alpha
{
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #without random variables
  mr0.06 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-1.2, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.low ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsalphalow <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsalphalow, file = "2_data/21, 23_r-K-alpha_Convergence/convergenceModelsAlphaLow.RData")
  bestalphalow <- mr0.08
  precis(bestalphalow, depth=2, prob = 0.95)
  
  #3.1) At low pH, compare distribution of r0 between groups
  mr0.01 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain] ,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.02 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bev*evolved + bevs[strain]*evolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bev ~ dnorm(0, 1),
      bevs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.03 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.04 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.05 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + as[strain]+ bh*highevolved + bhs[strain]*highevolved + bl*lowevolved + bls[strain]*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      as[strain] ~dnorm(0, 1),
      bh ~ dnorm(0, 1),
      bhs[strain] ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      bls[strain] ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  #Without random variables
  mr0.06 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a,
      a ~dlnorm(-1.2, 0.3),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.07 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bev*evolved,
      a ~dlnorm(-1.2, 0.3),
      bev ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.08 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  
  mr0.09 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  mr0.10 <- map2stan(
    alist(
      log_alphadiff ~ dnorm(mu,sigma),
      mu <- a + bh*highevolved + bl*lowevolved,
      a ~dlnorm(-1.2, 0.3),
      bh ~ dnorm(0, 1),
      bl ~ dnorm(0, 1),
      sigma ~dcauchy(0, 0.5)
    ) ,data=data.neutral ,
    iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  
  compare(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10)
  coeftab(mr0.01, mr0.02, mr0.03, mr0.04, mr0.05, mr0.06, mr0.07, mr0.08, mr0.09, mr0.10, digits = 2)
  modelsalphaneut <- list(m1 = mr0.01, m2 = mr0.02, m3 = mr0.03, m4 = mr0.04, m5 = mr0.05, m6 = mr0.06, m7 = mr0.07, m8 = mr0.08, m9 = mr0.09, m10 = mr0.10)
  save(modelsalphaneut, file = "2_data/21, 23_r-K-alpha_Convergence/convergenceModelsAlphaNeutral.RData")
  bestalphaneutral <- mr0.05
  precis(bestalphaneutral, depth=2, prob = 0.95)
}

#5) Save the selected models
{
  output_rKalphaconv <- list(r0low = bestr0low, r0neutral = bestr0neutral, Klow = bestKlow, Kneutral = bestKneutral, alphalow = bestalphalow, alphaneutral = bestalphaneutral)
  save(output_rKalphaconv, file="2_data/21, 23_r-K-alpha_Convergence/rKalphaConvergenceModels.RData")
  load("2_data/21, 23_r-K-alpha_Convergence/rKalphaConvergenceModels.RData")
}
typeof(output_rKalphaconv[[1]])

#5) Look at distribution of r0, K and alpha
{
  #3.2) At low pH
  {
    #r0
    f1 <- ggplot(data.low, aes(y=log_r0diff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean r0 (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #K
    f2 <- ggplot(data.low, aes(y=log_Kdiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3)+ylab("Diff. to mean K (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #alpha
    f3 <- ggplot(data.low, aes(y=log_alphadiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) +ylab("Diff. to mean alpha (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
  }
  
  #3.2) At neutral pH
  {
    #r0
    f4 <-ggplot(data.neutral, aes(y=log_r0diff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean r0 (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #K
    f5 <- ggplot(data.neutral, aes(y=log_Kdiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean K (log)") + 
      theme_light()+ theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                           axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    #alpha
    f6 <- ggplot(data.neutral, aes(y=log_alphadiff, colour=origin, fill=origin, x = origin)) + geom_boxplot(alpha=0.3) + ylab("Diff. to mean alpha (log)") + 
      theme_light() + theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), plot.title = element_text(size=20, hjust = 0.5), axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
  }
  
  #Show the confergence of the strains in figures
  {
    #1.1) r0
    f1.1<- ggplot(data.low, aes(y=log_r0mean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      ylab("r0 (log)") + xlab("origin of population")  + 
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
    f1.2 <-ggplot(data.neutral, aes(y=log_r0mean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      ylab("r0 (log)") + xlab("origin of population")  + 
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
    #1.2) K
    f2.1 <-ggplot(data.low, aes(y=log_Kmean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      xlab("origin of population")  + ylab("K (individuals/ml) (log)") +
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
    f2.2 <-ggplot(data.neutral, aes(y=log_Kmean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      xlab("origin of population")  + ylab("K (individuals/ml) (log)") +
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
    #1.3) alpha
    f3.1 <-ggplot(data.low, aes(y=log_alfamean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      ylab("alpha (log)") + xlab("origin of population")  + 
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
    f3.2 <- ggplot(data.neutral, aes(y=log_alfamean, x=paste(origin, strain), colour = origin, fill=origin)) + geom_boxplot(alpha=0.3) +
      ylab("alpha (log)") + xlab("origin of population")  + 
      theme_light() + 
      theme(axis.text=element_text(size=14), legend.text=element_text(size=14), strip.text.x=element_text(size=14),
            axis.title=element_text(size=14), plot.title = element_text(size=16, hjust = 0.5),
            axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank())+
      scale_color_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor")) + 
      scale_fill_manual(values=c("#2db3f9",  "#2d33f9", "#e41a1c"), breaks=unique(data.neutral$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestor"))
    
    
  }
  
flow <- ggarrange(f1.1, f1, f2.2, f2, f3.1, f3,  nrow=3, ncol=2, align = "hv", common.legend = T, legend =  "bottom")
flow <- annotate_figure(flow, top = text_grob("Low pH of medium", color = "black", face = "bold", size = 25))
fneut <- ggarrange(f1.2, f4, f2.2, f5, f3.2, f6, nrow=3, ncol=2, align = "hv", common.legend = T, legend =  "bottom")
fneut <- annotate_figure(fneut, top = text_grob("Neutral pH of medium", color = "black", face = "bold", size = 25))

ggarrange(flow, "blank",  fneut, ncol = 3, nrow = 1, align = "hv", common.legend = T, legend =  "bottom", widths = c(5, 1, 5))
}
