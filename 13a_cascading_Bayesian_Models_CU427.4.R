# clear memory
rm(list=ls())

#Set working directory
#setwd("D://Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(nlme)
library(rethinking)

#Load dataset
load("2_data/posteriorBH_means.RData")
sumdata <- filter(sumdata, timing != "nocg")
sumdata <- filter(sumdata, strain == "CU427.4")
sumdata$strainval <- sumdata$strain
sumdata$evolved <- ifelse(sumdata$timing=="cg", 1, 0)
sumdata$highevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==6.5, 1, 0)
sumdata$lowevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==5, 1, 0)
sumdata <- sumdata[complete.cases(sumdata), ]

#Prepare data as list
m1data <- list(
  log_r0_obs = sumdata$log_r0mean,
  log_r0_sd=sumdata$log_r0sd,
  pHmedium = ifelse(sumdata$pHmedium==6.5, 1, 0),
  evolved = sumdata$evolved,
  lowevolved = sumdata$lowevolved,
  highevolved = sumdata$highevolved,
  log_alfa_obs = sumdata$log_alfamean,
  log_alfa_sd=sumdata$log_alfasd,
  log_K_obs = sumdata$log_Kmean,
  log_K_sd=sumdata$log_Ksd,
  d_obs=sumdata$dmean,
  d_sd=sumdata$dsd
)

#1) Fit cascading Bayesian model for r0
#1.01) r0 ~ a
{
  mr0.01 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a ,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.01,digits=4, depth=2)
save(mr0.01, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.01.RData")


#1.02) r0 ~ a + pHmedium
{
  mr0.02 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium *pHmedium,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.02,digits=4, depth=2)
save(mr0.02, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.02.RData")


#1.03) r0 ~ a + evolved
{
  mr0.03 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bevolved *evolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.03,digits=4, depth=2)
save(mr0.03, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.03.RData")

#1.04) r0 ~ a + lowevolved
{
  mr0.04 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + blowevolved *lowevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.04,digits=4, depth=2)
save(mr0.04, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.04.RData")

#1.05) r0 ~ a + highevolved
{
  mr0.05 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.05,digits=4, depth=2)
save(mr0.05, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.05.RData")

#1.06) r0 ~ a + lowevolved + highevolved
{
  mr0.06 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved + blowevolved *lowevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.06,digits=4, depth=2)
save(mr0.06, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.06.RData")

#1.07) r0 ~ a + medium + evolved
{
  mr0.07 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.07,digits=4, depth=2)
save(mr0.07, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.07.RData")

#1.08) r0 ~ a + medium + evolved + medium*evolved
{
  mr0.08 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved + bint1*pHmedium*evolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.08,digits=4, depth=2)
save(mr0.08, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.08.RData")

#1.09) r0 ~ a + medium + lowevolved
{
  mr0.09 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.09,digits=4, depth=2)
save(mr0.09, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.09.RData")

#1.10) r0 ~ a + medium + lowevolved + medium*lowevolved
{
  mr0.10 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved + bint1*pHmedium*lowevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.10,digits=4, depth=2)
save(mr0.10, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.10.RData")

#1.11) r0 ~ a + medium + highevolved
{
  mr0.11 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.11,digits=4, depth=2)
save(mr0.11, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.11.RData")

#1.12) r0 ~ a + medium + highevolved + medium*highevolved
{
  mr0.12 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + bint1*pHmedium*highevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.12,digits=4, depth=2)
save(mr0.12, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.12.RData")

#1.13) r0 ~ a + medium + highevolved + lowevolved
{
  mr0.13 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.13,digits=4, depth=2)
save(mr0.13, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.13.RData")

#1.14) r0 ~ a + medium + highevolved + lowevolved + lowevolved*medium
{
  mr0.14 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*lowevolved*pHmedium,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.14,digits=4, depth=2)
save(mr0.14, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.14.RData")

#1.15) r0 ~ a + medium + highevolved + lowevolved + highevolved*medium
{
  mr0.15 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*highevolved*pHmedium,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.15,digits=4, depth=2)
save(mr0.15, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.15.RData")

#1.16) r0 ~ a + medium + highevolved + lowevolved + highevolved*medium + lowevolved*medium
{
  mr0.16 <- map2stan(
    alist(
      log_r0_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint2*highevolved*pHmedium + bint1*lowevolved*pHmedium,
      log_r0_obs ~ dnorm(log_r0_est,log_r0_sd),
      a ~dnorm(-2.4, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bint2  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_r0_est=m1data$log_r0_obs) ,
   iter=2e4, warmup=4e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mr0.16,digits=4, depth=2)
save(mr0.16, file="2_data/14,16,18_EvolutionModelsBayesianContinued/r0models/CU427.4/mr0.16.RData")