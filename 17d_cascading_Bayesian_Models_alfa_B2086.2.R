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
sumdata <- filter(sumdata, strain == "B2086.2")
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

#1) Fit cascading Bayesian model for α
#1.01) α ~ a
{
  malfa.01 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a ,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.01,digits=4, depth=2)
save(malfa.01, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.01.RData")


#1.02) α ~ a + pHmedium
{
  malfa.02 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium *pHmedium,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.02,digits=4, depth=2)
save(malfa.02, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.02.RData")


#1.03) α ~ a + evolved
{
  malfa.03 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bevolved *evolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.03,digits=4, depth=2)
save(malfa.03, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.03.RData")

#1.04) α ~ a + lowevolved
{
  malfa.04 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + blowevolved *lowevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.04,digits=4, depth=2)
save(malfa.04, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.04.RData")

#1.05) α ~ a + highevolved
{
  malfa.05 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.05,digits=4, depth=2)
save(malfa.05, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.05.RData")

#1.06) α ~ a + lowevolved + highevolved
{
  malfa.06 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved + blowevolved *lowevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.06,digits=4, depth=2)
save(malfa.06, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.06.RData")

#1.07) α ~ a + medium + evolved
{
  malfa.07 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.07,digits=4, depth=2)
save(malfa.07, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.07.RData")

#1.08) α ~ a + medium + evolved + medium*evolved
{
  malfa.08 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved + bint1*pHmedium*evolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.08,digits=4, depth=2)
save(malfa.08, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.08.RData")

#1.09) α ~ a + medium + lowevolved
{
  malfa.09 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.09,digits=4, depth=2)
save(malfa.09, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.09.RData")

#1.10) α ~ a + medium + lowevolved + medium*lowevolved
{
  malfa.10 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved + bint1*pHmedium*lowevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.10,digits=4, depth=2)
save(malfa.10, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.10.RData")

#1.11) α ~ a + medium + highevolved
{
  malfa.11 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.11,digits=4, depth=2)
save(malfa.11, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.11.RData")

#1.12) α ~ a + medium + highevolved + medium*highevolved
{
  malfa.12 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + bint1*pHmedium*highevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.12,digits=4, depth=2)
save(malfa.12, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.12.RData")

#1.13) α ~ a + medium + highevolved + lowevolved
{
  malfa.13 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.13,digits=4, depth=2)
save(malfa.13, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.13.RData")

#1.14) α ~ a + medium + highevolved + lowevolved + lowevolved*medium
{
  malfa.14 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*lowevolved*pHmedium,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.14,digits=4, depth=2)
save(malfa.14, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.14.RData")

#1.15) α ~ a + medium + highevolved + lowevolved + highevolved*medium
{
  malfa.15 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*highevolved*pHmedium,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.15,digits=4, depth=2)
save(malfa.15, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.15.RData")

#1.16) α ~ a + medium + highevolved + lowevolved + highevolved*medium + lowevolved*medium
{
  malfa.16 <- map2stan(
    alist(
      log_alfa_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint2*highevolved*pHmedium + bint1*lowevolved*pHmedium,
      log_alfa_obs ~ dnorm(log_alfa_est,log_alfa_sd),
      a ~dnorm(-14.75, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bint2  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_alfa_est=m1data$log_alfa_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(malfa.16,digits=4, depth=2)
save(malfa.16, file="2_data/14,16,18_EvolutionModelsBayesianContinued/alfamodels/B2086.2/malfa.16.RData")

