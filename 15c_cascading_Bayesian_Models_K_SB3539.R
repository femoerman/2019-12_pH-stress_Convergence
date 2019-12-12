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
sumdata <- filter(sumdata, strain == "SB3539")
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

#1) Fit cascading Bayesian model for K
#1.01) K ~ a
{
  mK.01 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a ,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.01,digits=4, depth=2)
save(mK.01, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.01.RData")


#1.02) K ~ a + pHmedium
{
  mK.02 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium *pHmedium,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.02,digits=4, depth=2)
save(mK.02, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.02.RData")


#1.03) K ~ a + evolved
{
  mK.03 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bevolved *evolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.03,digits=4, depth=2)
save(mK.03, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.03.RData")

#1.04) K ~ a + lowevolved
{
  mK.04 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + blowevolved *lowevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.04,digits=4, depth=2)
save(mK.04, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.04.RData")

#1.05) K ~ a + highevolved
{
  mK.05 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.05,digits=4, depth=2)
save(mK.05, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.05.RData")

#1.06) K ~ a + lowevolved + highevolved
{
  mK.06 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bhighevolved *highevolved + blowevolved *lowevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.06,digits=4, depth=2)
save(mK.06, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.06.RData")

#1.07) K ~ a + medium + evolved
{
  mK.07 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.07, depth=2)
save(mK.07, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.07.RData")

#1.08) K ~ a + medium + evolved + medium*evolved
{
  mK.08 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bevolved *evolved + bint1*pHmedium*evolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.08,digits=4, depth=2)
save(mK.08, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.08.RData")

#1.09) K ~ a + medium + lowevolved
{
  mK.09 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.09,digits=4, depth=2)
save(mK.09, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.09.RData")

#1.10) K ~ a + medium + lowevolved + medium*lowevolved
{
  mK.10 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + blowevolved *lowevolved + bint1*pHmedium*lowevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.10,digits=4, depth=2)
save(mK.10, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.10.RData")

#1.11) K ~ a + medium + highevolved
{
  mK.11 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.11,digits=4, depth=2)
save(mK.11, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.11.RData")

#1.12) K ~ a + medium + highevolved + medium*highevolved
{
  mK.12 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + bint1*pHmedium*highevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.12,digits=4, depth=2)
save(mK.12, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.12.RData")

#1.13) K ~ a + medium + highevolved + lowevolved
{
  mK.13 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.13,digits=4, depth=2)
save(mK.13, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.13.RData")

#1.14) K ~ a + medium + highevolved + lowevolved + lowevolved*medium
{
  mK.14 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*lowevolved*pHmedium,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.14,digits=4, depth=2)
save(mK.14, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.14.RData")

#1.15) K ~ a + medium + highevolved + lowevolved + highevolved*medium
{
  mK.15 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint1*highevolved*pHmedium,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.15,digits=4, depth=2)
save(mK.15, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.15.RData")

#1.16) K ~ a + medium + highevolved + lowevolved + highevolved*medium + lowevolved*medium
{
  mK.16 <- map2stan(
    alist(
      log_K_est ~ dnorm(mu,sigma),
      mu <- a + bmedium*pHmedium + bhighevolved *highevolved + blowevolved*lowevolved + 
        bint2*highevolved*pHmedium + bint1*lowevolved*pHmedium,
      log_K_obs ~ dnorm(log_K_est,log_K_sd),
      a ~dnorm(13, 1),
      bhighevolved  ~ dnorm(0, 1),
      blowevolved  ~ dnorm(0, 1),
      bint1  ~ dnorm(0, 1),
      bint2  ~ dnorm(0, 1),
      bmedium  ~ dnorm(0, 1),
      sigma ~dcauchy(0, 1)
    ) ,data=m1data ,
    start=list(log_K_est=m1data$log_K_obs) ,
    iter=2e5, warmup=4e4,
    control=list(adapt_delta=0.95), WAIC=FALSE )
}
precis(mK.16,digits=4, depth=2)
save(mK.16, file="2_data/14,16,18_EvolutionModelsBayesianContinued/Kmodels/SB3539/mK.16.RData")
