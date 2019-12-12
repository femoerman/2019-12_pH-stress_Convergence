# clear memory
rm(list=ls())

#Set working directory
# setwd("F:/Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#load packages
library(tidyverse)
library(rethinking)
library(nlme)
library(lme4)
library(MuMIn)
library(rstan)
library(ellipse)
library(mvtnorm)

#Load data and prepare for usage (filter out before common garden data; rename pH variables as numerical; make factors binary (needed for Bayesian model fitting))
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
data.all$Genotype <- as.factor(ifelse(data.all$strain=="B2086.2", "Genotype 1", ifelse(data.all$strain=="CU427.4", "Genotype 2", ifelse(data.all$strain == "CU428.2", "Genotype 3", "Genotype 4"))))
data.all <- arrange(data.all, strain)

data.all$pHmediumfact <- ifelse(data.all$pHmedium == 6.5, 1, 0)
data.all$r0mean <- exp(data.all$log_r0mean)
data.all$pHmname <- ifelse(data.all$pHmedium == 6.5 , "pH 6.5", "pH 4.5")

library(scales)
#Plot raw correlation data
fig1 <- ggplot(filter(data.all, pHmedium==4.5, (pHorigin != 6.5 | timing == "ancestral")), aes(y = log_r0mean, x = log_alfamean, shape=Genotype, colour = origin)) + geom_point(size=3) + theme_light() +
  stat_ellipse(inherit.aes = F, data = filter(data.all, pHmedium==4.5), mapping = aes(x = log_alfamean, y = log_r0mean, group=pHmname), type="norm") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) )  + xlab("α (mL/(h indiv)) (ln)") + ylab(expression("r"[0]*" (1/h) (ln)")) + 
  scale_color_manual(values=c("#2db3f9", "#e41a1c"), breaks=unique(data.all$origin)[2:3], name="Origin", labels=c("LpH", "ANC")) +
  theme(strip.background =element_rect(fill="white"))+ ylim(-3.3, -0.6) + xlim(-16.5, -13.3) +
  theme(strip.text = element_text(colour = 'black', size=16))

fig2 <- ggplot(filter(data.all, pHmedium==6.5, (pHorigin != 4.5 | timing == "ancestral")), aes(y = log_r0mean, x = log_alfamean, shape=Genotype, colour = origin)) + geom_point(size=3) + theme_light() +
  stat_ellipse(inherit.aes = F, data = filter(data.all, pHmedium==6.5), mapping = aes(x = log_alfamean, y = log_r0mean, group=pHmname), type="norm") + 
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) )  + xlab("α (mL/(h indiv)) (ln)") + ylab(expression("r"[0]*" (1/h) (ln)")) + 
  scale_color_manual(values=c("#2db3f9", "#2d33f9"), breaks=unique(data.all$origin)[c(1, 3)], name="Origin", labels=c("NpH", "ANC")) +
  theme(strip.background =element_rect(fill="white"))+ ylim(-3.3, -0.6) + xlim(-16.5, -13.3) +
  theme(strip.text = element_text(colour = 'black', size=16))

library(ggpubr)
ggarrange(fig1, fig2,  nrow=1, ncol=2, align = "hv", common.legend = T, legend =  "right",
          labels =list("A", "B"))
ggsave(filename = "Fig3.png", path = "4_results/Figures", width = 255, height = 204, units = "mm", device = "png")

#Separate data for model fitting (by pH of the assay medium)
data.low <- filter(data.all, pHmedium==4.5, (pHorigin != 4.5 | timing == "ancestral"))
data.neutral <- filter(data.all, pHmedium==6.5, (pHorigin != 4.5 | timing == "ancestral"))

#Assign the r0 and alpha variables to the new data set for low pH of the assay medium
xlow <- matrix(nrow=nrow(data.low), ncol=2)
xlow[, 1] <- data.low$log_r0mean
xlow[, 2] <- data.low$log_alfamean
xsdlow <- matrix(nrow=nrow(data.low), ncol=2)
xsdlow[, 1] <- data.low$log_r0sd
xsdlow[, 2] <- data.low$log_alfasd
nlow <- nrow(data.low)

#Assign the r0 and alpha variables to the new data set for neutral pH of the assay medium
xneutral <- matrix(nrow=nrow(data.neutral), ncol=2)
xneutral[, 1] <- data.neutral$log_r0mean
xneutral[, 2] <- data.neutral$log_alfamean
xsdneutral <- matrix(nrow=nrow(data.neutral), ncol=2)
xsdneutral[, 1] <- data.neutral$log_r0sd
xsdneutral[, 2] <- data.neutral$log_alfasd
nneutral <- nrow(data.neutral)


#Stan code for the correlation model
{
  #### Notes to Stan model #######################################################
  ## 1) Multivariate normal distribution in Stan uses covariance matrix instead of 
  ##    precision matrix.
  ## 2) Multivariate normal distribution can be (and is) also vectorized.
  ## 3) Warnings may occur during sampling, ignore them.
  ################################################################################
  model <- "
  // Pearson Correlation
  data { 
  int<lower=1> n;
  vector[2] x_obs[n];
  vector[2] x_sd[n];
  }
  parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1,upper=1> r;
  vector[2] x_est[n];
  
  } 
  transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;
  
  // Reparameterization
  sigma[1] = sqrt(lambda[1]);
  sigma[2] = sqrt(lambda[2]);
  T[1,1] = square(sigma[1]);
  T[1,2] = r * sigma[1] * sigma[2];
  T[2,1] = r * sigma[1] * sigma[2];
  T[2,2] = square(sigma[2]);
  }
  model {
  // Priors
  mu ~ normal(0, 10);
  lambda ~ normal(0,1);
  r ~ normal(0,1);
  
  //x_obs[1] ~ normal(x_est[1], x_sd[1]);
  //x_obs[2] ~ normal(x_est[2], x_sd[2]);
  
  // Data
  x_est ~ multi_normal(mu, T);
  
  for (i in 1:n){
  x_obs[i][1] ~ normal(x_est[i][1], x_sd[i][1]);
  x_obs[i][2] ~ normal(x_est[i][2], x_sd[i][2]);
  }
  
  }"
  }

#Fit Bayesian correlation model for r0 and alpha
{
  #Fit model at low pH
  {
  # complex model
  datalow <- list(x_obs=xlow, x_sd=xsdlow, n=nlow) # to be passed on to Stan
  myinits <- list(list(r=0.75, mu=c(-2, -15), lambda=c(0.5, 0.5), x_est = xlow))
  
  # parameters to be monitored: 
  parameters <- c("r", "mu", "sigma")
  
  # The following command calls Stan with specific options.
  # For a detailed description type "?rstan".
  sampleslow <- stan(model_code=model,   
                  data=datalow, 
                  init=myinits,  # If not specified, gives random inits
                  pars=parameters,
                  iter=10000, 
                  warmup=2000,
                  chains=1,
                  control = list(adapt_delta = 0.99)
                  # warmup = 100,  # Stands for burn-in; Default = iter/2
                  # seed = 123  # Setting seed; Default is random seed
  )
  r0alfphalow <- sampleslow
  extra
  rlow <- extract.samples(sampleslow)$r
  median(rlow)
  
  act_covlow <- extract.samples(sampleslow)$r * extract.samples(sampleslow)$sigma[,1] * extract.samples(sampleslow)$sigma[,2]
  median(act_covlow)
  
  
  }
  
  #Fit model at neutral pH
  {
    # complex model
    dataneutral <- list(x_obs=xneutral, x_sd=xsdneutral, n=nneutral) # to be passed on to Stan
    myinits <- list(list(r=0.75, mu=c(-1, -14), lambda=c(0.5, 0.5), x_est = xneutral))
    
    # parameters to be monitored: 
    parameters <- c("r", "mu", "sigma")
    
    # The following command calls Stan with specific options.
    # For a detailed description type "?rstan".
    samplesneutral <- stan(model_code=model,   
                       data=dataneutral, 
                       init=myinits,  # If not specified, gives random inits
                       pars=parameters,
                       iter=10000, 
                       warmup=2000,
                       chains=1,
                       control = list(adapt_delta = 0.99)
                       # warmup = 100,  # Stands for burn-in; Default = iter/2
                       # seed = 123  # Setting seed; Default is random seed
    )
    r0alphaneutral <- samplesneutral
    
    rneutral <- extract.samples(samplesneutral)$r
    median(rneutral)
    
    act_covneutral <- extract.samples(samplesneutral)$r * extract.samples(samplesneutral)$sigma[,1] * extract.samples(samplesneutral)$sigma[,2]
    median(act_covneutral)
  }
}

