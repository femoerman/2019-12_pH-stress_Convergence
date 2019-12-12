#1) Clear memory
rm(list=ls())

#2) Set working directory
#setwd("F:/Documenten/PhD/05_student_project_2017")
setwd("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017")

#3) Load datasets and packages
load("2_data/posteriorBH_means.RData")
load("2_data/DensitiesDuringEvolution.RData")
library(lme4)
library(MuMIn)
library(rethinking)
library(scales)
library("tidyverse")
library("ggpubr")
library("gridExtra")

#3.2) Remove extinct populations from analysis
evolve$survive <- ifelse(evolve$file %in% c("sample_00015", "sample_00005", "sample_00001"), " -extinct", " -surviving")
evolve <- filter(evolve, survive==" -surviving")
sumdata$dmean <- exp(sumdata$log_dmean)

#4) Plot histogram for densities over evolution experiment
########################################################
hist(evolve$indiv_per_ml[which(evolve$DaysSinceStart>max(evolve$DaysSinceStart)/2)])
hist(evolve$indiv_per_ml,breaks=1000)$mids


#5) Integrate fitness weighted over the densities experienced during evolution
#########################################################
# dens reg functions
drf <- function(N) (r0+d)/(1+(r0/(K*d))*N) -d

par(mfrow=c(2,2))

fit_sums <- numeric()
fit_sums_non_weighted <-numeric()
fit_sums_cat <- numeric()
fit_sums_strain <- numeric()
fit_sums_r0 <- numeric()

j <- 0

#Loop over evolved populations and ancestors, and do the integration of fitness (at neutral pH for neutral-evolved and ancestors)
for(act_strain in sort(unique(sumdata$strain))){
  j = j + 1
  act_subset <- which(sumdata$timing=="ancestral" & sumdata$pHmedium==6.5 & sumdata$strain==act_strain)
  
  r0 <- exp(sumdata$log_r0mean[act_subset[1]])
  K <- exp(sumdata$log_Kmean[act_subset[1]])
  d <- sumdata$dmean[act_subset[1]]
  
  par(bty="l")
  curve(drf,0,1000000,ylim=c(-0.1,max(exp(sumdata$log_r0mean))),col="lightblue",lwd=2)
  abline(h=0,lty=3)
  legend("topright",legend=act_strain,bty="n")
  fit_sums[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")],
                                 breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")]),
                                            len=1000000),plot=F)$mids,drf)*hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")],seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")]),len=1000000),plot=F)$counts)
  
  fit_sums_non_weighted[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")],breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="low")]),len=1000000),plot=F)$mids,drf))
  #fit_sums_non_weighted[j] <-sum(sapply(seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain)]),len=10000),drf))
  fit_sums_cat[j] <- "anc"
  fit_sums_strain[j] <- act_strain
  fit_sums_r0[j] <- r0
  
  #Do so for the evolved lines
  for(i in 2:length(act_subset)){
    r0 <- exp(sumdata$log_r0mean[act_subset[i]])
    K <- exp(sumdata$log_Kmean[act_subset[i]])
    d <- sumdata$dmean[act_subset[i]]
    
    curve(drf,0,1000000,col="lightblue",lwd=2,add=T)
    j <- j+1
    fit_sums[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],
                                   breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),
                                              len=1000000),plot=F)$mids,drf)*hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),len=1000000),plot=F)$counts)
    
    fit_sums_non_weighted[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),len=1000000),plot=F)$mids,drf))
    #fit_sums_non_weighted[j] <-sum(sapply(seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain)]),len=10000),drf))
    fit_sums_cat[j] <- "anc"
    fit_sums_strain[j] <- act_strain
    fit_sums_r0[j] <- r0
  }
  
  act_subset <- which(sumdata$timing=="cg" & sumdata$pHorigin==6.5 & sumdata$pHmedium==6.5 & sumdata$strain==act_strain)
  
  #Do so for the ancestor lines
  for(i in 1:length(act_subset)){
    r0 <- exp(sumdata$log_r0mean[act_subset[i]])
    K <- exp(sumdata$log_Kmean[act_subset[i]])
    d <- sumdata$dmean[act_subset[i]]
    
    curve(drf,0,1000000,col="darkblue",lwd=2,add=T)
    j <- j +1
    fit_sums[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),len=1000000),plot=F)$mids,drf)*hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),len=1000000),plot=F)$counts)
    fit_sums_non_weighted[j] <- sum(sapply(hist(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")],breaks=seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain & evolve$pH=="high")]),len=1000000),plot=F)$mids,drf))
    #fit_sums_non_weighted[j] <- sum(sapply(seq(0,max(evolve$indiv_per_ml[which(evolve$strain==act_strain)]),len=10000),drf))
    fit_sums_cat[j] <- "evol_npH"
    fit_sums_strain[j] <- act_strain
    fit_sums_r0[j] <- r0
  }
  
}

######

#6) Plot integrated fitness estimates, and fit lmer
########################################################
# test plot for fitness vs r0
fit_sums_r0_centered <- fit_sums_r0 - mean(fit_sums_r0)
plot(fit_sums ~ fit_sums_r0_centered, type="n")
points(fit_sums[which(fit_sums_cat=="anc")]~fit_sums_r0_centered[which(fit_sums_cat=="anc")], pch=21, col="white",bg="lightblue",cex=1.5)
points(fit_sums[which(fit_sums_cat=="evol_npH")]~fit_sums_r0_centered[which(fit_sums_cat=="evol_npH")], pch=21, col="white",bg="darkblue",cex=1.5)

#Fit lmer and compare
AICc(lmer(fit_sums ~ 1 + (1|fit_sums_strain)))
AICc(lmer(fit_sums ~ fit_sums_r0_centered + (1|fit_sums_strain)))
AICc(lmer(fit_sums ~ fit_sums_cat + (1|fit_sums_strain)))
AICc(lmer(fit_sums ~ fit_sums_cat+fit_sums_r0_centered + (1|fit_sums_strain)))
AICc(lmer(fit_sums ~ fit_sums_cat*fit_sums_r0_centered + (1|fit_sums_strain)))
summary(lmer(fit_sums ~ fit_sums_cat*fit_sums_r0_centered + (1|fit_sums_strain)))
m1 <- lmer(fit_sums ~ fit_sums_cat*fit_sums_r0_centered + (1|fit_sums_strain))
abline(fixef(m1)[[1]],fixef(m1)[[3]],col="blue")
abline(fixef(m1)[[1]]+fixef(m1)[[2]],fixef(m1)[[3]]+fixef(m1)[[4]],col="darkblue")

######7) Plot integrated fitness estimates, and fit Bayesian mixed models
{
  #7.0) Prepare data as list
  {
   dataBayes <- list(centeredr0 = fit_sums_r0_centered, fits = fit_sums, strains = fit_sums_strain, category = fit_sums_cat, evolved = ifelse(fit_sums_cat=="anc", 0, 1))
  }
  
  
  #7.1) Fit intercept model
  mB1 <- map2stan(
    alist(
      fits ~ dnorm(mu,sigma),
      mu <- a + agen[strains],
      a ~dnorm(3.8, 2),
      agen[strains] ~ dnorm(0, 2),
      sigma ~ dcauchy(0, 1)
    ) ,data=dataBayes,
    iter=4e4, warmup=8e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  precis(mB1, depth=2)
  
  #7.2 Fit model with r0-centered
  mB2 <- map2stan(
    alist(
      fits ~ dnorm(mu,sigma),
      mu <- a + agen[strains] + b1 * centeredr0,
      a ~dnorm(2.3, 2),
      agen[strains] ~ dnorm(0, 2),
      b1 ~ dnorm(23, 10),
      sigma ~ dcauchy(0, 1)
    ) ,data=dataBayes,
    iter=4e4, warmup=8e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  precis(mB2, depth=2)
  
  #7.3 Fit model with category (evolved vs ancestral)
  mB3 <- map2stan(
    alist(
      fits ~ dnorm(mu,sigma),
      mu <- a + agen[strains] + b2 * evolved,
      a ~dnorm(2.3, 2),
      agen[strains] ~ dnorm(0, 2),
      b2 ~ dnorm(0, 10),
      sigma ~ dcauchy(0, 1)
    ) ,data=dataBayes,
    iter=4e4, warmup=8e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  precis(mB3, depth=2)
  
  #7.4) Fit additive model with category and r0-centered
  mB4 <- map2stan(
    alist(
      fits ~ dnorm(mu,sigma),
      mu <- a + agen[strains] + b2*evolved + b1*centeredr0,
      a ~dnorm(2.3, 2),
      agen[strains] ~ dnorm(0, 2),
      b2 ~ dnorm(0, 10),
      b1 ~ dnorm(23, 10),
      sigma ~ dcauchy(0, 1)
    ) ,data=dataBayes,
    iter=4e4, warmup=8e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  precis(mB4, depth=2)
  
  #7.5) Fit multiplicative model with category and r0-centered
  mB5 <- map2stan(
    alist(
      fits ~ dnorm(mu,sigma),
      mu <- a + agen[strains] + b2*evolved + b1*centeredr0 + i1*evolved*centeredr0,
      a ~dnorm(2.3, 2),
      agen[strains] ~ dnorm(0, 2),
      b1 ~ dnorm(23, 10),
      b2 ~ dnorm(0, 10),
      i1 ~ dnorm(-30, 10),
      sigma ~ dcauchy(0, 1)
    ) ,data=dataBayes,
    iter=4e4, warmup=8e3,
    control=list(adapt_delta=0.95), WAIC=FALSE )
  precis(mB5, depth=2, prob = 0.95)
  
  #7.6) Compare the models 
  t <- compare(mB1, mB2, mB3, mB4, mB5)
}

######8) Average predictions
{
  #Create data for prediction
  dataBayes.df <- as.data.frame(dataBayes)
  data.low<- expand.grid(evolved=1, 
                         centeredr0=seq(min(filter(dataBayes.df, evolved==1)$centeredr0), max(filter(dataBayes.df, evolved==1)$centeredr0), length.out = 200))
  data.anc<- expand.grid(evolved=0, 
                         centeredr0=seq(min(filter(dataBayes.df, evolved==0)$centeredr0), max(filter(dataBayes.df, evolved==0)$centeredr0), length.out = 200))
  
  
  #Create averaged predictions (disregarding strain effects)
  {
    t
    #Define weights
    weight5 <- 1
    
    #Define sample lists for drawing parameters
    data5 <- extract.samples(mB5)
    
    #Create function to draw for specific values of centeredr0 and evolved
    weightedpred <- function(i, data, weight5, n=32000){

      centeredr0 <- data[i, "centeredr0"]
      evolved <- data[i, "evolved"]
      pred5 <- data5$a + data5$b1*centeredr0 + data5$b2*evolved + data5$i1 *evolved*centeredr0
      
      weighted.data <- c(sample(pred5, size = as.integer(weight5*n)))
      weightedmean <- mean(weighted.data)
      weightedsd <- sd(weighted.data)
      weightedupper <- weightedmean + 1.96*weightedsd
      weightedlower <- weightedmean - 1.96*weightedsd
      
      return(c(weightedmean, weightedupper, weightedlower))
    }
  }
  
  #Generate predictions for low-evolved and ancestor data
  data.low.pred <- matrix(unlist(lapply(1:nrow(data.low), weightedpred, data.low, weight5)), ncol=3, byrow=T)
  data.anc.pred <- matrix(unlist(lapply(1:nrow(data.anc), weightedpred, data.anc, weight5)), ncol=3, byrow=T)
  data.low$mean <- data.low.pred[, 1]
  data.low$upper <- data.low.pred[, 2]
  data.low$lower <- data.low.pred[, 3]
  
  data.anc$mean <- data.anc.pred[, 1]
  data.anc$upper <- data.anc.pred[, 2]
  data.anc$lower <- data.anc.pred[, 3]
  
  #Concatenate predictions together
  data <- rbind(data.low, data.anc)
  data$category <- ifelse(data$evolved==1, "evol_npH", "anc")
  
  #Create figure showing data + predictions
  ggplot(data, aes(x= centeredr0, y = mean, colour=evolved)) + geom_point()
  category <- unique(data$category)
  dataBayes$Strain <- ifelse(dataBayes$strain=="B2086.2", "Genotype 1", ifelse(dataBayes$strain=="CU427.4", "Genotype 2", 
                                                                             ifelse(dataBayes$strain=="CU428.2", "Genotype 3", "Genotype 4")))
  
  #Prepare data for plotting points that have been used in other figures
  temp <- as.data.frame(dataBayes)
  data.selected <- temp[c(1, 5, 21, 23), ]
  
  figC <- ggplot(as.data.frame(dataBayes), aes(y=fits,  x=centeredr0, colour=category, shape=Strain)) + 
    #geom_point(data = data.selected, inherit.aes=F, mapping = aes(y = fits, x = centeredr0), colour="darkgrey", size=7) +
    geom_point(size=4) + 
    geom_line(data=data, inherit.aes = F, aes(x=centeredr0, y=mean, colour=category, group=category), size=2) +
    geom_ribbon(data=data, inherit.aes = F, aes(ymax=upper, ymin=lower, x=centeredr0, fill=category, group=category),alpha=0.3) +
    theme_light() + xlab(expression("Centered r"[0]*" (1/h)")) + ylab("Density-dependent fitness") + 
    theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
           axis.text.y=element_text(size=16), axis.text.x=element_text(size=16),
           plot.title = element_text(size=20, hjust = 0.5)) + 
    scale_color_manual(values=c("#2db3f9", "#2d33f9"), breaks=category, name="Origin", labels=c( "NpH", "ANC")) + 
    scale_fill_manual(values=c("#2db3f9", "#2d33f9"), breaks=category, name="Origin", labels=c("NpH", "ANC")) +
    scale_shape_discrete(name = "Genotype")

}

######9) Calculate density trgulation functions and plot
{
  #Select good example and calculate drf
  {
    sumdata$strainval <- sumdata$strain
    sumdata$evolved <- ifelse(sumdata$timing=="cg", 1, 0)
    sumdata$pHmedium <- ifelse(sumdata$pHmedium==6.5, 6.5, 4.5)
    sumdata$pHorigin <- ifelse(sumdata$pHorigin==6.5, 6.5, 4.5)
    sumdata$highevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==6.5, 1, 0)
    sumdata$lowevolved <- ifelse(sumdata$timing=="cg"&sumdata$pHorigin==4.5, 1, 0)
    sumdata$pHmedium <- as.factor(sumdata$pHmedium)
    data.all <- sumdata %>% filter(timing != "nocg")
    data.all <- mutate(data.all, origin=ifelse(timing!="cg", "ancestor", ifelse(pHorigin==6.5, "high-evolved", "low-evolved")))
    data.all$strain <- as.factor(ifelse(data.all$strain=="B2086.2", "Genotype 1", ifelse(data.all$strain=="CU427.4", "Genotype 2", ifelse(data.all$strain == "CU428.2", "Genotype 3", "Genotype 4"))))
    data.all <- arrange(data.all, strain)
    data.all <- data.all[which(data.all$name %in% c("ancestral 6.5 5 CU428.2 ph_anc_1", "cg 6.5 CU428.2 ph_high_2",
                                                    "cg 5 5 SB3539 ph_low_2", "ancestral 6.5 5 SB3539 ph_anc_3")), ]
    data.all$log_dmean <- log(data.all$dmean)
    
    {
      data.all <- mutate(data.all, fitness = ((exp(log_r0mean) + exp(log_dmean)))/(1+(exp(log_r0mean)/(exp(log_Kmean)*exp(log_dmean)))) - exp(log_dmean))
      data.all <- mutate(data.all, fitnesshalfK = ((exp(log_r0mean) + exp(log_dmean))/(1 + ((exp(log_r0mean)/(exp(log_Kmean)*exp(log_dmean))) * exp(log_Kmean)/999)) - exp(log_dmean))*exp(log_Kmean)/2)
    }
    
    predict.data <- expand.grid(name = unique(data.all$name), percK = 1:120)
    predict.data$fitness = NA
    predict.data$strain <- NA
    predict.data$pHmedium <- NA
    predict.data$origin <- NA
    predict.data$dens <- NA
    for (i in unique(data.all$name)){
      tempdata <- filter(data.all, name==i)
      for (j in 1:120){
        predict.data[which(predict.data$name==i & predict.data$percK==j), "fitness"] <- ((exp(tempdata$log_r0mean) + exp(tempdata$log_dmean))/(1 + ((exp(tempdata$log_r0mean)/(exp(tempdata$log_Kmean)*exp(tempdata$log_dmean))) * (j/100*exp(tempdata$log_Kmean)))) - exp(tempdata$log_dmean))
        predict.data[which(predict.data$name==i & predict.data$percK==j), "dens"] <- (j/100)*exp(tempdata$log_Kmean)
        predict.data[which(predict.data$name==i & predict.data$percK==j), "strain"] <- as.character(tempdata$strain)
        predict.data[which(predict.data$name==i & predict.data$percK==j), "pHmedium"] <- as.character(tempdata$pHmedium)
        predict.data[which(predict.data$name==i & predict.data$percK==j), "origin"] <- tempdata$origin
      }
    }
  }
  
  #Get density histogram for chosen example
  data.ex1 <- filter(evolve, strain=="B2086.2" & pH=="high")
  data.ex2 <- filter(evolve, strain=="SB3539" & pH=="high")
  case1 <- predict.data[which(predict.data$name %in% c("ancestral 6.5 6.5 SB3539 ph_anc_1", "cg 6.5 6.5 SB3539 ph_low_2")), ]
  case2 <- predict.data[which(predict.data$name %in% c("cg 6.5 6.5 SB3539 ph_low_2", "ancestral 6.5 6.5 SB3539 ph_anc_3")), ]
  

  #Plot the figure for the low r0 case
  figA <- ggplot(case1, aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) + 
    geom_hline(yintercept=0, size=1.5) +
    geom_histogram(data = data.ex1, inherit.aes = F, mapping = aes(x = indiv_per_ml, y=50000*..density..), binwidth = 5e4, alpha = 0.4) +
    theme_light()  + xlab("Population density, N (indiv/mL)") + ylab("Pop. growth rate, r (indiv/(mL h))") +
    theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
           axis.text.y=element_text(size=16), axis.text.x=element_text(size=16),
           plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific, limits = c(0, 1.5e6)) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c"), breaks=unique(predict.data$origin), name="Origin", labels=c("LpH", "ANC"))
  
  
  #Plot the figure fr the high r0 case
  figB <- ggplot(case2, aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) + 
    geom_hline(yintercept=0, size=1.5) +
    geom_histogram(data = data.ex2, inherit.aes = F, mapping = aes(x = indiv_per_ml, y=50000*..density..), binwidth = 5e4, alpha = 0.4) +
    theme_light()  + xlab("Population density, N (indiv/mL)") + ylab("Pop. growth rate, r (indiv/(mL h))") + 
    theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
           axis.text.y=element_text(size=16), axis.text.x=element_text(size=16),
           plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific, limits = c(0, 1.5e6)) +
    scale_color_manual(values=c("#2db3f9", "#e41a1c"), breaks=unique(predict.data$origin), name="Origin", labels=c("LpH", "ANC"))
  
}

#10) Arrange figures
{
  figAB <- ggarrange(figA, figB, ncol = 1, nrow=2, common.legend = T, legend = "none",  labels =list("A", "B"), align = "hv")
  
  fig1 <- ggarrange(figAB, figC, nrow=1, ncol=2, widths = c(0.6, 1), labels = list("", "C"))
  ggsave(plot = figC,filename = "Fig4-neut.png", path = "4_results/Figures", width = 350, height = 250, units = "mm", device = "png")
}
