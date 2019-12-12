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
data.all$strain <- as.factor(ifelse(data.all$strain=="B2086.2", "Strain 1", ifelse(data.all$strain=="CU427.4", "Strain 2", ifelse(data.all$strain == "CU428.2", "Strain 3", "Strain 4"))))
data.all <- arrange(data.all, strain)


#Try fitness function idea / density regulation function
{
    data.all <- mutate(data.all, fitness = ((exp(log_r0mean) + exp(log_dmean)))/(1+(exp(log_r0mean)/(exp(log_Kmean)*exp(log_dmean)))) - exp(log_dmean))
    data.all <- mutate(data.all, fitnesshalfK = ((exp(log_r0mean) + exp(log_dmean))/(1 + ((exp(log_r0mean)/(exp(log_Kmean)*exp(log_dmean))) * exp(log_Kmean)/999)) - exp(log_dmean))*exp(log_Kmean)/2)
}

ggplot(data.all, aes(x = paste(origin, pHmedium, strain), y = fitnesshalfK, colour = pHmedium)) + geom_boxplot()
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
predict.data$pHmedium <- ifelse(predict.data$pHmedium==6.5, "pH 6.5", "pH 4.5")
library(scales)
data.lowonly <- filter(predict.data, origin != "high-evolved" & pHmedium == 1)
fig1 <- ggplot(filter(predict.data, pHmedium=="pH 4.5" & strain == "Strain 1"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light() +  ylab("Population growth rate r (indiv/(mL hours)")  +ylim(-0.05, 0.4) + ggtitle("Genotype 1") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig2 <- ggplot(filter(predict.data, pHmedium=="pH 4.5" & strain == "Strain 2"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light() + ggtitle("Genotype 2") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific)  +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig3 <- ggplot(filter(predict.data, pHmedium=="pH 4.5" & strain == "Strain 3"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light()  + ggtitle("Genotype 3") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig4 <- ggplot(filter(predict.data, pHmedium=="pH 4.5" & strain == "Strain 4"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light() + ggtitle("Genotype 4") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig5 <- ggplot(filter(predict.data, pHmedium=="pH 6.5" & strain == "Strain 1"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light() +  ylab("Population growth rate r (indiv/(mL hours)") + xlab("Population density (indiv/mL)") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig6 <- ggplot(filter(predict.data, pHmedium=="pH 6.5" & strain == "Strain 2"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light() + xlab("Population density (indiv/mL)") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig7 <- ggplot(filter(predict.data, pHmedium=="pH 6.5" & strain == "Strain 3"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light()  + xlab("Population density (indiv/mL)") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) ) + scale_x_continuous(labels=scientific) +ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))

fig8 <- ggplot(filter(predict.data, pHmedium=="pH 6.5" & strain == "Strain 4"), aes(x = dens, y = fitness, group = name, colour = origin)) + geom_line(size=1.5) +
  theme_light()  + xlab("Population density (indiv/mL)") +
  theme( axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=20), strip.text.x=element_text(20),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
         plot.title = element_text(size=20, hjust = 0.5),axis.text.x = element_text(angle=90, hjust=1) ) + scale_x_continuous(labels=scientific) + ylim(-0.05, 0.4) +
  scale_color_manual(values=c("#2db3f9", "#2d33f9", "#e41a1c"), breaks=unique(predict.data$origin), name="origin", labels=c("neutral-evolved", "low-evolved", "ancestral"))


#Combine figures
library(ggpubr)
ggarrange(fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8,  nrow=2, ncol=4, align = "hv", common.legend = T, legend =  "bottom",
          labels =list("A", "B", "C", "D", "E", "F", "G", "H"))
