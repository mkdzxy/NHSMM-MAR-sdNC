
# Top line ----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(ggsignif)


# Simulation1 -------------------------------------------------------------

ACC_HMM <- read.csv("Simulation1/ACC_HMM_MAR.txt", header=FALSE)
ACC_NHSMM <- read.csv("Simulation1/ACC_NHSMM_MAR.txt", header=FALSE)
VN_HMM <- read.csv("Simulation1/VN_HMM_MAR.txt", header=FALSE)
VN_NHSMM <- read.csv("Simulation1/VN_NHSMM_MAR.txt", header=FALSE)

apply(ACC_HMM, 2, mean)
apply(ACC_HMM, 2, sd)
apply(ACC_NHSMM, 2, mean)
apply(ACC_NHSMM, 2, sd)

meanHMM <- 1:8
sdHMM <- 1:8
meanNHSMM <- 1:8
sdNHSMM <- 1:8

Abnormalvaluedelete <- function(df1,df2){
  r <- length(df1)
  ddf1 <- c()
  ddf2 <- c()
  mean1 <- mean(df1)
  mean2 <- mean(df2)
  sd1 <- sd(df1)
  sd2 <- sd(df2)
  i <- 1
  sdtimes <- 1.5
  while (i<=r) {
    if (df1[i]>=(mean1-sdtimes*sd1) && df1[i]<=(mean1+sdtimes*sd1) && df2[i]>=(mean2-sdtimes*sd2) && df2[i]<=(mean2+sdtimes*sd2)) {
      ddf1 <- c(ddf1,df1[i])
      ddf2 <- c(ddf2,df2[i])
      i <- i+1
    } else{
      i <- i+1
    }
  } 
  res <- list()
  res[[1]] <- ddf1
  res[[2]] <- ddf2
  res
}

for (i in 1:8) {
  l <- Abnormalvaluedelete(ACC_HMM[,i],ACC_NHSMM[,i])
  meanHMM[i] <- mean(l[[1]])
  sdHMM[i] <- sd(l[[1]])
  meanNHSMM[i] <- mean(l[[2]])
  sdNHSMM[i] <- sd(l[[2]])
}

ACC_HMM_matrix <- matrix(data.matrix(ACC_HMM),50*8,1)
VN_HMM_matrix <- matrix(data.matrix(VN_HMM),50*8,1)
ACC_NHSMM_matrix <- matrix(data.matrix(ACC_NHSMM),50*8,1)
VN_NHSMM_matrix <- matrix(data.matrix(VN_NHSMM),50*8,1)

dataAcc_merge <- data.frame(kbar=rep(rep(3:10,each=50),2),
                            Acc=c(ACC_HMM_matrix,ACC_NHSMM_matrix),
                            group=rep(c("HMM-MAR","NHSMM-MAR"),each=400))
dataAcc_merge$kbar <- as.factor(dataAcc_merge$kbar)

sig <- c()

for (i in 3:10){
  a <- t.test(ACC_HMM[,i-2],ACC_NHSMM[,i-2])
  a$p.value
  sig <- c(sig,a$p.value)
}

p_Acc <- dataAcc_merge %>% ggplot(aes(x=kbar,y=Acc,fill=group))+
  # stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+
  stat_boxplot(geom = "errorbar",width=0.5,aes(x=kbar,y=Acc))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(dodge.width = 0.5,jitter.width = 0.06),
             pch=21,
             alpha=0.4,
             size=2
  )+
  theme_bw()+
  scale_fill_manual(values = c("HMM-MAR"="#baadd0", "NHSMM-MAR"="#fbbf87"))+
  geom_signif(y_position = rep(0.87,8),
              xmin=c(0.87,1.87,2.87,3.87,4.87,5.87,6.87,7.87),
              xmax=c(1.13,2.13,3.13,4.13,5.13,6.13,7.13,8.13),
              annotation = c(rep("**",8)),
              size=0.6,
              textsize = 7,
              family = "serif",
              vjust = 0.7,
              tip_length = 0.008,
  )+
  theme(aspect.ratio = 6/10)+
  scale_y_continuous(limits = c(0.3,0.9),breaks = seq(0.3,0.9,0.1),labels =seq(0.3,0.9,0.1) )+
  theme(aspect.ratio = 6/8,
        strip.text.x = element_text(
          size = 15
        ), 
        strip.text.y = element_text(
          size = 15
        ),
        legend.key.size = unit(50, "pt")
  )+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=28),
        axis.title.x =element_text(size=28), 
        axis.title.y=element_text(size=28),
        axis.text.y=element_text(size=28),
        legend.key = element_rect(
          color = "black",
          fill = "white"
        ),
        legend.key.size = unit(0.5, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(
          size =  28,
          color = "black",
        ),
        legend.title = element_text(
          size =  28,
          face = "bold",
          hjust = 0,
          color = "black",
        ),
        legend.text.align = NULL,
        legend.title.align = NULL,
        legend.direction = "vertical",
        legend.box = NULL,
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),)

# Simulation2 -------------------------------------------------------------

HMM_geo_ACC <- read.csv("Simulation2/ACC_HMMMAR_geo.txt", header=FALSE)
HMM_non_ACC <- read.csv("Simulation2/ACC_HMMMAR_non.txt", header=FALSE)
NHSMM_geo_ACC <- read.csv("Simulation2/ACC_NHSMMMAR_geo.txt", header=FALSE)
NHSMM_non_ACC <- read.csv("Simulation2/ACC_NHSMMMAR_non.txt", header=FALSE)
HMM_geo_VN <- read.csv("Simulation2/VN_HMMMAR_geo.txt", header=FALSE)
HMM_non_VN <- read.csv("Simulation2/VN_HMMMAR_non.txt", header=FALSE)
NHSMM_geo_VN <- read.csv("Simulation2/VN_NHSMMMAR_geo.txt", header=FALSE)
NHSMM_non_VN <- read.csv("Simulation2/VN_NHSMMMAR_non.txt", header=FALSE)

apply(HMM_geo_VN, 2, mean)
apply(HMM_geo_VN, 2, sd)
apply(HMM_non_VN, 2, mean)
apply(HMM_non_VN, 2, sd)
apply(NHSMM_geo_VN, 2, mean)
apply(NHSMM_geo_VN, 2, sd)
apply(NHSMM_non_VN, 2, mean)
apply(NHSMM_non_VN, 2, sd)

meanHMMgeo <- 1:5
sdHMMgeo <- 1:5
meanHMMnon <- 1:5
sdHMMnon <- 1:5
meanNHSMMgeo <- 1:5
sdNHSMMgeo <- 1:5
meanNHSMMnon <- 1:5
sdNHSMMnon <- 1:5
for (i in 1:5) {
  l <- Abnormalvaluedelete(HMM_geo_ACC[,i],HMM_non_ACC[,i])
  meanHMMgeo[i] <- mean(l[[1]])
  sdHMMgeo[i] <- sd(l[[1]])
  meanHMMnon[i] <- mean(l[[2]])
  sdHMMnon[i] <- sd(l[[2]])
}

for (i in 1:5) {
  l <- Abnormalvaluedelete(NHSMM_geo_ACC[,i],NHSMM_non_ACC[,i])
  meanNHSMMgeo[i] <- mean(l[[1]])
  sdNHSMMgeo[i] <- sd(l[[1]])
  meanNHSMMnon[i] <- mean(l[[2]])
  sdNHSMMnon[i] <- sd(l[[2]])
}


HMM_geo_ACC_matrix <- matrix(data.matrix(HMM_geo_ACC),50*5,1)
HMM_non_ACC_matrix <- matrix(data.matrix(HMM_non_ACC),50*5,1)
NHSMM_geo_ACC_matrix <- matrix(data.matrix(NHSMM_geo_ACC),50*5,1)
NHSMM_non_ACC_matrix <- matrix(data.matrix(NHSMM_non_ACC),50*5,1)

geodf <- data.frame(ACC=c(HMM_geo_ACC_matrix,NHSMM_geo_ACC_matrix),
                    model=rep(c("HMM-MAR","NHSMM-MAR"),each=250),
                    N=rep(rep(2:6,each=50),2))
geodf$N <- as.factor(geodf$N)

nondf <- data.frame(ACC=c(HMM_non_ACC_matrix,NHSMM_non_ACC_matrix),
                    model=rep(c("HMM-MAR","NHSMM-MAR"),each=250),
                    N=rep(rep(2:6,each=50),2))
nondf$N <- as.factor(nondf$N)


HMMdf <- data.frame(ACC=c(HMM_geo_ACC_matrix,HMM_non_ACC_matrix),
                    duration=rep(c("geo","non"),each=250),
                    N=rep(rep(2:6,each=50),2))
HMMdf$N <- as.factor(HMMdf$N)

NHSMMdf <- data.frame(ACC=c(NHSMM_geo_ACC_matrix,NHSMM_non_ACC_matrix),
                      duration=rep(c("geo","non"),each=250),
                      N=rep(rep(2:6,each=50),2))
NHSMMdf$N <- as.factor(NHSMMdf$N)

# geop
signif <- c()

for (i in 2:6){
  a <- t.test(HMM_geo_ACC[,i-1],NHSMM_geo_ACC[,i-1])
  a$p.value
  signif <- c(signif,a$p.value)
}

geop <- geodf %>% 
  ggplot(aes(x=N,y=ACC,fill=model))+
  # stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+
  stat_boxplot(geom = "errorbar",width=0.5,aes(x=N,y=ACC))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(dodge.width = 0.5,jitter.width = 0.06),
             pch=21,
             alpha=0.4,
             size=2
  )+
  theme_bw()+
  scale_fill_manual(values = c("HMM-MAR"="#baadd0", "NHSMM-MAR"="#fbbf87"))+
  scale_y_continuous(limits = c(0.1,0.9),breaks = seq(0.1,0.9,0.1),labels =seq(0.1,0.9,0.1) )+
  geom_signif(y_position = rep(0.87,5),
              xmin=c(0.87,1.87,2.87,3.87,4.87),
              xmax=c(1.13,2.13,3.13,4.13,5.13),
              annotation = c(rep("**",5)),
              size=0.8,
              textsize = 7,
              family = "serif",
              vjust = 0.7,
              tip_length = 0.008,
  )+
  theme(aspect.ratio = 6/10)+
  theme(aspect.ratio = 6/8,
        strip.text.x = element_text(
          size = 15
        ), 
        strip.text.y = element_text(
          size = 15
        ),
        legend.key.size = unit(50, "pt")
  )+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=28,color = "black",),
        axis.title.x =element_text(size=28), 
        axis.title.y=element_text(size=28),
        axis.text.y=element_text(size=28,color = "black",),
        legend.key = element_rect(
          color = "black",
          fill = "white"
        ),
        legend.key.size = unit(0.5, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(
          size =  28,
          color = "black",
        ),
        legend.title = element_text(
          size =  28,
          face = "bold",
          hjust = 0,
          color = "black",
        ),
        legend.text.align = NULL,
        legend.title.align = NULL,
        legend.direction = "vertical",
        legend.box = NULL,
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),)

# nonp
signif <- c()

for (i in 2:6){
  a <- t.test(HMM_non_ACC[,i-1],NHSMM_non_ACC[,i-1])
  a$p.value
  signif <- c(signif,a$p.value)
}

nonp <- nondf %>% 
  ggplot(aes(x=N,y=ACC,fill=model))+
  # stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+
  stat_boxplot(geom = "errorbar",width=0.5,aes(x=N,y=ACC))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(dodge.width = 0.5,jitter.width = 0.06),
             pch=21,
             alpha=0.4,
             size=2
  )+
  theme_bw()+
  scale_fill_manual(values = c("HMM-MAR"="#baadd0", "NHSMM-MAR"="#fbbf87"))+
  scale_y_continuous(limits = c(0.1,0.9),breaks = seq(0.1,0.9,0.1),labels =seq(0.1,0.9,0.1) )+
  geom_signif(y_position = rep(0.9,5),
              xmin=c(0.87,1.87,2.87,3.87,4.87),
              xmax=c(1.13,2.13,3.13,4.13,5.13),
              annotation = c(rep("**",5)),
              size=0.8,
              textsize = 7,
              family = "serif",
              vjust = 0.7,
              tip_length = 0.008,
  )+
  theme(aspect.ratio = 6/10)+
  theme(aspect.ratio = 6/8,
        strip.text.x = element_text(
          size = 15
        ), 
        strip.text.y = element_text(
          size = 15
        ),
        legend.key.size = unit(50, "pt")
  )+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=28,color = "black",),
        axis.title.x =element_text(size=28), 
        axis.title.y=element_text(size=28),
        axis.text.y=element_text(size=28,color = "black",),
        legend.key = element_rect(
          color = "black",
          fill = "white"
        ),
        legend.key.size = unit(0.5, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(
          size =  28,
          color = "black",
        ),
        legend.title = element_text(
          size =  28,
          face = "bold",
          hjust = 0,
          color = "black",
        ),
        legend.text.align = NULL,
        legend.title.align = NULL,
        legend.direction = "vertical",
        legend.box = NULL,
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),)

# HMMp
signif <- c()

for (i in 2:6){
  l <- Abnormalvaluedelete(HMM_geo_ACC[,i-1],HMM_non_ACC[,i-1])
  a <- t.test(l[[1]],l[[2]])
  a$p.value
  signif <- c(signif,a$p.value)
}

HMMp <- HMMdf %>% 
  ggplot(aes(x=N,y=ACC,fill=duration))+
  # stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+
  stat_boxplot(geom = "errorbar",width=0.5,aes(x=N,y=ACC))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(dodge.width = 0.5,jitter.width = 0.06),
             pch=21,
             alpha=0.4,
             size=2
  )+
  theme_bw()+
  scale_fill_manual(values = c("geo"="#65cbcd", "non"="#ffbfbf"))+
  scale_y_continuous(limits = c(0.1,0.9),breaks = seq(0.1,0.9,0.1),labels =seq(0.1,0.9,0.1) )+
  geom_signif(y_position = rep(0.9,4),
              xmin=c(1.87,2.87,3.87,4.87),
              xmax=c(2.13,3.13,4.13,5.13),
              annotation = c("*","**","**","**"),
              size=0.8,
              textsize = 7,
              family = "serif",
              vjust = 0.7,
              tip_length = 0.008,
  )+
  theme(aspect.ratio = 6/10)+
  theme(aspect.ratio = 6/8,
        strip.text.x = element_text(
          size = 15
        ), 
        strip.text.y = element_text(
          size = 15
        ),
        legend.key.size = unit(50, "pt")
  )+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=28,color = "black",),
        axis.title.x =element_text(size=28), 
        axis.title.y=element_text(size=28),
        axis.text.y=element_text(size=28,color = "black",),
        legend.key = element_rect(
          color = "black",
          fill = "white"
        ),
        legend.key.size = unit(0.5, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(
          size =  28,
          color = "black",
        ),
        legend.title = element_text(
          size =  28,
          face = "bold",
          hjust = 0,
          color = "black",
        ),
        legend.text.align = NULL,
        legend.title.align = NULL,
        legend.direction = "vertical",
        legend.box = NULL,
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),)

# NHSMMp
signif <- c()

for (i in 2:6){
  l <- Abnormalvaluedelete(NHSMM_geo_ACC[,i-1],NHSMM_non_ACC[,i-1])
  a <- t.test(l[[1]],l[[2]])
  a$p.value
  signif <- c(signif,a$p.value)
}

NHSMMp <- NHSMMdf %>% 
  ggplot(aes(x=N,y=ACC,fill=duration))+
  # stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+
  stat_boxplot(geom = "errorbar",width=0.5,aes(x=N,y=ACC))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(dodge.width = 0.5,jitter.width = 0.15),
             pch=21,
             alpha=0.4,
             size=2
  )+
  theme_bw()+
  scale_fill_manual(values = c("geo"="#65cbcd", "non"="#ffbfbf"))+
  scale_y_continuous(limits = c(0.1,0.9),breaks = seq(0.1,0.9,0.1),labels =seq(0.1,0.9,0.1) )+
  theme(aspect.ratio = 6/10)+
  theme(aspect.ratio = 6/8,
        strip.text.x = element_text(
          size = 15
        ), 
        strip.text.y = element_text(
          size = 15
        ),
        legend.key.size = unit(50, "pt")
  )+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(size=28,color = "black",),
        axis.title.x =element_text(size=28), 
        axis.title.y=element_text(size=28),
        axis.text.y=element_text(size=28,color = "black",),
        legend.key = element_rect(
          color = "black",
          fill = "white"
        ),
        legend.key.size = unit(0.5, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(
          size =  28,
          color = "black",
        ),
        legend.title = element_text(
          size =  28,
          face = "bold",
          hjust = 0,
          color = "black",
        ),
        legend.text.align = NULL,
        legend.title.align = NULL,
        legend.direction = "vertical",
        legend.box = NULL,
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),)

