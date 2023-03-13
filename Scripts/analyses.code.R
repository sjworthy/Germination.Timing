library(ggplot2)
library(lme4)
library(lmerTest)
library(phytools)
library(ape)
library(dplyr)
library(MuMIn)
library(piecewiseSEM)
library(lmtest)
library(tidyr)
library(boot)
library(emmeans)
library(nlme)
library(chisq.posthoc.test)


#### Calculating Germination Rate ####
# 1/days2germ2
germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")

germ.pheno[,21]=1/germ.pheno$days2germ2
colnames(germ.pheno)[21]="germ.rate" 

#### Linear model of Germination Rate~mean temperature ####
ggplot(germ.pheno, aes(x=mean.Temp, y=germ.rate, color=as.factor(Cohort), group=Pop))+
  geom_point()+
  geom_smooth(method="lm",SE=FALSE)+
  facet_wrap(~Pop)

caam.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                        data=caam.2) # significant
caan1.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=caan1.2) # significant
caan2.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=caan2.2) # significant
caco.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                        data=caco.2) # significant
cain3.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=cain3.2) # significant
cain4.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=cain4.2) # significant
stbr3.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=stbr3.2) # significant
stdi.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                        data=stdi.2) # significant
stdr2.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=stdr2.2) # significant
stgl1.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=stgl1.2) # significant
stin.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                        data=stin.2) # significant
stpo1.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                         data=stpo1.2) # significant
stto.mod.mean.Temp=lmer(germ.rate~mean.Temp+(1|Cohort), 
                        data=stto.2) # significant

#### Global LM model with of germination rate~mean.temperature ####

germ.pheno.2=germ.pheno[!is.na(germ.pheno$mean.Temp),] # remove NAs

global.rate.mod=lmer(germ.rate~mean.Temp*Pop + Cohort + (1|Block), 
                     data=germ.pheno.2)

germ.rate.emm=emtrends(global.rate.mod, "Pop", var="mean.Temp")
#write.csv(germ.rate.emm, file="germrate.meantemp.slopes.csv")
germ.rate.emm.pairs=emtrends(global.rate.mod, pairwise~Pop, var="mean.Temp")
#write.csv(germ.rate.emm.pairs$contrasts, file="germrate.meantemp.contrasts.csv")
emmip(global.rate.mod, Pop ~ mean.Temp, cov.reduce = range)

#### Linear model of germ rate~cohort ####
# need to test if linear or quadratic is better fit for each species

# make a separate dataframe for each species
# remove NAs so poly() will work
caam=germ.pheno[germ.pheno$Pop=="CAAM-GB",]
caam.2=caam[!is.na(caam$mean.Temp),]
caan1=germ.pheno[germ.pheno$Pop =="CAAN1",]
caan1.2=caan1[!is.na(caan1$mean.Temp),]
caan2=germ.pheno[germ.pheno$Pop =="CAAN2",]
caan2.2=caan2[!is.na(caan2$mean.Temp),]
caco=germ.pheno[germ.pheno$Pop =="CACO1",]
caco.2=caco[!is.na(caco$mean.Temp),]
cain3=germ.pheno[germ.pheno$Pop =="CAIN3",]
cain3.2=cain3[!is.na(cain3$mean.Temp),]
cain4=germ.pheno[germ.pheno$Pop =="CAIN4",]
cain4.2=cain4[!is.na(cain4$mean.Temp),]
stbr3=germ.pheno[germ.pheno$Pop =="STBR3",]
stbr3.2=stbr3[!is.na(stbr3$mean.Temp),]
stdi=germ.pheno[germ.pheno$Pop =="STDI",]
stdi.2=stdi[!is.na(stdi$mean.Temp),]
stdr2=germ.pheno[germ.pheno$Pop =="STDR2",]
stdr2.2=stdr2[!is.na(stdr2$mean.Temp),]
stgl1=germ.pheno[germ.pheno$Pop =="STGL1",]
stgl1.2=stgl1[!is.na(stgl1$mean.Temp),]
stin=germ.pheno[germ.pheno$Pop =="STIN",]
stin.2=stin[!is.na(stin$mean.Temp),]
stpo1=germ.pheno[germ.pheno$Pop =="STPO1",]
stpo1.2=stpo1[!is.na(stpo1$mean.Temp),]
stto=germ.pheno[germ.pheno$Pop =="STTO-BH",]
stto.2=stto[!is.na(stto$mean.Temp),]

# cohort is not treated as factor
caam.2.mod.cohort=lm(germ.rate~Cohort,
                     data=caam.2)
caan1.mod.cohort=lm(germ.rate~Cohort, 
                    data=caan1.2)
caan2.mod.cohort=lm(germ.rate~Cohort, 
                    data=caan2.2)
caco.mod.cohort=lm(germ.rate~Cohort, 
                   data=caco.2)
cain3.mod.cohort=lm(germ.rate~Cohort, 
                    data=cain3.2)
cain4.mod.cohort=lm(germ.rate~Cohort, 
                    data=cain4.2)
stbr3.mod.cohort=lm(germ.rate~Cohort, 
                    data=stbr3.2)
stdi.mod.cohort=lm(germ.rate~Cohort, 
                   data=stdi.2)
stdr2.mod.cohort=lm(germ.rate~Cohort, 
                    data=stdr2.2)
stgl1.mod.cohort=lm(germ.rate~Cohort, 
                    data=stgl1.2)
stin.mod.cohort=lm(germ.rate~Cohort, 
                   data=stin.2)
stpo1.mod.cohort=lm(germ.rate~Cohort, 
                    data=stpo1.2)
stto.mod.cohort=lm(germ.rate~Cohort, 
                   data=stto.2)

#### Quadratic model of germ rate~cohort ####

# cohort is not treated as factor
caam.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                      data=caam.2)
caan1.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=caan1.2)
caan2.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=caan2.2)
caco.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                   data=caco.2)
cain3.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=cain3.2)
cain4.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=cain4.2)
stbr3.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=stbr3.2)
stdi.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                   data=stdi.2)
stdr2.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=stdr2.2)
stgl1.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=stgl1.2)
stin.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                   data=stin.2)
stpo1.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                    data=stpo1.2)
stto.mod.cohort.q=lm(germ.rate~poly(Cohort,2), 
                   data=stto.2)

#### Seeing which model is best fit ####
# use AIC
# use Likelihood ratio test
# use Rsquare

AIC(caam.2.mod.cohort, caam.mod.cohort.q)
lrtest(caam.mod.cohort.q, caam.2.mod.cohort) # no difference in fit
rsquared(caam.2.mod.cohort)
rsquared(caam.mod.cohort.q)

AIC(caan1.mod.cohort, caan1.mod.cohort.q)
lrtest(caan1.mod.cohort.q, caan1.mod.cohort) # no difference in fit
rsquared(caan1.mod.cohort)
rsquared(caan1.mod.cohort.q)

AIC(caan2.mod.cohort, caan2.mod.cohort.q)
lrtest(caan2.mod.cohort.q, caan2.mod.cohort) # no difference in fit
rsquared(caan2.mod.cohort)
rsquared(caan2.mod.cohort.q)

AIC(caco.mod.cohort, caco.mod.cohort.q)
lrtest(caco.mod.cohort.q, caco.mod.cohort)
rsquared(caco.mod.cohort)
rsquared(caco.mod.cohort.q)

AIC(cain3.mod.cohort, cain3.mod.cohort.q)
lrtest(cain3.mod.cohort.q, cain3.mod.cohort)
rsquared(cain3.mod.cohort)
rsquared(cain3.mod.cohort.q)

AIC(cain4.mod.cohort, cain4.mod.cohort.q)
lrtest(cain4.mod.cohort.q, cain4.mod.cohort) # no difference in fit
rsquared(cain4.mod.cohort)
rsquared(cain4.mod.cohort.q)

AIC(stbr3.mod.cohort, stbr3.mod.cohort.q)
lrtest(stbr3.mod.cohort.q, stbr3.mod.cohort)
rsquared(stbr3.mod.cohort)
rsquared(stbr3.mod.cohort.q)

AIC(stdi.mod.cohort, stdi.mod.cohort.q)
lrtest(stdi.mod.cohort.q, stdi.mod.cohort)
rsquared(stdi.mod.cohort)
rsquared(stdi.mod.cohort.q)

AIC(stdr2.mod.cohort, stdr2.mod.cohort.q)
lrtest(stdr2.mod.cohort.q, stdr2.mod.cohort)
rsquared(stdr2.mod.cohort)
rsquared(stdr2.mod.cohort.q)

AIC(stgl1.mod.cohort, stgl1.mod.cohort.q)
lrtest(stgl1.mod.cohort.q, stgl1.mod.cohort) # no difference in fit
rsquared(stgl1.mod.cohort)
rsquared(stgl1.mod.cohort.q)

AIC(stin.mod.cohort, stin.mod.cohort.q)
lrtest(stin.mod.cohort.q, stin.mod.cohort)
rsquared(stin.mod.cohort)
rsquared(stin.mod.cohort.q)

AIC(stpo1.mod.cohort, stpo1.mod.cohort.q)
lrtest(stpo1.mod.cohort.q, stpo1.mod.cohort) # no difference in fit
rsquared(stpo1.mod.cohort)
rsquared(stpo1.mod.cohort.q)

AIC(stto.mod.cohort, stto.mod.cohort.q)
lrtest(stto.mod.cohort.q, stto.mod.cohort)
rsquared(stto.mod.cohort)
rsquared(stto.mod.cohort.q)

#### Testing for difference between cohort 2 and 4 ####

# cohort has to be factor to compare with emmeans

caam.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=caam.2)
caam.m1=emmeans(caam.mod.cohort,pairwise~Cohort)
# p = 0.9226

caan1.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=caan1.2)
caan1.m1=emmeans(caan1.mod.cohort,pairwise~Cohort)
# p = 0.8522

caan2.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=caan2.2)
caan2.m1=emmeans(caan2.mod.cohort,pairwise~Cohort)
# p = 0.3676

caco.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                   data=caco.2)
caco.m1=emmeans(caco.mod.cohort,pairwise~Cohort)
# best fit was quadratic
# no cohort 4 for comparison

cain3.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=cain3.2)
cain3.m1=emmeans(cain3.mod.cohort,pairwise~Cohort)
# p = <.0001

cain4.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=cain4.2)
cain4.m1=emmeans(cain4.mod.cohort,pairwise~Cohort)
# p = 0.9975

stbr3.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=stbr3.2)
stbr3.m1=emmeans(stbr3.mod.cohort,pairwise~Cohort)
# p = 0.2244

stdi.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                   data=stdi.2)
stdi.m1=emmeans(stdi.mod.cohort,pairwise~Cohort)
# p = <0.001

stdr2.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=stdr2.2)
stdr2.m1=emmeans(stdr2.mod.cohort,pairwise~Cohort)
# p = 0.9996

stgl1.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=stgl1.2)
stgl1.m1=emmeans(stgl1.mod.cohort,pairwise~Cohort)
# p = 0.1995

stin.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                   data=stin.2)
stin.m1=emmeans(stin.mod.cohort,pairwise~Cohort)
# p = <0.001

stpo1.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                    data=stpo1.2)
stpo1.m1=emmeans(stpo1.mod.cohort,pairwise~Cohort)
# p = <0.001

stto.mod.cohort=lm(germ.rate~as.factor(Cohort), 
                   data=stto.2)
stto.m1=emmeans(stto.mod.cohort,pairwise~Cohort)
# p = <0.001

#### Calculating Germination Proportions ####
# calculation for species not considering each block separately
germ.proport=matrix(data=NA, ncol=13, nrow=7)
colnames(germ.proport)=c("CAAM_GB","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3","STDI","STDR2","STGL1",
                         "STIN","STPO1","STTO_BH")
row.names(germ.proport)=c("Cohort.1","Cohort.2","Cohort.3","Cohort.4","Cohort.5","Cohort.6","Cohort.7")

# repeat for loop below for each pop, also change column number in germ.per
for(i in 1:7){
  pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
  cohort=subset(pop, pop$Cohort==i)
  germ.proport[i,13]=(sum(cohort$germinated)/sum(cohort$planted))
}
germ.proport=as.data.frame(germ.proport)
germ.proport.2=t(germ.proport) # transpose the dataframe

write.csv(germ.proport, file="germ.proportion.csv")

# calculating germination proportion for each block of each species/pop
germ.proport.block=matrix(data=NA, ncol=13, nrow=21)
colnames(germ.proport.block)=c("CAAM_GB","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3","STDI","STDR2","STGL1",
                               "STIN","STPO1","STTO_BH")
row.names(germ.proport.block)=c("Cohort.1.Block.1","Cohort.1.Block.2","Cohort.1.Block.3",
                                "Cohort.2.Block.1","Cohort.2.Block.2","Cohort.2.Block.3",
                                "Cohort.3.Block.1","Cohort.3.Block.2","Cohort.3.Block.3",
                                "Cohort.4.Block.1","Cohort.4.Block.2","Cohort.4.Block.3",
                                "Cohort.5.Block.1","Cohort.5.Block.2","Cohort.5.Block.3",
                                "Cohort.6.Block.1","Cohort.6.Block.2","Cohort.6.Block.3",
                                "Cohort.7.Block.1","Cohort.7.Block.2","Cohort.7.Block.3")
# code repeated for each species/pop
pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==1)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[1,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==1)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[2,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==1)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[3,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==2)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[4,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==2)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[5,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==2)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[6,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==3)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[7,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==3)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[8,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==3)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[9,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==4)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[10,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==4)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[11,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==4)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[12,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==5)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[13,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==5)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[14,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==5)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[15,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==6)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[16,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==6)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[17,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==6)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[18,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==7)
block=subset(cohort, cohort$Block == 1)
germ.proport.block[19,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==7)
block=subset(cohort, cohort$Block == 2)
germ.proport.block[20,13]=(sum(block$germinated)/sum(block$planted))

pop=subset(germ.pheno, germ.pheno$Pop=="STTO-BH")
cohort=subset(pop, pop$Cohort==7)
block=subset(cohort, cohort$Block == 3)
germ.proport.block[21,13]=(sum(block$germinated)/sum(block$planted))

germ.proport.block=as.data.frame(germ.proport.block)
write.csv(germ.proport.block, file="germ.proportion.blocks.csv")

# adding count of germination vs total seeds to germ.proportion.tempdiff.csv
# code repeated for each species/pop

# gives number that germinated
block.germ=germ.pheno  %>%
  filter(Pop=="CAAM-GB")%>%
  group_by(Cohort,Block) %>% 
  summarize(germinated.block = sum(germinated))

write.csv(block.plant, file="block.germ.csv")

# gives total number planted
block.plant=germ.pheno%>%
  filter(Pop=="CAAM-GB")%>%
  group_by(Cohort,Block) %>% 
  summarize(planted.block = sum(planted))

# gives mean Temp.Differences for each block within each Cohort
block.temp.diff=stto   %>%
  group_by(Cohort,Block) %>% 
  summarize(average.temp.diff = mean(Temp.Differance, na.rm=TRUE))

#### Plotting Germination Proportions ####
germ.proport.block.long=read.csv("Germination.Timing/Formatted.Data/germ.proport.block.long.csv", header=T, row.names = 1)

germ.proport.means <- germ.proport.block.long %>% 
  group_by(Pop, Cohort) %>% 
  summarize(mean.germ.proport=mean(Germ.Proport))
germ.proport.means

germ.proport.sd <- germ.proport.block.long %>% 
  group_by(Pop, Cohort) %>% 
  summarize(sd.germ.proport=sd(Germ.Proport))
germ.proport.sd

germ.proport.mean.sd=merge(germ.proport.means, germ.proport.sd)
write.csv(germ.proport.mean.sd, file="germ.proport.mean.sd.csv")

#### GLMER model of Germination Proportion~mean temperature ####

# gives mean of mean Temp for each block within each Cohort
# output for each pop and add to germ.proportion.block.tempdiff.csv
block.mean.temp=germ.pheno %>%
  filter(Pop=="CAAM-GB")%>%
  group_by(Cohort,Block) %>% 
  summarize(mean.temp= mean(mean.Temp, na.rm=TRUE))
write.csv(block.mean.temp, file="block.mean.temp.csv")

germ.proport.temp.block=read.csv("Germination.Timing/Formatted.Data/germ.proportion.block.tempdiff.csv")

# remove NAs so poly() works 
caam.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caam.mean.temp),]
caan1.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caan1.mean.temp),]
caan2.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caan2.mean.temp),]
caco.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caco.mean.temp),]
cain3.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$cain3.mean.temp),]
stdr2.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$stdr2.mean.temp),]
stgl1.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$stgl1.mean.temp),]

# models
y=cbind(caam.germ.proport$caam.count, (caam.germ.proport$caam.sample-caam.germ.proport$caam.count))
caam.mod.glmer=glmer(y~caam.mean.temp + (1|Block), data=caam.germ.proport, family = binomial)

y=cbind(caan1.germ.proport$caan1.count, (caan1.germ.proport$caan1.sample-caan1.germ.proport$caan1.count))
caan1.mod.glmer=glmer(y~caan1.mean.temp + (1|Block), data=caan1.germ.proport, family = binomial)

y=cbind(caan2.germ.proport$caan2.count, (caan2.germ.proport$caan2.sample-caan2.germ.proport$caan2.count))
caan2.mod.glmer=glmer(y~caan2.mean.temp + (1|Block), data=caan2.germ.proport, family = binomial)

y=cbind(caco.germ.proport$caco1.count, (caco.germ.proport$caco1.sample-caco.germ.proport$caco1.count))
caco1.mod.glmer=glmer(y~caco.mean.temp + (1|Block), data=caco.germ.proport, family = binomial)

y=cbind(cain3.germ.proport$cain3.count, (cain3.germ.proport$cain3.sample-cain3.germ.proport$cain3.count))
cain3.mod.glmer=glmer(y~cain3.mean.temp + (1|Block), data=cain3.germ.proport, family = binomial) 
# significant

y=cbind(germ.proport.temp.block$cain4.count, (germ.proport.temp.block$cain4.sample-germ.proport.temp.block$cain4.count))
cain4.mod.glmer=glmer(y~cain4.mean.temp + (1|Block), data=germ.proport.temp.block, family = binomial)
# significant

y=cbind(germ.proport.temp.block$stbr3.count, (germ.proport.temp.block$stbr3.sample-germ.proport.temp.block$stbr3.count))
stbr3.mod.glmer=glmer(y~stbr3.mean.temp + (1|Block), data=germ.proport.temp.block, family = binomial)
# significant

y=cbind(germ.proport.temp.block$stdi.count, (germ.proport.temp.block$stdi.sample-germ.proport.temp.block$stdi.count))
stdi.mod.glmer=glmer(y~stdi.mean.temp + (1|Block), data=germ.proport.temp.block, family = binomial)
# significant

y=cbind(stdr2.germ.proport$stdr2.count, (stdr2.germ.proport$stdr2.sample-stdr2.germ.proport$stdr2.count))
stdr2.mod.glmer=glmer(y~stdr2.mean.temp + (1|Block), data=stdr2.germ.proport, family = binomial)
# significant

y=cbind(stgl1.germ.proport$stgl1.count, (stgl1.germ.proport$stgl1.sample-stgl1.germ.proport$stgl1.count))
stgl1.mod.glmer=glmer(y~stgl1.mean.temp + (1|Block), data=stgl1.germ.proport, family = binomial)
# significant

y=cbind(germ.proport.temp.block$stin.count, (germ.proport.temp.block$stin.sample-germ.proport.temp.block$stin.count))
stin.mod.glmer=glmer(y~stin.mean.temp + (1|Block), data=germ.proport.temp.block, family = binomial) 
# significant

y=cbind(germ.proport.temp.block$stpo1.count, (germ.proport.temp.block$stpo1.sample-germ.proport.temp.block$stpo1.count))
stpo1.mod.glmer=glmer(y~stpo1.mean.temp + (1|Block), data=germ.proport.temp.block, family = binomial)
# significant

y=cbind(germ.proport.temp.block$stto.count, (germ.proport.temp.block$stto.sample-germ.proport.temp.block$stto.count))
stto.mod.glmer=glmer(y~stto.mean.temp + (1|Block),data=germ.proport.temp.block, family = binomial)
# significant

#### Global GLMER model with of germination proportion~mean.temperature ####
global.germ.proport.temp.block=read.csv("Germination.timing/Formatted.Data/global.germ.proportion.block.tempdiff.csv")
global.germ.proport.temp.block$Pop=as.factor(global.germ.proport.temp.block$Pop)

y=cbind(global.germ.proport.temp.block$count, (global.germ.proport.temp.block$sample-global.germ.proport.temp.block$count))
global.germ.proport.mod=glmer(y~mean.temp*Pop + (1|Block), data=global.germ.proport.temp.block, family = binomial)

germ.proport.emm=emtrends(global.germ.proport.mod, "Pop", var="mean.temp")
write.csv(germ.proport.emm, file="germproport.meantemp.slopes.csv")
germ.proport.emm.pairs=emtrends(global.germ.proport.mod, pairwise~Pop, var="mean.temp")
write.csv(germ.proport.emm.pairs$contrasts, file="germproport.meantemp.contrasts.csv")
emmip(global.germ.proport.mod, Pop ~ mean.temp, cov.reduce = range)

global.germ.proport.temp.block$proport=global.germ.proport.temp.block$count/global.germ.proport.temp.block$sample

#### GLM of germination proportion~Cohort ####
germ.proport.temp.block=read.csv("Germination.Timing/Formatted.Data/germ.proportion.block.tempdiff.csv")

# remove NAs so poly() works 
caam.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caam.mean.temp),]
caan1.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caan1.mean.temp),]
caan2.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caan2.mean.temp),]
caco.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$caco.mean.temp),]
cain3.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$cain3.mean.temp),]
stdr2.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$stdr2.mean.temp),]
stgl1.germ.proport=germ.proport.temp.block[!is.na(germ.proport.temp.block$stgl1.mean.temp),]

# models
y=cbind(caam.germ.proport$caam.count, (caam.germ.proport$caam.sample-caam.germ.proport$caam.count))
caam.mod.glm.cohort=glm(y~Cohort, data=caam.germ.proport, family = binomial)

y=cbind(caan1.germ.proport$caan1.count, (caan1.germ.proport$caan1.sample-caan1.germ.proport$caan1.count))
caan1.mod.glm.cohort=glm(y~Cohort, data=caan1.germ.proport, family = binomial)

y=cbind(caan2.germ.proport$caan2.count, (caan2.germ.proport$caan2.sample-caan2.germ.proport$caan2.count))
caan2.mod.glm.cohort=glm(y~Cohort , data=caan2.germ.proport, family = binomial)

y=cbind(caco.germ.proport$caco1.count, (caco.germ.proport$caco1.sample-caco.germ.proport$caco1.count))
caco1.mod.glm.cohort=glm(y~Cohort , data=caco.germ.proport, family = binomial)

y=cbind(cain3.germ.proport$cain3.count, (cain3.germ.proport$cain3.sample-cain3.germ.proport$cain3.count))
cain3.mod.glm.cohort=glm(y~Cohort , data=cain3.germ.proport, family = binomial) 

y=cbind(germ.proport.temp.block$cain4.count, (germ.proport.temp.block$cain4.sample-germ.proport.temp.block$cain4.count))
cain4.mod.glm.cohort=glm(y~Cohort , data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stbr3.count, (germ.proport.temp.block$stbr3.sample-germ.proport.temp.block$stbr3.count))
stbr3.mod.glm.cohort=glm(y~Cohort , data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stdi.count, (germ.proport.temp.block$stdi.sample-germ.proport.temp.block$stdi.count))
stdi.mod.glm.cohort=glm(y~Cohort , data=germ.proport.temp.block, family = binomial)

y=cbind(stdr2.germ.proport$stdr2.count, (stdr2.germ.proport$stdr2.sample-stdr2.germ.proport$stdr2.count))
stdr2.mod.glm.cohort=glm(y~Cohort , data=stdr2.germ.proport, family = binomial)

y=cbind(stgl1.germ.proport$stgl1.count, (stgl1.germ.proport$stgl1.sample-stgl1.germ.proport$stgl1.count))
stgl1.mod.glm.cohort=glm(y~Cohort, data=stgl1.germ.proport, family = binomial)

y=cbind(germ.proport.temp.block$stin.count, (germ.proport.temp.block$stin.sample-germ.proport.temp.block$stin.count))
stin.mod.glm.cohort=glm(y~Cohort , data=germ.proport.temp.block, family = binomial) 

y=cbind(germ.proport.temp.block$stpo1.count, (germ.proport.temp.block$stpo1.sample-germ.proport.temp.block$stpo1.count))
stpo1.mod.glm.cohort=glm(y~Cohort, data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stto.count, (germ.proport.temp.block$stto.sample-germ.proport.temp.block$stto.count))
stto.mod.glm.cohort=glm(y~Cohort,data=germ.proport.temp.block, family = binomial)

#### Quadratic germination proportion~Cohort ####

y=cbind(caam.germ.proport$caam.count, (caam.germ.proport$caam.sample-caam.germ.proport$caam.count))
caam.mod.glm.cohort.q=glm(y~poly(Cohort,2), data=caam.germ.proport, family = binomial)

y=cbind(caan1.germ.proport$caan1.count, (caan1.germ.proport$caan1.sample-caan1.germ.proport$caan1.count))
caan1.mod.glm.cohort.q=glm(y~poly(Cohort,2), data=caan1.germ.proport, family = binomial)

y=cbind(caan2.germ.proport$caan2.count, (caan2.germ.proport$caan2.sample-caan2.germ.proport$caan2.count))
caan2.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=caan2.germ.proport, family = binomial)

y=cbind(caco.germ.proport$caco1.count, (caco.germ.proport$caco1.sample-caco.germ.proport$caco1.count))
caco1.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=caco.germ.proport, family = binomial)

y=cbind(cain3.germ.proport$cain3.count, (cain3.germ.proport$cain3.sample-cain3.germ.proport$cain3.count))
cain3.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=cain3.germ.proport, family = binomial) 

y=cbind(germ.proport.temp.block$cain4.count, (germ.proport.temp.block$cain4.sample-germ.proport.temp.block$cain4.count))
cain4.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stbr3.count, (germ.proport.temp.block$stbr3.sample-germ.proport.temp.block$stbr3.count))
stbr3.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stdi.count, (germ.proport.temp.block$stdi.sample-germ.proport.temp.block$stdi.count))
stdi.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=germ.proport.temp.block, family = binomial)

y=cbind(stdr2.germ.proport$stdr2.count, (stdr2.germ.proport$stdr2.sample-stdr2.germ.proport$stdr2.count))
stdr2.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=stdr2.germ.proport, family = binomial)

y=cbind(stgl1.germ.proport$stgl1.count, (stgl1.germ.proport$stgl1.sample-stgl1.germ.proport$stgl1.count))
stgl1.mod.glm.cohort.q=glm(y~poly(Cohort,2), data=stgl1.germ.proport, family = binomial)

y=cbind(germ.proport.temp.block$stin.count, (germ.proport.temp.block$stin.sample-germ.proport.temp.block$stin.count))
stin.mod.glm.cohort.q=glm(y~poly(Cohort,2) , data=germ.proport.temp.block, family = binomial) 

y=cbind(germ.proport.temp.block$stpo1.count, (germ.proport.temp.block$stpo1.sample-germ.proport.temp.block$stpo1.count))
stpo1.mod.glm.cohort.q=glm(y~poly(Cohort,2), data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stto.count, (germ.proport.temp.block$stto.sample-germ.proport.temp.block$stto.count))
stto.mod.glm.cohort.q=glm(y~poly(Cohort,2),data=germ.proport.temp.block, family = binomial)

#### Seeing which model is best fit ####
# use AIC
# use Likelihood ratio test
# use Rsquare
AIC(caam.mod.glm.cohort, caam.mod.glm.cohort.q)
lrtest(caam.mod.glm.cohort.q, caam.mod.glm.cohort)
rsquared(caam.mod.glm.cohort)
rsquared(caam.mod.glm.cohort.q)

AIC(caan1.mod.glm.cohort, caan1.mod.glm.cohort.q)
lrtest(caan1.mod.glm.cohort.q, caan1.mod.glm.cohort) # linear best fit
rsquared(caan1.mod.glm.cohort)
rsquared(caan1.mod.glm.cohort.q)

AIC(caan2.mod.glm.cohort, caan2.mod.glm.cohort.q)
lrtest(caan2.mod.glm.cohort.q, caan2.mod.glm.cohort) # linear best fit
rsquared(caan2.mod.glm.cohort)
rsquared(caan2.mod.glm.cohort.q)

AIC(caco1.mod.glm.cohort, caco1.mod.glm.cohort.q)
lrtest(caco1.mod.glm.cohort.q, caco1.mod.glm.cohort) # linear best fit
rsquared(caco.mod.glm.cohort)
rsquared(caco.mod.glm.cohort.q)

AIC(cain3.mod.glm.cohort, cain3.mod.glm.cohort.q)
lrtest(cain3.mod.glm.cohort.q, cain3.mod.glm.cohort) # linear best fit
rsquared(cain3.mod.glm.cohort)
rsquared(cain3.mod.glm.cohort.q)

AIC(cain4.mod.glm.cohort, cain4.mod.glm.cohort.q)
lrtest(cain4.mod.glm.cohort.q, cain4.mod.glm.cohort) # linear best fit
rsquared(cain4.mod.glm.cohort)
rsquared(cain4.mod.glm.cohort.q)

AIC(stbr3.mod.glm.cohort, stbr3.mod.glm.cohort.q)
lrtest(stbr3.mod.glm.cohort.q, stbr3.mod.glm.cohort) # linear best fit
rsquared(stbr3.mod.glm.cohort)
rsquared(stbr3.mod.glm.cohort.q)

AIC(stdi.mod.glm.cohort, stdi.mod.glm.cohort.q)
lrtest(stdi.mod.glm.cohort.q, stdi.mod.glm.cohort)
rsquared(stdi.mod.glm.cohort)
rsquared(stdi.mod.glm.cohort.q)

AIC(stdr2.mod.glm.cohort, stdr2.mod.glm.cohort.q)
lrtest(stdr2.mod.glm.cohort.q, stdr2.mod.glm.cohort) # linear best fit
rsquared(stdr2.mod.glm.cohort)
rsquared(stdr2.mod.glm.cohort.q)

AIC(stgl1.mod.glm.cohort, stgl1.mod.glm.cohort.q)
lrtest(stgl1.mod.glm.cohort.q, stgl1.mod.glm.cohort)
rsquared(stgl1.mod.glm.cohort)
rsquared(stgl1.mod.glm.cohort.q)

AIC(stin.mod.glm.cohort, stin.mod.glm.cohort.q)
lrtest(stin.mod.glm.cohort.q, stin.mod.glm.cohort) # linear best fit
rsquared(stin.mod.glm.cohort)
rsquared(stin.mod.glm.cohort.q)

AIC(stpo1.mod.glm.cohort, stpo1.mod.glm.cohort.q)
lrtest(stpo1.mod.glm.cohort.q, stpo1.mod.glm.cohort) # linear best fit
rsquared(stpo1.mod.glm.cohort)
rsquared(stpo1.mod.glm.cohort.q)

AIC(stto.mod.glm.cohort, stto.mod.glm.cohort.q)
lrtest(stto.mod.glm.cohort.q, stto.mod.glm.cohort) # linear best fit
rsquared(stto.mod.glm.cohort)
rsquared(stto.mod.glm.cohort.q)

#### Testing for difference between cohort 2 and 4 ####
# cohort has to be factor to compare with emmeans

y=cbind(caam.germ.proport$caam.count, (caam.germ.proport$caam.sample-caam.germ.proport$caam.count))
caam.mod.glm=glm(y~as.factor(Cohort), data=caam.germ.proport, family = binomial)
caam.m1.glm=emmeans(caam.mod.glm,pairwise~Cohort)
# p = 0.9988

y=cbind(caan1.germ.proport$caan1.count, (caan1.germ.proport$caan1.sample-caan1.germ.proport$caan1.count))
caan1.mod.glm=glm(y~as.factor(Cohort), data=caan1.germ.proport, family = binomial)
caan1.m1.glm=emmeans(caan1.mod.glm,pairwise~Cohort)
# p = 0.9192

y=cbind(caan2.germ.proport$caan2.count, (caan2.germ.proport$caan2.sample-caan2.germ.proport$caan2.count))
caan2.mod.glm=glm(y~as.factor(Cohort) , data=caan2.germ.proport, family = binomial)
caan2.m1.glm=emmeans(caan2.mod.glm,pairwise~Cohort)
# p = 0.9999

y=cbind(caco.germ.proport$caco1.count, (caco.germ.proport$caco1.sample-caco.germ.proport$caco1.count))
caco1.mod.glm=glm(y~as.factor(Cohort) , data=caco.germ.proport, family = binomial)
caco1.m1.glm=emmeans(caco1.mod.glm,pairwise~Cohort)
# no cohort 4

y=cbind(cain3.germ.proport$cain3.count, (cain3.germ.proport$cain3.sample-cain3.germ.proport$cain3.count))
cain3.mod.glm=glm(y~as.factor(Cohort) , data=cain3.germ.proport, family = binomial) 
cain3.m1.glm=emmeans(cain3.mod.glm,pairwise~Cohort)
# p = 1.0000

y=cbind(germ.proport.temp.block$cain4.count, (germ.proport.temp.block$cain4.sample-germ.proport.temp.block$cain4.count))
cain4.mod.glm=glm(y~as.factor(Cohort) , data=germ.proport.temp.block, family = binomial)
cain4.m1.glm=emmeans(cain4.mod.glm,pairwise~Cohort)
# p = 0.0013

y=cbind(germ.proport.temp.block$stbr3.count, (germ.proport.temp.block$stbr3.sample-germ.proport.temp.block$stbr3.count))
stbr3.mod.glm=glm(y~as.factor(Cohort) , data=germ.proport.temp.block, family = binomial)
stbr3.m1.glm=emmeans(stbr3.mod.glm,pairwise~Cohort)
# p = 0.2458

y=cbind(germ.proport.temp.block$stdi.count, (germ.proport.temp.block$stdi.sample-germ.proport.temp.block$stdi.count))
stdi.mod.glm=glm(y~as.factor(Cohort) , data=germ.proport.temp.block, family = binomial)
stdi.m1.glm=emmeans(stdi.mod.glm,pairwise~Cohort)
# p = 

y=cbind(stdr2.germ.proport$stdr2.count, (stdr2.germ.proport$stdr2.sample-stdr2.germ.proport$stdr2.count))
stdr2.mod.glm=glm(y~as.factor(Cohort) , data=stdr2.germ.proport, family = binomial)
stdr2.m1.glm=emmeans(stdr2.mod.glm,pairwise~Cohort)
# p = 0.9486

y=cbind(stgl1.germ.proport$stgl1.count, (stgl1.germ.proport$stgl1.sample-stgl1.germ.proport$stgl1.count))
stgl1.mod.glm=glm(y~as.factor(Cohort) , data=stgl1.germ.proport, family = binomial)
stgl1.m1.glm=emmeans(stgl1.mod.glm,pairwise~Cohort)
# p = 

y=cbind(germ.proport.temp.block$stin.count, (germ.proport.temp.block$stin.sample-germ.proport.temp.block$stin.count))
stin.mod.glm=glm(y~as.factor(Cohort) , data=germ.proport.temp.block, family = binomial) 
stin.m1.glm=emmeans(stin.mod.glm,pairwise~Cohort)
# p = 0.7477

y=cbind(germ.proport.temp.block$stpo1.count, (germ.proport.temp.block$stpo1.sample-germ.proport.temp.block$stpo1.count))
stpo1.mod.glm=glm(y~as.factor(Cohort), data=germ.proport.temp.block, family = binomial)
stpo1.m1.glm=emmeans(stpo1.mod.glm,pairwise~Cohort)
# p = 0.9985

y=cbind(germ.proport.temp.block$stto.count, (germ.proport.temp.block$stto.sample-germ.proport.temp.block$stto.count))
stto.mod.glm=glm(y~as.factor(Cohort),data=germ.proport.temp.block, family = binomial)
stto.m1.glm=emmeans(stto.mod.glm,pairwise~Cohort)
# p = 0.0029

#### Reading in Phylogeny ####
all.phylo <- read.tree("Germination.Timing/Raw.Data/tree_pruned.new")

sp.list=c("Caulanthus_amplexicaulis","Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)

# artificially add branch for CAAN2 and CAIN4

plotTree(phylo)
nodelabels()
edgelabels()

which(phylo$tip.label=="Caulanthus_anceps") # 3,edge 5

tree.1=bind.tip(phylo, tip.label = "Caulanthus_anceps.2", edge.length=0.001497901, where = 3, position = 0)
plot(tree.1)
nodelabels()
edgelabels()

which(tree.1$tip.label=="Caulanthus_inflatus") # 1, edge 2

tree.2=bind.tip(tree.1, tip.label = "Caulanthus_inflatus.2", edge.length=0.001158696, where = 1, position = 0)
plot(tree.2)

#### Phylogenetic Signal in Slopes ####
# testing for phylogenetic signal in slopes of germ.rate ~ temp.diff

# read in slopes  dataframe
model.slopes=read.csv("Germination.Timing/Formatted.Data/slopes.csv", row.names = 1)
colnames(model.slopes)[1]="sp"
phylosig(tree.2, model.slopes$germrate.slope, method = "K", test = TRUE, nsim=1000)
# Lapack routine dgesv: system is exactly singular: U[13,13] = 0
# matrix is singular for tree.2 but not for original phylo
# can't have edge length = 0 so adding small amount.

# populations have edge lengths of 0
phylosig(tree.2, model.slopes$germrate.slope, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.990482, P = 0.0107472
phylosig(tree.2, model.slopes$germproport.slope, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.844796, P = 0.0783106
phylosig(tree.2, model.slopes$germrate.cohort, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.996092, P = 0.00373091
phylosig(tree.2, model.slopes$germproport.cohort, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.867766, P = 0.0518399

# make populations have edge length not equal to 0, this makes K work
tree.2$edge.length <-  tree.2$edge.length + 0.001

phylosig(tree.2, model.slopes$germrate.slope, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.999934, P = 0.656032
phylosig(tree.2, model.slopes$germproport.slope, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.999934, P = 0.118127
phylosig(tree.2, model.slopes$germrate.cohort, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.999934, P = 0.504427
phylosig(tree.2, model.slopes$germproport.cohort, method = "lambda", test = TRUE, nsim=1000)
# lambda = 0.999934, P = 0.0731766

#### Differences in proportion between Round 1 and Round 2 ####

# Round 1 data
global.germ.proport.temp.block=read.csv("Germination.timing/Formatted.Data/global.germ.proportion.block.tempdiff.csv")
global.germ.proport.temp.block$Pop=as.factor(global.germ.proport.temp.block$Pop)
global.germ.proport.temp.block$Round=1

global.germ.proport.temp.block = global.germ.proport.temp.block %>%
mutate(Pop = recode(Pop, "caan1" = "CAAN1", "cain3" = "CAIN3", "caam" = "CAAM", "caco" = "CACO", 
                        "stto" = "STTO", "stdr2" = "STDR", "stpo1" = "STPO", "stbr3" = "STBR", 
                        "stdi" = "STDI", "stgl1" = "STGL", "stin" = "STIN", "caan2"="CAAN2","cain4"="CAIN4"))

# Round 2 data
germ.pheno.R2.all=read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.formatted.csv", row.names = 1)
germ.pheno.R2=subset(germ.pheno.R2.all, germ.pheno.R2.all$Population !="blank")

R2.germinated=germ.pheno.R2 %>%
  group_by(Population,Cohort,Old.Block) %>%
  summarise(count=sum(germinated),
            sample=sum(planted))

R2.germinated$Round=2
colnames(R2.germinated)[3]="Block"
colnames(R2.germinated)[1]="Pop"
R2.germinated=R2.germinated %>%
mutate(Pop = recode(Pop, "CAAN1" = "CAAN1", "CAIN3" = "CAIN3", "CAAM-GB" = "CAAM", "CACO1" = "CACO", 
                    "STTO-BH" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                    "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN", "CAAN2"="CAAN2","CAIN4"="CAIN4"))

all.proport.data=full_join(global.germ.proport.temp.block, R2.germinated)
all.proport.data$Cohort=as.factor(all.proport.data$Cohort)
all.proport.data$Round=as.factor(all.proport.data$Round)
all.proport.data$Pop=as.factor(all.proport.data$Pop)

write.csv(all.proport.data, file="all.proport.data.long.csv")

#### Contingency analyses ####

germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")

germ.pheno.R2.all=read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.formatted.csv", row.names = 1)
germ.pheno.R2=subset(germ.pheno.R2.all, germ.pheno.R2.all$Population !="blank")

total.planted=germ.pheno %>%
  group_by(Pop) %>%
  summarise(total.count=sum(planted),
            germ.R1=sum(germinated))

total.germ.R2=germ.pheno.R2 %>%
  group_by(Population)%>%
  summarise(germ.R2=sum(germinated),
            planted.R2=sum(planted))

pop.table=read.csv("Germination.Timing/Formatted.Data/Pop.contingeny.table.csv", header=T, row.names=1)

pop.chisq=chisq.test(pop.table)
# significant X = 1006.1, df = 24, p < 2.2e-16
pop.chisq$residuals
# CAAN1, CAAN2 positive association with year 2, negative association with year 1
# CAAM, CACO, CAIN4, STIN, negative association with both years, positive association with Never germinating
pop.posthoc=chisq.posthoc.test(pop.table) # bonferroni corrected p values
# Significant
# CAAM, CAAN2, STTO STBR all three
# CAAN1 year 1 and year 2
# CAIN4, STGL year 2 and never
# STDR year 2
# STIN, STPO year 1 and never

# write.csv(pop.posthoc, file="pop.chisq.posthoc.results.csv")

# Pop and Cohort
germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")

germ.pheno.R2.all=read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.formatted.csv", row.names = 1)
germ.pheno.R2=subset(germ.pheno.R2.all, germ.pheno.R2.all$Population !="blank")

total.planted=germ.pheno %>%
  group_by(Pop,Cohort) %>%
  summarise(total.count=sum(planted),
            germ.R1=sum(germinated))
# write.csv(total.planted, file="test.output.csv")

total.germ.R2=germ.pheno.R2 %>%
  group_by(Population,Cohort)%>%
  summarise(germ.R2=sum(germinated),
            planted.R2=sum(planted))
#write.csv(total.germ.R2, file="test.output.csv")

cohort.table=read.csv("Germination.Timing/Formatted.Data/cohort.contingency.table.csv", header=T, row.names=1)

cohort.chisq=chisq.test(cohort.table)
# In chisq.test(cohort.table) : Chi-squared approximation may be incorrect because some cells have < 5
# significant X = 1919, df = 180, p < 2.2e-16

cohort.chisq$residuals
# CAAM negative association for all cohort for both years, positive association for all cohorts with Never
# CAAN1 negative association for all cohort for Year 1, positive association for all cohorts Year 2
# CAAN2 negative association for all cohort for Year 1, positive association for all cohorts Year 2 except cohort 7
# 11 of 13 species have negative association between cohort 7 and Year 1, all except STPO and STTO

cohort.posthoc=chisq.posthoc.test(cohort.table) # bonferroni corrected p values

