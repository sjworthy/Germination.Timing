# Calculate average minimum, mean, and maximum temperature for each individual
# between rainfall onset data and germ date
#rm(list=ls(all=T)) #code to clear console if needed

library(dplyr)
library(lubridate)
library(ggplot2)

#### read in data sets ####
ibutton_all <- read.csv("Germination.Timing/Formatted.Data/ibutton_allcohorts.csv")
germ.pheno <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.formatted.csv", row.names = 1)

# Change character dates to actual Date
germ.pheno$germdate=as.Date(germ.pheno$germdate)
germ.pheno$plantdate=as.Date(germ.pheno$plantdate)
ibutton_all$Date=as.Date(ibutton_all$Date)

#### calculate mean temperature between plant date and germ date ####
  for(i in 1:dim(germ.pheno)[1]){ #JRG: made this more general in case germ pheno dataframe changes
  data=subset(ibutton_all, between(ibutton_all$Date, germ.pheno[i,8], germ.pheno[i,6])) 
  mean=mean(data$temp, na.rm = TRUE)
  germ.pheno[i,17]=mean
  }

summary(germ.pheno)
colnames(germ.pheno)[17]="mean.Temp"

for(i in 1:dim(germ.pheno)[1]){
  data=subset(ibutton_all, between(ibutton_all$Date, germ.pheno[i,8], germ.pheno[i,6]))
  if(nrow(data)>0){
  min=aggregate(data$temp, by=list(data$Date), min, na.rm=TRUE)
  mean=mean(min$x, na.rm = TRUE)
  germ.pheno[i,18]=mean
  }
}

colnames(germ.pheno)[18]="mean.Min.Temp"

for(i in 1:dim(germ.pheno)[1]){
  data=subset(ibutton_all, between(ibutton_all$Date, germ.pheno[i,8], germ.pheno[i,6]))
  if(nrow(data)>0){
    max=aggregate(data$temp, by=list(data$Date), max, na.rm=TRUE)
    mean=mean(max$x, na.rm = TRUE)
    germ.pheno[i,19]=mean
  }
}

colnames(germ.pheno)[19]="mean.Max.Temp"

write.csv(germ.pheno, file = "Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv")


#look at histograms
ggplot(data= germ.pheno, aes(x= mean.Temp)) + geom_histogram(binwidth = 1) + facet_wrap(~ Pop)

### getting frequency of temperatures experienced by each seed
#### read in data sets ####
ibutton_all <- read.csv("Germination.Timing/Formatted.Data/ibutton_allcohorts.csv")
germ.pheno <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.formatted.csv", row.names = 1)

# Change character dates to actual Date
germ.pheno$germdate=as.Date(germ.pheno$germdate)
germ.pheno$plantdate=as.Date(germ.pheno$plantdate)
ibutton_all$Date=as.Date(ibutton_all$Date)

#### take difference between hourly temp and Topt per seed ####
# get the average of all these differences
# get the sum of all these differences

#### read in data sets ####
ibutton_all <- read.csv("Germination.Timing/Formatted.Data/ibutton_allcohorts.csv")
ibutton_all$Date=as.Date(ibutton_all$Date)

#### calculate mean temperature between plant date and germ date ####
for(i in 1:dim(germ.pheno2)[1]){
  data=subset(ibutton_all, between(ibutton_all$Date, germ.pheno2[i,8], germ.pheno2[i,6]))
  diff=data$Value-germ.pheno2[i,20]
  mean.hourly.diff=mean(diff, na.rm = TRUE)
  sum.hourly.diff=sum(diff, na.rm = TRUE)
  germ.pheno2[i,23]=mean.hourly.diff
  germ.pheno2[i,24]=sum.hourly.diff
}

# gives the same value if you just take the Differences between Topt and mean Temp.

#### get hourly temperature between plant date and germ date ####
# need to figure out how to get the data for all seeds in a cohort
# then histogram of all hourly values within a cohort or some other summary

# subset by each pop/species
caam=germ.pheno[germ.pheno$Pop =="CAAM-GB",]
caan1=germ.pheno[germ.pheno$Pop =="CAAN1",]
caan2=germ.pheno[germ.pheno$Pop =="CAAN2",]
caco=germ.pheno[germ.pheno$Pop =="CACO1",]
cain3=germ.pheno[germ.pheno$Pop =="CAIN3",]
cain4=germ.pheno[germ.pheno$Pop =="CAIN4",]
stbr3=germ.pheno[germ.pheno$Pop =="STBR3",]
stdi=germ.pheno[germ.pheno$Pop =="STDI",]
stdr2=germ.pheno[germ.pheno$Pop =="STDR2",]
stgl1=germ.pheno[germ.pheno$Pop =="STGL1",]
stin=germ.pheno[germ.pheno$Pop =="STIN",]
stpo1=germ.pheno[germ.pheno$Pop =="STPO1",]
stto=germ.pheno[germ.pheno$Pop =="STTO-BH",]

# extract hourly data from first sow date to last germination date for each pop/species, then plot
# This would give the same values for all species/pop
range(caam$germdate, na.rm=NA)
# last germ date = "2021-01-04"
range(caam$plantdate)
# first plant date = "2020-09-17"
caam.hourly=subset(ibutton_all, between(ibutton_all$Date, as.Date("2020-09-17"), as.Date("2021-01-04"))) 
hist(caam.hourly$Value)

# Getting 5% and 95% quantile of mean temp for each Pop
germ.pheno.temp <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
germ.pheno.temp=subset(germ.pheno.temp, germ.pheno.temp$Pop !="blank")


my_quantile <- function(x, probs) {
  tibble(x = quantile(x, probs, na.rm=TRUE), probs = probs)
}

quant.temps=germ.pheno.temp  %>%
  group_by(Pop) %>% 
  summarize(my_quantile(mean.Temp, c(0.05,0.95)))

#### Get mean temperature between re-water date and germ date for Round 2 ####
ibutton_all <- read.csv("Germination.Timing/Formatted.Data/ibutton_round_2.csv")
germ.pheno <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.formatted.csv", row.names = 1)

# Change character dates to actual Date
germ.pheno$germdate=as.Date(germ.pheno$germdate)
germ.pheno$water.date=as.Date(germ.pheno$water.date)
ibutton_all$Date=as.Date(ibutton_all$Date)

for(i in 1:dim(germ.pheno)[1]){ 
  data=subset(ibutton_all, between(ibutton_all$Date, germ.pheno[i,18], germ.pheno[i,7])) 
  mean=mean(data$temp, na.rm = TRUE)
  germ.pheno[i,21]=mean
}

summary(germ.pheno)
colnames(germ.pheno)[21]="mean.Temp"

write.csv(germ.pheno, file = "Germination.Timing/Formatted.Data/germ.pheno.round.2.temps.ranges.csv")
# a couple of blanks have germ dates. Is an error.

#look at histograms
ggplot(data= germ.pheno, aes(x= mean.Temp)) + geom_histogram(binwidth = 1) + facet_wrap(~ Population)






