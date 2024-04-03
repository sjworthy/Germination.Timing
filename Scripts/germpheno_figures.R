# Code to reproduce figures

#rm(list = ls())
#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges) #for density ridges function
library(cowplot) #if we want to manually panel/arrange plots
library(lubridate) #for formatting date time data (e.g. ibutton data)
library(ape) # for phylogeny
library(phytools) # for phylogeny and dotTree
library(emmeans) # for marginal means
library(boot)
library(jtools)
library(ggprism)
library(lme4)
library(ggbiplot)
library(ggrepel)
library(nationalparkcolors)
library(lmerTest)
library(weathermetrics)
library(gridExtra)

#### Figure 1A: Map ####

sp.points=read.csv("Germination.Timing/Formatted.Data/lat.long.csv")

sp.points = sp.points %>% 
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(species %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(justsp = substr(species, 1,4)) %>%
  mutate(genus = substr(species, 1,1)) %>%
  mutate(fullgenus = recode(genus, "C" = "Caulanthus", "S" = "Streptanthus")) %>%
  mutate(altsp = ifelse(justsp %in% c("CAAN", "CAIN"), paste0(justsp, altpop), justsp))

states <- map_data("state") %>% filter(region == "california")

map <- ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = "gray", color="black") +
  coord_quickmap(xlim = c(-124, -114.5), ylim = c(32.5, 42))+
  geom_point(data = sp.points,
             aes(x = long, y = lat),
             size = 3.5) +
  geom_label_repel(data = sp.points, 
                   aes(x = long, y = lat, label = species)) +
  theme_classic(base_size = 15 )+
  labs(x="Longitude",y="Latitude") + theme(legend.position = "none")

map

#ggsave("Germination.Timing/Plots/map_no color.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/map_no color.png", height = 10, width = 12)

#### Figure 1B: PC biplot of climate ####
# Get climate data for pops and species from Flint and Flint 2014
# data here: https://ucdavis.app.box.com/file/949830830865?s=gmpbjzprzbswqukbal1gf2e8bxtbykwg
# subset for year 1991-2016, 25 years

# do 1991 - 2015 because 2016 only has up until month 9 (25 years of values)
flint.data.mothly = read.csv("~/Library/CloudStorage/Box-Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv") %>% 
 #read.csv("C:/Users/Jenny/Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv") %>% 
  filter(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                 "STGL1","STIN","STPO1","STTO-BH")) %>%
 filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>%
  filter(clim_month %in% c(9,10,11,12)) %>%
  mutate(tmean=(tmin+tmax)/2) %>%
  group_by(id, clim_month) %>%
  dplyr::summarize(CWD = mean(cwd),PPT_CV =(sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax), Tmean_SD=sd(tmean),
                   Tmin = mean(tmin), Tmax = mean(tmax),Tmean=mean(tmean)) %>%
  mutate(genus = "S") %>%
  mutate(genus= ifelse(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4"), "C", genus))

flint.data.seasonal = read.csv("~/Library/CloudStorage/Box-Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv") %>% 
  #read.csv("C:/Users/Jenny/Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv") %>%
  filter(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                   "STGL1","STIN","STPO1","STTO-BH")) %>%
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>% 
  filter(clim_month %in% c(9,10,11,12)) %>%
  mutate(tmean=(tmin+tmax)/2) %>%
  group_by(id) %>%
  dplyr::summarize(CWD = mean(cwd),PPT_CV =(sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax), Tmean_SD=sd(tmean),
                   Tmin = mean(tmin), Tmax = mean(tmax),Tmean=mean(tmean))%>%
  mutate(genus = "S") %>%
  mutate(genus= ifelse(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4"), "C", genus))

# seasonal correlations
cor(flint.data.seasonal[,c(2:10)])
# cwd and ppt r = -0.93

#monthly cor
flint.data.mothly %>%
  ungroup() %>% #for some reason this drops ID column, maybe from grouping in pipe above
  select(CWD, PPT, Tmin, Tmax, Tmean) %>%
  cor(.)
#all of the monthly variables are pretty highly correlated- 0.65 < |r|< 0.99

# pca space
# all variables
pc.monthly = prcomp(flint.data.mothly[c(3:11)], scale  = TRUE, center = TRUE)
ggbiplot(pc.monthly, labels = flint.data.mothly$clim_month)
ggbiplot(pc.monthly, labels = flint.data.mothly$id, groups = flint.data.mothly$genus, ellipse = TRUE)

pc.seasonal = prcomp(flint.data.seasonal[c(2:10)], scale  = TRUE, center = TRUE)
ggbiplot(pc.seasonal, labels = flint.data.seasonal$id, groups = flint.data.seasonal$genus, ellipse = TRUE)+
  theme_classic()
pc.seasonal
#first pc is driven by mean ppt and standard dev of temp.  2nd pc is really driven by sd.tmin.  

#May remove tmean, since it is calc'd by averaging Tmax/tmin (and my spot check is consistent with this)
#flint.data.seasonal %>% mutate(tmean2 = (tmin + Tmax)/2)

# remove tmean from seasonal because of correlation
pc.seasonal.2a = prcomp(flint.data.seasonal[c(3:6, 8:9 )], scale  = TRUE, center = TRUE)
ggbiplot(pc.seasonal.2a, labels = flint.data.seasonal$id, groups = flint.data.seasonal$genus, ellipse = TRUE)+
  theme_classic()
pc.seasonal.2a #first PC looks like CWD, ppt, and sdtmax, 2nd is var in ppt and var in tmin, which is cool, because CACO and CAAM and CAIN4 have highest

#screeplot
plot(pc.seasonal.2a, type = "lines")

pc.seasonal.2a$sdev #may want to plot 3rd PC too

pcvarexplained = tibble(var_explained = (pc.seasonal.2a$sdev^2) / (sum(pc.seasonal.2a$sdev^2))) %>%
                mutate(pc = seq(1,length(var_explained), 1)) %>%
                mutate(cumvarexplained = cumsum(var_explained))
pcvarexplained

#Clean up 2D PCA
pc_data = data.frame(pc.seasonal.2a$x)

pops_clim_pc = cbind(flint.data.seasonal, pc_data)

pops_clim_pc$ID.2=c("CAAM","CAAN1","CAAN2","CACO","CAIN3","CAIN4","STBR","STDI","STDR",
                    "STGL","STIN","STPO","STTO")

 pops_clim_pc = pops_clim_pc %>%
               mutate(justsp = substr(id, 1,4)) %>% 
               mutate(altpop = 1) %>%
               mutate(altpop= ifelse(id %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
               mutate(altsp= justsp) %>%
               mutate(altsp = ifelse(justsp %in% c("CAAN", "CAIN"), paste0(justsp, altpop), justsp)) %>%
               mutate(fullgenus = recode(genus, "C" = "Caulanthus", "S" = "Streptanthus")) 

#try biplot
#define annotations to axis labels for PCA axes
 text_x_low = textGrob("Hot & dry", gp=gpar(fontsize = 9,col="gray35"))
 text_x_high = textGrob("Cool & wet", gp=gpar(fontsize = 9,col="gray35"))
 text_y_low = textGrob("Low CV(PPT)", gp=gpar(fontsize = 9,col="gray35"), rot=90)
 text_y_high = textGrob("High CV(PPT)", gp=gpar(fontsize = 9,col="gray35"), rot=90)
 
seasonalbiplot = ggbiplot(pc.seasonal.2a,  varname.adjust = 1.1) +
  geom_point() + 
  geom_text(vjust = 1.2, hjust = 1.1, label =  pops_clim_pc$ID.2) +
  theme_classic() + 
  labs(x = paste0("Standardized PC1\n (", round(pcvarexplained$var_explained[1]*100,1), "% explained var.)"), 
      y = paste0("Standardized PC2\n (", round(pcvarexplained$var_explained[2]*100,1), "% explained var.)") #,
     # color = "Genus" #add back for colored points
     ) + xlim(-2,2) + ylim(-1.1,2.3)+
 annotation_custom(text_x_high,xmin=1.8,xmax=2,ymin=-1.6,ymax=-1.6) + #here I just fiddled with coordinates to get titles where we want them
  annotation_custom(text_x_low,xmin=-2,xmax=-1.9,ymin=-1.6,ymax=-1.6) + 
  annotation_custom(text_y_high,xmin=-2.5,xmax=-2.5,ymin=1.8,ymax=2) + 
  annotation_custom(text_y_low,xmin=-2.5,xmax=-2.5,ymin=-1.3,ymax=-1.1) + 
  coord_cartesian( clip="off") #this keeps it from clipping off the stuff outside the plot

seasonalbiplot

#some code to customize biplot
#get indexes for what to change
seg <- which(sapply(seasonalbiplot$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(seasonalbiplot$layers, function(x) class(x$geom)[1] == 'GeomText'))

#change segment and text labels for loading arrows to gray
seasonalbiplot$layers[[seg]]$aes_params$colour <- 'gray60'
seasonalbiplot$layers[[txt[1]]]$aes_params$colour <- 'gray60'

fig1b = seasonalbiplot 
fig1b

#ggsave("Germination.Timing/Plots/climatePCA_fig1b.pdf", width = 5, height = 5)
#ggsave("Germination.Timing/Plots/climatePCA_fig1b.jpg", width = 5, height = 5)

#put panel together

fig1_draft = plot_grid(map, fig1b, labels = c("A.", "B."),ncol= 2, rel_heights = c(1,0.3),
                       rel_heights = c(1,0.3))

plot_grid(map, fig1b, labels = c("A.", "B."), rel_widths = c(0.5,0.5), rel_heights = c(1, 0.3))
fig1_draft

#ggsave("Germination.Timing/Plots/fig1_mapandPCA_no color.jpg", height = 8, width = 12)
# ggsave("Germination.Timing/Plots/fig1_mapandPCA_no color.pdf", height = 8, width = 12)
#cowplot::save_plot( "Germination.Timing/Plots/fig1_mapandPCA_no color.jpg", fig1_draft, base_height = 8, base_width = 12)

#### Figure 2:  Cohort timing and temperatures ####
#see prep iButton data script for manipulations to this dataset

ibuttondat = read.csv("Germination.Timing/Formatted.Data/ibutton_allcohorts.csv")
head(ibuttondat)
str(ibuttondat) #didn't read in dates in date format
summary(ibuttondat)

ibuttondat = ibuttondat %>%
              mutate(Date = as.Date(Date), Date.Time = as.Date(Date.Time))
str(ibuttondat) 
summary(ibuttondat)

ibutton.mean=ibuttondat %>%
  group_by(Date.Time) %>%
  dplyr::summarise(mean.temp=mean(temp),
            sd.temp=sd(temp))

mean(ibutton.mean$mean.temp) # 14.79415
mean(ibutton.mean$sd.temp) # 5.474449
# october 2 temp: 21.33333, sd = 3.383889
# october 30 temp: 13.69271, sd = 4.253652

ggplot(ibutton.mean, aes(x=Date.Time, y=mean.temp)) +
  geom_line()+
  geom_ribbon(aes(ymin=mean.temp-sd.temp, ymax=mean.temp+sd.temp), alpha=0.2)+
  xlab("Date")+
  ylab("Mean Daily Temperature")+
  theme_classic(base_size = 15)+
  geom_vline(xintercept=as.numeric(ibutton.mean[1,1]), linetype="dashed")+ # 9/17/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[16,1]), linetype="dashed")+ # 10/2/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[30,1]), linetype="dashed")+ # 10/16/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[44,1]), linetype="dashed")+ # 10/30/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[58,1]), linetype="dashed")+ # 11/13/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[72,1]), linetype="dashed")+ # 11/27/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[86,1]), linetype="dashed") # 12/11/2020

# subset the temperature data to only go as high as the last date any seed germinated (1/25/2021)
ibutton.mean.2=ibutton.mean[c(1:131),]

Fig2=ggplot(ibutton.mean.2, aes(x=Date.Time, y=mean.temp)) +
  geom_line()+
  geom_ribbon(aes(ymin=mean.temp-sd.temp, ymax=mean.temp+sd.temp), alpha=0.2)+
  xlab("Date")+
  ylab("Mean Daily Temperature  (°C)\n Year One")+
  ylim(0,35)+
  theme_classic(base_size = 15)+
  geom_vline(xintercept=as.numeric(ibutton.mean[1,1]), linetype="dashed", color="#5495CF",linewidth=1)+ # 9/17/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[16,1]), linetype="dashed", color="#847CA3",linewidth=1)+ # 10/2/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[30,1]), linetype="dashed", color="#E45A5A",linewidth=1)+ # 10/16/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[44,1]), linetype="dashed", color="#F4A65E",linewidth=1)+ # 10/30/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[58,1]), linetype="dashed",color="#80792B",linewidth=1)+ # 11/13/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[72,1]), linetype="dashed", color="#F2D56F",linewidth=1)+ # 11/27/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[86,1]), linetype="dashed", color="#359F8B",linewidth=1)+ # 12/11/2020
  geom_vline(xintercept=as.numeric(ibutton.mean[15,1]), color="black",linewidth=1)+ # 10/1/2020 (historical)
  geom_vline(xintercept=as.numeric(ibutton.mean[42,1]), color="gray",linewidth=1) # 10/28/2020 (contemporary)
Fig2

#ggsave("Germination.Timing/Plots/cohort.timing.and.temps.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/cohort.timing.and.temps.png", height = 10, width = 12)

# adding year 2 ibutton data to plot

ibutton.y2 = read.csv("Germination.Timing/Formatted.Data/ibutton_round_2.csv")

ibutton.y2.2 = ibutton.y2 %>%
  mutate(Date = as.Date(Date, "%m/%d/%y"))

ibutton.y2.mean=ibutton.y2.2 %>%
  group_by(Date) %>%
  dplyr::summarise(mean.temp=mean(temp),
                   sd.temp=sd(temp))

# Oct 1 2021: 20.370833, sd = 6.3040923
# Oct 30 2021: 16.645833, sd = 3.7238702

# rain event on Sept 15, 2021
# last day Dec 27, 2021

ibutton.y2.mean.2 = ibutton.y2.mean[c(34:137),]

mean(ibutton.y2.mean.2$mean.temp) # 14.07
sd(ibutton.y2.mean.2$mean.temp) # 4.76

Fig2.1=ggplot(ibutton.y2.mean.2, aes(x=Date, y=mean.temp)) +
  geom_line()+
  geom_ribbon(aes(ymin=mean.temp-sd.temp, ymax=mean.temp+sd.temp), alpha=0.2)+
  xlab("Date")+
  ylab("Mean Daily Temperature  (°C)\n Year Two")+
  theme_classic(base_size = 15)+
  ylim(0,35)+
  geom_vline(xintercept=as.numeric(ibutton.y2.mean.2[1,1]), linetype="dashed", color="#325731",linewidth=1) # 9/15/2020
Fig2.1

# figure with both year 1 and year 2 temperatures (Figure S2A)

ibutton.y2.mean.2$Date.2 = ibutton.y2.mean.2$Date-years(1)

all.years.temps = left_join(ibutton.mean.2,ibutton.y2.mean.2, by = join_by("Date.Time"=="Date.2"))

#output to organize


#write.csv(all.years.temps, file = "./Germination.Timing/Formatted.Data/all.years.temp.csv")

all.years.temps = read.csv("./Germination.Timing/Formatted.Data/all.years.temp.csv") # didn't work

cols <- c("YEAR1"="black","YEAR2"= "orange")
year1.temp = ggplot(ibutton.mean.2, aes(x=Date.Time, y=mean.temp)) +
  geom_line(aes(color = "YEAR1"),linewidth=1)+
  geom_line(data = ibutton.y2.mean.2, aes(x=Date.2, y=mean.temp, color = "YEAR2"), linewidth=1)+
  xlab("Date")+
  ylab("Mean Daily Temperature (°C)")+
  ylim(0,35)+
  theme_classic(base_size = 15)+
  scale_colour_manual(name="Years",values=cols)
year1.temp

#ggsave("Germination.Timing/Plots/all.year.temp.fig.pdf", height = 8, width = 12)
#ggsave("Germination.Timing/Plots/all.year.temp.fig.jpg", height = 8, width = 12)

year1.temp.nolegend = ggplot(ibutton.mean.2, aes(x=Date.Time, y=mean.temp)) +
  geom_line(aes(color = "YEAR1"),linewidth=1)+
  geom_line(data = ibutton.y2.mean.2, aes(x=Date.2, y=mean.temp, color = "YEAR2"), linewidth=1)+
  xlab("Date")+
  ylab("Mean Daily Temperature (°C)")+
  ylim(0,35)+
  theme_classic(base_size = 15)+
  scale_colour_manual(name="Years",values=cols)+
  theme(legend.position = "None")
year1.temp.nolegend


#### Figure 3A: Temperatures that seeds germinated under ####
#density ridge plots that Jenny started

germ.pheno = read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1) # Sam reading in file
summary(germ.pheno)
dim(germ.pheno)

germpheno_2plot = germ.pheno %>%
  filter(Pop != "blank")
dim(germpheno_2plot)
summary(as.factor(germpheno_2plot$Pop))
summary(germpheno_2plot)

# Adding species name to the dataframe
germpheno_2plot$Species=germpheno_2plot$Pop
germpheno_2plot$Species <- gsub("CAAN1", "CAAN", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("CAAN2", "CAAN", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("CAIN3", "CAIN", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("CAIN4", "CAIN", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("CAAM-GB", "CAAM", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("STTO-BH", "STTO", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("STPO1", "STPO", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("STGL1", "STGL", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("STDR2", "STDR", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("STBR3", "STBR", germpheno_2plot$Species)
germpheno_2plot$Species <- gsub("CACO1", "CACO", germpheno_2plot$Species)

#Jenny adding pop
germpheno_2plot = germpheno_2plot %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop))

table(germpheno_2plot$Pop, germpheno_2plot$altpop)

# read in Round 2 data and add as additional cohort
germ.pheno.R2 = read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.temps.ranges.csv", row.names = 1)
summary(germ.pheno.R2)
dim(germ.pheno.R2)

germphenoR2_2plot = germ.pheno.R2 %>%
  filter(Population != "blank")
dim(germphenoR2_2plot)
summary(as.factor(germphenoR2_2plot$Population))
summary(germphenoR2_2plot)

# Adding species name to the dataframe
germphenoR2_2plot$Species=germphenoR2_2plot$Population
germphenoR2_2plot$Species <- gsub("CAAN1", "CAAN", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("CAAN2", "CAAN", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("CAIN3", "CAIN", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("CAIN4", "CAIN", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("CAAM-GB", "CAAM", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("STTO-BH", "STTO", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("STPO1", "STPO", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("STGL1", "STGL", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("STDR2", "STDR", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("STBR3", "STBR", germphenoR2_2plot$Species)
germphenoR2_2plot$Species <- gsub("CACO1", "CACO", germphenoR2_2plot$Species)

# add cohort number to match other df
colnames(germphenoR2_2plot)[2]="Old.Cohort"
germphenoR2_2plot$Cohort=8
colnames(germphenoR2_2plot)[1]="Pop"

# make a new data frame with all data

germ.pheno.all.rounds=full_join(germpheno_2plot, germphenoR2_2plot)

cols=c("#5495CF","#847CA3", "#E45A5A", "#F4A65E", "#80792B", "#F2D56F", "#359F8B", "gray")

#can panel by cohorts
# using geom_density_ridges2 to draw the complete line for each cohort
fig3a=ggplot(germ.pheno.all.rounds, aes(x = mean.Temp, y = as.factor(Cohort), group = as.factor(Cohort), fill = as.factor(Cohort))) + 
  ylab("Rainfall Onset Date")  + xlab("Mean temperature (°C) during germination") +
  geom_density_ridges2(aes(y = reorder(as.factor(Cohort), desc(as.factor(Cohort)))), show.legend = FALSE)+
  scale_y_discrete(breaks=c("1","2","3","4","5","6","7","8"),labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec","15-Sept \n (Year Two)"))+
  theme_classic(base_size = 15)+
  scale_fill_manual(values = cols)
fig3a

#ggsave("Germination.Timing/Plots/temps.experienced.cohort.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/temps.experienced.cohort.png", height = 10, width = 12)

FigureS2=plot_grid(year1.temp.nolegend,fig3a, labels = c("A","B"))

#ggsave("Germination.Timing/Plots/FigureS2.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/FigureS2.png", height = 10, width = 12)

#### Figure 3B: Germination Fraction for Rounds 1 and 2 #####
# calculating germination proportion

# Round 1
germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")

germ.pheno.R2.all=read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.formatted.csv", row.names = 1)
germ.pheno.R2=subset(germ.pheno.R2.all, germ.pheno.R2.all$Population !="blank")

total.planted=germ.pheno %>%
  group_by(Pop,Cohort) %>%
  summarise(total.count=sum(planted))

total.germ.R1=germ.pheno %>%
  group_by(Pop,Cohort)%>%
  summarise(germ.R1=sum(germinated))

total.germ.R2=germ.pheno.R2 %>%
  group_by(Population,Cohort)%>%
  summarise(germ.R2=sum(germinated))

total.planted.R2=germ.pheno.R2 %>%
  group_by(Population,Cohort)%>%
  summarise(planted.R2=sum(planted))

total.planted=as.data.frame(total.planted)
total.planted[,4]=total.germ.R1$germ.R1
colnames(total.planted)[4]="germ.R1"
total.planted[,5]=total.germ.R2$germ.R2
colnames(total.planted)[5]="germ.R2"
total.planted[,6]=total.planted.R2$planted.R2
colnames(total.planted)[6]="planted.R2"
total.planted[,7]=total.planted$total.count-total.planted$germ.R1
colnames(total.planted)[7]="math"
total.planted[,8]=total.planted$germ.R1/total.planted$total.count
colnames(total.planted)[8]="proportion.R1"
total.planted[,9]=total.planted$germ.R2/total.planted$total.count
colnames(total.planted)[9]="proportion.R2"
total.planted[,10]=(total.planted$germ.R1+total.planted$germ.R2)/total.planted$total.count
colnames(total.planted)[10]="total.proportion"

# made a new data table in long format

cohort.proportion.data=read.csv("Germination.Timing/Formatted.Data/proportion.rounds.cohort.plot.csv")
cohort.proportion.data$Pop=cohort.proportion.data$Species # adding Pop column

cohort.proportion.data.2 = cohort.proportion.data %>%
  mutate(Species = as.factor(substr(Pop, 1, 4))) %>%
  mutate(Species = recode(Species, "caan" = "CAAN", "cain" = "CAIN", "caam" = "CAAM", "caco" = "CACO", 
                          "stto" = "STTO", "stdr" = "STDR", "stpo" = "STPO", "stbr" = "STBR", 
                          "stdi" = "STDI", "stgl" = "STGL", "stin" = "STIN")) %>% 
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# plot with pops separate
fig3b=ggplot(cohort.proportion.data.2, aes(fill=Year, y=Germination.Proportion, x=factor(Cohort))) + 
  geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
  theme_classic(base_size = 15)+scale_fill_manual(values = c(16,1))+
  facet_wrap(.~phy_order, ncol=4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))+
  labs(x="Rainfall Onset Date", y="Germination Fraction")
fig3b

#ggsave("Germination.Timing/Plots/fig3b_legend.pdf", height = 8, width = 12)

#Fig3_draft = plot_grid(fig3a,fig3b, labels = c("A.", "B."), rel_widths = c(0.9, 1),ncol = 2)
#fig3_draft
# ggsave("Germination.Timing/Plots/fig3.2.pdf", height = 8, width = 12)

Fig2_draft = plot_grid(Fig2,fig3b, labels = c("A.", "B."), rel_widths = c(0.9, 1),ncol = 2)
Fig2_draft
#ggsave("Germination.Timing/Plots/fig2.revised.pdf", height = 8, width = 12)
#ggsave("Germination.Timing/Plots/fig2.revised.jpg", height = 8, width = 12)

fig.2.revised = grid.arrange(arrangeGrob(Fig2,Fig2.1), fig3b, ncol = 2)
fig.2.revised

Fig2.1
#ggsave("Germination.Timing/Plots/year2temp.pdf", height = 8, width = 12)
#ggsave("Germination.Timing/Plots/year2temp.jpg", height = 8, width = 12)


# ggsave("Germination.Timing/Plots/fig3_new.pdf", height = 8, width = 12)

#ggsave("Germination.Timing/Plots/germproport.rounds.separate.pops.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/germproport.rounds.separate.pops.png", height = 10, width = 12)

#### Figure 4: Germination Fraction~Cohort plots ####
#GLM of germination proportion~Cohort
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

# get predicted values
caam.df=caam.germ.proport[,c(1,55)]
colnames(caam.df)[2]="mean.temp"
caam.df=caam.df %>%
  mutate(Pop="CAAM")%>%
  mutate(logitproport = predict(caam.mod.glm.cohort, re.form=NA))%>%
  mutate(proportion = inv.logit(logitproport))

caan1.df=caan1.germ.proport[,c(1,43)]
colnames(caan1.df)[2]="mean.temp"
caan1.df=caan1.df %>%
  mutate(Pop="CAAN1")%>%
  mutate(logitproport = predict(caan1.mod.glm.cohort, re.form=NA))%>%
  mutate(proportion = inv.logit(logitproport))

caan2.df=caan2.germ.proport[,c(1,44)]
colnames(caan2.df)[2]="mean.temp"
caan2.df=caan2.df %>%
  mutate(Pop="CAAN2")%>%
  mutate(logitproport = predict(caan2.mod.glm.cohort, re.form=NA))%>%
  mutate(proportion = inv.logit(logitproport))

caco.df=caco.germ.proport[,c(1,45)]
colnames(caco.df)[2]="mean.temp"
caco.df=caco.df %>%
  mutate(Pop="CACO")%>%
  mutate(logitproport = predict(caco1.mod.glm.cohort, re.form=NA))%>%
  mutate(proportion = inv.logit(logitproport))

cain3.df=cain3.germ.proport[,c(1,46)]
colnames(cain3.df)[2]="mean.temp"
cain3.test = predict(cain3.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(cain3.df$mean.temp), max(cain3.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
cain3.df$logitproport = cain3.test$fit
cain3.df$logit.se = cain3.test$se.fit
cain3.df=cain3.df %>%
  mutate(Pop="CAIN3")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

cain4.df=germ.proport.temp.block[,c(1,47)]
colnames(cain4.df)[2]="mean.temp"
cain4.test = predict(cain4.mod.glm.cohort, re.form = NA,
                     new.data = seq(min(cain4.df$mean.temp), max(cain4.df$mean.temp), length.out = 50),
                     type = "link", se.fit = TRUE)
cain4.df$logitproport = cain4.test$fit
cain4.df$logit.se = cain4.test$se.fit
cain4.df=cain4.df %>%
  mutate(Pop="CAIN4")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stbr.df=germ.proport.temp.block[,c(1,48)]
colnames(stbr.df)[2]="mean.temp"
stbr.test = predict(stbr3.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stbr.df$mean.temp), max(stbr.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stbr.df$logitproport = stbr.test$fit
stbr.df$logit.se = stbr.test$se.fit
stbr.df=stbr.df %>%
  mutate(Pop="STBR")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stdi.df=germ.proport.temp.block[,c(1,49)]
colnames(stdi.df)[2]="mean.temp"
stdi.test = predict(stdi.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stdi.df$mean.temp), max(stdi.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stdi.df$logitproport = stdi.test$fit
stdi.df$logit.se = stdi.test$se.fit
stdi.df=stdi.df %>%
  mutate(Pop="STDI")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stdr.df=stdr2.germ.proport[,c(1,50)]
colnames(stdr.df)[2]="mean.temp"
stdr.test = predict(stdr2.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stdr.df$mean.temp), max(stdr.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stdr.df$logitproport = stdr.test$fit
stdr.df$logit.se = stdr.test$se.fit
stdr.df=stdr.df %>%
  mutate(Pop="STDR")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stgl.df=stgl1.germ.proport[,c(1,51)]
colnames(stgl.df)[2]="mean.temp"
stgl.test = predict(stgl1.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stgl.df$mean.temp), max(stgl.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stgl.df$logitproport = stgl.test$fit
stgl.df$logit.se = stgl.test$se.fit
stgl.df=stgl.df %>%
  mutate(Pop="STGL")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stin.df=germ.proport.temp.block[,c(1,52)]
colnames(stin.df)[2]="mean.temp"
stin.test = predict(stin.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stin.df$mean.temp), max(stin.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stin.df$logitproport = stin.test$fit
stin.df$logit.se = stin.test$se.fit
stin.df=stin.df %>%
  mutate(Pop="STIN")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stpo.df=germ.proport.temp.block[,c(1,53)]
colnames(stpo.df)[2]="mean.temp"
stpo.test = predict(stpo1.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stpo.df$mean.temp), max(stpo.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stpo.df$logitproport = stpo.test$fit
stpo.df$logit.se = stpo.test$se.fit
stpo.df=stpo.df %>%
  mutate(Pop="STPO")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

stto.df=germ.proport.temp.block[,c(1,54)]
colnames(stto.df)[2]="mean.temp"
stto.test = predict(stto.mod.glm.cohort, re.form = NA,
                    new.data = seq(min(stto.df$mean.temp), max(stto.df$mean.temp), length.out = 50),
                    type = "link", se.fit = TRUE)
stto.df$logitproport = stto.test$fit
stto.df$logit.se = stto.test$se.fit
stto.df=stto.df %>%
  mutate(Pop="STTO")%>%
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(upper = inv.logit(logitproport + logit.se)) %>%
  mutate(lower = inv.logit(logitproport - logit.se))

all.contain.cohort=list(caam.df,caan1.df,caan2.df,caco.df,cain3.df,cain4.df,stbr.df,stdi.df,stdr.df,stgl.df,stin.df,
                        stpo.df,stto.df) %>% reduce(full_join)

all.contain.cohort = all.contain.cohort %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# Testing for difference between cohort 2 and 4
# cohort has to be factor to compare with emmeans

y=cbind(caam.germ.proport$caam.count, (caam.germ.proport$caam.sample-caam.germ.proport$caam.count))
caam.germ.proport$cohort.factor=as.factor(caam.germ.proport$Cohort)
caam.mod.glm.cohort.fact=glm(y~cohort.factor, data=caam.germ.proport, family = binomial)

y=cbind(caan1.germ.proport$caan1.count, (caan1.germ.proport$caan1.sample-caan1.germ.proport$caan1.count))
caan1.germ.proport$cohort.factor=as.factor(caan1.germ.proport$Cohort)
caan1.mod.glm.cohort.fact=glm(y~cohort.factor, data=caan1.germ.proport, family = binomial)

y=cbind(caan2.germ.proport$caan2.count, (caan2.germ.proport$caan2.sample-caan2.germ.proport$caan2.count))
caan2.germ.proport$cohort.factor=as.factor(caan2.germ.proport$Cohort)
caan2.mod.glm.cohort.fact=glm(y~cohort.factor, data=caan2.germ.proport, family = binomial)

y=cbind(caco.germ.proport$caco1.count, (caco.germ.proport$caco1.sample-caco.germ.proport$caco1.count))
caco.germ.proport$cohort.factor=as.factor(caco.germ.proport$Cohort)
caco1.mod.glm.cohort.fact=glm(y~cohort.factor, data=caco.germ.proport, family = binomial)

y=cbind(cain3.germ.proport$cain3.count, (cain3.germ.proport$cain3.sample-cain3.germ.proport$cain3.count))
cain3.germ.proport$cohort.factor=as.factor(cain3.germ.proport$Cohort)
cain3.mod.glm.cohort.fact=glm(y~cohort.factor, data=cain3.germ.proport, family = binomial) 

y=cbind(germ.proport.temp.block$cain4.count, (germ.proport.temp.block$cain4.sample-germ.proport.temp.block$cain4.count))
germ.proport.temp.block$cohort.factor=as.factor(germ.proport.temp.block$Cohort)
cain4.mod.glm.cohort.fact=glm(y~cohort.factor, data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stbr3.count, (germ.proport.temp.block$stbr3.sample-germ.proport.temp.block$stbr3.count))
stbr3.mod.glm.cohort.fact=glm(y~cohort.factor, data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stdi.count, (germ.proport.temp.block$stdi.sample-germ.proport.temp.block$stdi.count))
stdi.mod.glm.cohort.fact=glm(y~cohort.factor, data=germ.proport.temp.block, family = binomial)

y=cbind(stdr2.germ.proport$stdr2.count, (stdr2.germ.proport$stdr2.sample-stdr2.germ.proport$stdr2.count))
stdr2.germ.proport$cohort.factor=as.factor(stdr2.germ.proport$Cohort)
stdr2.mod.glm.cohort.fact=glm(y~cohort.factor, data=stdr2.germ.proport, family = binomial)

y=cbind(stgl1.germ.proport$stgl1.count, (stgl1.germ.proport$stgl1.sample-stgl1.germ.proport$stgl1.count))
stgl1.germ.proport$cohort.factor=as.factor(stgl1.germ.proport$Cohort)
stgl1.mod.glm.cohort.fact=glm(y~cohort.factor, data=stgl1.germ.proport, family = binomial)

y=cbind(germ.proport.temp.block$stin.count, (germ.proport.temp.block$stin.sample-germ.proport.temp.block$stin.count))
stin.mod.glm.cohort.fact=glm(y~cohort.factor, data=germ.proport.temp.block, family = binomial) 

y=cbind(germ.proport.temp.block$stpo1.count, (germ.proport.temp.block$stpo1.sample-germ.proport.temp.block$stpo1.count))
stpo1.mod.glm.cohort.fact=glm(y~cohort.factor, data=germ.proport.temp.block, family = binomial)

y=cbind(germ.proport.temp.block$stto.count, (germ.proport.temp.block$stto.sample-germ.proport.temp.block$stto.count))
stto.mod.glm.cohort.fact=glm(y~cohort.factor,data=germ.proport.temp.block, family = binomial)

# predicting values
caam.df.fact=caam.germ.proport[,c(5,6,56)]
caam.df.fact$y=caam.df.fact$caam.count/caam.df.fact$caam.sample
caam.proport.factor.predict=as.data.frame(effect_plot(caam.mod.glm.cohort.fact, pred = cohort.factor, data=caam.df.fact, interval = TRUE)$data)
caam.proport.factor.predict = caam.proport.factor.predict%>%
  mutate(Pop="CAAM")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

caan1.df.fact=caan1.germ.proport[,c(7,8,56)]
caan1.df.fact$y=caan1.df.fact$caan1.count/caan1.df.fact$caan1.sample
caan1.proport.factor.predict=as.data.frame(effect_plot(caan1.mod.glm.cohort.fact, pred = cohort.factor, data=caan1.df.fact, interval = TRUE)$data)
caan1.proport.factor.predict = caan1.proport.factor.predict%>%
  mutate(Pop="CAAN1")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

caan2.df.fact=caan2.germ.proport[,c(10,11,56)]
caan2.df.fact$y=caan2.df.fact$caan2.count/caan2.df.fact$caan2.sample
caan2.proport.factor.predict=as.data.frame(effect_plot(caan2.mod.glm.cohort.fact, pred = cohort.factor, data=caan2.df.fact, interval = TRUE)$data)
caan2.proport.factor.predict = caan2.proport.factor.predict%>%
  mutate(Pop="CAAN2")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

caco.df.fact=caco.germ.proport[,c(13,14,56)]
caco.df.fact$y=caco.df.fact$caco1.count/caco.df.fact$caco1.sample
caco.proport.factor.predict=as.data.frame(effect_plot(caco1.mod.glm.cohort.fact, pred = cohort.factor, data=caco.df.fact, interval = TRUE)$data)
caco.proport.factor.predict = caco.proport.factor.predict%>%
  mutate(Pop="CACO")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

cain3.df.fact=cain3.germ.proport[,c(16,17,56)]
cain3.df.fact$y=cain3.df.fact$cain3.count/cain3.df.fact$cain3.sample
cain3.proport.factor.predict=as.data.frame(effect_plot(cain3.mod.glm.cohort.fact, pred = cohort.factor, data=cain3.df.fact, interval = TRUE)$data)
cain3.proport.factor.predict = cain3.proport.factor.predict%>%
  mutate(Pop="CAIN3")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

cain4.df.fact=germ.proport.temp.block[,c(19,20,56)]
cain4.df.fact$y=cain4.df.fact$cain4.count/cain4.df.fact$cain4.sample
cain4.proport.factor.predict=as.data.frame(effect_plot(cain4.mod.glm.cohort.fact, pred = cohort.factor, data=cain4.df.fact, interval = TRUE)$data)
cain4.proport.factor.predict = cain4.proport.factor.predict%>%
  mutate(Pop="CAIN4")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stbr.df.fact=germ.proport.temp.block[,c(22,23,56)]
stbr.df.fact$y=stbr.df.fact$stbr3.count/stbr.df.fact$stbr3.sample
stbr.proport.factor.predict=as.data.frame(effect_plot(stbr3.mod.glm.cohort.fact, pred = cohort.factor, data=stbr.df.fact, interval = TRUE)$data)
stbr.proport.factor.predict = stbr.proport.factor.predict%>%
  mutate(Pop="STBR")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stdi.df.fact=germ.proport.temp.block[,c(25,26,56)]
stdi.df.fact$y=stdi.df.fact$stdi.count/stdi.df.fact$stdi.sample
stdi.proport.factor.predict=as.data.frame(effect_plot(stdi.mod.glm.cohort.fact, pred = cohort.factor, data=stdi.df.fact, interval = TRUE)$data)
stdi.proport.factor.predict = stdi.proport.factor.predict%>%
  mutate(Pop="STDI")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stdr.df.fact=stdr2.germ.proport[,c(28,29,56)]
stdr.df.fact$y=stdr.df.fact$stdr2.count/stdr.df.fact$stdr2.sample
stdr.proport.factor.predict=as.data.frame(effect_plot(stdr2.mod.glm.cohort.fact, pred = cohort.factor, data=stdr.df.fact, interval = TRUE)$data)
stdr.proport.factor.predict = stdr.proport.factor.predict%>%
  mutate(Pop="STDR")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stgl.df.fact=stgl1.germ.proport[,c(31,32,56)]
stgl.df.fact$y=stgl.df.fact$stgl1.count/stgl.df.fact$stgl1.sample
stgl.proport.factor.predict=as.data.frame(effect_plot(stgl1.mod.glm.cohort.fact, pred = cohort.factor, data=stgl.df.fact, interval = TRUE)$data)
stgl.proport.factor.predict = stgl.proport.factor.predict%>%
  mutate(Pop="STGL")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stin.df.fact=germ.proport.temp.block[,c(34,35,56)]
stin.df.fact$y=stin.df.fact$stin.count/stin.df.fact$stin.sample
stin.proport.factor.predict=as.data.frame(effect_plot(stin.mod.glm.cohort.fact, pred = cohort.factor, data=stin.df.fact, interval = TRUE)$data)
stin.proport.factor.predict = stin.proport.factor.predict%>%
  mutate(Pop="STIN")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stpo.df.fact=germ.proport.temp.block[,c(37,38,56)]
stpo.df.fact$y=stpo.df.fact$stpo1.count/stpo.df.fact$stpo1.sample
stpo.proport.factor.predict=as.data.frame(effect_plot(stpo1.mod.glm.cohort.fact, pred = cohort.factor, data=stpo.df.fact, interval = TRUE)$data)
stpo.proport.factor.predict = stpo.proport.factor.predict%>%
  mutate(Pop="STPO")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

stto.df.fact=germ.proport.temp.block[,c(40,41,56)]
stto.df.fact$y=stto.df.fact$stto.count/stto.df.fact$stto.sample
stto.proport.factor.predict=as.data.frame(effect_plot(stto.mod.glm.cohort.fact, pred = cohort.factor, data=stto.df.fact, interval = TRUE)$data)
stto.proport.factor.predict = stto.proport.factor.predict%>%
  mutate(Pop="STTO")%>%
  mutate(proportion = inv.logit(y))%>%
  mutate(max.proportion = inv.logit(ymax))%>%
  mutate(min.proportion = inv.logit(ymin))

all.proport.factor.predict=list(caam.proport.factor.predict,caan1.proport.factor.predict,caan2.proport.factor.predict,
                                caco.proport.factor.predict,cain3.proport.factor.predict,cain4.proport.factor.predict,
                                stbr.proport.factor.predict,stdi.proport.factor.predict,stdr.proport.factor.predict,
                                stgl.proport.factor.predict,stin.proport.factor.predict,stpo.proport.factor.predict,
                                stto.proport.factor.predict) %>% reduce(full_join)

all.proport.factor.predict = all.proport.factor.predict %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(Species = recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CAAM" = "CAAM", "CACO" = "CACO", 
                          "STTO" = "STTO", "STDR" = "STDR", "STPO" = "STPO", "STBR" = "STBR", 
                          "STDI" = "STDI", "STGL" = "STGL", "STIN" = "STIN")) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# Making each plot individually
# CAAM
caam.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="CAAM")
caam.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="CAAM")

caam.plot=ggplot(caam.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CAAM")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

# CAAN
caan.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Species=="CAAN")
caan.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Species=="CAAN")

caan.plot=ggplot(caan.all.proport.factor.predict, aes(x=cohort.factor, y=y, group=altpop, color=altpop))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75,show.legend = FALSE)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","gray","black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CAAN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

caan.plot.2=ggplot(caan.all.proport.factor.predict, aes(x=cohort.factor, y=y, group=altpop, color=altpop))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","gray","black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CAAN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))
caan.plot.2

# CACO
caco.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="CACO")
caco.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="CACO")

# artificially insert a 0 for cohorts 4 and 5 to make graph match
caco.all.proport.factor.predict[6,1]=NA
caco.all.proport.factor.predict[7,1]=NA
caco.all.proport.factor.predict[6,2]=4
caco.all.proport.factor.predict[7,2]=5

caco.plot=ggplot(caco.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CACO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

# CAIN
cain.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Species=="CAIN")
cain.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Species=="CAIN")

# subset out the two cohorts that are significant
sig.cohort=cain.all.proport.factor.predict[c(9,11),]

cain.plot=ggplot(cain.all.proport.factor.predict, aes(x=cohort.factor, y=y, group=altpop,color=altpop))+
  geom_point(show.legend = FALSE)+
  geom_point(data=sig.cohort, color="red")+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  geom_errorbar(data=sig.cohort, aes(ymin=ymin, ymax=ymax),width=0.75, color="red")+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CAIN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

cain.plot.2=cain.plot+
  geom_line(data=cain.all.contain.cohort, aes(x=Cohort,y=proportion, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=cain.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=cain.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

cain.all.proport.factor.predict.2=cain.all.proport.factor.predict
cain.all.proport.factor.predict.2$altpop=c(3,3,3,3,3,3,3,4,4,4,4,4,4,4)
cain.all.proport.factor.predict.2$altpop=as.factor(cain.all.proport.factor.predict.2$altpop)
cain.all.contain.cohort.2=cain.all.contain.cohort
cain.all.contain.cohort.2$altpop=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
cain.all.contain.cohort.2$altpop=as.factor(cain.all.contain.cohort.2$altpop)

cain.plot.3=ggplot(cain.all.proport.factor.predict.2, aes(x=cohort.factor, y=y, group=altpop,color=altpop))+
  geom_point()+
  geom_point(data=sig.cohort, color="red")+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  geom_errorbar(data=sig.cohort, aes(ymin=ymin, ymax=ymax),width=0.75, color="red")+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("CAIN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

cain.plot.4=cain.plot.3+
  geom_line(data=cain.all.contain.cohort.2, aes(x=Cohort,y=proportion, group=altpop, linetype=altpop), show.legend = FALSE)
cain.plot.4

# STDR
stdr.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STDR")
stdr.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STDR")

stdr.plot=ggplot(stdr.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STDR")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stdr.plot.2=stdr.plot+
  geom_line(data=stdr.all.contain.cohort, aes(x=Cohort,y=proportion, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=stdr.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stdr.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STBR
stbr.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STBR")
stbr.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STBR")

stbr.plot=ggplot(stbr.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STBR")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stbr.plot.2=stbr.plot+
  geom_line(data=stbr.all.contain.cohort, aes(x=Cohort,y=proportion, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=stbr.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stbr.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

# STTO
stto.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STTO")
stto.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STTO")
stto.all.proport.factor.predict$color=c("black","red","black","red","black","black","black")

stto.plot=ggplot(stto.all.proport.factor.predict, aes(x=cohort.factor, y=y,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STTO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stto.plot.2=stto.plot+
  geom_line(data=stto.all.contain.cohort, aes(x=Cohort,y=proportion, color="black"), show.legend = FALSE)+
  geom_line(data=stto.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stto.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STDI
stdi.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STDI")
stdi.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STDI")
stdi.all.proport.factor.predict$color=c("black","red","black","red","black","black","black")

stdi.plot=ggplot(stdi.all.proport.factor.predict, aes(x=cohort.factor, y=y,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STDI")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stdi.plot.2=stdi.plot+
  geom_line(data=stdi.all.contain.cohort, aes(x=Cohort,y=proportion, color="black"), show.legend = FALSE)+
  geom_line(data=stdi.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stdi.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

# STPO
stpo.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STPO")
stpo.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STPO")

stpo.plot=ggplot(stpo.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STPO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stpo.plot.2=stpo.plot+
  geom_line(data=stpo.all.contain.cohort, aes(x=Cohort,y=proportion, color="black"), show.legend = FALSE)+
  geom_line(data=stpo.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stpo.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

# STIN
stin.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STIN")
stin.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STIN")

stin.plot=ggplot(stin.all.proport.factor.predict, aes(x=cohort.factor, y=y))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STIN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stin.plot.2=stin.plot+
  geom_line(data=stin.all.contain.cohort, aes(x=Cohort,y=proportion, color="black"), show.legend = FALSE)+
  geom_line(data=stin.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stin.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

# STGL
stgl.all.proport.factor.predict=subset(all.proport.factor.predict, all.proport.factor.predict$Pop=="STGL")
stgl.all.contain.cohort=subset(all.contain.cohort, all.contain.cohort$Pop=="STGL")
stgl.all.proport.factor.predict$color=c("black","red","black","red","black","black","black")

stgl.plot=ggplot(stgl.all.proport.factor.predict, aes(x=cohort.factor, y=y,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  ylim(0.0,1.0)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Fraction", color="Population")+
  ggtitle("STGL")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stgl.plot.2=stgl.plot+
  geom_line(data=stgl.all.contain.cohort, aes(x=Cohort,y=proportion, color="black"), show.legend = FALSE)+
  geom_line(data=stgl.all.contain.cohort, aes(x=Cohort,y=upper, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stgl.all.contain.cohort, aes(x=Cohort,y=lower, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


plot_grid(stdr.plot.2,stbr.plot.2,stto.plot.2,stdi.plot.2,stpo.plot.2,caam.plot,stin.plot.2,
          stgl.plot.2,caan.plot,caco.plot,cain.plot.2,
          nrow=3, ncol=4)

#ggsave("Germination.Timing/Plots/stdr.plot.pdf", height = 10, width = 12)

#ggsave("Germination.Timing/Plots/Germproport_cohort_panelbyspecies.CI.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/Germproport_cohort_panelbyspecies.CI.png", height = 10, width = 12)

#### Figure 5: Germination Rate~Cohort plots ####

# 1/days2germ2
germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")

germ.pheno[,21]=1/germ.pheno$days2germ2
colnames(germ.pheno)[21]="germ.rate" 

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

# making cohort factor
caam.2$cohort.factor=as.factor(caam.2$Cohort)
caan1.2$cohort.factor=as.factor(caan1.2$Cohort)
caan2.2$cohort.factor=as.factor(caan2.2$Cohort)
caco.2$cohort.factor=as.factor(caco.2$Cohort)
cain3.2$cohort.factor=as.factor(cain3.2$Cohort)
cain4.2$cohort.factor=as.factor(cain4.2$Cohort)
stbr3.2$cohort.factor=as.factor(stbr3.2$Cohort)
stdi.2$cohort.factor=as.factor(stdi.2$Cohort)
stdr2.2$cohort.factor=as.factor(stdr2.2$Cohort)
stgl1.2$cohort.factor=as.factor(stgl1.2$Cohort)
stin.2$cohort.factor=as.factor(stin.2$Cohort)
stpo1.2$cohort.factor=as.factor(stpo1.2$Cohort)
stto.2$cohort.factor=as.factor(stto.2$Cohort)

rate.test=lmer(germ.rate~Cohort*Pop + (1|Block), 
               data=germ.pheno.2)

# Linear models with cohort as continuous to look at trend and as factor to look for differences between cohort 2 and 4
caam.mod.cohort=lm(germ.rate~cohort.factor,
                   data=caam.2)
caam.mod=lm(germ.rate~Cohort,
            data=caam.2)
caan1.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=caan1.2)
caan1.mod=lm(germ.rate~Cohort, 
             data=caan1.2)
caan2.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=caan2.2)
caan2.mod=lm(germ.rate~Cohort, 
             data=caan2.2)
caco.mod.cohort=lm(germ.rate~cohort.factor, 
                   data=caco.2)
caco.mod=lm(germ.rate~Cohort, 
            data=caco.2)
cain3.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=cain3.2)
cain3.mod=lm(germ.rate~Cohort,
             data=cain3.2)
cain4.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=cain4.2)
cain4.mod=lm(germ.rate~Cohort,
             data=cain4.2)
stbr3.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=stbr3.2)
stbr3.mod=lm(germ.rate~Cohort, 
             data=stbr3.2)
stdi.mod.cohort=lm(germ.rate~cohort.factor, 
                   data=stdi.2)
stdi.mod=lm(germ.rate~Cohort, 
            data=stdi.2)
stdr2.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=stdr2.2)
stdr2.mod=lm(germ.rate~Cohort, 
             data=stdr2.2)
stgl1.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=stgl1.2)
stgl1.mod=lm(germ.rate~Cohort, 
             data=stgl1.2)
stin.mod.cohort=lm(germ.rate~cohort.factor, 
                   data=stin.2)
stin.mod=lm(germ.rate~Cohort, 
            data=stin.2)
stpo1.mod.cohort=lm(germ.rate~cohort.factor, 
                    data=stpo1.2)
stpo1.mod=lm(germ.rate~Cohort, 
             data=stpo1.2)
stto.mod.cohort=lm(germ.rate~cohort.factor, 
                   data=stto.2)
stto.mod=lm(germ.rate~Cohort, 
            data=stto.2)

# predict values for all species to get one dataframe to use facet_wrap instead of individual plots

# predict linear relationship rate~ continuous value cohort
caam.new.data=predict(caam.mod, interval = "confidence")
caam.2 = cbind(caam.2, caam.new.data)
caan1.new.data=predict(caan1.mod, interval = "confidence")
caan1.2 = cbind(caan1.2, caan1.new.data)
caan2.new.data=predict(caan2.mod, interval = "confidence")
caan2.2 = cbind(caan2.2, caan2.new.data)
caco.new.data=predict(caco.mod, interval = "confidence")
caco.2 = cbind(caco.2, caco.new.data)
cain3.new.data=predict(cain3.mod, interval = "confidence")
cain3.2 = cbind(cain3.2, cain3.new.data)
cain4.new.data=predict(cain4.mod, interval = "confidence")
cain4.2 = cbind(cain4.2, cain4.new.data)
stbr3.new.data=predict(stbr3.mod, interval = "confidence")
stbr3.2 = cbind(stbr3.2, stbr3.new.data)
stdi.new.data=predict(stdi.mod, interval = "confidence")
stdi.2 = cbind(stdi.2, stdi.new.data)
stdr2.new.data=predict(stdr2.mod, interval = "confidence")
stdr2.2 = cbind(stdr2.2, stdr2.new.data)
stgl1.new.data=predict(stgl1.mod, interval = "confidence")
stgl1.2 = cbind(stgl1.2, stgl1.new.data)
stin.new.data=predict(stin.mod, interval = "confidence")
stin.2 = cbind(stin.2, stin.new.data)
stpo1.new.data=predict(stpo1.mod, interval = "confidence")
stpo1.2 = cbind(stpo1.2, stpo1.new.data)
stto.new.data=predict(stto.mod, interval = "confidence")
stto.2 = cbind(stto.2, stto.new.data)

all.contin.data=list(caam.2,caan1.2,caan2.2,caco.2,cain3.2,cain4.2,stbr3.2,stdi.2,stdr2.2,stgl1.2,stin.2,stpo1.2,stto.2) %>% reduce(full_join)

all.contin.data = all.contin.data %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(Species = recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CAAM-GB" = "CAAM", "CACO1" = "CACO", 
                          "STTO-BH" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                          "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# predict linear relationship rate~ factor value cohort
caam.cohort.factor.predict=as.data.frame(effect_plot(caam.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
caam.cohort.factor.predict$Pop="CAAM"
caan1.cohort.factor.predict=as.data.frame(effect_plot(caan1.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
caan1.cohort.factor.predict$Pop="CAAN1"
caan2.cohort.factor.predict=as.data.frame(effect_plot(caan2.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
caan2.cohort.factor.predict$Pop="CAAN2"
caco.cohort.factor.predict=as.data.frame(effect_plot(caco.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
caco.cohort.factor.predict$Pop="CACO"
cain3.cohort.factor.predict=as.data.frame(effect_plot(cain3.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
cain3.cohort.factor.predict$Pop="CAIN3"
cain4.cohort.factor.predict=as.data.frame(effect_plot(cain4.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
cain4.cohort.factor.predict$Pop="CAIN4"
stbr3.cohort.factor.predict=as.data.frame(effect_plot(stbr3.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stbr3.cohort.factor.predict$Pop="STBR"
stdi.cohort.factor.predict=as.data.frame(effect_plot(stdi.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stdi.cohort.factor.predict$Pop="STDI"
stdr2.cohort.factor.predict=as.data.frame(effect_plot(stdr2.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stdr2.cohort.factor.predict$Pop="STDR"
stgl1.cohort.factor.predict=as.data.frame(effect_plot(stgl1.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stgl1.cohort.factor.predict$Pop="STGL"
stin.cohort.factor.predict=as.data.frame(effect_plot(stin.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stin.cohort.factor.predict$Pop="STIN"
stpo1.cohort.factor.predict=as.data.frame(effect_plot(stpo1.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stpo1.cohort.factor.predict$Pop="STPO"
stto.cohort.factor.predict=as.data.frame(effect_plot(stto.mod.cohort, pred = cohort.factor, interval = TRUE)$data)
stto.cohort.factor.predict$Pop="STTO"

all.cohort.factor.predict=list(caam.cohort.factor.predict,caan1.cohort.factor.predict,caan2.cohort.factor.predict,
                               caco.cohort.factor.predict,cain3.cohort.factor.predict,cain4.cohort.factor.predict,
                               stbr3.cohort.factor.predict,stdi.cohort.factor.predict,stdr2.cohort.factor.predict,
                               stgl1.cohort.factor.predict,stin.cohort.factor.predict,stpo1.cohort.factor.predict,
                               stto.cohort.factor.predict) %>% reduce(full_join)

all.cohort.factor.predict = all.cohort.factor.predict %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(Species = recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CAAM-GB" = "CAAM", "CACO1" = "CACO", 
                          "STTO-BH" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                          "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# plotting all species together
rate.cohort.plots=ggplot(all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate, group=altpop, color=altpop))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  scale_color_manual(values=c("black","gray","black","gray"))+
  theme_classic(base_size=15)+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  facet_wrap(.~ phy_order)+
  scale_x_discrete(breaks = c(2,4,7), labels = c("2-Oct","30-Oct","11-Dec"))

rate.cohort.plots.2=rate.cohort.plots+
  geom_line(data=all.contin.data, aes(x=Cohort,y=fit, group=altpop, linetype=altpop), show.legend = FALSE)

# Making each plot individually
# CAAM
caam.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="CAAM")
caam.all.contin.data=subset(all.contin.data, all.contin.data$Species=="CAAM")

caam.plot=ggplot(caam.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate")+
  ggtitle("CAAM")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

caam.plot.2=caam.plot+
  geom_line(data=caam.all.contin.data, aes(x=Cohort,y=fit), show.legend = FALSE)+
  geom_line(data=caam.all.contin.data, aes(x=Cohort,y=lwr), color = "gray74",show.legend = FALSE)+
  geom_line(data=caam.all.contin.data, aes(x=Cohort,y=upr), color = "gray74",show.legend = FALSE)
  

# CAAN
caan.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Species=="CAAN")
caan.all.contin.data=subset(all.contin.data, all.contin.data$Species=="CAAN")

caan.plot=ggplot(caan.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate, group=altpop, color=altpop))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75,show.legend = FALSE)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","gray","black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("CAAN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

caan.plot.2=caan.plot+
  geom_line(data=caan.all.contin.data, aes(x=Cohort,y=fit,group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=caan.all.contin.data, aes(x=Cohort,y=lwr,group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=caan.all.contin.data, aes(x=Cohort,y=upr,group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)

# CACO
caco.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="CACO")
caco.all.contin.data=subset(all.contin.data, all.contin.data$Species=="CACO")

# artificially insert a 0 for cohorts 4 and 5 to make graph match
caco.all.cohort.factor.predict[6,1]=NA
caco.all.cohort.factor.predict[7,1]=NA
caco.all.cohort.factor.predict[6,2]=4
caco.all.cohort.factor.predict[7,2]=5

caco.plot=ggplot(caco.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("CACO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

caco.plot.2=caco.plot+
  geom_line(data=caco.all.contin.data, aes(x=Cohort,y=fit), show.legend = FALSE)+
  geom_line(data=caco.all.contin.data, aes(x=Cohort,y=lwr), color = "gray74",show.legend = FALSE)+
  geom_line(data=caco.all.contin.data, aes(x=Cohort,y=upr), color = "gray74",show.legend = FALSE)


# CAIN
cain.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Species=="CAIN")
cain.all.contin.data=subset(all.contin.data, all.contin.data$Species=="CAIN")

# subset out the two cohorts that are significant
sig.cohort=cain.all.cohort.factor.predict[c(2,4),]

cain.plot=ggplot(cain.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate, group=altpop,color=altpop))+
  geom_point(show.legend = FALSE)+
  geom_point(data=sig.cohort, color="red")+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  geom_errorbar(data=sig.cohort, aes(ymin=ymin, ymax=ymax),width=0.75, color="red")+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","gray"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("CAIN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

cain.plot.2=cain.plot+
  geom_line(data=cain.all.contin.data, aes(x=Cohort,y=fit, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=cain.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=cain.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STDR
stdr.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STDR")
stdr.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STDR")

stdr.plot=ggplot(stdr.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("STDR")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stdr.plot.2=stdr.plot+
  geom_line(data=stdr.all.contin.data, aes(x=Cohort,y=fit, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=stdr.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stdr.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STBR
stbr.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STBR")
stbr.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STBR")

stbr.plot=ggplot(stbr.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("STBR")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stbr.plot.2=stbr.plot+
  geom_line(data=stbr.all.contin.data, aes(x=Cohort,y=fit, group=altpop, linetype=altpop), show.legend = FALSE)+
  geom_line(data=stbr.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stbr.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STTO
stto.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STTO")
stto.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STTO")
stto.all.cohort.factor.predict$color=c("black","red","black","red","black","black","black")

stto.plot=ggplot(stto.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("STTO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stto.plot.2=stto.plot+
  geom_line(data=stto.all.contin.data, aes(x=Cohort,y=fit, color="black"), show.legend = FALSE)+
  geom_line(data=stto.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stto.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STDI
stdi.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STDI")
stdi.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STDI")
stdi.all.cohort.factor.predict$color=c("black","red","black","red","black","black","black")

stdi.plot=ggplot(stdi.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75, show.legend = FALSE)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("STDI")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stdi.plot.2=stdi.plot+
  geom_line(data=stdi.all.contin.data, aes(x=Cohort,y=fit, color="black"), show.legend = FALSE)+
  geom_line(data=stdi.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stdi.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STPO
stpo.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STPO")
stpo.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STPO")
stpo.all.cohort.factor.predict$color=c("black","red","black","red","black","black","black")

stpo.plot=ggplot(stpo.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate,color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75,show.legend = FALSE)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate")+
  ggtitle("STPO")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stpo.plot.2=stpo.plot+
  geom_line(data=stpo.all.contin.data, aes(x=Cohort,y=fit, color="black"), show.legend = FALSE)+
  geom_line(data=stpo.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stpo.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STIN
stin.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STIN")
stin.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STIN")
stin.all.cohort.factor.predict$color=c("black","red","black","red","black","black","black")

stin.plot=ggplot(stin.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate, color=color))+
  geom_point(show.legend = FALSE)+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75,show.legend = FALSE)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black","red","black","red","black","black","black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate", color="Population")+
  ggtitle("STIN")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stin.plot.2=stin.plot+
  geom_line(data=stin.all.contin.data, aes(x=Cohort,y=fit, color="black"), show.legend = FALSE)+
  geom_line(data=stin.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stin.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)


# STGL
stgl.all.cohort.factor.predict=subset(all.cohort.factor.predict, all.cohort.factor.predict$Pop=="STGL")
stgl.all.contin.data=subset(all.contin.data, all.contin.data$Species=="STGL")

stgl.plot=ggplot(stgl.all.cohort.factor.predict, aes(x=cohort.factor, y=germ.rate))+
  geom_point()+
  geom_errorbar(aes(ymin=ymin, ymax=ymax),width=0.75)+
  ylim(-0.07,0.34)+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  labs(x="Rainfall Onset Date",y="Germination Rate")+
  ggtitle("STGL")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7), labels = c("17-Sept","2-Oct","16-Oct","30-Oct","13-Nov","27-Nov","11-Dec"))

stgl.plot.2=stgl.plot+
  geom_line(data=stgl.all.contin.data, aes(x=Cohort,y=fit, color="black"), show.legend = FALSE)+
  geom_line(data=stgl.all.contin.data, aes(x=Cohort,y=lwr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)+
  geom_line(data=stgl.all.contin.data, aes(x=Cohort,y=upr, group=altpop, linetype=altpop), color = "gray74",show.legend = FALSE)



plot_grid(stdr.plot.2,stbr.plot.2,stto.plot.2,stdi.plot.2,stpo.plot.2,caam.plot.2,stin.plot.2,
          stgl.plot.2,caan.plot.2,caco.plot.2,cain.plot.2,
          nrow=3, ncol=4)

#ggsave("Germination.Timing/Plots/Germrate_cohort_panelbyspecies.CI.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/Germrate_cohort_panelbyspecies.CI.png", height = 10, width = 12)

# make individual effects plots if need be
# https://cran.r-project.org/web/packages/jtools/vignettes/effect_plot.html

#### Figure 6:  Germination Fraction as a function of cohort and temperature ####

global.germ.proport.temp.block=read.csv("Germination.timing/Formatted.Data/global.germ.proportion.block.tempdiff.csv")
global.germ.proport.temp.block$Pop=as.factor(global.germ.proport.temp.block$Pop)

#Jenny adding species name and altpop to germ proportion dataframe
global.germ.proport.temp.block = global.germ.proport.temp.block%>%
                                  mutate(Species = as.factor(substr(Pop, 1, 4))) %>%
                                  mutate(Species = recode(Species, "caan" = "CAAN", "cain" = "CAIN", "caam" = "CAAM", "caco" = "CACO", 
                                                          "stto" = "STTO", "stdr" = "STDR", "stpo" = "STPO", "stbr" = "STBR", 
                                                          "stdi" = "STDI", "stgl" = "STGL", "stin" = "STIN")) %>% 
                                  mutate(altpop = 3) %>%
                                  mutate(altpop= ifelse(Pop %in% c("caan2", "cain4"), 4, altpop)) %>%
                                  mutate(altpop = as.factor(altpop)) %>%
                                  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                                                 "STGL", "CAAN", "CACO", "CAIN"))

table(global.germ.proport.temp.block$Pop, global.germ.proport.temp.block$altpop)
table(global.germ.proport.temp.block$Species, global.germ.proport.temp.block$altpop)

# global generalized linear mixed effects model
y=cbind(global.germ.proport.temp.block$count, (global.germ.proport.temp.block$sample-global.germ.proport.temp.block$count))
global.germ.proport.mod=glmer(y~mean.temp*Pop + (1|Block), data=global.germ.proport.temp.block, family = binomial)

# getting population specific slopes and CIs
germ.proport.emm=emtrends(global.germ.proport.mod, "Pop", var="mean.temp")
germ.proport.emm.df=as.data.frame(germ.proport.emm)

# plotted based on Jenny's code: https://github.com/jrgremer/GRAMPS_IPMs/blob/main/scripts/GRAMPS_demog_survival.R
Pops = unique(global.germ.proport.temp.block$Pop)
Temp = seq(min(global.germ.proport.temp.block$mean.temp, na.rm=T), max(global.germ.proport.temp.block$mean.temp, na.rm=T), length=100)
altpop = unique(global.germ.proport.temp.block$altpop)
Species = unique(global.germ.proport.temp.block$Species)
Block = unique(global.germ.proport.temp.block$Block)

proportpred = expand_grid(mean.temp =  Temp, Pop = Pops) %>%
  mutate(logitproport =  predict(global.germ.proport.mod, ., re.form = NA)) %>% #re.form should predict at higher level (not subject), similar to old level = 0 call
  mutate(proportion = inv.logit(logitproport)) %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(Species = recode(Species, "caan" = "CAAN", "cain" = "CAIN", "caam" = "CAAM", "caco" = "CACO", 
                          "stto" = "STTO", "stdr" = "STDR", "stpo" = "STPO", "stbr" = "STBR", 
                          "stdi" = "STDI", "stgl" = "STGL", "stin" = "STIN")) %>% 
  mutate(altpop = 3) %>%
  mutate(altpop= ifelse(Pop %in% c("caan2", "cain4"), 4, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

global.germ.proport.temp.block$proportion=global.germ.proport.temp.block$count/global.germ.proport.temp.block$sample

germproport_bin = ggplot(global.germ.proport.temp.block, aes(x = mean.temp, y = proportion, shape = altpop, group = altpop, linetype = altpop)) + 
  geom_point(show.legend = FALSE) + 
  theme_classic(base_size=15) + scale_shape_manual(values = c(16,1))  + #scale_color_manual(values = c("black", "gray")) +
  labs(x = "Mean Temperature (°C)", y= "Germination Fraction", linetype = "Population", shape = "Population") +
  facet_wrap(.~ phy_order)

# remove data that isn't significant so it doesn't plot the line
listy=c("CAAM","CAAN","CACO")
proportpred.2 = proportpred %>%
  group_by(Species) %>%
  filter(!any(Species %in% listy))
                             
germproport_wpred = germproport_bin + geom_line(data = proportpred.2,aes(x= mean.temp, y = proportion, group = altpop, 
                                                                       linetype = altpop), size=1, show.legend = FALSE) #+
germproport_wpred

#ggsave("Germination.Timing/Plots/Germproport_meantemp_panelbyspecies.test.2.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/Germproport_meantemp_panelbyspecies.png", height = 10, width = 12)

# get proportion for each block for re-watering in year 2

year.2.germ = read.csv("Germination.Timing/Formatted.Data/germ.pheno.round.2.temps.ranges.csv", row.names = 1)
year.2.germ.2=subset(year.2.germ, year.2.germ$Population !="blank")

total.germ.R2=year.2.germ.2 %>%
  group_by(Population,Cohort,Old.Block)%>%
  dplyr::summarise(count=sum(germinated),
            sample=sum(planted))

#write.csv(total.germ.R2, file = "./Germination.Timing/Formatted.Data/global.germ.proportion.block.round.2.csv")

# gives mean of mean Temp for each block within each Cohort
# output for each pop and add to global.germ.proportion.block.round.2.csv
block.mean.temp=year.2.germ.2 %>%
  filter(Population=="STTO-BH")%>%
  group_by(Cohort,Old.Block) %>% 
  dplyr::summarize(mean.temp= mean(mean.Temp, na.rm=TRUE))
# write.csv(block.mean.temp, file="block.mean.temp.csv")

global.germ.proport.block.R2=read.csv("Germination.timing/Formatted.Data/global.germ.proportion.block.round.2.csv", row.names = 1)

global.germ.proport.block.R2$cohort.proportion = global.germ.proport.block.R2$count/global.germ.proport.block.R2$sample

global.germ.proport.block.R2.4 = global.germ.proport.block.R2 %>%
  mutate(Species = substr(Population, 1, 4)) %>%
  mutate(altpop = 3) %>%
  mutate(altpop= ifelse(Population %in% c("CAAN2", "CAIN4"), 4, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))


germproport_wpred + geom_point(data = global.germ.proport.block.R2.4,aes(x= mean.temp, y = cohort.proportion, shape = altpop, 
                                                                         group = altpop), color = "red",show.legend = FALSE)

# taking mean across all cohorts
global.germ.proport.block.R2.2 = global.germ.proport.block.R2 %>%
  group_by(Population,Old.Block) %>%
  dplyr::summarise(count.2 = sum(count),
                   sample.2 = sum(sample),
                   mean.temp.2 = mean(mean.temp, na.rm = TRUE))

global.germ.proport.block.R2.2$proportion=global.germ.proport.block.R2.2$count.2/global.germ.proport.block.R2.2$sample.2


global.germ.proport.block.R2.3 = global.germ.proport.block.R2.2 %>%
  mutate(Species = substr(Population, 1, 4)) %>%
  mutate(altpop = 3) %>%
  mutate(altpop= ifelse(Population %in% c("CAAN2", "CAIN4"), 4, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))


germproport_wpred + geom_point(data = global.germ.proport.block.R2.3,aes(x= mean.temp.2, y = proportion, shape = altpop, 
                                                                         group = altpop), color = "red",show.legend = FALSE)


#### Figure 7: Germination rate as a function of cohort and temperature ####

germ.pheno.all <- read.csv("Germination.Timing/Formatted.Data/germ.pheno.temps.ranges.csv", row.names = 1)
# remove blanks from germ.pheno.all
germ.pheno=subset(germ.pheno.all, germ.pheno.all$Pop !="blank")
# calculate germination rate for each seed
germ.pheno[,21]=1/germ.pheno$days2germ2
colnames(germ.pheno)[21]="germ.rate" 
# remove NAs
germ.pheno.2=germ.pheno[!is.na(germ.pheno$mean.Temp),] # remove NAs

# global linear mixed effects model

global.rate.mod=lmer(germ.rate~mean.Temp*Pop + Cohort + (1|Block), 
                   data=germ.pheno.2)

#check model
plot(global.rate.mod)
qqnorm(residuals(global.rate.mod))

#checking linearity in variables
ggplot(data.frame(x1=germ.pheno.2$mean.Temp,pearson=residuals(global.rate.mod,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw() #more variance at warmer temps...

ggplot(data.frame(x1=germ.pheno.2$Pop,pearson=residuals(global.rate.mod,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw() #seems fine

# getting population specific slopes and CIs
germ.rate.emm=emtrends(global.rate.mod, "Pop", var="mean.Temp")
germ.rate.emm.df=as.data.frame(germ.rate.emm)

# Predicting values across range of temps for each Pop and plotting these.
# website with code to help with this: https://stats.oarc.ucla.edu/r/seminars/interactions-r/

mylist <- list(mean.Temp=seq(7,22,by=1), Pop=c("CAAM-GB","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3",
                                               "STDI","STDR2","STGL1","STIN","STPO1","STTO-BH"))

dat=emmip(global.rate.mod,Pop~mean.Temp,at=mylist, CIs=TRUE, plotit = FALSE)

# trying new way to predict and get exact same graphs

Pops = unique(germ.pheno.2$Pop)
Temp = seq(min(germ.pheno.2$mean.Temp, na.rm=T), max(germ.pheno.2$mean.Temp, na.rm=T), length=100)

pred.rate.temp = expand_grid(mean.Temp =  Temp, Pop = Pops) %>%
  mutate(yvar=predict(global.rate.mod, ., re.form = NA)) %>% #re.form should predict at higher level (not subject), similar to old level = 0 call
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(Species = recode(Species, "CAAN" = "CAAN", "CAIN" = "CAIN", "CAAM-GB" = "CAAM", "CACO1" = "CACO", 
                          "STTO-BH" = "STTO", "STDR2" = "STDR", "STPO1" = "STPO", "STBR3" = "STBR", 
                          "STDI" = "STDI", "STGL1" = "STGL", "STIN" = "STIN")) 

# add species column to dat
#dat$Species=c("CAAM","CAAN","CAAN","CACO","CAIN","CAIN","STBR","STDI","STDR","STGL","STIN","STPO","STTO")

plot.pop.names=c("CAAM","CAAN1","CAAN2","CACO","CAIN3","CAIN4","STBR","STDI","STDR","STGL","STIN","STPO","STTO")

#observed values
germ.pheno.2 = germ.pheno.2 %>%
                mutate(Species = substr(Pop, 1, 4)) %>%
                mutate(altpop = 1) %>%
                mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
                mutate(altpop = as.factor(altpop)) %>%
                mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN"))

# predicted values
dat = dat %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) %>%
  rename(germ.rate= yvar)

pred.rate.temp.2 = pred.rate.temp %>%
  mutate(Species = substr(Pop, 1, 4)) %>%
  mutate(altpop = 1) %>%
  mutate(altpop= ifelse(Pop %in% c("CAAN2", "CAIN4"), 2, altpop)) %>%
  mutate(altpop = as.factor(altpop)) %>%
  mutate(phy_order = fct_relevel(Species, "STDR", "STBR", "STTO", "STDI", "STPO", "CAAM", "STIN", 
                                 "STGL", "CAAN", "CACO", "CAIN")) %>%
  rename(germ.rate= yvar)

table(dat$Species, dat$altpop)

germrate_obs = ggplot(germ.pheno.2, 
                      aes(x = mean.Temp, y = germ.rate, shape = altpop, group = altpop, linetype = altpop)) + 
  geom_point() + theme_classic(base_size=15) + scale_shape_manual(values = c(16,1))  + #scale_color_manual(values = c("black", "gray")) +
  labs(x = "Mean Temperature (°C)", y= "Germination Rate", linetype = "Population", shape = "Population") +
  facet_wrap(.~ phy_order)
  
germrate_wpred = germrate_obs + geom_line(data = dat,aes(x= mean.Temp, y = germ.rate, group = altpop, 
                                                                       linetype = altpop), size=1)

# new predict way, same graph
germrate_wpred = germrate_obs + geom_line(data = pred.rate.temp.4,aes(x= mean.Temp, y = germ.rate, group = altpop, 
                                                         linetype = altpop), size=1)

germrate_wpred

#ggsave("Germination.Timing/Plots/Germrate_meantemp_panelbyspecies.cohort.fixed.pdf", height = 10, width = 12)
#ggsave("Germination.Timing/Plots/Germrate_meantemp_panelbyspecies.cohort.fixed.png", height = 10, width = 12)

#just checking points/model fit...
germrate_obs_all = ggplot(germ.pheno.2, 
       aes(x = mean.Temp, y = germ.rate, shape = altpop, group = Species, linetype = Species, color = Species)) + 
  geom_point(size = 4) + theme_classic(base_size=15) + scale_shape_manual(values = c(16,1))  + #scale_color_manual(values = c("black", "gray")) +
  labs(x = "Mean Temperature (degC)", y= "Germination rate\n (1/Days to germination)")

germrate_all_wpred = germrate_obs_all + geom_line(data = dat,aes(x= mean.Temp, y = germ.rate,  
                                                                 shape = altpop, group = Species, linetype = Species, color = Species), size=1) #+
germrate_all_wpred

#### Figure 8: comparing species slopes #####
# read in phylogeny

all.phylo <- read.tree("Germination.Timing/Raw.Data/tree_pruned.new")

sp.list=c("Caulanthus_amplexicaulis","Caulanthus_anceps","Caulanthus_coulteri","Caulanthus_inflatus",
          "Streptanthus_breweri","Streptanthus_diversifolius","Streptanthus_drepanoides",
          "Streptanthus_glandulosus","Streptanthus_insignis","Streptanthus_polygaloides",
          "Streptanthus_tortuosus")

# prune phylogeny to only includes our species
phylo=keep.tip(all.phylo, sp.list)
phylo$tip.label=c("Caulanthus_inflatus.3", "Caulanthus_coulteri","Caulanthus_anceps.1",         
                  "Streptanthus_glandulosus","Streptanthus_insignis","Caulanthus_amplexicaulis", 
                  "Streptanthus_polygaloides", "Streptanthus_diversifolius", "Streptanthus_tortuosus",    
                  "Streptanthus_breweri","Streptanthus_drepanoides")

# artificially add branch for CAAN2 and CAIN4

plotTree(phylo)
nodelabels()
edgelabels()

which(phylo$tip.label=="Caulanthus_anceps.1") # 3,edge 5

tree.1=bind.tip(phylo, tip.label = "Caulanthus_anceps.2", edge.length=0.001497901, where = 3, position = 0)
plot(tree.1)
nodelabels()
edgelabels()

which(tree.1$tip.label=="Caulanthus_inflatus.3") # 1, edge 2

tree.2=bind.tip(tree.1, tip.label = "Caulanthus_inflatus.4", edge.length=0.001158696, where = 1, position = 0)
plot(tree.2)

# read in the slopes of the relationships
slopes <- read.csv("Germination.Timing/Formatted.Data/slopes.csv", row.names = 1)

# making slopes named vectors
germrate.slope=slopes$germrate.slope
names(germrate.slope)=slopes$Binomial

germrate.cohort=slopes$germrate.cohort
names(germrate.cohort)=slopes$Binomial

germproport.slope=slopes$germproport.slope
names(germproport.slope)=slopes$Binomial

germproport.cohort=slopes$germproport.cohort
names(germproport.cohort)=slopes$Binomial

germrate.cohort.absolute=slopes$absolute.germrate.cohort
names(germrate.cohort.absolute)=slopes$Binomial

germproport.cohort.absolute=slopes$absolute.germproport.cohort
names(germproport.cohort.absolute)=slopes$Binomial

# best bet may be to make each of these separately and manually put them together in word
# need to output each phylo without the tip names to make merging easier

png(filename = "Germination.Timing/Plots/phylo.germrate.slopes.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germrate.slope, ftype="off", colors="black")
dev.off()

pdf(file = "Germination.Timing/Plots/phylo.germrate.slopes.pdf")
dotTree(tree.2,germrate.slope, ftype="off", colors="black")
dev.off()

png(filename = "Germination.Timing/Plots/phylo.germproport.slopes.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germproport.slope, colors="gray", ftype="off")
dev.off()

pdf(file = "Germination.Timing/Plots/phylo.germproport.slopes.pdf")
dotTree(tree.2,germproport.slope, colors="gray", ftype="off")
dev.off()

png(filename = "Germination.Timing/Plots/phylo.germrate.cohort.slopes.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germrate.cohort.absolute,ftype="off",colors="#00B089") # use absolute to get negative dots larger since that is the direction of the slope
dev.off()

png(filename = "Germination.Timing/Plots/phylo.germrate.cohort.slopes.legend.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germrate.cohort,ftype="off",colors="#00B089") # use absolute to get negative dots larger since that is the direction of the slope
dev.off()

pdf(file = "Germination.Timing/Plots/phylo.germrate.cohort.slopes.pdf")
dotTree(tree.2,germrate.cohort.absolute,ftype="off",colors="#00B089") # use absolute to get negative dots larger since that is the direction of the slope
dev.off()

png(filename = "Germination.Timing/Plots/phylo.germproport.cohort.slopes.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germproport.cohort.absolute,ftype="off", colors="#DDC70F")
dev.off()

png(filename = "Germination.Timing/Plots/phylo.germproport.cohort.slopes.legend.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germproport.cohort,ftype="off", colors="#DDC70F")
dev.off()

pdf(file = "Germination.Timing/Plots/phylo.germproport.cohort.slopes.pdf")
dotTree(tree.2,germproport.cohort.absolute,ftype="off", colors="#DDC70F")
dev.off()

# plot of tree with tip labels
png(filename = "Germination.Timing/Plots/phylo.tips.png",  width = 20, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germrate.slope, ftype="i", colors="black")
dev.off()

pdf(file = "Germination.Timing/Plots/phylo.tips.slopes.pdf")
dotTree(tree.2,germrate.slope, ftype="i", colors="black")
dev.off()

# plot of tree with tip labels
png(filename = "Germination.Timing/Plots/test.tips.png",  width = 15, 
    height = 16, units = "cm", res = 300)
dotTree(tree.2,germrate.slope, ftype="i", colors="black")
dev.off()

# making rownames the Binomial
slopes.2=slopes
rownames(slopes.2)=slopes.2$Binomial
slopes.2=slopes.2[,c(2:5,7,8)]

# temperature slopes
temp.slopes=slopes.2[,c(1:2)]
mean.temp.dot.tree=dotTree(tree.2,temp.slopes[c(1:2)],legend = TRUE,colors="black")
# export as pdf device size 5 x 5

# cohort slopes
dotTree(tree.2,slopes.2[,c(3:4)], colors="black")
# absolute value cohort slopes
dotTree(tree.2,slopes.2[,c(5:6)], legend = TRUE, colors="black")
# export as pdf device size 5 x 5
# change the legend manually for positive and negative

#### historical and contemporary temps ####
# do up to 2015 because 2016 only has up until month 9 (25 years of values)

flint.data.mothly = read.csv("~/Library/CloudStorage/Box-Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv") %>% 
  filter(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                   "STGL1","STIN","STPO1","STTO-BH")) %>%
  mutate(tmean=(tmin+tmax)/2) %>%
  filter(clim_year != 2016) %>%
  group_by(id, clim_year, clim_month) %>%
  dplyr::summarize(Tmin = mean(tmin), Tmax = mean(tmax),Tmean=mean(tmean)) %>%
  mutate(genus = "S") %>%
  mutate(genus= ifelse(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4"), "C", genus))

ggplot(flint.data.mothly, aes(x = clim_month, y = Tmean, color = as.factor(clim_year))) +
  geom_line()+
  facet_wrap(~id)

flint.data.contemporary = flint.data.mothly %>%
  group_by(id) %>%
  filter(clim_year > 1965) %>%
  dplyr::summarise(historical.tmean = mean(Tmean),
                   historial.sd = sd(Tmean))
  
flint.data.historical = flint.data.mothly %>%
  group_by(id) %>%
  filter(clim_year < 1966) %>%
  filter(clim_year > 1914) %>%
  dplyr::summarise(historical.tmean = mean(Tmean),
                   historial.sd = sd(Tmean))

#### average temp on Oct 2 and Oct 30 ####

# compile the daily PRISM data

# Make a list of all .csv files in prism directory
files <- dir("./Germination.Timing/PRISM.daily/",
             pattern="PRISM_tmean.*csv",
             full.names=TRUE)
files

# Create a tibble to hold the imported data

dat <- tibble(path=files, filename=basename(path))
dat

# import the data
dat.2 <- dat %>%
  mutate(sheets=purrr::map(path, 
                           read_csv, skip = 10)) %>%
  dplyr::select(-path)

# combine the data

dat.3 <- dat.2 %>% unnest(sheets) %>%
  dplyr::select(-filename)

colnames(dat.3)[2] = "tmean_F"

# convert F to C

dat.3$tmean_C = fahrenheit.to.celsius(dat.3$tmean_F, round = 2)

dat.4 = dat.3 %>% 
  separate(Date, into = c("Year", "Month", "Day"), sep = "-")

Oct.2 = dat.4 %>%
  filter(Day == "02") %>%
  dplyr::summarise(oct.2.tmean = mean(tmean_C),
                   oct.2.sd = sd(tmean_C))
Oct.2

Oct.30 = dat.4 %>%
  filter(Day == "30") %>%
  dplyr::summarise(oct.30.tmean = mean(tmean_C),
                   oct.20.sd = sd(tmean_C))
Oct.30 

