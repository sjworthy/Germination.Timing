library(tidyverse)
library(rgdal)
library(ncdf4)
library(magrittr)
library(raster)
library(ggfortify) #for autoplot
library(ggrepel) #for labels
library(cowplot)

#### species info ####

species = read.csv("./Germination.Timing/Raw.Data/species.info.csv")
# Figure 1b uses years 1991 - 2016 and months September - December

#### Determining species climate space #####
## code and data here from Herbarium Study: https://ucdavis.app.box.com/folder/126809089300

setwd("~/Library/CloudStorage/Box-Box/StreptanthusDimensions/HerbariumStudy/merged_data")

locs = read.csv("georeferencing_clean.csv")
climate = read.csv("all_herbarium_climate.csv") %>%
  mutate(type = "herbarium")

climsummaries = climate %>% 
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>% 
  group_by(specimen,type) %>% 
  dplyr::summarize(Cwd = mean(cwd), PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  left_join(., locs) %>%
  rename(id = specimen)
# data is summarized over years 1991-2016

climsummaries.month = climate %>% 
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>%
  filter(clim_month %in% c(9,10,11,12)) %>%
  group_by(specimen,type) %>% 
  dplyr::summarize(Cwd = mean(cwd), PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  left_join(., locs) %>%
  rename(id = specimen)
# data is summarized over years 1991-2016 and months 9-12

#### Get climate variables of populations from Flint ####

flint.data = read.csv("~/Library/CloudStorage/Box-Box/StreptanthusDimensions/FlintBCM/HTG_climate_data.csv")
  
flint.data.2 = flint.data %>%
  filter(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                   "STGL1","STIN","STPO1","STTO-BH")) %>%
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>%
  group_by(id) %>%
  dplyr::summarize(Cwd = mean(cwd),PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  mutate(type= "seedpop")

flint.data.months = flint.data %>%
  filter(id %in% c("CAAM","CAAN1","CAAN2","CACO1","CAIN3","CAIN4","STBR3", "STDI","STDR2",
                   "STGL1","STIN","STPO1","STTO-BH")) %>%
  filter(clim_year > 1990) %>% 
  filter(clim_year < 2016) %>%
  filter(clim_month %in% c(9,10,11,12)) %>%
  group_by(id) %>%
  dplyr::summarize(Cwd = mean(cwd),PPT_CV = (sd(ppt_mm)/mean(ppt_mm)),PPT = mean(ppt_mm), 
                   Tmin_SD=sd(tmin), Tmax_SD=sd(tmax),Tmin = mean(tmin), Tmax = mean(tmax)) %>%
  mutate(type= "seedpop")

##### temp,precip, variance months #####
climsummaries.month.2 = subset(climsummaries.month, climsummaries.month$folder %in% c("c_amplexicaulis","c_anceps",
                                                                                      "c_coulteri","c_inflatus",
                                                                                      "s_breweri","s_diversifolius",
                                                                                      "s_drepanoides","s_glandulosus",
                                                                                      "s_insignis","s_polygaloides",
                                                                                      "s_tortuosus"))

all.data.months = bind_rows(climsummaries.month.2, flint.data.months)


all.data.4pc.months = all.data.months %>%
  ungroup() %>%
  dplyr:: select(PPT,PPT_CV,Tmin,Tmax,Tmax_SD,Tmin_SD)


all.data.pc = prcomp(all.data.4pc.months, scale  = TRUE, center = TRUE)
all.pc.dat = data.frame(all.data.pc$x)
all_locs_pc = cbind(all.data.months, all.pc.dat)
all_loadings = data.frame(varnames=rownames(all.data.pc$rotation), all.data.pc$rotation)
autoplot(all.data.pc, all_loadings = TRUE, loadings.label = TRUE, loadings.colour = "grey", loadings.label.colour = "black")

tibble(var_explained = (all.data.pc$sdev^2) / (sum(all.data.pc$sdev^2))) %>%
  mutate(pc = seq(1,length(var_explained), 1)) %>%
  mutate(cumvarexplained = cumsum(var_explained))

ggplot() +
  geom_point(data = filter(all_locs_pc, type == "herbarium"), aes(x = PC1, y = PC2, color = folder), alpha = 0.3, size = 1) +
  geom_point(data = filter(all_locs_pc, type == "seedpop"), aes(x = PC1, y = PC2, color = folder)) +
  geom_text_repel(data = filter(all_locs_pc, type == "seedpop"), aes(x = PC1, y = PC2, label = id))+
  theme_classic()+
  labs(x = paste0("Standardized PC1\n (39.2% explained var.)"), 
       y = paste0("Standardized PC2\n (22.0% explained var.)"))+
  ggtitle("ppt,temp,variance,month")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

all_locs_pc[c(1803:1815),10] = c("c_amplexicaulis","c_anceps","c_anceps",
                                 "c_coulteri","c_inflatus","c_inflatus",
                                 "s_breweri","s_diversifolius",
                                 "s_drepanoides","s_glandulosus",
                                 "s_insignis","s_polygaloides",
                                 "s_tortuosus")

caam = subset(all_locs_pc, all_locs_pc$folder == "c_amplexicaulis")
caan = subset(all_locs_pc, all_locs_pc$folder == "c_anceps")
caco = subset(all_locs_pc, all_locs_pc$folder == "c_coulteri")
cain = subset(all_locs_pc, all_locs_pc$folder == "c_inflatus")
stbr = subset(all_locs_pc, all_locs_pc$folder == "s_breweri")
stdi = subset(all_locs_pc, all_locs_pc$folder == "s_diversifolius")
stdr = subset(all_locs_pc, all_locs_pc$folder == "s_drepanoides")
stgl = subset(all_locs_pc, all_locs_pc$folder == "s_glandulosus")
stin = subset(all_locs_pc, all_locs_pc$folder == "s_insignis")
stpo = subset(all_locs_pc, all_locs_pc$folder == "s_polygaloides")
stto = subset(all_locs_pc, all_locs_pc$folder == "s_tortuosus")

caam.monthly.climate = ggplot() +
  geom_point(data = filter(caam, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(caam, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CAAM")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

caan.monthly.climate = ggplot() +
  geom_point(data = filter(caan, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(caan, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4,color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CAAN")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

caco.monthly.climate = ggplot() +
  geom_point(data = filter(caco, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(caco, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CACO")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cain.monthly.climate = ggplot() +
  geom_point(data = filter(cain, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(cain, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("CAIN")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stbr.monthly.climate = ggplot() +
  geom_point(data = filter(stbr, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stbr, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STBR")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stdi.monthly.climate = ggplot() +
  geom_point(data = filter(stdi, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stdi, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STDI")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stdr.monthly.climate = ggplot() +
  geom_point(data = filter(stdr, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stdr, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STDR")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stgl.monthly.climate = ggplot() +
  geom_point(data = filter(stgl, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stgl, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STGL")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stin.monthly.climate = ggplot() +
  geom_point(data = filter(stin, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stin, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STIN")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stpo.monthly.climate = ggplot() +
  geom_point(data = filter(stpo, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(stpo, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STPO")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

stto.monthly.climate = ggplot() +
  geom_point(data = filter(stto, type == "herbarium"), aes(x = PC1, y = PC2, color = minimumElevationInMeters), size = 1) +
  scale_color_gradient(low = "black", high = "gray", "Elevation (m)")+
  geom_point(data = filter(stto, type == "seedpop"), aes(x = PC1, y = PC2), shape = 17, size = 4, color = "black") +
  theme_classic()+
  labs(x = paste0("Standardized PC1"), 
       y = paste0("Standardized PC2"))+
  ggtitle("STTO")+
  xlim(-3.5,7)+
  ylim(-2.5,6)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

germ.pheno.climate.plots = plot_grid(caam.monthly.climate,caan.monthly.climate,caco.monthly.climate,
                                    cain.monthly.climate,stbr.monthly.climate,stdi.monthly.climate,
                                    stdr.monthly.climate,stgl.monthly.climate,stin.monthly.climate,
                                    stpo.monthly.climate,stto.monthly.climate)

# ggsave("Germination.Timing/Plots/germ.pheno.climate.plots.legend.pdf", height = 10, width = 12)
# ggsave("Germination.Timing/Plots/germ.pheno.climate.plots.png", height = 10, width = 12)
