####################################################
#################### SITE MAP ######################
####################################################

source("Master_modelling.R")

################################
#### LOAD REQUIRED PACKAGES ####
################################

library(phyloseq)
library(tidyverse)
library(emmeans)
library(car)
library(maps)
library(sf)
library(lme4)
library(lmerTest)
library(vegan)
library(DHARMa)
library(viridis)
library(RColorBrewer)
library(glmmTMB)
library(cowplot)

################################
########### LOAD DATA ##########
################################

phyloseq_all_tree <- readRDS("data/phyloseq_all_tree_AK.BR.KS2019.RDS") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_all_blast <- readRDS("data/phyloseq_all_blast_AK.BR.KS2019.RDS") %>% prune_taxa(taxa_sums(.) > 0, .)

# removes bunknown versions
smd_all_tree_fam <- readRDS("data/smd_all_tree_fam_AK.BR.KS2019.RDS") %>%
  filter(!(grepl("_bunknown", AMF_family))) 
smd_all_blast_fam <- readRDS("data/smd_all_blast_fam_AK.BR.KS2019.RDS")
smd_blasttotree <- readRDS("data/smd_blasttotree.RDS") 

smd_all_tree_fam_bunknown <- readRDS("data/smd_all_tree_fam_AK.BR.KS2019.RDS") %>% 
  filter((grepl("_bunknown", AMF_family)))

################################
######### CREATE COLORS ########
################################

colScaleR <- scale_colour_manual(name = "Land use", values = c("orange3","olivedrab4"))
colFillR <- scale_fill_manual(name = "Land use", values = c("orange3","olivedrab4")) 

colScaleR3 <- scale_colour_manual(name = "Land use", values = c("orange3","goldenrod1","olivedrab4"))
colFillR3 <- scale_fill_manual(name = "Land use", values = c("orange3","goldenrod1","olivedrab4")) 

##############################
######## MAP OF SITES ########
##############################

mapdat <- smd_all_tree_fam %>% 
  select(c("Country", "State", "Site", "Remnant", "Remnant_2","Lat","Long","Biome"))%>%
  mutate(latlon = paste(Lat,Long,sep="_")) %>% 
  distinct(latlon, .keep_all = TRUE) %>%
  select(-latlon)

# read world data
world <- map_data("world")

# write out
png("figures/MAP_REMNANT.jpg", width = 4, height = 6, units = 'in', res = 300)
ggplot()+
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray70", alpha = 0.5) +
  geom_point(data = mapdat, aes(x = Long, y = Lat,color=factor(Remnant_2), fill=factor(Remnant_2)), pch = 19, size = 3, alpha = 0.8,position=position_jitter(h=2, w=2)) +
  colScaleR + colFillR +
  xlab("") + ylab ("") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 30)) +
  coord_sf(ylim = c(-65, 85), xlim = c(-180, -30), expand = FALSE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# write out, 3 cat
png("figures/MAP_REMNANT3.jpg", width = 4, height = 6, units = 'in', res = 300)
ggplot()+
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray70", alpha = 0.5) +
  geom_point(data = mapdat, aes(x = Long, y = Lat,color=factor(Remnant), fill=factor(Remnant)), pch = 19, size = 3, alpha = 0.8,position=position_jitter(h=2, w=2)) +
  colScaleR3 + colFillR3 +
  xlab("") + ylab ("") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 30)) +
  coord_sf(ylim = c(-65, 85), xlim = c(-180, -30), expand = FALSE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
