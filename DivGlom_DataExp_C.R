####################################################
################ DATA EXPLORATION ##################
####################################################

################################
#### LOAD REQUIRED FUNCTIONS ###
################################

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
####### DATA EXPLORATION #######
########### TREE DATA ##########
################################

# conclusions here:
# MAT, MAP and PH seem pretty useless (but are correlated)
# hard to tell what's going on with boxplots categorical vars
# no obvious outliers in dotplots

names(smd_all_tree_fam)
cont.var <- c("MAT","MAP","PH")
cat.var <- c("Biome","AMF_family","Remnant_2")
d.dat <- smd_all_tree_fam
Mydotplot(d.dat[,cont.var]) # dotplot for outliers; looks ok
Mypairs(d.dat[,cont.var]) # all correlated; choose one. 

# observed
# cont var v response
p1 <- ggplot(data = d.dat, aes(x = PH, y = observed)) + geom_point() + geom_jitter() + geom_smooth() 
p2 <- ggplot(data = d.dat, aes(x = MAT, y = observed)) + geom_point() + geom_jitter()+ geom_smooth() 
p3 <- ggplot(data = d.dat, aes(x = MAP, y = observed)) + geom_point() + geom_jitter()+ geom_smooth() 
plot_grid(p1, p2, p3, ncol = 3)
# cat var v response
p1 <- ggplot(data = d.dat, aes(x = Biome, y = observed)) + geom_point() + geom_jitter() #+ geom_boxplot() 
p2 <- ggplot(data = d.dat, aes(x = Remnant_2, y = observed)) + geom_point() + geom_jitter()#+ geom_boxplot() 
p3 <- ggplot(data = d.dat, aes(x = AMF_family, y = observed)) + geom_point() + geom_jitter()#+ geom_boxplot() 
plot_grid(p1, p2, p3, ncol = 3)

ggplot(data = d.dat, aes(x = Biome, y = observed)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Remnant_2, y = observed)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Biome, y = observed)) + geom_point() + geom_jitter() + facet_grid(~Remnant_2)#+ geom_boxplot() 

# shannon
# cont var v response
p1 <- ggplot(data = d.dat, aes(x = PH, y = shannon)) + geom_point() + geom_jitter() + geom_smooth() 
p2 <- ggplot(data = d.dat, aes(x = MAT, y = shannon)) + geom_point() + geom_jitter()+ geom_smooth() 
p3 <- ggplot(data = d.dat, aes(x = MAP, y = shannon)) + geom_point() + geom_jitter()+ geom_smooth() 
plot_grid(p1, p2, p3, ncol = 3)
# cat var v response
p1 <- ggplot(data = d.dat, aes(x = Biome, y = shannon)) + geom_point() + geom_jitter() #+ geom_boxplot() 
p2 <- ggplot(data = d.dat, aes(x = Remnant_2, y = shannon)) + geom_point() + geom_jitter()#+ geom_boxplot() 
p3 <- ggplot(data = d.dat, aes(x = AMF_family, y = shannon)) + geom_point() + geom_jitter()#+ geom_boxplot() 
plot_grid(p1, p2, p3, ncol = 3)

ggplot(data = d.dat, aes(x = Biome, y = shannon)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Remnant_2, y = shannon)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Biome, y = shannon)) + geom_point() + geom_jitter() + facet_grid(~Remnant_2)#+ geom_boxplot() 

# chao1
# cont var v response
p1 <- ggplot(data = d.dat, aes(x = PH, y = chao1)) + geom_point() + geom_jitter() + geom_smooth() 
p2 <- ggplot(data = d.dat, aes(x = MAT, y = chao1)) + geom_point() + geom_jitter()+ geom_smooth() 
p3 <- ggplot(data = d.dat, aes(x = MAP, y = chao1)) + geom_point() + geom_jitter()+ geom_smooth() 
plot_grid(p1, p2, p3, ncol = 3)
# cat var v response
p1 <- ggplot(data = d.dat, aes(x = Biome, y = chao1)) + geom_point() + geom_jitter() #+ geom_boxplot() 
p2 <- ggplot(data = d.dat, aes(x = Remnant_2, y = chao1)) + geom_point() + geom_jitter()#+ geom_boxplot() 
p3 <- ggplot(data = d.dat, aes(x = AMF_family, y = chao1)) + geom_point() + geom_jitter()#+ geom_boxplot() 
plot_grid(p1, p2, p3, ncol = 3)

ggplot(data = d.dat, aes(x = Biome, y = chao1)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Remnant_2, y = chao1)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = d.dat, aes(x = Biome, y = chao1)) + geom_point() + geom_jitter() + facet_grid(~Remnant_2)#+ geom_boxplot() 

# prop
site.dat <- smd_all_tree_fam %>% 
  select(c("Site","MAT","MAP","PH","Biome", "State")) %>%
  distinct(Site, .keep_all=TRUE)

sitep_smd_all_tree_fam <- smd_all_tree_fam %>% 
  group_by(Site, Replicate, Remnant_2, AMF_family) %>%
  #ASV
  #summarize(N = sum(observed)) %>%
  summarize(N = sum(reads_sample)) %>%
  group_by(Site, Replicate, Remnant_2) %>%
  mutate(prop = N / sum(N)) %>%
  mutate(prop = prop +1) %>%
  drop_na() %>%
  left_join(site.dat, by = "Site") 

p.dat <- sitep_smd_all_tree_fam 

# cont var v response
p1 <- ggplot(data = p.dat, aes(x = PH, y = prop)) + geom_point() + geom_smooth() 
p2 <- ggplot(data = p.dat, aes(x = MAT, y = prop)) + geom_point() + geom_smooth() 
p3 <- ggplot(data = p.dat, aes(x = MAP, y = prop)) + geom_point() + geom_smooth() 
plot_grid(p1, p2, p3, ncol = 3)
# cat var v response
p1 <- ggplot(data = p.dat, aes(x = Biome, y = prop)) + geom_boxplot() 
p2 <- ggplot(data = p.dat, aes(x = Remnant_2, y = prop)) + geom_boxplot() 
p3 <- ggplot(data = p.dat, aes(x = AMF_family, y = prop)) + geom_boxplot() 
plot_grid(p1, p2, p3, ncol = 3)

ggplot(data = p.dat, aes(x = Biome, y = prop)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = p.dat, aes(x = Remnant_2, y = prop)) + geom_point() + geom_jitter() + facet_grid(~AMF_family)#+ geom_boxplot() 
ggplot(data = p.dat, aes(x = Biome, y = prop)) + geom_point() + geom_jitter() + facet_grid(~Remnant_2)#+ geom_boxplot() 
