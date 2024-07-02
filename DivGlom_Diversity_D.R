####################################################
######### MODELS DIVERSITY FAMILY * BIOME ##########
######## MODELS DIVERSITY FAMILY * REMNANT #########
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
library(scales)

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

#jitter y values ASV 3 diversity metrics
smd_all_tree_fam <- smd_all_tree_fam %>%
  mutate(observed_jitter = abs(jitter(observed)), shannon_jitter = abs(jitter(shannon)), chao1_jitter = abs(jitter(chao1)))
smd_all_blast_fam <- smd_all_blast_fam %>%
  mutate(observed_jitter = abs(jitter(observed)), shannon_jitter = abs(jitter(shannon)), chao1_jitter = abs(jitter(chao1)))

################################
######### CREATE COLORS ########
################################

colScaleB <- scale_colour_manual(name = "Biome", values = c("skyblue3","darkgreen","darkorchid4"))
colFillB <- scale_fill_manual(name = "Biome", values = c("skyblue3","darkgreen","darkorchid4")) 

colScaleR <- scale_colour_manual(name = "Land use", values = c("orange3","olivedrab4"))
colFillR <- scale_fill_manual(name = "Land use", values = c("orange3","olivedrab4")) 

colScaleR3 <- scale_colour_manual(name = "Land use", values = c("orange3","goldenrod1","olivedrab4"))
colFillR3 <- scale_fill_manual(name = "Land use", values = c("orange3","goldenrod1","olivedrab4")) 

b_order <- c('Boreal', 'Temperate', 'Tropical')
r_order <- c('Remnant','Disturbed')
r_order3 <- c('Remnant','Post-ag', 'Disturbed')

#########################################
################## TREE #################
#########################################

dat.full <- smd_all_tree_fam
dat.sub <- smd_all_tree_fam %>%
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(Biome!= "Boreal") 

#########################################
############### CONTRASTS ###############
#########################################

tree.contrasts.biome <- list(
  #biome
  "Temp v Trop all" = c(-1,-1,-1,-1,1,1,1,1),
  #family by biome
  "Cla by Temp v Trop" = c(-1,0,0,0, 1,0,0,0),
  "Div by Temp v Trop" = c(0,-1,0,0, 0,1,0,0),
  "Glo by Temp v Trop" = c(0,0,-1,0, 0,0,1,0),
  "Unk by Temp v Trop" = c(0,0,0,-1, 0,0,0,1)
)

tree.contrasts.remn <- list(
  #biome
  "Remnant v Disturbed all" = c(1,1,1,1,-1,-1,-1,-1),
  #family by remnant
  "Cla by Remnant v Disturbed" = c(1,0,0,0, -1,0,0,0),
  "Div by Remnant v Disturbed" = c(0,1,0,0, 0,-1,0,0),
  "Glo by Remnant v Disturbed" = c(0,0,1,0, 0,0,-1,0),
  "Unk by Remnant v Disturbed" = c(0,0,0,1, 0,0,0,-1)
)

#########################################
########### OBSERVED RICHNESS ###########
#########################################

#keep all data: observed poisson, too many zeros
div.tree.obs.full <- glmer(observed ~ Biome*AMF_family*Remnant_2 + reads_sample + (1|Biome:State:Site:Replicate), family = "poisson", data = dat.full)

#drop data & validated
dat <- dat.sub %>%
  filter(., !rownames(.) == 'Brazil_Tocantins_Na1_S48_L0015' & !rownames(.) == 'Brazil_Tocantins_Na4_S51_L0015'& !rownames(.) == 'Brazil_Tocantins_Na5_S52_L0015' &
           !rownames(.) == 'Brazil_Tocantins_Na2_S49_L0015' & !rownames(.) == 'Brazil_Roraima_RR.FL5.D_S79_L0017'& !rownames(.) == 'Brazil_Roraima_RR.Nat1.C_S60_L0017' &
           !rownames(.) == '2019_AMF_20_S20_L0015')

div.tree.obs.b <- glmmTMB(observed ~ AMF_family*Biome + reads_sample + (1|Biome:State:Site:Replicate), data = dat, family = "nbinom2")
div.tree.obs.r <- glmmTMB(observed ~ AMF_family*Remnant_2 + reads_sample + (1|Biome:State:Site:Replicate), data = dat, family = "nbinom2")

#model validation
mod <- div.tree.obs.b
mod <- div.tree.obs.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 40
#Brazil_Tocantins_Na1_S48_L0015
#Brazil_Tocantins_Na4_S51_L0015
#Brazil_Tocantins_Na5_S52_L0015

#Brazil_Tocantins_Na2_S49_L0015
#Brazil_Roraima_RR.FL5.D_S79_L0017
#Brazil_Roraima_RR.Nat1.C_S60_L0017

#2019_AMF_20_S20_L0015

#run contrasts 
means <- emmeans(div.tree.obs.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,tree.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.tree.obs.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,tree.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.tree.obs.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.tree.obs.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_OBS_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = smd_all_tree_fam, aes(Biome, observed_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0), size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha = 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("ASV richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,3) +
  theme(strip.text.x = element_text(angle = 60)) 
dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.obs.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(div.tree.obs.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_OBS_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Remnant_2, observed_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("ASV richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,2) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
########### SHANNON DIVERSITY ###########
#########################################

#keep all data: lmer, too many zeros
div.tree.shan.full <- lmer(shannon ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
#glmer tweedie
#div.tree.shan <- glmmTMB(shannon ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

div.tree.shan.b <- glmmTMB(shannon ~ AMF_family*Biome + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())
div.tree.shan.r <- glmmTMB(shannon ~ AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

#model validation
mod <- div.tree.shan.b
mod <- div.tree.shan.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(div.tree.shan.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,tree.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.tree.shan.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,tree.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.tree.shan.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.tree.shan.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_SHAN_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Biome, shannon_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.shan.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(div.tree.shan.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_SHAN_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Remnant_2, shannon_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
################## CHAO1 ################
#########################################

#keep all data: lmer, too many zeros
div.tree.chao.full <- lmer(chao1 ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
dat <- dat.sub %>%
  filter(., !rownames(.) == 'Brazil_Tocantins_Na1_S48_L0015' & !rownames(.) == 'Brazil_Tocantins_Na4_S51_L0015'& !rownames(.) == 'Brazil_Tocantins_Na5_S52_L0015' &
           !rownames(.) == 'Brazil_Tocantins_Na2_S49_L0015')
#glmer tweedie
div.tree.chao.b <- glmmTMB(chao1 ~ AMF_family*Biome  + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())
div.tree.chao.r <- glmmTMB(chao1 ~ AMF_family*Remnant_2  + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

#model validation
mod <- div.tree.chao.b
mod <- div.tree.chao.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 30

#run contrasts 
means <- emmeans(div.tree.chao.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,tree.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.tree.chao.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,tree.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.tree.chao.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.tree.chao.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_CHAO_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Biome, chao1_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Chao diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,5) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.chao.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(div.tree.chao.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_CHAO_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Remnant_2, chao1_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2,height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Chao diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,6) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
########## TREE - 3 CATEGORY ############
#########################################

dat.full <- smd_all_tree_fam
dat.sub <- smd_all_tree_fam %>%
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(Biome!= "Boreal")

#########################################
############### CONTRASTS ###############
#########################################

#only across remn
tree.contrasts.remn.3 <- list(
  #biome
  "Remnant v Disturbed all" = c(1,1,1,1,0,0,0,0,-1,-1,-1,-1),
  "Remnant v Post-ag all" = c(0,0,0,0,1,1,1,1,-1,-1,-1,-1),
  "Post-ag v Disturbed all" = c(1,1,1,1,-1,-1,-1,-1,0,0,0,0),
  #family by remnant
  "Cla by Remnant v Disturbed" = c(1,0,0,0, 0,0,0,0, -1,0,0,0),
  "Div by Remnant v Disturbed" = c(0,1,0,0, 0,0,0,0, 0,-1,0,0),
  "Glo by Remnant v Disturbed" = c(0,0,1,0, 0,0,0,0, 0,0,-1,0),
  "Unk by Remnant v Disturbed" = c(0,0,0,1, 0,0,0,0, 0,0,0,-1),
  "Cla by Remnant v Post-ag" = c(0,0,0,0, 1,0,0,0, -1,0,0,0),
  "Div by Remnant v Post-ag" = c(0,0,0,0, 0,1,0,0, 0,-1,0,0),
  "Glo by Remnant v Post-ag" = c(0,0,0,0, 0,0,1,0, 0,0,-1,0),
  "Unk by Remnant v Post-ag" = c(0,0,0,0, 0,0,0,1, 0,0,0,-1),
  "Cla by Post-ag v Disturbed" = c(1,0,0,0, -1,0,0,0, 0,0,0,0),
  "Div by Post-ag v Disturbed" = c(0,1,0,0, 0,-1,0,0, 0,0,0,0),
  "Glo by Post-ag v Disturbed" = c(0,0,1,0, 0,0,-1,0, 0,0,0,0),
  "Unk by Post-ag v Disturbed" = c(0,0,0,1, 0,0,0,-1, 0,0,0,0)
  
)

#########################################
########### OBSERVED RICHNESS ###########
#########################################

#keep all data: observed poisson, too many zeros
div.tree.obs.full <- glmer(observed ~ Biome*AMF_family*Remnant + reads_sample + (1|Biome:State:Site:Replicate), family = "poisson", data = dat.full)

#drop data & validated
dat <- dat.sub %>%
  filter(., !rownames(.) == 'Brazil_Tocantins_Na1_S48_L0015' & !rownames(.) == 'Brazil_Tocantins_Na4_S51_L0015'& !rownames(.) == 'Brazil_Tocantins_Na5_S52_L0015' &
           !rownames(.) == 'Brazil_Tocantins_Na2_S49_L0015' & !rownames(.) == 'Brazil_Roraima_RR.FL5.D_S79_L0017'& !rownames(.) == 'Brazil_Roraima_RR.Nat1.C_S60_L0017' &
           !rownames(.) == '2019_AMF_10_S10_L0015' & !rownames(.) == '2019_AMF_20_S20_L0015')

div.tree.obs.r <- glmmTMB(observed ~ AMF_family*Remnant + reads_sample + (1|Biome:State:Site:Replicate), data = dat, family = "nbinom2")

#model validation
mod <- div.tree.obs.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 50
#Brazil_Tocantins_Na1_S48_L0015
#Brazil_Tocantins_Na4_S51_L0015
#Brazil_Tocantins_Na5_S52_L0015

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 20
#Brazil_Tocantins_Na2_S49_L0015
#Brazil_Roraima_RR.FL5.D_S79_L0017
#Brazil_Roraima_RR.Nat1.C_S60_L0017

#2019_AMF_10_S10_L0015
#2019_AMF_20_S20_L0015

#run contrasts 
means <- emmeans(div.tree.obs.r, ~AMF_family*Remnant)
means

results <- lsmeans::contrast(means,tree.contrasts.remn.3)
results.df <- as.data.frame(results)
results.df

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.obs.r, pairwise ~ AMF_family*Remnant, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant, sep="_"))

f.ref <- lsmeans(div.tree.obs.full, pairwise ~ AMF_family*Remnant, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_OBS_AMFXRemn3.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant, lsmean, color = Remnant, fill = Remnant)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Remnant, observed_jitter, color = Remnant), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("ASV richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR3 + colFillR3 +
  scale_x_discrete(limits = r_order3) +
  ylim(-0.1,2) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
########### SHANNON DIVERSITY ###########
#########################################

dat <- dat.sub %>%
  filter(., !rownames(.) == 'Brazil_Tocantins_Na1_S48_L0015' & !rownames(.) == 'Brazil_Tocantins_Na4_S51_L0015'& !rownames(.) == 'Brazil_Tocantins_Na5_S52_L0015')
           
#keep all data: lmer, too many zeros
div.tree.shan.full <- lmer(shannon ~ Biome*AMF_family*Remnant + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
#glmer tweedie
div.tree.shan.r <- glmmTMB(shannon ~ AMF_family*Remnant  + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

#model validation
mod <- div.tree.shan.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 2
#Brazil_Tocantins_Na1_S48_L0015
#Brazil_Tocantins_Na4_S51_L0015
#Brazil_Tocantins_Na5_S52_L0015

#run contrasts 
means <- emmeans(div.tree.shan.r, ~AMF_family*Remnant)
means

results <- lsmeans::contrast(means,tree.contrasts.remn.3)
results.df <- as.data.frame(results)
results.df

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.shan.r, pairwise ~ AMF_family*Remnant, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant, sep="_"))

f.ref <- lsmeans(div.tree.shan.full, pairwise ~ AMF_family*Remnant, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_SHAN_AMFXRemn3.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant, lsmean, color = Remnant, fill = Remnant)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_tree_fam, aes(Remnant, shannon_jitter, color = Remnant), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR3 + colFillR3 +
  scale_x_discrete(limits = r_order3) +
  ylim(-0.1,1.2) +
    theme(strip.text.x = element_text(angle = 60))
dev.off()









#########################################
################## CHAO1 ################
#########################################

#keep all data: lmer, too many zeros
div.tree.chao.full <- lmer(chao1 ~ Biome*AMF_family*Remnant + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
dat <- dat %>%
  filter(., !rownames(.) == 'Brazil_Tocantins_Na2_S49_L0015' & !rownames(.) == 'Brazil_Tocantins_Na1_S48_L0015' & !rownames(.) == 'Brazil_Tocantins_Na4_S51_L0015'& !rownames(.) == 'Brazil_Tocantins_Na5_S52_L0015')

#glmer tweedie
div.tree.chao.r <- glmmTMB(chao1 ~ AMF_family*Remnant  + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

#model validation
mod <- div.tree.chao.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 30
#Brazil_Tocantins_Na2_S49_L0015
#Brazil_Tocantins_Na1_S48_L0015
#Brazil_Tocantins_Na4_S51_L0015
#Brazil_Tocantins_Na5_S52_L0015

#run contrasts 
means <- emmeans(div.tree.chao.r, ~AMF_family*Remnant)
means

results <- lsmeans::contrast(means,tree.contrasts.remn.3)
results.df <- as.data.frame(results)
results.df

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.tree.chao.r, pairwise ~ AMF_family*Remnant, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant, sep="_"))

f.ref <- lsmeans(div.tree.chao.full, pairwise ~ AMF_family*Remnant, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  mutate(lsmean = emmean)

png("figures/TREE_CHAO_AMFXRemn3.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant, lsmean, color = Remnant, fill = Remnant)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = dat, aes(Remnant, chao1_jitter, color = Remnant), 
             position = position_jitter(width = 0.2,height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Chao diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR3 + colFillR3 +
  scale_x_discrete(limits = r_order3) +
  ylim(-0.1,4) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
################# BLAST #################
#########################################

dat.full <- smd_all_blast_fam
dat.sub <- smd_all_blast_fam %>%
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(AMF_family!="Sacculosporaceae") %>%
  filter(AMF_family!="Pervestustaceae") %>%
  filter(Biome!= "Boreal") 

#########################################
############### CONTRASTS ###############
#########################################

blast.contrasts.biome <- list(
  #biome
  "Temp v Trop all" = c(-1,-1,-1,-1,-1,1,1,1,1,1),
  #family by biome
  "Arc by Temp v Trop" = c(-1,0,0,0,0, 1,0,0,0,0),
  "Cla by Temp v Trop" = c(0,-1,0,0,0, 0,1,0,0,0),
  "Div by Temp v Trop" = c(0,0,-1,0,0, 0,0,1,0,0),
  "Glo by Temp v Trop" = c(0,0,0,-1,0, 0,0,0,1,0),
  "Unk by Temp v Trop" = c(0,0,0,0,-1, 0,0,0,0,1)
)

blast.contrasts.remn <- list(
  #biome
  "Remnant v Disturbed all" = c(-1,-1,-1,-1,-1,1,1,1,1,1),
  #family by remnant
  "Arc by Remnant v Disturbed" = c(1,0,0,0,0, -1,0,0,0,0),
  "Cla by Remnant v Disturbed" = c(0,1,0,0,0, 0,-1,0,0,0),
  "Div by Remnant v Disturbed" = c(0,0,1,0,0, 0,0,-1,0,0),
  "Glo by Remnant v Disturbed" = c(0,0,0,1,0, 0,0,0,-1,0),
  "Unk by Remnant v Disturbed" = c(0,0,0,0,1, 0,0,0,0,-1)
)

#########################################
########### OBSERVED RICHNESS ###########
#########################################

#keep all data: observed poisson, too many zeros
div.blast.obs.full <- glmer(observed ~ Biome*AMF_family*Remnant_2 + reads_sample + (1|Biome:State:Site:Replicate), family = "poisson", data = dat.full)

#drop data & validated
dat <- dat.sub %>%
  filter(., !rownames(.) == '2019_AMF_19_S19_L0016' & !rownames(.) == '2019_AMF_28_S28_L0016' & !rownames(.) == '2019_AMF_25_S25_L00110')

div.blast.obs.b <- glmmTMB(observed ~ AMF_family*Biome + reads_sample + (1|Biome:State:Site:Replicate), data = dat, family = "nbinom2")
div.blast.obs.r <- glmmTMB(observed ~ AMF_family*Remnant_2 + reads_sample + (1|Biome:State:Site:Replicate), data = dat, family = "nbinom2")

#model validation
mod <- div.blast.obs.b
mod <- div.blast.obs.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > 50
#2019_AMF_19_S19_L0016
#2019_AMF_28_S28_L0016
#2019_AMF_25_S25_L00110

#run contrasts 
means <- emmeans(div.blast.obs.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,blast.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.blast.obs.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,blast.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.blast.obs.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.blast.obs.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

#png("figures/BLAST_OBS_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_blast_fam, aes(Biome, observed_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("ASV richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,10) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.blast.obs.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) 

f.ref <- lsmeans(div.blast.obs.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() 

#png("figures/BLAST_OBS_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_blast_fam, aes(Remnant_2, observed_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2,height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("ASV richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,2) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

#########################################
########### SHANNON DIVERSITY ###########
#########################################

#keep all data: lmer, too many zeros
div.blast.shan.full <- lmer(shannon ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
#glmer tweedie
div.blast.shan.b <- glmmTMB(shannon ~ AMF_family*Biome + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())
#this one won't converge:
div.blast.shan.r <- glmmTMB(shannon ~ AMF_family*Remnant_2 , data = dat, family = tweedie())

#model validation
mod <- div.blast.shan.b
mod <- div.blast.shan.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(div.blast.shan.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,blast.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.blast.shan.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,blast.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.blast.shan.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.blast.shan.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  rename(lsmean = emmean)

#png("figures/BLAST_SHAN_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Biome, shannon_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.blast.shan.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(div.blast.shan.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  rename(lsmean = emmean)

#png("figures/BLAST_SHAN_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_blast_fam, aes(Remnant_2, shannon_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Shannon diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

#########################################
################## CHAO1 ################
#########################################

#keep all data: lmer, too many zeros
div.blast.chao.full <- lmer(chao1 ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat.full)

#drop data & validated
dat <- dat.sub 

#glmer tweedie
div.blast.chao.b <- glmmTMB(chao1 ~ AMF_family*Biome + (1|Biome:State:Site:Replicate) , data = dat, family = tweedie())
div.blast.chao.r <- glmmTMB(chao1 ~ AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), data = dat, family = tweedie())

#model validation
mod <- div.blast.chao.b
mod <- div.blast.chao.r
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(div.blast.chao.b, ~AMF_family*Biome)
means

results <- lsmeans::contrast(means,blast.contrasts.biome)
results.df <- as.data.frame(results)
results.df

means <- emmeans(div.blast.chao.r, ~AMF_family*Remnant_2)
means

results <- lsmeans::contrast(means,blast.contrasts.remn)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(div.blast.chao.b, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(div.blast.chao.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  rename(lsmean = emmean)

#png("figures/BLAST_CHAO_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_blast_fam, aes(Biome, chao1_jitter, color = Biome), 
             position = position_jitter(width = 0.2,height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Chao diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(div.blast.chao.r, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(div.blast.chao.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na() %>%
  rename(lsmean = emmean)

#png("figures/BLAST_CHAO_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  smd_all_blast_fam, aes(Remnant_2, chao1_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) +  
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Chao diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1, 0.6) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()
