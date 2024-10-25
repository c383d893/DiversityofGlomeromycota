####################################################
######### MODELS PROPORTION FAMILY * BIOME #########
######## MODELS PROPORTION FAMILY * REMNANT ########
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

#################################
########## BLAST TO TREE ########
#################################

p_smd_blasttotree_unmatched <- smd_blasttotree %>%
  separate(blasttotree_AMF_family,c("blast_AMF_family","tree_AMF_family"), remove=FALSE) %>%
  mutate(fam_matched = ifelse(blast_AMF_family==tree_AMF_family,"match","nomatch")) %>%
  filter(fam_matched == "nomatch" | blasttotree_AMF_family=="Unknown_Unknown") %>%
  mutate(blasttotree_AMF_family = ifelse(tree_AMF_family=="Unknown", "Tree_Unknown", blasttotree_AMF_family)) #%>%
  filter(!tree_AMF_family=="Unknown") 

################ OVERALL ###############
# misassigned v unknown
p_smd_btot_o <- p_smd_blasttotree_unmatched %>%
  group_by(blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = ifelse(blast_AMF_family=="Unknown", "Unknown","Misassigned")) %>%
  mutate(Biome = "All")

# unknown
p_smd_btot_unknown_o <- p_smd_blasttotree_unmatched %>%
  filter(blast_AMF_family == "Unknown") %>%
  group_by(blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = "Unknown") %>%
  mutate(Biome = "All")

#misassigned 
p_smd_btot_misassigned_o <- p_smd_blasttotree_unmatched %>%
  filter(!blast_AMF_family == "Unknown") %>%
  group_by(blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = "Misassigned") %>%
  mutate(Biome = "All")

################# BIOME ################

# misassigned v unknown
p_smd_btot <- p_smd_blasttotree_unmatched %>%
  group_by(Biome, blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  group_by(Biome) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = ifelse(blast_AMF_family=="Unknown", "Unknown","Misassigned")) %>%
  rbind(p_smd_btot_o)

#relevel
p_smd_btot$Biome = factor(p_smd_btot$Biome, levels = c("Tropical", "Temperate",  "Boreal", "All"), ordered = TRUE)

# unknown
p_smd_btot_unknown <- p_smd_blasttotree_unmatched %>%
  filter(blast_AMF_family == "Unknown") %>%
  group_by(Biome, blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  group_by(Biome) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = "Unknown") %>%
  rbind(p_smd_btot_unknown_o)

#misassigned 
p_smd_btot_misassigned <- p_smd_blasttotree_unmatched %>%
  filter(!blast_AMF_family == "Unknown") %>%
  group_by(Biome, blasttotree_AMF_family) %>%
  summarize(N = sum(observed)) %>%
  group_by(Biome) %>%
  mutate(prop = N / sum(N)) %>%
  drop_na() %>%
  separate(blasttotree_AMF_family, c("blast_AMF_family", "tree_AMF_family"), remove = FALSE) %>%
  mutate(status = "Misassigned") %>%
  rbind(p_smd_btot_misassigned_o)

#relevel
p_smd_btot_misassigned$Biome = factor(p_smd_btot_misassigned$Biome, levels = c("Tropical", "Temperate",  "Boreal", "All"), ordered = TRUE)

########## PLOT BIOME ##########
##### UNKNOWN V MISASSIGNED ####

png("figures/BLAST_TREE_A.jpg", width = 6, height = 6, units ='in', res = 300)
ggplot(p_smd_btot, aes(x = Biome, y = prop, fill = status)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  colScaleUM + colFillUM +
  #scale_fill_viridis(discrete = TRUE) +
  #scale_color_brewer(palette = "Dark2")
  ylab("Proportion ASVs") +
  xlab("") +
  coord_flip() +
  guides(fill=guide_legend(title="BLAST_TREE")) +
  theme(legend.position="none")
dev.off()

##### UNKNOWN V MISASSIGNED ####

png("figures/BLAST_TREE_B.jpg", width = 6, height = 6, units ='in', res = 300)
ggplot(p_smd_btot_misassigned, aes(x = Biome, y = prop, fill = blasttotree_AMF_family)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  #colScale14 + colFill14 +
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_brewer(palette = "Dark2")
  ylab("Proportion ASVs") +
  xlab("") +
  coord_flip() +
  guides(fill=guide_legend(title="BLAST_TREE")) +
  theme(legend.position="none")
dev.off()

##### UNKNOWN V MISASSIGNED ####
png("figures/BLAST_TREE_C.jpg", width = 6, height = 6, units ='in', res = 300)
ggplot(p_smd_btot_unknown, aes(x = Biome, y = prop, fill = blasttotree_AMF_family)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  #colScale14 + colFill14 +
  scale_fill_viridis(discrete = TRUE) +
  #scale_color_brewer(palette = "Dark2")
  ylab("Proportion ASVs") +
  xlab("") +
  coord_flip() +
  guides(fill=guide_legend(title="BLAST_TREE"))
dev.off()

################################
######### CREATE COLORS ########
################################

colScaleUM <- scale_colour_manual(name = "Land use", values = c("grey40","black"))
colFillUM <- scale_fill_manual(name = "Land use", values = c("grey40", "black")) 

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

totalreads <- smd_all_tree_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  group_by(sampleid) %>%                       
  summarise(totalreads = sum(reads_sample))

sitep_smd_all_tree_fam <- smd_all_tree_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  left_join(totalreads, by = "sampleid") %>%
  mutate(prop = reads_sample / totalreads) %>%
  mutate(rest = totalreads - reads_sample) %>%
  drop_na() 

dat.full <- sitep_smd_all_tree_fam 

dat <- sitep_smd_all_tree_fam %>% 
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(Biome!= "Boreal")

#jitter y values proportions
dat <- dat %>%
  mutate(prop_jitter = abs(jitter(prop)))

#keep all data, binomial too many zeros
prop.fam.tree.full <- glmer(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), family = "binomial", data = dat)

#drop data & validated
#prop.fam.tree <- glmmTMB(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), ziformula =~ 1, family = "binomial", data = dat)
prop.fam.tree <- glmmTMB(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant_2 , ziformula =~ 1, family = "binomial", data = dat)

#model validation
mod <- prop.fam.tree
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(prop.fam.tree, ~Remnant_2*AMF_family*Biome)
means

tree.contrasts <- list(
  #biome
  "Temp v Trop all" = c(-1,-1,-1,-1,-1,-1,-1,-1, 1,1,1,1,1,1,1,1),
  #remnant
  "Remnant v Disturbed" = c(1,-1,1,-1,1,-1,1,-1, 1,-1,1,-1,1,-1,1,-1),
  #family by biome
  "Cla by Temp v Trop" = c(-1,-1,0,0, 0,0,0,0, 1,1,0,0, 0,0,0,0),
  "Div by Temp v Trop" = c(0,0,-1,-1, 0,0,0,0, 0,0,1,1, 0,0,0,0),
  "Glo by Temp v Trop" = c(0,0,0,0, -1,-1,0,0, 0,0,0,0, 1,1,0,0),
  "Unk by Temp v Trop" = c(0,0,0,0, 0,0,-1,-1, 0,0,0,0, 0,0,1,1),
  #family by remnant
  "Cla by Remnant v Disturbed" = c(1,-1,0,0, 0,0,0,0, 1,-1,0,0, 0,0,0,0),
  "Div by Remnant v Disturbed" = c(0,0,1,-1, 0,0,0,0, 0,0,1,-1, 0,0,0,0),
  "Glo by Remnant v Disturbed" = c(0,0,0,0, 1,-1,0,0, 0,0,0,0, 1,-1,0,0),
  "Unk by Remnant v Disturbed" = c(0,0,0,0, 0,0,1,-1, 0,0,0,0, 0,0,1,-1),
  #family by remnant within biome
  "Cla by Land use Temp" = c(1,-1,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0),
  "Div by Land use Temp" = c(0,0,1,-1, 0,0,0,0, 0,0,0,0, 0,0,0,0),
  "Glo by Land use Temp" = c(0,0,0,0, 1,-1,0,0, 0,0,0,0, 0,0,0,0),
  "Unk by Land use Temp" = c(0,0,0,0, 0,0,1,-1, 0,0,0,0, 0,0,0,0),
  "Cla by Land use Trop" = c(0,0,0,0, 0,0,0,0, 1,-1,0,0, 0,0,0,0),
  "Div by Land use Trop" = c(0,0,0,0, 0,0,0,0, 0,0,1,-1, 0,0,0,0),
  "Glo by Land use Trop" = c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,-1,0,0),
  "Unk by Land use Trop" = c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,-1)
)

results <- lsmeans::contrast(means,tree.contrasts)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = dat, aes(Biome, prop_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,0.7) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Remnant_2, prop_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,0.9) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY BY BIOME ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Remnant_2*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Remnant_2*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXRemnXBiome.jpg", width = 17, height = 10, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Remnant_2, prop_jitter, color = Remnant_2), 
             position = position_jitter(width = .2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(rows = vars(Biome), cols = vars(AMF_family)) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
########## TREE - 3 CATEGORY ############
#########################################

totalreads <- smd_all_tree_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  group_by(sampleid) %>%                       
  summarise(totalreads = sum(reads_sample))

sitep_smd_all_tree_fam <- smd_all_tree_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  left_join(totalreads, by = "sampleid") %>%
  mutate(prop = reads_sample / totalreads) %>%
  mutate(rest = totalreads - reads_sample) %>%
  drop_na() 

dat.full <- sitep_smd_all_tree_fam 

dat <- sitep_smd_all_tree_fam %>% 
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(Biome!= "Boreal")

#jitter y values proportions
dat <- dat %>%
  mutate(prop_jitter = abs(jitter(prop)))

#keep all data, binomial too many zeros
prop.fam.tree.full <- glmer(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant + (1|Biome:State:Site:Replicate), family = "binomial", data = dat)

#drop data & validated
#prop.fam.tree <- glmmTMB(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant + (1|Biome:State:Site:Replicate), ziformula =~ 1, family = "binomial", data = dat)
prop.fam.tree <- glmmTMB(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant , ziformula =~ 1, family = "binomial", data = dat)

#model validation
mod <- prop.fam.tree
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(prop.fam.tree, ~Remnant*AMF_family*Biome)
means

tree.contrasts.3 <- list(
  #biome
  "Temp v Trop all" = c(-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1,  1,1,1, 1,1,1, 1,1,1, 1,1,1),
  #remnant v ag
  "Remnant v Ag all" = c(1,0,-1, 1,0,-1, 1,0,-1, 1,0,-1, 1,0,-1, 1,0,-1, 1,0,-1, 1,0,-1),
  #remnant v post-ag
  "Remnant v Post-ag all" = c(0,1,-1, 0,1,-1, 0,1,-1, 0,1,-1, 0,1,-1, 0,1,-1, 0,1,-1, 0,1,-1),
  #ag v post-ag
  "Ag v Post-ag" = c(-1,1,0, -1,1,0, -1,1,0, -1,1,0, -1,1,0, -1,1,0, -1,1,0, -1,1,0),
  
  #family by biome
  "Cla by Temp v Trop" = c(-1,-1,-1, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0),
  "Div by Temp v Trop" = c(0,0,0, -1,-1,-1, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0),
  "Glo by Temp v Trop" = c(0,0,0, 0,0,0, -1,-1,-1, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0),
  "Unk by Temp v Trop" = c(0,0,0, 0,0,0, 0,0,0, -1,-1,-1, 0,0,0, 0,0,0, 0,0,0, 1,1,1),
  
  #family by remnant, remn v ag
  "Cla by Remnant v Ag" = c(1,0,-1, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Ag" = c(0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0),
  "Glo by Remnant v Ag" = c(0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0),
  "Unk by Remnant v Ag" = c(0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0, 1,0,-1),
  #family by remnant, remn v post-ag
  "Cla by Remnant v Post-ag" = c(0,1,-1, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Post-ag" = c(0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0),
  "Glo by Remnant v Post-ag" = c(0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0),
  "Unk by Remnant v Post-ag" = c(0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0, 0,1,-1),
  #family by remnant, ag v post-ag
  "Cla by Ag v Post-ag" = c(-1,1,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0),
  "Div by Ag v Post-ag" = c(0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0),
  "Glo by Ag v Post-ag" = c(0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0),
  "Unk by Ag v Post-ag" = c(0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0),
  
  #family by remnant within biome
  #temperate
  #family by remnant, remn v ag
  "Cla by Remnant v Ag Temp" = c(1,0,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Ag Temp" = c(0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0),
  "Glo by Remnant v Ag Temp" = c(0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0),
  "Unk by Remnant v Ag Temp" = c(0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0),
  #family by remnant, remn v post-ag
  "Cla by Remnant v Post-ag Temp" = c(0,1,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Post-ag Temp" = c(0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0),
  "Glo by Remnant v Post-ag Temp" = c(0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0),
  "Unk by Remnant v Post-ag Temp" = c(0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0,  0,0,0),
  #family by remnant, ag v post-ag
  "Cla by Ag v Post-ag Temp" = c(-1,1,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0),
  "Div by Ag v Post-ag Temp" = c(0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0),
  "Glo by Ag v Post-ag Temp" = c(0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0),
  "Unk by Ag v Post-ag Temp" = c(0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0),
  #tropical
  #family by remnant, remn v ag
  "Cla by Remnant v Ag Trop" = c( 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Ag Trop" = c(0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0, 0,0,0),
  "Glo by Remnant v Ag Trop" = c(0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,-1, 0,0,0),
  "Unk by Remnant v Ag Trop" = c(0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,-1),
  #family by remnant, remn v post-ag
  "Cla by Remnant v Post-ag Trop" = c( 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0, 0,0,0),
  "Div by Remnant v Post-ag Trop" = c(0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0, 0,0,0),
  "Glo by Remnant v Post-ag Trop" = c(0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,-1, 0,0,0),
  "Unk by Remnant v Post-ag Trop" = c(0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,-1),
  #family by remnant, ag v post-ag
  "Cla by Ag v Post-ag Trop" = c( 0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0, 0,0,0),
  "Div by Ag v Post-ag Trop" = c(0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0, 0,0,0),
  "Glo by Ag v Post-ag Trop" = c(0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 0,0,0),
  "Unk by Ag v Post-ag Trop" = c(0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0)
  
)

results <- lsmeans::contrast(means,tree.contrasts.3)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXBiome3.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Biome, prop_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,0.7) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Remnant, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Remnant, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXRemn3.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant, lsmean, color = Remnant, fill = Remnant)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Remnant, prop_jitter, color = Remnant), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR3 + colFillR3 +
  scale_x_discrete(limits = r_order3) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

########### REMNANT BY FAMILY BY BIOME ###########
ref <- lsmeans(prop.fam.tree, pairwise ~ AMF_family*Remnant*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.tree.full, pairwise ~ AMF_family*Remnant*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

png("figures/TREE_PROP_AMFXRemnXBiome3.jpg", width = 17, height = 12, units ='in', res = 300)
ggplot(ref.table, aes(Remnant, lsmean, color = Remnant, fill = Remnant)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data =  dat, aes(Remnant, prop_jitter, color = Remnant), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(rows = vars(Biome), cols = vars(AMF_family)) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR3 + colFillR3 +
  scale_x_discrete(limits = r_order3) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
dev.off()

#########################################
################## BLAST ################
#########################################

totalreads <- smd_all_blast_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  group_by(sampleid) %>%                       
  summarise(totalreads = sum(reads_sample))

sitep_smd_all_blast_fam <- smd_all_blast_fam %>% 
  rownames_to_column("sampleid") %>%
  mutate(sampleid = gsub("\\_L0.*", "_L001", sampleid)) %>%
  left_join(totalreads, by = "sampleid") %>%
  mutate(prop = reads_sample / totalreads) %>%
  mutate(rest = totalreads - reads_sample) %>%
  drop_na() 

dat.full <- sitep_smd_all_blast_fam 

dat <- sitep_smd_all_blast_fam %>% 
  filter(AMF_family!="Acaulosporaceae") %>%
  filter(AMF_family!="Ambisporaceae") %>%
  filter(AMF_family!="Paraglomeraceae") %>%
  filter(AMF_family!="Gigasporaceae") %>%
  filter(AMF_family!="Sacculosporaceae") %>%
  filter(AMF_family!="Pervestustaceae") %>%
  filter(Biome!= "Boreal")

#jitter y values proportions
dat <- dat %>%
  mutate(prop_jitter = abs(jitter(prop)))

#keep all data, binomial too many zeros
prop.fam.blast.full <- glmer(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), family = "binomial", data = dat)

#drop data & validated
prop.fam.blast <- glmmTMB(cbind(reads_sample,rest) ~ Biome*AMF_family*Remnant_2 + (1|Biome:State:Site:Replicate), ziformula =~ 1, family = "binomial", data = dat)

#model validation
mod <- prop.fam.blast
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$Biome)
indep_cat_plot(mod,dat,dat$Remnant_2)
indep_cat_plot(mod,dat,dat$AMF_family)
testZeroInflation(mod) 

#run contrasts 
means <- emmeans(prop.fam.blast, ~Remnant_2*AMF_family*Biome)
means

blast.contrasts <- list(
  #biome
  "Temp v Trop all" = c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,1,1,1,1,1,1,1,1,1),
  #biome
  "Remnant v Disturbed" = c(1,-1,1,-1,1,-1,1,-1,1,-1, 1,-1,1,-1,1,-1,1,-1,1,-1),
  #family by biome
  "Arc by Temp v Trop" = c(-1,-1,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0,0),
  "Cla by Temp v Trop" = c(0,0,-1,-1,0,0,0,0,0,0, 0,0,1,1,0,0,0,0,0,0),
  "Div by Temp v Trop" = c(0,0,0,0,-1,-1,0,0,0,0, 0,0,0,0,1,1,0,0,0,0),
  "Glo by Temp v Trop" = c(0,0,0,0,0,0,-1,-1,0,0, 0,0,0,0,0,0,1,1,0,0),
  "Unk by Temp v Trop" = c(0,0,0,0,0,0,0,0,-1,-1, 0,0,0,0,0,0,0,0,1,1),
  #family by remnant
  "Arc by Remnant v Disturbed" = c(1,-1,0,0,0,0,0,0,0,0, 1,-1,0,0,0,0,0,0,0,0),
  "Cla by Remnant v Disturbed" = c(0,0,1,-1,0,0,0,0,0,0, 0,0,1,-1,0,0,0,0,0,0),
  "Div by Remnant v Disturbed" = c(0,0,0,0,1,-1,0,0,0,0, 0,0,0,0,1,-1,0,0,0,0),
  "Glo by Remnant v Disturbed" = c(0,0,0,0,0,0,1,-1,0,0, 0,0,0,0,0,0,1,-1,0,0),
  "Unk by Remnant v Disturbed" = c(0,0,0,0,0,0,0,0,1,-1, 0,0,0,0,0,0,0,0,1,-1),
  #family by remnant within biome
  "Arc by Land use Temp" = c(1,-1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0),
  "Cla by Land use Temp" = c(0,0,1,-1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0),
  "Div by Land use Temp" = c(0,0,0,0,1,-1,0,0,0,0, 0,0,0,0,0,0,0,0,0,0),
  "Glo by Land use Temp" = c(0,0,0,0,0,0,1,-1,0,0, 0,0,0,0,0,0,0,0,0,0),
  "Unk by Land use Temp" = c(0,0,0,0,0,0,0,0,1,-1, 0,0,0,0,0,0,0,0,0,0),
  "Arc by Land use Trop" = c(0,0,0,0,0,0,0,0,0,0, 1,-1,0,0,0,0,0,0,0,0),
  "Cla by Land use Trop" = c(0,0,0,0,0,0,0,0,0,0, 0,0,1,-1,0,0,0,0,0,0),
  "Div by Land use Trop" = c(0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,-1,0,0,0,0),
  "Glo by Land use Trop" = c(0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,1,-1,0,0),
  "Unk by Land use Trop" = c(0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,-1)
)

results <- lsmeans::contrast(means,blast.contrasts)
results.df <- as.data.frame(results)
results.df

########### BIOME BY FAMILY ###########
ref <- lsmeans(prop.fam.blast, pairwise ~ AMF_family*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.blast.full, pairwise ~ AMF_family*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

#png("figures/blast_PROP_AMFXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Biome, lsmean, color = Biome, fill = Biome)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = dat, aes(Biome, prop_jitter, color = Biome), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleB + colFillB +
  scale_x_discrete(limits = b_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

########### REMNANT BY FAMILY ###########
ref <- lsmeans(prop.fam.blast, pairwise ~ AMF_family*Remnant_2, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, sep="_"))

f.ref <- lsmeans(prop.fam.blast.full, pairwise ~ AMF_family*Remnant_2, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

#png("figures/blast_PROP_AMFXRemn.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = dat, aes(Remnant_2, prop_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(~AMF_family) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()

########### REMNANT BY FAMILY BY BIOME ###########
ref <- lsmeans(prop.fam.blast, pairwise ~ AMF_family*Remnant_2*Biome, data = dat, type = "response")
ref.table <- as.data.frame(ref[1]$lsmeans)  %>%
  mutate(comb = paste(AMF_family, Remnant_2, Biome, sep="_"))

f.ref <- lsmeans(prop.fam.blast.full, pairwise ~ AMF_family*Remnant_2*Biome, data = dat.full, type = "response")
f.ref.table <- as.data.frame(f.ref[1]$lsmeans) 
f.ref.table <- f.ref.table %>% 
  mutate(comb = paste(AMF_family, Remnant_2, Biome, sep="_")) %>%
  filter(!comb %in% ref.table$comb) %>%
  drop_na()

#png("figures/blast_PROP_AMFXRemnXBiome.jpg", width = 17, height = 8, units ='in', res = 300)
ggplot(ref.table, aes(Remnant_2, lsmean, color = Remnant_2, fill = Remnant_2)) + 
  geom_bar(stat="identity", alpha = 0.7) +
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width = 0.5,size=1, position=position_dodge(1)) +
  geom_point(data = dat, aes(Remnant_2, prop_jitter, color = Remnant_2), 
             position = position_jitter(width = 0.2, height = 0),size=2, alpha=0.3) +
  geom_bar(data = f.ref.table, stat="identity", alpha = 0.1) +
  geom_point(data = f.ref.table, position=position_dodge(1), size =4, alpha= 0.7) + 
  facet_grid(rows = vars(Biome), cols = vars(AMF_family)) +
  theme_minimal(base_size = 25) +
  ylab("Proportion reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  colScaleR + colFillR +
  scale_x_discrete(limits = r_order) +
  ylim(-0.1,1) +
  theme(strip.text.x = element_text(angle = 60))
#dev.off()
