####################################################
###### MODELS DIVERSITY ASV * REMNANT * BIOME ######
####################################################

source("Master_modelling.R")

################################
#### LOAD REQUIRED PACKAGES ####
################################

library(phyloseq)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(DESeq2)

################################
######### CREATE COLORS ########
################################

colScaleR <- scale_colour_manual(name = "Land use", values = c("orange3","olivedrab4"))
colFillR <- scale_fill_manual(name = "Land use", values = c("orange3","olivedrab4")) 

################################
########### LOAD DATA ##########
################################

# read phyloseq data
phyloseq_all_tree <- readRDS("data/phyloseq_all_tree_AK.BR.KS2019.RDS") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_all_blast <- readRDS("data/phyloseq_all_blast_AK.BR.KS2019.RDS") %>% prune_taxa(taxa_sums(.) > 0, .)

# add pseudocount of 1
otu_table(phyloseq_all_tree) <- otu_table(phyloseq_all_tree) + 1
otu_table(phyloseq_all_blast) <- otu_table(phyloseq_all_blast) + 1

# family data
taxa_tree <- read.table('data/ASV_family.tsv',header=TRUE) %>% mutate(ASV = X.ASV.ID) %>% select(-X.ASV.ID) %>% mutate(across('ASV', str_replace, 'ASV_', ''))
taxa_blast <- readRDS('data/taxa_blast.RDS') #%>% rownames_to_column(., "ASV") %>% select(-c("NCBI_ID","Percentage")) %>% mutate(across('ASV', str_replace, 'ASV_', ''))

# look at metadata
smd.tree <- as(sample_data(phyloseq_all_tree), 'data.frame')
smd.blast <- as(sample_data(phyloseq_all_blast), 'data.frame')

# subset Biomes (no postag)
phyloseq_all_tree_boreal <- phyloseq_all_tree %>% subset_samples(Biome == "Boreal") %>% subset_samples(!Remnant == "Post-ag")
phyloseq_all_tree_temp <- phyloseq_all_tree %>% subset_samples(Biome == "Temperate") %>% subset_samples(!Remnant == "Post-ag")
phyloseq_all_tree_trop <- phyloseq_all_tree %>% subset_samples(Biome == "Tropical") %>% subset_samples(!Remnant == "Post-ag")

# subset no postag
phyloseq_all_tree_remdist <- phyloseq_all_tree %>% subset_samples(!Remnant == "Post-ag")

# look at metadata
smd.tree.boreal <- as(sample_data(phyloseq_all_tree_boreal), 'data.frame')
smd.tree.temp <- as(sample_data(phyloseq_all_tree_temp), 'data.frame')
smd.tree.trop <- as(sample_data(phyloseq_all_tree_trop), 'data.frame')

# check whether site nested within land use:
smd.tree.boreal %>% group_by(Site, Remnant) %>% tally() # no
smd.tree.temp %>% group_by(Site, Remnant) %>% tally() # no
smd.tree.trop %>% group_by(Site, Remnant) %>% tally() # no

#########################################
################## TREE #################
#########################################

#########################################
############# LAND USE TEST #############
#########################################

#choose one
dat.deseq <- phyloseq_to_deseq2(phyloseq_all_tree_boreal, ~ Remnant)
dat.deseq <- phyloseq_to_deseq2(phyloseq_all_tree_temp, ~ Remnant)
dat.deseq <- phyloseq_to_deseq2(phyloseq_all_tree_trop, ~ Remnant)

#including all biomes 
dat.deseq <- phyloseq_to_deseq2(phyloseq_all_tree_remdist, ~ Remnant + Biome)

#check reference level
levels(dat.deseq$Remnant)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU 
# (+ val = higher in remnant; - val = higher in disturbed sites), along with p-values for the log2 fold change.
deseq.out <- DESeq(dat.deseq, fitType = 'local')

# access the raw results
res <- results(deseq.out)

# sort results by p-value 
resOrdered <- res[order(res$pvalue),]
# more conservative: use for TROP, all biomes combined
resSig <- subset(resOrdered, padj < 0.1) 
# more lenient: use for Boreal and temp
#resSig <- subset(resOrdered, pvalue < 0.1)

# convert to dataframe
res.df <- as.data.frame(resSig) %>%
  rownames_to_column(., "ASV") %>%
  mutate(across('ASV', str_replace, 'ASV_', '')) %>%
  mutate(Remnant = ifelse(log2FoldChange > 0, "disturbed","remnant")) %>%
  left_join(taxa_tree, by = "ASV") %>% 
  drop_na(padj)

# write table
write_csv(res.df, "data/DA_boreal_landuse.csv")
#write_csv(res.df, "data/DA_temp_landuse.csv")
write_csv(res.df, "data/DA_trop_landuse.csv")

write_csv(res.df, "data/DA_landuse.csv")

# across families
# choose one.
png("figures/TREE_DA_Boreal.jpg", width = 17, height = 8, units ='in', res = 300)
png("figures/TREE_DA_Temp.jpg", width = 17, height = 8, units ='in', res = 300)
png("figures/TREE_DA_Trop.jpg", width = 17, height = 8, units ='in', res = 300)

png("figures/TREE_DA_AllBiomes.jpg", width = 17, height = 8, units ='in', res = 300)

ggplot(data = res.df, aes(x = fct_reorder(ASV, AMF_family), y = log2FoldChange, fill = AMF_family)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Differential Abundance \n (log2fold change)") +
  theme(axis.text.x = element_text(angle = 45))
  
dev.off()

#########################################
############## BIOME TEST ###############
#########################################

#choose one
dat.deseq <- phyloseq_to_deseq2(phyloseq_all_tree_tt, ~ Biome)

#check reference level
levels(dat.deseq$Biome)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU 
# (+ val = higher in remnant; - val = higher in disturbed sites), along with p-values for the log2 fold change.
deseq.out <- DESeq(dat.deseq, fitType = 'local')

# access the raw results
res <- results(deseq.out)

# sort results by p-value 
resOrdered <- res[order(res$pvalue),]
# more conservative: use for TEMPTROP
resSig <- subset(resOrdered, pvalue < 0.05) 

# convert to dataframe
res.df <- as.data.frame(resSig) %>%
  rownames_to_column(., "ASV") %>%
  mutate(across('ASV', str_replace, 'ASV_', '')) %>%
  left_join(taxa_tree, by = "ASV") %>% 
  drop_na(padj)

# write table
write_csv(res.df, "data/biome.csv")

# across families
# choose one.
png("figures/TREE_DA_TempvTrop.jpg", width = 17, height = 8, units ='in', res = 300)

ggplot(data = res.df, aes(x = fct_reorder(ASV, AMF_family), y = log2FoldChange, fill = AMF_family)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Differential Abundance \n (log2fold change)") +
  theme(axis.text.x = element_text(angle = 45))

dev.off()
