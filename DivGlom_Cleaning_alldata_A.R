## Script for :
## 1. double unknowns, 
## 2. tree unknowns,
## 3. blast to tree mismatches
## 4. blast unknown per tree fam: what is family assigned in tree but not classified to family in BLAST?

################################
#### LOAD REQUIRED PACKAGES ####
################################

library(phyloseq)
library(tidyverse)

################################
####### CLEAN SEQ & META #######
################################

#import asv table: 535 vars, 20263 ASVs
asvs <- read.table("data/ASVtable_clean_BLASTonly.tsv", header = TRUE, sep = '\t') %>%
   column_to_rownames('X.ASV.ID') %>%
   rename_with(~gsub("X2019", "2019", .))

#import meta data: MAT, MAP, PH
meta_comp <- read.table("data/composite_data.csv", header = TRUE, sep = ',') %>%
  select(c("CHELSA_BIO_Annual_Mean_Temperature", "CHELSA_BIO_Annual_Precipitation", "SG_Soil_pH_H2O_000cm", "Sample")) %>%
  rename(MAT = CHELSA_BIO_Annual_Mean_Temperature, MAP = CHELSA_BIO_Annual_Precipitation, PH = SG_Soil_pH_H2O_000cm) 

#import meta data & merge with meta_comp : sample data (country, depths, sample name, remnant...)
#keep kansas years and depths
meta <- read.table("data/MetaDataJointforR_divglom_CDMar13.3.23.txt", header = TRUE, sep = '\t') %>%
  select(-checked) %>%
  #filter out Ecuador; poor quality
  filter(!str_detect(Sample, "MMIAMF")) %>%
  #keep only Kansas with sample name containing 2019 (remove other years)
  #filter(!str_detect(Sample, "Spr18|Fall18")) %>%
  #remove KS depths other than 15
  #or not if want to keep for unknowns
  #filter(!(Depth > 20 & State == "KS")) %>% filter(!(Depth < 10 & State == "KS")) %>%
  filter(!Remnant=="No") %>% 
  #create 2-cat remnant column
  mutate(Remnant_2 = ifelse(Remnant %in% c("Post-ag", "Disturbed"), "Disturbed", "Remnant")) %>%
  #add biome col
  mutate(Biome = case_when(Country == "Brazil" ~ "Tropical",
                           State == "KS" ~ "Temperate",
                           State =="Alaska" ~ "Boreal")) %>%
  mutate(Year = case_when(str_detect(Sample, "Spr18") ~ "2018",
                          str_detect(Sample, "Fall18") ~ "2018",
                          str_detect(Sample, "2019") ~ "2019",
                          Country == "Brazil" ~ "2017",
                          State == "Alaska" ~ "2017")) %>%
  drop_na() %>%
  left_join(meta_comp, by = "Sample") %>%
  column_to_rownames(., "Sample")

write.csv(meta, "data/cleaned_meta_dat.csv")

meta.site <- meta %>% group_by(Site, Remnant) %>% slice(1)
write.csv(meta.site, "data/cleaned_site_meta_dat.csv")

################################
####### CLEAN BLAST DATA #######
################################

#https://docs.google.com/spreadsheets/d/1vHHj9XtfIYj9B60L7UP4rUci8bJeolUP/edit#gid=1484552018

#read BLAST information (to extract family info)
#ID_list contains only ID scores (all % match) filtered in the "Fetch_NCBI_ID_comas.R" script
ID_list <- read.table("data/ncbi_species_output.txt", header = FALSE, sep = '\t', fill = TRUE) %>%
  rename(NCBI_ID = V1) %>%
  mutate(AMF_family = "Undetermined") %>%
  mutate(across(where(is.character), tolower)) 

#read BLAST match (ASV to NCBI_ID, percentage)
blast_doc_all <- read.table("data/blast_allglom_tophit_BLAST_clean07.06.23.txt", header = FALSE, sep = '\t', fill = TRUE) %>%
  select(c("V1","V2","V3")) %>%
  rename(Sample_ID = V1, NCBI_ID = V2, Percentage = V3) %>%
  mutate(NCBI_ID = tolower(NCBI_ID)) 

length(unique(blast_doc_all$NCBI_ID)) # 765 (the same input # ASVs Alexis put into the fetch)

################################
####### BLAST FAM ASSIGN #######
################################

#Remove rows that are carried over from below: aka make one row per ID.
ID_list <- ID_list %>%
# If "gene" or "subunit" in first column, remove:
  filter(!str_detect(NCBI_ID, "gene|subunit")) %>%
# If "gene" is second column, remove:
  filter(!str_detect(V2, "gene")) 

####### UNKNOWN FAMILY #######

ID_list.unknown <- ID_list %>%
  filter_all(any_vars(str_detect(., "glomeromycotina") |
                        str_detect(., "glomerales") | 
                        str_detect(., "glomeromycete") |
                        str_detect(., "glomeromycetes") | 
                        str_detect(., "glomeromycota") |
                        str_detect(., "fungus") | 
                        str_detect(., "ascomycete"))) %>%
  mutate(AMF_family = "Unknown")

######### GLOMERACEAE ########

ID_list.Glo <- ID_list %>%
  filter_all(any_vars(str_detect(., "^glomus$") |
                        str_detect(., "dominikia") | 
                        str_detect(., "rhizoglomus") |
                        str_detect(., "rhizophagus") | 
                        str_detect(., "sclerocarpum") |
                        str_detect(., "redeckera") | 
                        str_detect(., "glomeraceae") |
                        str_detect(., "kamienskia") | 
                        str_detect(., "septoglomus") | 
                        str_detect(., "microdominikia") |
                        str_detect(., "funneliformis") |
                        str_detect(., "nanoglomus") | 
                        str_detect(., "glomeraceae") |
                        str_detect(., "glmous"))) %>%
  mutate(AMF_family = "Glomeraceae")

######### CLAROIDEOGLOMERACEAE ######## 
ID_list.Cla <- ID_list %>%
  filter_all(any_vars(str_detect(., "claroid"))) %>%
  mutate(AMF_family = "Claroideoglomeraceae")

######### ACAULOSPORACEAE ######## 
ID_list.Aca <- ID_list %>%
  filter_all(any_vars(str_detect(., "acau"))) %>%
  mutate(AMF_family = "Acaulosporaceae")

######### PERVESTUSTACEAE ######## 
ID_list.Per <- ID_list %>%
  filter_all(any_vars(str_detect(., "perv"))) %>%
  mutate(AMF_family = "Pervestustaceae")

######### AMBISPORACEAE ######## 
ID_list.Amb <- ID_list %>%
  filter_all(any_vars(str_detect(., "ambi"))) %>%
  mutate(AMF_family = "Ambisporaceae")

######### DIVERSIPORACEAE ######## 
ID_list.Div <- ID_list %>%
  filter_all(any_vars(str_detect(., "dive")| str_detect(., "corymb"))) %>%
  mutate(AMF_family = "Diversisporaceae")

######### GIGASPORACEAE ######## 
ID_list.Gig <- ID_list %>%
  filter_all(any_vars(str_detect(., "gig")| 
                      str_detect(., "cetraspor")| 
                        str_detect(., "orbis")| 
                        str_detect(., "scutell")| 
                        str_detect(., "dentisc"))) %>%
  mutate(AMF_family = "Gigasporaceae")

######### PARAGLOMERACEAE ######## 
ID_list.Par <- ID_list %>%
  filter_all(any_vars(str_detect(., "para"))) %>%
  mutate(AMF_family = "Paraglomeraceae")

######### ARCHAEOSPORACEAE ######## 
ID_list.Arc <- ID_list %>%
  filter_all(any_vars(str_detect(., "archae"))) %>%
  mutate(AMF_family = "Archaesporaceae")

######### PACISPORACEAE ######## 
ID_list.Pac <- ID_list %>%
  filter_all(any_vars(str_detect(., "paci"))) %>%
  mutate(AMF_family = "Pacisporaceae")

######### SACCULOSPORACEAE ######## 
ID_list.Sac <- ID_list %>%
  filter_all(any_vars(str_detect(., "sac")|
                      str_detect(., "entrophospora")| 
                      str_detect(., "baltica"))) %>%
  mutate(AMF_family = "Sacculosporaceae")

###### JOIN ######

ID_list_undetermined <- ID_list %>% filter(!(NCBI_ID %in% ID_list.Aca$NCBI_ID |
                                           NCBI_ID %in% ID_list.Amb$NCBI_ID |
                                           NCBI_ID %in% ID_list.Arc$NCBI_ID |
                                           NCBI_ID %in% ID_list.Cla$NCBI_ID |
                                           NCBI_ID %in% ID_list.Div$NCBI_ID |
                                           NCBI_ID %in% ID_list.Gig$NCBI_ID |
                                           NCBI_ID %in% ID_list.Glo$NCBI_ID |
                                           NCBI_ID %in% ID_list.Pac$NCBI_ID |
                                           NCBI_ID %in% ID_list.Par$NCBI_ID |
                                           NCBI_ID %in% ID_list.Per$NCBI_ID |
                                           NCBI_ID %in% ID_list.Sac$NCBI_ID |
                                           NCBI_ID %in% ID_list.unknown$NCBI_ID))

ID_list_undetermined_solved <- ID_list_undetermined %>%
  mutate(AMF_family = case_when(V2 == "otospora" ~ "Diversisporaceae",
                           V2 == "kuklospora" ~ "Acaulosporaceae",
                           V2 =="palaeospora" ~ "Archaeosporaceae",
                           V2 =="geosiphon" ~ "Unknown",
                           V2 == "polonospora" ~ "Unknown" )) %>%
  #remove [uncultured hymenoscyphus]
  filter(!V2== "uncultured") 
  
ID_list_full <- rbind(ID_list.Aca, ID_list.Amb, ID_list.Arc, ID_list.Cla,
                          ID_list.Div, ID_list.Gig, ID_list.Glo, ID_list.Pac,
                          ID_list.Par, ID_list.Per, ID_list.Sac, ID_list.unknown,ID_list_undetermined_solved)

ID_list_duplicates <- ID_list_full %>% 
  group_by(NCBI_ID) %>% 
  filter(n()>1)

ID_list_no_duplicates <- ID_list_full %>% 
  group_by(NCBI_ID) %>% 
  filter(n()==1)

# glomus claroideum = Claroideoglomeraceae
# entrophospora = Claroideoglomeraceae
# paradentiscutata = Gigasporaceae
# paraglomeraceae = Paraglomeraceae
# claroideoglomeraceae = Claroideoglomeraceae

ID_list_duplicates_clean <- ID_list_duplicates %>%
  filter(V2 == "glomus" & AMF_family =="Claroideoglomeraceae" |
         V2 == "entrophospora" & AMF_family == "Claroideoglomeraceae" |
         V2 == "paradentiscutata" & AMF_family == "Gigasporaceae" |
         V2 == "paraglomeraceae" & AMF_family == "Paraglomeraceae"|
         V3 == "paraglomerales" & AMF_family == "Paraglomeraceae"|
         V3 == "claroideoglomeraceae" & AMF_family == "Claroideoglomeraceae") 

ID_list_final <- rbind(ID_list_no_duplicates, ID_list_duplicates_clean) %>%
  select(c(NCBI_ID, AMF_family))

################################
###### PHYLOSEQ OBJECTS ########
################################

# 20263 (including all levels of matching)
# taxa_blast <- blast_doc_all %>%
#  left_join(ID_list_final, by = "NCBI_ID") %>%
#  column_to_rownames("Sample_ID") 

taxa_blast <- blast_doc_all %>%
  left_join(ID_list_final, by = "NCBI_ID") %>%
  column_to_rownames("Sample_ID") %>%
  filter(Percentage >=98)

saveRDS(taxa_blast, "data/taxa_blast.RDS")

length(unique(taxa_blast$NCBI_ID)) #418

#shared between blast and tree
#asvs from pre_blast input. same for blast and tree
ph.asvs <- otu_table(asvs, taxa_are_rows = TRUE) 
#meta same
ph.meta <- sample_data(meta) # 535 observations

#bring in tree taxa data
taxa_tree <- read.table('data/ASV_family.tsv',header=TRUE) %>% column_to_rownames("X.ASV.ID") # 2540 asvs
ph.taxa.tree <- tax_table(as.matrix(taxa_tree))
ph.taxa.blast <- tax_table(as.matrix(taxa_blast))

#import as phyloseq object
phyloseq_all_tree <- phyloseq(ph.asvs, ph.meta, ph.taxa.tree)
phyloseq_all_blast <- phyloseq(ph.asvs, ph.meta, ph.taxa.blast)

# remove zero sum OTUS in this cleaned ps
phyloseq_all_tree <- prune_taxa(taxa_sums(phyloseq_all_tree) > 0, phyloseq_all_tree) 
phyloseq_all_blast <- prune_taxa(taxa_sums(phyloseq_all_blast) > 0, phyloseq_all_blast) 

ntaxa(phyloseq_all_tree) #2529
ntaxa(phyloseq_all_blast) #4817

nsamples(phyloseq_all_tree) #529
nsamples(phyloseq_all_blast) #529

# extract taxa
taxa_tree.nozero <- tax_table(phyloseq_all_tree) %>% as.data.frame()
taxa_blast.nozero <- tax_table(phyloseq_all_blast) %>% as.data.frame()

################################
###### JOIN TREE & BLAST #######
################################

taxa_tree_merge <- taxa_tree.nozero %>% rownames_to_column("ASV") %>% rename(tree_AMF_family = AMF_family)
taxa_blast_merge <- taxa_blast.nozero %>% rownames_to_column("ASV") %>% rename(blast_AMF_family = AMF_family, NCBI_Per = Percentage)
taxa_merge <- taxa_tree_merge %>% 
  full_join(taxa_blast_merge, by = "ASV")

#how many BLAST are unknown to family?
blastbyfam <- taxa_blast_merge %>% group_by(blast_AMF_family) %>% tally()
blastbyfam[11,2]/ sum(blastbyfam$n)

#how many Tree are unknown to family?
treebyfam <- taxa_tree_merge %>% group_by(tree_AMF_family) %>% tally()
treebyfam[9,2]/ sum(treebyfam$n)

#how many are found in tree, not in blast? 747
taxa_merge_treenoblast <- taxa_merge %>% filter(is.na(blast_AMF_family) & !is.na(tree_AMF_family))

#how many are found in blast, not in tree? 3035
taxa_merge_blastnotree <- taxa_merge %>% filter(is.na(tree_AMF_family) & !is.na(blast_AMF_family))

#how many are found in both as AMF? 1782
taxa_merge_overlap <- taxa_tree_merge %>% 
  inner_join(taxa_blast_merge, by = "ASV") %>%
  mutate(blasttotree_AMF_family = paste(blast_AMF_family, tree_AMF_family, sep="_"))

#what is the proportion unknown of these overlaps blast v tree
taxa_merge_overlap %>% 
  group_by(tree_AMF_family)%>%
  tally()

taxa_merge_overlap %>% 
  group_by(blast_AMF_family)%>%
  tally()

#what are the blast unknowns in the tree?
taxa_merge_overlap_unknown <- taxa_merge_overlap %>% filter(blast_AMF_family=="Unknown")
taxa_merge_overlap_unknown %>%
  group_by(tree_AMF_family) %>%
  tally()

#what are the tree versions for each family in blast?
familyassignments <- taxa_merge_overlap %>%
  group_by(blast_AMF_family,tree_AMF_family) %>%
  tally()

write.csv(familyassignments, "data/blast-tree.familyassignments.csv")

#what proportion matches?
prop.match <- familyassignments %>%
  group_by(blast_AMF_family, tree_AMF_family) %>%
  summarize(N = sum(n)) %>%
  mutate(prop = N / sum(N)) 

prop.match.true <- prop.match %>% filter(blast_AMF_family==tree_AMF_family)
write.csv(prop.match.true, "data/blast-tree.familyassignments.csv")

count.mismatch <- familyassignments %>%
  group_by(blast_AMF_family, tree_AMF_family) %>%
  summarize(N = sum(n)) %>%
  filter(!blast_AMF_family==tree_AMF_family)

#proportion unknown
count.mismatch.bu <- count.mismatch %>% filter(blast_AMF_family=="Unknown")
sum(count.mismatch.bu$N)/sum(familyassignments$n) #34%

#proportion misassigned
count.mismatch.ma <- count.mismatch %>% filter(!blast_AMF_family=="Unknown")
sum(count.mismatch.ma$N)/sum(familyassignments$n) #6%

#save list of ASVs that are double unknown
double.unknown <- taxa_merge_overlap %>% 
  filter(tree_AMF_family=="Unknown" & blast_AMF_family=="Unknown")
write.csv(double.unknown, "data/blast-tree.double.unknown.csv")

#save list of ASVs that are TREE unknown
TREE.unknown <- taxa_merge_overlap %>% 
  filter(tree_AMF_family=="Unknown")
write.csv(TREE.unknown, "data/blast-tree.TREE.unknown.csv")

#extract from fasta file
#library(Biostrings)
#library(stringr)

#DOUBLE UNKNOWNS
#read in fasta and format seqs based on study cutoffs
#AMFseqs <- readDNAStringSet('data/ASVrepseqs_clean_AMFonly.fasta')
#doubleunknowns <- AMFseqs[double.unknown$ASV]
#writeXStringSet(doubleunknowns, "data/ASVrepseqs_clean_doubleunknowns.fasta",width=10000) #write out.

#TREE UNKNOWNS
#read in fasta and format seqs based on study cutoffs
#AMFseqs <- readDNAStringSet('data/ASVrepseqs_clean_AMFonly.fasta')
#TREEunknowns <- AMFseqs[TREE.unknown$ASV]
#writeXStringSet(TREEunknowns, "data/ASVrepseqs_clean_TREEunknowns.fasta",width=10000) #write out.

################################
###### PROP TREE NO BLAST ######
################################

# take all tree identified ASVS, see what BLAST has identified.
taxa_merge_TB <- taxa_tree_merge %>% 
  left_join(taxa_blast_merge, by = "ASV") 

saveRDS(taxa_merge_TB, "data/taxa_merge_TB.rds")

# count total tree ASVS
taxa_merge_Ttotal <- taxa_merge_TB %>%
  group_by(tree_AMF_family) %>%
  tally() %>%
  rename(ttot = n)

# count unknown blast but exist in tree  
taxa_merge_ToutB <- taxa_merge_TB %>%
  filter(is.na(blast_AMF_family)) %>%
  group_by(tree_AMF_family) %>%
  tally() %>%
  rename(toutb = n)

# get proportion
prop_noBinT <- taxa_merge_Ttotal %>%
  full_join(taxa_merge_ToutB, by = "tree_AMF_family") %>%
  mutate(prop_noBT = toutb/ttot)

################################
####### PHY TREE & BLAST #######
################################

# bring in tree taxa data
taxa.blastottree <- taxa_merge_overlap %>% column_to_rownames("ASV") 
ph.taxa.blastottree <- tax_table(as.matrix(taxa.blastottree))

# import as phyloseq object
phyloseq_blasttotree <- phyloseq(ph.asvs, ph.meta, ph.taxa.blastottree)

# remove zero sum OTUS in this cleaned ps
pruned.phyloseq_blasttotree <- prune_taxa(taxa_sums(phyloseq_blasttotree) > 0, phyloseq_blasttotree)
# reextract taxa
pruned.taxa.blastottree <- tax_table(pruned.phyloseq_blasttotree) %>% as.data.frame()

# get observed for each blasttotree combination per sample
family <- unique(pruned.taxa.blastottree$blasttotree_AMF_family)
smd_blasttotree <- data.frame()
for (f in family) { 
  print(f)
  phyloseq_subset <- subset_taxa(phyloseq_blasttotree, blasttotree_AMF_family==f)
  sample_data(phyloseq_subset)$observed <- estimate_richness(phyloseq_subset, measures=c('Observed'))$Observed
  smd.subset <- as(sample_data(phyloseq_subset), 'data.frame')
  smd.subset$blasttotree_AMF_family = f
  smd_blasttotree <- rbind(smd_blasttotree, smd.subset)
}

# save 
saveRDS(smd_blasttotree, "data/smd_blasttotree.RDS")

#############################
#### BIOME SUBSETS ASVS #####
#############################

#CHOOSE DOUBLE UNKNOWN OR TREE
#DOUBLE UNKNOWN
# bring in tree taxa data; double unknown starts with 48 taxa
taxa.blastottree.dunk <- double.unknown %>% column_to_rownames("ASV") 
ph.taxa.blastottree.dunk <- tax_table(as.matrix(taxa.blastottree.dunk))

#import as phyloseq object; also 48 taxa
phyloseq_blasttotree.dunk <- phyloseq(ph.asvs, ph.meta, ph.taxa.blastottree.dunk)

# extract biome subsets
#phyloseq_blasttotree_boreal <- phyloseq_blasttotree.dunk %>% subset_samples(., Biome =="Boreal") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_trop <- phyloseq_blasttotree.dunk %>% subset_samples(., Biome =="Tropical") %>% prune_taxa(taxa_sums(.) > 0, .) 
phyloseq_blasttotree_temp <- phyloseq_blasttotree.dunk %>% subset_samples(., Biome =="Temperate") %>% prune_taxa(taxa_sums(.) > 0, .)

# extract temp land use
phyloseq_blasttotree_temp.ag <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant == "Disturbed") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant == "Post-ag") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant == "Remnant") %>% prune_taxa(taxa_sums(.) > 0, .)

# extract temp at depth by land use
phyloseq_blasttotree_temp.ag.5 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.15 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.30 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="30") %>% prune_taxa(taxa_sums(.) > 0, .)

phyloseq_blasttotree_temp.postag.5 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.15 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.30 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="30") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.75 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="75") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.120 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="120") %>% prune_taxa(taxa_sums(.) > 0, .)

phyloseq_blasttotree_temp.remn.5 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.15 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)

#taxa.blasttotree_boreal <- tax_table(phyloseq_blasttotree_boreal) %>% as.data.frame()
taxa.blasttotree_trop <- tax_table(phyloseq_blasttotree_trop) %>% as.data.frame()
taxa.blasttotree_temp <- tax_table(phyloseq_blasttotree_temp) %>% as.data.frame()

taxa.blasttotree_temp.ag.5 <- tax_table(phyloseq_blasttotree_temp.ag.5) %>% as.data.frame()
taxa.blasttotree_temp.ag.15 <- tax_table(phyloseq_blasttotree_temp.ag.15) %>% as.data.frame()
taxa.blasttotree_temp.ag.30 <- tax_table(phyloseq_blasttotree_temp.ag.30) %>% as.data.frame()

taxa.blasttotree_temp.postag.5 <- tax_table(phyloseq_blasttotree_temp.postag.5) %>% as.data.frame()
taxa.blasttotree_temp.postag.15 <- tax_table(phyloseq_blasttotree_temp.postag.15) %>% as.data.frame()
taxa.blasttotree_temp.postag.30 <- tax_table(phyloseq_blasttotree_temp.postag.30) %>% as.data.frame()
taxa.blasttotree_temp.postag.75 <- tax_table(phyloseq_blasttotree_temp.postag.75) %>% as.data.frame()
taxa.blasttotree_temp.postag.120 <- tax_table(phyloseq_blasttotree_temp.postag.120) %>% as.data.frame()

taxa.blasttotree_temp.remn.5 <- tax_table(phyloseq_blasttotree_temp.remn.5) %>% as.data.frame()
taxa.blasttotree_temp.remn.15 <- tax_table(phyloseq_blasttotree_temp.remn.15) %>% as.data.frame()

#rownames(taxa.blasttotree_boreal)
rnamestrop <- rownames(taxa.blasttotree_trop) %>% as.list
#"ASV_1786"  "ASV_2373"  "ASV_3734"  "ASV_3755"  "ASV_4173"  "ASV_6320"  "ASV_10244" "ASV_11042" "ASV_17461" "ASV_18754" "ASV_23553"
#"ASV_25219" "ASV_26206" "ASV_30814" "ASV_32361" "ASV_37976" "ASV_38043" "ASV_40327" "ASV_53895"
rnamestemp <- rownames(taxa.blasttotree_temp) %>% as.list
# "ASV_6379"  "ASV_6941"  "ASV_11663" "ASV_11809" "ASV_13092" "ASV_13347" "ASV_14555" "ASV_14651" "ASV_15108" "ASV_22229" "ASV_23061"
# "ASV_24456" "ASV_26322" "ASV_30961" "ASV_36559" "ASV_38495" "ASV_41126" "ASV_44264" "ASV_46933" "ASV_48033" "ASV_48750" "ASV_49440"
# "ASV_51364" "ASV_52723" "ASV_54375" "ASV_56133" "ASV_60333" "ASV_60377" "ASV_61072"

#which overlap?
common_rnames <- intersect(rnamestrop, rnamestemp)

rownames(taxa.blasttotree_temp.ag.5) # ASV_13347
rownames(taxa.blasttotree_temp.ag.15) # ASV_6941; ASV_24456
rownames(taxa.blasttotree_temp.ag.30) # ASV_61072

rownames(taxa.blasttotree_temp.postag.5) # ASV_6379; ASV_6941; ASV_11663; ASV_11809;ASV_13092; ASV_13347; ASV_14555; 
#ASV_14651; ASV_15108; ASV_22229; ASV_23061; ASV_24456; ASV_26322; ASV_44264; ASV_46933; ASV_54375; ASV_60333
rownames(taxa.blasttotree_temp.postag.15) # ASV_6379; ASV_6941; ASV_11663; ASV_13092; ASV_13347; ASV_14651; 
#ASV_15108; ASV_24456; ASV_26322; ASV_36559; ASV_41126; ASV_51364
rownames(taxa.blasttotree_temp.postag.30) # ASV_6379; ASV_6941; ASV_11663; ASV_11809; ASV_13092; ASV_13347;
#ASV_14555; ASV_14651; ASV_15108; ASV_22229; ASV_23061; ASV_24456; ASV_26322; ASV_30961; ASV_49440; ASV_60377

rownames(taxa.blasttotree_temp.postag.75) # ASV_11663
rownames(taxa.blasttotree_temp.postag.120) # ASV_15108; ASV_48033

rownames(taxa.blasttotree_temp.remn.5) #ASV_38495; ASV_48750; ASV_52723
rownames(taxa.blasttotree_temp.remn.15) #ASV_38495; ASV_56133

#TREE UNKNOWN
# bring in tree taxa data; tree unknown starts with 37 taxa
taxa.blastottree.treeunk <- TREE.unknown %>% column_to_rownames("ASV") 
ph.taxa.blastottree.treeunk <- tax_table(as.matrix(taxa.blastottree.treeunk))

#import as phyloseq object; also 17 taxa
phyloseq_blasttotree.treeunk <- phyloseq(ph.asvs, ph.meta, ph.taxa.blastottree.treeunk)

# extract biome subsets
phyloseq_blasttotree_boreal <- phyloseq_blasttotree.treeunk %>% subset_samples(., Biome =="Boreal") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_trop <- phyloseq_blasttotree.treeunk %>% subset_samples(., Biome =="Tropical") %>% prune_taxa(taxa_sums(.) > 0, .) 
phyloseq_blasttotree_temp <- phyloseq_blasttotree.treeunk %>% subset_samples(., Biome =="Temperate") %>% prune_taxa(taxa_sums(.) > 0, .)

# extract temp land use
phyloseq_blasttotree_temp.ag <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant = "Disturbed") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant = "Post-ag") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn <- phyloseq_blasttotree_temp %>% subset_samples(., Remnant = "Remnant") %>% prune_taxa(taxa_sums(.) > 0, .)

# extract temp at depth by land use
phyloseq_blasttotree_temp.ag.5 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.15 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.30 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="30") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.75 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="75") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.120 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="120") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.ag.150 <- phyloseq_blasttotree_temp.ag %>% subset_samples(., Depth =="150") %>% prune_taxa(taxa_sums(.) > 0, .)

phyloseq_blasttotree_temp.postag.5 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.15 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.30 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="30") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.75 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="75") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.120 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="120") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.postag.150 <- phyloseq_blasttotree_temp.postag %>% subset_samples(., Depth =="150") %>% prune_taxa(taxa_sums(.) > 0, .)

phyloseq_blasttotree_temp.remn.5 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="5") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.15 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="15") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.30 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="30") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.75 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="75") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.120 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="120") %>% prune_taxa(taxa_sums(.) > 0, .)
phyloseq_blasttotree_temp.remn.150 <- phyloseq_blasttotree_temp.remn %>% subset_samples(., Depth =="150") %>% prune_taxa(taxa_sums(.) > 0, .)

taxa.blasttotree_boreal <- tax_table(phyloseq_blasttotree_boreal) %>% as.data.frame()
taxa.blasttotree_trop <- tax_table(phyloseq_blasttotree_trop) %>% as.data.frame()

taxa.blasttotree_temp <- tax_table(phyloseq_blasttotree_temp) %>% as.data.frame()

rownames(taxa.blasttotree_boreal)
#"ASV_1569"
rownames(taxa.blasttotree_trop)
# "ASV_1786"  "ASV_2373"  "ASV_3734"  "ASV_3755"  "ASV_4173"  "ASV_6320"  "ASV_10244" "ASV_11042" "ASV_17461" "ASV_18754" "ASV_19207"
# "ASV_23553" "ASV_24306" "ASV_25219" "ASV_26206" "ASV_30814" "ASV_32361" "ASV_37976" "ASV_38043" "ASV_40327" "ASV_53732" "ASV_53895"
# "ASV_58530"
rownames(taxa.blasttotree_temp)
#"ASV_1300"  "ASV_2124"  "ASV_3179"  "ASV_3811"  "ASV_4586"  "ASV_4590"  "ASV_5005"  "ASV_6379"  "ASV_6630"  "ASV_6941"  "ASV_9201" 
#"ASV_9202"  "ASV_10502" "ASV_10796" "ASV_11663" "ASV_11809" "ASV_12717" "ASV_13092" "ASV_13347" "ASV_13604" "ASV_14555" "ASV_14651"
# "ASV_15108" "ASV_15953" "ASV_20663" "ASV_22229" "ASV_23061" "ASV_24438" "ASV_24456" "ASV_26094" "ASV_26322" "ASV_27924" "ASV_29961"
# "ASV_30961" "ASV_33085" "ASV_34746" "ASV_36558" "ASV_36559" "ASV_38495" "ASV_39160" "ASV_41126" "ASV_41130" "ASV_41737" "ASV_44264"
# "ASV_46658" "ASV_46703" "ASV_46933" "ASV_48033" "ASV_48561" "ASV_48750" "ASV_48894" "ASV_49440" "ASV_51121" "ASV_51364" "ASV_52723"
# "ASV_52724" "ASV_52855" "ASV_54375" "ASV_56133" "ASV_58992" "ASV_60333" "ASV_60377" "ASV_60688" "ASV_61072" "ASV_61143"

rownames(taxa.blasttotreetemp.ag.5)
rownames(taxa.blasttotreetemp.ag.15)
rownames(taxa.blasttotreetemp.ag.30)
rownames(taxa.blasttotreetemp.ag.75)
rownames(taxa.blasttotreetemp.ag.120)
rownames(taxa.blasttotreetemp.ag.150)

rownames(taxa.blasttotreetemp.postag.5)
rownames(taxa.blasttotreetemp.postag.15)
rownames(taxa.blasttotreetemp.postag.30)
rownames(taxa.blasttotreetemp.postag.75)
rownames(taxa.blasttotreetemp.postag.120)
rownames(taxa.blasttotreetemp.postag.150)

rownames(taxa.blasttotreetemp.remn.5)
rownames(taxa.blasttotreetemp.remn.15)
rownames(taxa.blasttotreetemp.remn.30)
rownames(taxa.blasttotreetemp.remn.75)
rownames(taxa.blasttotreetemp.remn.120)
rownames(taxa.blasttotreetemp.remn.150)

