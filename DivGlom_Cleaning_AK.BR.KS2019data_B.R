## Script for :
## 1. KS depth 15, 2019
## 2. tree v blast diversity calculations
## 3. tree v blast prop calculations

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
  filter(!str_detect(Sample, "Spr18|Fall18")) %>%
  #remove KS depths other than 15
  #or not if want to keep for unknowns
  filter(!(Depth > 20 & State == "KS")) %>% filter(!(Depth < 10 & State == "KS")) %>%
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

write.csv(meta, "data/cleaned_meta_dat_AK.BR.KS2019.csv")

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

length(unique(taxa_blast$NCBI_ID)) #418

#shared between blast and tree
#asvs from pre_blast input. same for blast and tree
ph.asvs <- otu_table(asvs, taxa_are_rows = TRUE) 
#meta same
ph.meta <- sample_data(meta) # 102 observations

#bring in tree taxa data
#merge taxa tree with blast AMF family

taxa_blast_famknown <- taxa_blast %>% 
  mutate(blast_fam_known = ifelse(AMF_family == "Unknown", 'no','yes')) %>%
  rownames_to_column(., "X.ASV.ID") %>%
  select(c(X.ASV.ID, blast_fam_known)) # 20263 asvs
taxa_tree <- read.table('data/ASV_family.tsv',header=TRUE) %>%
  left_join(taxa_blast_famknown, by = "X.ASV.ID") %>% 
  column_to_rownames("X.ASV.ID") # 2540 asvs

ph.taxa.tree <- tax_table(as.matrix(taxa_tree))
ph.taxa.blast <- tax_table(as.matrix(taxa_blast))

#import as phyloseq object
phyloseq_all_tree <- phyloseq(ph.asvs, ph.meta, ph.taxa.tree)
phyloseq_all_blast <- phyloseq(ph.asvs, ph.meta, ph.taxa.blast)

# remove zero sum OTUS in this cleaned ps
phyloseq_all_tree <- prune_taxa(taxa_sums(phyloseq_all_tree) > 0, phyloseq_all_tree) 
phyloseq_all_blast <- prune_taxa(taxa_sums(phyloseq_all_blast) > 0, phyloseq_all_blast) 

ntaxa(phyloseq_all_tree) #604
ntaxa(phyloseq_all_blast) #948

nsamples(phyloseq_all_tree) #102
nsamples(phyloseq_all_blast) #102

# extract taxa
taxa_tree.nozero <- tax_table(phyloseq_all_tree) %>% as.data.frame()
taxa_blast.nozero <- tax_table(phyloseq_all_blast) %>% as.data.frame()

################################
###### DIVERSITY METRICS #######
################################

sample_data(phyloseq_all_tree)$shannon <- estimate_richness(phyloseq_all_tree, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree)$observed <- estimate_richness(phyloseq_all_tree, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree)$chao1 <- estimate_richness(phyloseq_all_tree, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree)$invsimpson <- estimate_richness(phyloseq_all_tree, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree)$reads_sample <- sample_sums(phyloseq_all_tree)

sample_data(phyloseq_all_blast)$shannon <- estimate_richness(phyloseq_all_blast, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast)$observed <- estimate_richness(phyloseq_all_blast, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast)$chao1 <- estimate_richness(phyloseq_all_blast, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast)$invsimpson <- estimate_richness(phyloseq_all_blast, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast)$reads_sample <- sample_sums(phyloseq_all_blast)

# save 
saveRDS(phyloseq_all_tree, "data/phyloseq_all_tree_AK.BR.KS2019.RDS")
saveRDS(phyloseq_all_blast, "data/phyloseq_all_blast_AK.BR.KS2019.RDS")

################################
###### DIVERSITY METRICS #######
########## BY FAMILY ###########
################################

################################
########## TREE DATA ###########
################################

####### UNKNOWN FAMILY #######
phyloseq_all_tree.unknown <- subset_taxa(phyloseq_all_tree, AMF_family=="Unknown")

sample_data(phyloseq_all_tree.unknown)$shannon <- estimate_richness(phyloseq_all_tree.unknown, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.unknown)$observed <- estimate_richness(phyloseq_all_tree.unknown, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.unknown)$chao1 <- estimate_richness(phyloseq_all_tree.unknown, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.unknown)$invsimpson <- estimate_richness(phyloseq_all_tree.unknown, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.unknown)$reads_sample <- sample_sums(phyloseq_all_tree.unknown)

smd.unknown <- as(sample_data(phyloseq_all_tree.unknown), 'data.frame')
smd.unknown$AMF_family = "Unknown"

######### GLOMERACEAE ########
phyloseq_all_tree.glo <- subset_taxa(phyloseq_all_tree, AMF_family=="Glomeraceae")

sample_data(phyloseq_all_tree.glo)$shannon <- estimate_richness(phyloseq_all_tree.glo, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.glo)$observed <- estimate_richness(phyloseq_all_tree.glo, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.glo)$chao1 <- estimate_richness(phyloseq_all_tree.glo, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.glo)$invsimpson <- estimate_richness(phyloseq_all_tree.glo, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.glo)$reads_sample <- sample_sums(phyloseq_all_tree.glo)

smd.glo <- as(sample_data(phyloseq_all_tree.glo), 'data.frame')
smd.glo$AMF_family = "Glomeraceae"

######### CLAROIDEOGLOMERACEAE ######## 
phyloseq_all_tree.cla <- subset_taxa(phyloseq_all_tree, AMF_family=="Claroideoglomeraceae")

sample_data(phyloseq_all_tree.cla)$shannon <- estimate_richness(phyloseq_all_tree.cla, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.cla)$observed <- estimate_richness(phyloseq_all_tree.cla, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.cla)$chao1 <- estimate_richness(phyloseq_all_tree.cla, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.cla)$invsimpson <- estimate_richness(phyloseq_all_tree.cla, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.cla)$reads_sample <- sample_sums(phyloseq_all_tree.cla)

smd.cla <- as(sample_data(phyloseq_all_tree.cla), 'data.frame')
smd.cla$AMF_family = "Claroideoglomeraceae"

######### ACAULOSPORACEAE ######## 
phyloseq_all_tree.aca <- subset_taxa(phyloseq_all_tree, AMF_family=="Acaulosporaceae")

sample_data(phyloseq_all_tree.aca)$shannon <- estimate_richness(phyloseq_all_tree.aca, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.aca)$observed <- estimate_richness(phyloseq_all_tree.aca, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.aca)$chao1 <- estimate_richness(phyloseq_all_tree.aca, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.aca)$invsimpson <- estimate_richness(phyloseq_all_tree.aca, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.aca)$reads_sample <- sample_sums(phyloseq_all_tree.aca)

smd.aca <- as(sample_data(phyloseq_all_tree.aca), 'data.frame')
smd.aca$AMF_family = "Acaulosporaceae"

######### PERVESTUSTACEAE ######## 
#phyloseq_all_tree.per <- subset_taxa(phyloseq_all_tree, AMF_family == "Pervestustaceae")

#sample_data(phyloseq_all_tree.per)$shannon <- estimate_richness(phyloseq_all_tree.per, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree.per)$observed <- estimate_richness(phyloseq_all_tree.per, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree.per)$chao1 <- estimate_richness(phyloseq_all_tree.per, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree.per)$invsimpson <- estimate_richness(phyloseq_all_tree.per, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree.per)$reads_sample <- sample_sums(phyloseq_all_tree.per)

#smd.per <- as(sample_data(phyloseq_all_tree.per), 'data.frame')
#smd.per$AMF_family = "Pervestustaceae"

######### AMBISPORACEAE ######## 
phyloseq_all_tree.amb <- subset_taxa(phyloseq_all_tree, AMF_family=="Ambisporaceae")

sample_data(phyloseq_all_tree.amb)$shannon <- estimate_richness(phyloseq_all_tree.amb, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.amb)$observed <- estimate_richness(phyloseq_all_tree.amb, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.amb)$chao1 <- estimate_richness(phyloseq_all_tree.amb, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.amb)$invsimpson <- estimate_richness(phyloseq_all_tree.amb, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.amb)$reads_sample <- sample_sums(phyloseq_all_tree.amb)

smd.amb <- as(sample_data(phyloseq_all_tree.amb), 'data.frame')
smd.amb$AMF_family = "Ambisporaceae"

######### DIVERSIPORACEAE ######## 
phyloseq_all_tree.div <- subset_taxa(phyloseq_all_tree, AMF_family=="Diversisporaceae")

sample_data(phyloseq_all_tree.div)$shannon <- estimate_richness(phyloseq_all_tree.div, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.div)$observed <- estimate_richness(phyloseq_all_tree.div, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.div)$chao1 <- estimate_richness(phyloseq_all_tree.div, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.div)$invsimpson <- estimate_richness(phyloseq_all_tree.div, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.div)$reads_sample <- sample_sums(phyloseq_all_tree.div)

smd.div <- as(sample_data(phyloseq_all_tree.div), 'data.frame')
smd.div$AMF_family = "Diversisporaceae"

######### GIGASPORACEAE ######## 
phyloseq_all_tree.gig <- subset_taxa(phyloseq_all_tree, AMF_family=="Gigasporaceae")

sample_data(phyloseq_all_tree.gig)$shannon <- estimate_richness(phyloseq_all_tree.gig, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.gig)$observed <- estimate_richness(phyloseq_all_tree.gig, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.gig)$chao1 <- estimate_richness(phyloseq_all_tree.gig, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.gig)$invsimpson <- estimate_richness(phyloseq_all_tree.gig, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.gig)$reads_sample <- sample_sums(phyloseq_all_tree.gig)

smd.gig <- as(sample_data(phyloseq_all_tree.gig), 'data.frame')
smd.gig$AMF_family = "Gigasporaceae"

######### PARAGLOMERACEAE ######## 
phyloseq_all_tree.par <- subset_taxa(phyloseq_all_tree, AMF_family=="Paraglomeraceae")

sample_data(phyloseq_all_tree.par)$shannon <- estimate_richness(phyloseq_all_tree.par, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree.par)$observed <- estimate_richness(phyloseq_all_tree.par, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree.par)$chao1 <- estimate_richness(phyloseq_all_tree.par, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree.par)$invsimpson <- estimate_richness(phyloseq_all_tree.par, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree.par)$reads_sample <- sample_sums(phyloseq_all_tree.par)

smd.par <- as(sample_data(phyloseq_all_tree.par), 'data.frame')
smd.par$AMF_family = "Paraglomeraceae"

######### ARCHAEOSPORACEAE ######## 
#phyloseq_all_tree.arc <- subset_taxa(phyloseq_all_tree, AMF_family=="Archaesporaceae")

#sample_data(phyloseq_all_tree.arc)$shannon <- estimate_richness(phyloseq_all_tree.arc, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree.arc)$observed <- estimate_richness(phyloseq_all_tree.arc, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree.arc)$chao1 <- estimate_richness(phyloseq_all_tree.arc, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree.arc)$invsimpson <- estimate_richness(phyloseq_all_tree.arc, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree.arc)$reads_sample <- sample_sums(phyloseq_all_tree.arc)

#smd.arc <- as(sample_data(phyloseq_all_tree.arc), 'data.frame')
#smd.arc$AMF_family = "Archaesporaceae"

######### PACISPORACEAE ######## 
#phyloseq_all_tree.paca <- subset_taxa(phyloseq_all_tree, AMF_family=="Pacisporaceae")

#sample_data(phyloseq_all_tree.paca)$shannon <- estimate_richness(phyloseq_all_tree.paca, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree.paca)$observed <- estimate_richness(phyloseq_all_tree.paca, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree.paca)$chao1 <- estimate_richness(phyloseq_all_tree.paca, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree.paca)$invsimpson <- estimate_richness(phyloseq_all_tree.paca, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree.paca)$reads_sample <- sample_sums(phyloseq_all_tree.paca)

#smd.paca <- as(sample_data(phyloseq_all_tree.paca), 'data.frame')
#smd.paca$AMF_family = "Pacisporaceae"

######### SACCULOSPORACEAE ######## 
#phyloseq_all_tree.sac <- subset_taxa(phyloseq_all_tree, AMF_family=="Sacculosporaceae")

#sample_data(phyloseq_all_tree.sac)$shannon <- estimate_richness(phyloseq_all_tree.sac, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree.sac)$observed <- estimate_richness(phyloseq_all_tree.sac, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree.sac)$chao1 <- estimate_richness(phyloseq_all_tree.sac, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree.sac)$invsimpson <- estimate_richness(phyloseq_all_tree.sac, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree.sac)$reads_sample <- sample_sums(phyloseq_all_tree.sac)

#smd.sac <- as(sample_data(phyloseq_all_tree.sac), 'data.frame')
#smd.sac$AMF_family = "Sacculosporaceae"

################################
########## TREE DATA ###########
####### UNKNOWN TO BLAST #######
######### KNOWN TREE ###########
################################

######### GLOMERACEAE ########
phyloseq_all_tree_bunknown.glo <- subset_taxa(phyloseq_all_tree, AMF_family=="Glomeraceae" & blast_fam_known == "no")

sample_data(phyloseq_all_tree_bunknown.glo)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.glo, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree_bunknown.glo)$observed <- estimate_richness(phyloseq_all_tree_bunknown.glo, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree_bunknown.glo)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.glo, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree_bunknown.glo)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.glo, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree_bunknown.glo)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.glo)

smd.glo.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.glo), 'data.frame')
smd.glo.bunknown$AMF_family = "Glomeraceae_bunknown"

######### CLAROIDEOGLOMERACEAE ######## 
phyloseq_all_tree_bunknown.cla <- subset_taxa(phyloseq_all_tree, AMF_family=="Claroideoglomeraceae" & blast_fam_known == "no")

sample_data(phyloseq_all_tree_bunknown.cla)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.cla, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree_bunknown.cla)$observed <- estimate_richness(phyloseq_all_tree_bunknown.cla, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree_bunknown.cla)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.cla, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree_bunknown.cla)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.cla, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree_bunknown.cla)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.cla)

smd.cla.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.cla), 'data.frame')
smd.cla.bunknown$AMF_family = "Claroideoglomeraceae_bunknown"

######### ACAULOSPORACEAE ######## 
phyloseq_all_tree_bunknown.aca <- subset_taxa(phyloseq_all_tree, AMF_family=="Acaulosporaceae" & blast_fam_known == "no")

sample_data(phyloseq_all_tree_bunknown.aca)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.aca, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree_bunknown.aca)$observed <- estimate_richness(phyloseq_all_tree_bunknown.aca, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree_bunknown.aca)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.aca, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree_bunknown.aca)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.aca, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree_bunknown.aca)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.aca)

smd.aca.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.aca), 'data.frame')
smd.aca.bunknown$AMF_family = "Acaulosporaceae_bunknown"

######### PERVESTUSTACEAE ######## 
#phyloseq_all_tree_bunknown.per <- subset_taxa(phyloseq_all_tree_bunknown, AMF_family == "Pervestustaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.per)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.per, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.per)$observed <- estimate_richness(phyloseq_all_tree_bunknown.per, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.per)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.per, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.per)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.per, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.per)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.per)

#smd.per.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.per), 'data.frame')
#smd.per.bunknown$AMF_family = "Pervestustaceae_bunknown"

######### AMBISPORACEAE ######## 
#phyloseq_all_tree_bunknown.amb <- subset_taxa(phyloseq_all_tree, AMF_family=="Ambisporaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.amb)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.amb, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.amb)$observed <- estimate_richness(phyloseq_all_tree_bunknown.amb, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.amb)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.amb, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.amb)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.amb, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.amb)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.amb)

#smd.amb.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.amb), 'data.frame')
#smd.amb.bunknown$AMF_family = "Ambisporaceae_bunknown"

######### DIVERSIPORACEAE ######## 
phyloseq_all_tree_bunknown.div <- subset_taxa(phyloseq_all_tree, AMF_family=="Diversisporaceae" & blast_fam_known == "no")

sample_data(phyloseq_all_tree_bunknown.div)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.div, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree_bunknown.div)$observed <- estimate_richness(phyloseq_all_tree_bunknown.div, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree_bunknown.div)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.div, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree_bunknown.div)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.div, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree_bunknown.div)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.div)

smd.div.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.div), 'data.frame')
smd.div.bunknown$AMF_family = "Diversisporaceae_bunknown"

######### GIGASPORACEAE ######## 
#phyloseq_all_tree_bunknown.gig <- subset_taxa(phyloseq_all_tree, AMF_family=="Gigasporaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.gig)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.gig, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.gig)$observed <- estimate_richness(phyloseq_all_tree_bunknown.gig, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.gig)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.gig, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.gig)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.gig, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.gig)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.gig)

#smd.gig.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.gig), 'data.frame')
#smd.gig.bunknown$AMF_family = "Gigasporaceae_bunknown"

######### PARAGLOMERACEAE ######## 
phyloseq_all_tree_bunknown.par <- subset_taxa(phyloseq_all_tree, AMF_family=="Paraglomeraceae" & blast_fam_known == "no")

sample_data(phyloseq_all_tree_bunknown.par)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.par, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_tree_bunknown.par)$observed <- estimate_richness(phyloseq_all_tree_bunknown.par, measures=c('Observed'))$Observed
sample_data(phyloseq_all_tree_bunknown.par)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.par, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_tree_bunknown.par)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.par, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_tree_bunknown.par)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.par)

smd.par.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.par), 'data.frame')
smd.par.bunknown$AMF_family = "Paraglomeraceae_bunknown"

######### ARCHAEOSPORACEAE ######## 
#phyloseq_all_tree_bunknown.arc <- subset_taxa(phyloseq_all_tree_bunknown, AMF_family=="Archaesporaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.arc)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.arc, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.arc)$observed <- estimate_richness(phyloseq_all_tree_bunknown.arc, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.arc)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.arc, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.arc)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.arc, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.arc)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.arc)

#smd.arc.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.arc), 'data.frame')
#smd.arc.bunknown$AMF_family = "Archaesporaceae_bunknown"

######### PACISPORACEAE ######## 
#phyloseq_all_tree_bunknown.paca <- subset_taxa(phyloseq_all_tree_bunknown, AMF_family=="Pacisporaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.paca)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.paca, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.paca)$observed <- estimate_richness(phyloseq_all_tree_bunknown.paca, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.paca)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.paca, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.paca)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.paca, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.paca)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.paca)

#smd.paca.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.paca), 'data.frame')
#smd.paca.bunknown$AMF_family = "Pacisporaceae_bunknown"

######### SACCULOSPORACEAE ######## 
#phyloseq_all_tree_bunknown.sac <- subset_taxa(phyloseq_all_tree_bunknown, AMF_family=="Sacculosporaceae" & blast_fam_known == "no")

#sample_data(phyloseq_all_tree_bunknown.sac)$shannon <- estimate_richness(phyloseq_all_tree_bunknown.sac, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_tree_bunknown.sac)$observed <- estimate_richness(phyloseq_all_tree_bunknown.sac, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_tree_bunknown.sac)$chao1 <- estimate_richness(phyloseq_all_tree_bunknown.sac, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_tree_bunknown.sac)$invsimpson <- estimate_richness(phyloseq_all_tree_bunknown.sac, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_tree_bunknown.sac)$reads_sample <- sample_sums(phyloseq_all_tree_bunknown.sac)

#smd.sac.bunknown <- as(sample_data(phyloseq_all_tree_bunknown.sac), 'data.frame')
#smd.sac.bunknown$AMF_family = "Sacculosporaceae_bunknown"

###### JOIN ######

smd_all_tree_fam <- rbind(smd.aca, smd.amb, smd.cla, smd.div, smd.gig, smd.glo, smd.par, smd.unknown,
                          smd.aca.bunknown,  smd.cla.bunknown, smd.div.bunknown, smd.glo.bunknown, smd.par.bunknown)

# save 
saveRDS(smd_all_tree_fam, "data/smd_all_tree_fam_AK.BR.KS2019.RDS")

################################
########## BLAST DATA ##########
################################

####### UNKNOWN FAMILY #######
phyloseq_all_blast.unknown <- subset_taxa(phyloseq_all_blast, AMF_family == "Unknown")

sample_data(phyloseq_all_blast.unknown)$shannon <- estimate_richness(phyloseq_all_blast.unknown, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.unknown)$observed <- estimate_richness(phyloseq_all_blast.unknown, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.unknown)$chao1 <- estimate_richness(phyloseq_all_blast.unknown, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.unknown)$invsimpson <- estimate_richness(phyloseq_all_blast.unknown, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.unknown)$reads_sample <- sample_sums(phyloseq_all_blast.unknown)

smd.unknown <- as(sample_data(phyloseq_all_blast.unknown), 'data.frame')
smd.unknown$AMF_family = "Unknown"

######### GLOMERACEAE ########
phyloseq_all_blast.glo <- subset_taxa(phyloseq_all_blast, AMF_family=="Glomeraceae")

sample_data(phyloseq_all_blast.glo)$shannon <- estimate_richness(phyloseq_all_blast.glo, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.glo)$observed <- estimate_richness(phyloseq_all_blast.glo, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.glo)$chao1 <- estimate_richness(phyloseq_all_blast.glo, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.glo)$invsimpson <- estimate_richness(phyloseq_all_blast.glo, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.glo)$reads_sample <- sample_sums(phyloseq_all_blast.glo)

smd.glo <- as(sample_data(phyloseq_all_blast.glo), 'data.frame')
smd.glo$AMF_family = "Glomeraceae"

######### CLAROIDEOGLOMERACEAE ######## 
phyloseq_all_blast.cla <- subset_taxa(phyloseq_all_blast, AMF_family=="Claroideoglomeraceae")

sample_data(phyloseq_all_blast.cla)$shannon <- estimate_richness(phyloseq_all_blast.cla, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.cla)$observed <- estimate_richness(phyloseq_all_blast.cla, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.cla)$chao1 <- estimate_richness(phyloseq_all_blast.cla, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.cla)$invsimpson <- estimate_richness(phyloseq_all_blast.cla, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.cla)$reads_sample <- sample_sums(phyloseq_all_blast.cla)

smd.cla <- as(sample_data(phyloseq_all_blast.cla), 'data.frame')
smd.cla$AMF_family = "Claroideoglomeraceae"

######### ACAULOSPORACEAE ######## 
phyloseq_all_blast.aca <- subset_taxa(phyloseq_all_blast, AMF_family=="Acaulosporaceae")

sample_data(phyloseq_all_blast.aca)$shannon <- estimate_richness(phyloseq_all_blast.aca, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.aca)$observed <- estimate_richness(phyloseq_all_blast.aca, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.aca)$chao1 <- estimate_richness(phyloseq_all_blast.aca, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.aca)$invsimpson <- estimate_richness(phyloseq_all_blast.aca, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.aca)$reads_sample <- sample_sums(phyloseq_all_blast.aca)

smd.aca <- as(sample_data(phyloseq_all_blast.aca), 'data.frame')
smd.aca$AMF_family = "Acaulosporaceae"

######### PERVESTUSTACEAE ######## 
phyloseq_all_blast.per <- subset_taxa(phyloseq_all_blast, AMF_family == "Pervestustaceae")

sample_data(phyloseq_all_blast.per)$shannon <- estimate_richness(phyloseq_all_blast.per, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.per)$observed <- estimate_richness(phyloseq_all_blast.per, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.per)$chao1 <- estimate_richness(phyloseq_all_blast.per, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.per)$invsimpson <- estimate_richness(phyloseq_all_blast.per, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.per)$reads_sample <- sample_sums(phyloseq_all_blast.per)

smd.per <- as(sample_data(phyloseq_all_blast.per), 'data.frame')
smd.per$AMF_family = "Pervestustaceae"

######### AMBISPORACEAE ######## 
phyloseq_all_blast.amb <- subset_taxa(phyloseq_all_blast, AMF_family=="Ambisporaceae")

sample_data(phyloseq_all_blast.amb)$shannon <- estimate_richness(phyloseq_all_blast.amb, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.amb)$observed <- estimate_richness(phyloseq_all_blast.amb, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.amb)$chao1 <- estimate_richness(phyloseq_all_blast.amb, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.amb)$invsimpson <- estimate_richness(phyloseq_all_blast.amb, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.amb)$reads_sample <- sample_sums(phyloseq_all_blast.amb)

smd.amb <- as(sample_data(phyloseq_all_blast.amb), 'data.frame')
smd.amb$AMF_family = "Ambisporaceae"

######### DIVERSIPORACEAE ######## 
phyloseq_all_blast.div <- subset_taxa(phyloseq_all_blast, AMF_family=="Diversisporaceae")

sample_data(phyloseq_all_blast.div)$shannon <- estimate_richness(phyloseq_all_blast.div, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.div)$observed <- estimate_richness(phyloseq_all_blast.div, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.div)$chao1 <- estimate_richness(phyloseq_all_blast.div, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.div)$invsimpson <- estimate_richness(phyloseq_all_blast.div, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.div)$reads_sample <- sample_sums(phyloseq_all_blast.div)

smd.div <- as(sample_data(phyloseq_all_blast.div), 'data.frame')
smd.div$AMF_family = "Diversisporaceae"

######### GIGASPORACEAE ######## 
phyloseq_all_blast.gig <- subset_taxa(phyloseq_all_blast, AMF_family=="Gigasporaceae")

sample_data(phyloseq_all_blast.gig)$shannon <- estimate_richness(phyloseq_all_blast.gig, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.gig)$observed <- estimate_richness(phyloseq_all_blast.gig, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.gig)$chao1 <- estimate_richness(phyloseq_all_blast.gig, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.gig)$invsimpson <- estimate_richness(phyloseq_all_blast.gig, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.gig)$reads_sample <- sample_sums(phyloseq_all_blast.gig)

smd.gig <- as(sample_data(phyloseq_all_blast.gig), 'data.frame')
smd.gig$AMF_family = "Gigasporaceae"

######### PARAGLOMERACEAE ######## 
phyloseq_all_blast.par <- subset_taxa(phyloseq_all_blast, AMF_family=="Paraglomeraceae")

sample_data(phyloseq_all_blast.par)$shannon <- estimate_richness(phyloseq_all_blast.par, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.par)$observed <- estimate_richness(phyloseq_all_blast.par, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.par)$chao1 <- estimate_richness(phyloseq_all_blast.par, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.par)$invsimpson <- estimate_richness(phyloseq_all_blast.par, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.par)$reads_sample <- sample_sums(phyloseq_all_blast.par)

smd.par <- as(sample_data(phyloseq_all_blast.par), 'data.frame')
smd.par$AMF_family = "Paraglomeraceae"

######### ARCHAEOSPORACEAE ######## 
phyloseq_all_blast.arc <- subset_taxa(phyloseq_all_blast, AMF_family=="Archaesporaceae")

sample_data(phyloseq_all_blast.arc)$shannon <- estimate_richness(phyloseq_all_blast.arc, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.arc)$observed <- estimate_richness(phyloseq_all_blast.arc, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.arc)$chao1 <- estimate_richness(phyloseq_all_blast.arc, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.arc)$invsimpson <- estimate_richness(phyloseq_all_blast.arc, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.arc)$reads_sample <- sample_sums(phyloseq_all_blast.arc)

smd.arc <- as(sample_data(phyloseq_all_blast.arc), 'data.frame')
smd.arc$AMF_family = "Archaesporaceae"

######### PACISPORACEAE ######## 
#phyloseq_all_blast.paca <- subset_taxa(phyloseq_all_blast, AMF_family=="Pacisporaceae")

#sample_data(phyloseq_all_blast.paca)$shannon <- estimate_richness(phyloseq_all_blast.paca, measures=c('Shannon'))$Shannon
#sample_data(phyloseq_all_blast.paca)$observed <- estimate_richness(phyloseq_all_blast.paca, measures=c('Observed'))$Observed
#sample_data(phyloseq_all_blast.paca)$chao1 <- estimate_richness(phyloseq_all_blast.paca, measures=c('Chao1'))$Chao1
#sample_data(phyloseq_all_blast.paca)$invsimpson <- estimate_richness(phyloseq_all_blast.paca, measures=c('InvSimpson'))$InvSimpson
#sample_data(phyloseq_all_blast.paca)$reads_sample <- sample_sums(phyloseq_all_blast.paca)

#smd.paca <- as(sample_data(phyloseq_all_blast.paca), 'data.frame')
#smd.paca$AMF_family = "Pacisporaceae"

######### SACCULOSPORACEAE ######## 
phyloseq_all_blast.sac <- subset_taxa(phyloseq_all_blast, AMF_family=="Sacculosporaceae")

sample_data(phyloseq_all_blast.sac)$shannon <- estimate_richness(phyloseq_all_blast.sac, measures=c('Shannon'))$Shannon
sample_data(phyloseq_all_blast.sac)$observed <- estimate_richness(phyloseq_all_blast.sac, measures=c('Observed'))$Observed
sample_data(phyloseq_all_blast.sac)$chao1 <- estimate_richness(phyloseq_all_blast.sac, measures=c('Chao1'))$Chao1
sample_data(phyloseq_all_blast.sac)$invsimpson <- estimate_richness(phyloseq_all_blast.sac, measures=c('InvSimpson'))$InvSimpson
sample_data(phyloseq_all_blast.sac)$reads_sample <- sample_sums(phyloseq_all_blast.sac)

smd.sac <- as(sample_data(phyloseq_all_blast.sac), 'data.frame')
smd.sac$AMF_family = "Sacculosporaceae"

###### JOIN ######
smd_all_blast_fam <- rbind(smd.aca, smd.amb, smd.arc, smd.cla, smd.div, smd.gig, smd.glo, smd.par, smd.per, smd.sac, smd.unknown)
# save 
saveRDS(smd_all_blast_fam, "data/smd_all_blast_fam_AK.BR.KS2019.RDS")