
blast_tophit <- read.table("data/blast_allglom_tophit_BLAST_clean07.06.23.txt", header = FALSE, sep = "\t")

#take out all hits under 98%
blast_tophit <- blast_tophit[blast_tophit$V3 >= 98,]
list_ncbi_ID <- as.character(blast_tophit$V2)

# Combine the ID's into a single string separated by a comma
ID_string <- paste(unique(list_ncbi_ID), collapse = ",")


# Write the string to a file
write(ID_string, "ID_list.txt")

Alexis_orig_list <- read.table("data/ID_list.txt", sep=",")
