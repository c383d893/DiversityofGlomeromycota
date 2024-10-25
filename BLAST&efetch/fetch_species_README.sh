#!/bin/bash
#SBATCH -n 24                     # 24 cores
#SBATCH --time 08:00:00                   # 8-hour run-time
#SBATCH --mem-per-cpu=4000     # 4000 MB per core
#SBATCH -J analysis1
#SBATCH -o analysis1.out
#SBATCH -e analysis1.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alexis.aellen@bluewin.ch


python3 Fetch_NCBI.py

#READ-ME 

#GOAL : This folder is to extract species name from a list of NCBI_IDs in coma separated file. Place a coma after the last occurence in the file, 
#this will be useful later. 

# 1) To be able to retrieve NCBI data, you first have to install Entrez direct software  : to do this, run the following : 
#  sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
#answer "y" to the questions
#  export PATH=${HOME}/edirect:${PATH}

# 2) TEST : You should then be able to retrieve information from ncbi from your terminal. To test this, run the following : 
#efetch -db nuccore -id HG798895.1 -format fasta | awk '/^>/ {sub(/.*>/,"",$0); sub(/,.*/,"",$0); print}'
#If this doesn't work, a problem likely occured in the installation process

# 3) To run the following code Python3 is necessary and make sure you have installed the Biopython library : to do this run : 
#python3 -m pip install biopython
#if a problem persist, add "export PATH="/usr/local/bin:$PATH" at the end of your .bashrc file in your home directory. 

#throughout this code, if a problem persist and the answer is not here, feel free of asking chatgpt, he is a very good friend to have ! 

# 4) Don't forget to change your email adress in the Fetch_NCBI.py

