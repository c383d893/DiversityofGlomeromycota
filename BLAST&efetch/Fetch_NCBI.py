
import subprocess
from Bio import Entrez


# Set your email address (required by NCBI)
Entrez.email = "alexis.aellen@bluewin.ch"

# Fetch the data for each ID
with open("ncbi_species_output.txt", "w") as output_file:
    # Write the header row
    #output_file.write("NCBI ID\tGenus\t...\n")

    with open("ID_list.txt", "r") as id_file:
        # Read the contents of the ID file and split by commas
        id_list = id_file.read().split(",")
        #  WARNING !!!!!!!  : at the end of the ID_list_copie.txt file there is sometimes '\n' : to check this execute : print(id_list)
        # make sure the last ID is followed by a ",". The for loop with the continue is used to ignore this '\n'. 

        for ncbi_id in id_list:
            if ncbi_id == "\n" : continue
            command = f"efetch -db nuccore -id {ncbi_id} -format fasta | awk '/^>/ {{sub(/.*>/,\"\",$0); sub(/,.*/,\"\",$0); print}}'"

            # Execute the command in the terminal and capture the output
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            output, _ = process.communicate()

            # Decode the output and print the fetched data
            fasta_data = output.decode("utf-8").strip()
            data_list = fasta_data.split()

            #print(fasta_data)
            # Write the data to the output file
            for i in range(len(data_list)):
                output_file.write(f"{data_list[i]}\t")
            output_file.write(f"\n")
