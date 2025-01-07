#!/usr/bin/env python3

import sys
#NB: this has to be run in my Broccoli env to work. Otherwise, decomment the following line (which is ugly)
#sys.path.extend(['/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python36.zip', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/lib-dynload', '/users/mirimia/fmantica/software/anaconda3/envs/env-broccoli/lib/python3.6/site-packages'])

import argparse
from Bio import SeqIO
import pandas as pd
import re
import os
import glob
import subprocess

parser = argparse.ArgumentParser(description="Generate multiple alignments among all the proteins in an orthogroup containing chimeric genes")
parser.add_argument("--input_fastas", "-f", required=True, metavar="input_fastas", help="directory containing all proteome fasta files of species in orthogroups")
parser.add_argument("--input_OG", "-i", required=True, metavar="input_orthogroups", help="File containing all orthogroups to subset")
parser.add_argument("--selected_orthogroups", "-g", required=True, metavar="selected_orthogroups", help="comma-separated list of orthogroup IDs for which to perform the alignment")
parser.add_argument("--mafft", "-m", required=True, metavar="mafft_executable", help="path to mafft execution file")
parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output dir")


###### Read arguments
args = parser.parse_args()
my_input_fastas = args.input_fastas
my_input_OG = args.input_OG
my_selected_orthogroups = args.selected_orthogroups.split(",") #split the original batch
my_mafft = args.mafft
my_output = args.output

print(my_selected_orthogroups)
#No header
orthogroups_df = pd.read_table(str(my_input_OG), sep="\t", header=None, names=["OG_ID", "Species", "GeneID"], index_col=False)

#Read and join fasta files all together
all_fastas_files = glob.glob(my_input_fastas+"/*") #list all the fasta files in the input directory
all_fastas_entries = []
for my_file in all_fastas_files: #cycle on all the fasta files
  print(my_file) #This is just for the output
  my_fasta = list(list(SeqIO.parse(my_file, "fasta")))
  for element in my_fasta: ### Match the geneID to Broccoli output 
    element.id = re.sub(".*\|", "", element.id)
    element.description = ""
  all_fastas_entries = all_fastas_entries+my_fasta

##################################################
################### MAIN #########################
##################################################

for my_OG_ID in  my_selected_orthogroups:
  #### Select all genes in orhotgroups and relative fasta entries
  OG_genes = list(orthogroups_df[orthogroups_df["OG_ID"]==my_OG_ID]["GeneID"])
  fastas_entries = [element for element in all_fastas_entries if element.id in OG_genes]

  #### Save fastas entries to temporary file
  input_temp = my_OG_ID+"-input_fasta_tmp.fa"
  with open(input_temp, "w") as input_to_aln:
    SeqIO.write(fastas_entries, input_to_aln, "fasta")
  input_to_aln.close()

  #### Print shell command dunning the mafft alignment to logs (just to show what has been run)
  my_command = "%s --quiet --retree 2 --localpair --maxiterate 500 input_fasta_tmp.fa > %s/%s-multiple_aln" % (my_mafft, my_output, my_OG_ID)
  print(my_command) #this is just to know what has been run.

  #### Define output file
  output_file =  "%s/%s-multiple_aln" % (my_output, my_OG_ID)

  #### Run alignments
  with open(output_file, "w") as output:
    p = subprocess.Popen([my_mafft, "--quiet", "--retree", "2", "--localpair", "--maxiterate", "500", input_temp], stdout=output)
    p.wait()

  #### Remove temporary input
  remove_command = "rm %s" % (input_temp) 
  os.system(remove_command)    
