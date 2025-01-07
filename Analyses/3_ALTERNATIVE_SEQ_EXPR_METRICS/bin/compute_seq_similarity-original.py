#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import re
import os
import glob
import itertools
from Bio import SeqIO
import subprocess
import time
import numpy as np

parser = argparse.ArgumentParser(description="Script to recluster the gene orthogroups based on the orthopairs connections between a subset of species")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--score_matrix", "-a", required=True, metavar="score_matrix", help="Score matrix in the long format: col1=AA1, col2=AA2, col3=Score")
parser.add_argument("--OG_ID", "-id", required=True, metavar="OG_ID", help="OrthogroupID of the orthogroup of interest")
parser.add_argument("--alignments", "-aln", required=True, metavar="input_alignments", help="Alignments of all pairwise combinations of genes within a given orthogroup")
parser.add_argument("--scoring_system", "-ss", required=True, metavar="scoring_system", help="Scoring system to use when computing sequence similarity: it could include NoEndGaps")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output directory")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
score_file = args.score_matrix
my_orthogroup = args.OG_ID
aln_file = args.alignments
scoring_system = args.scoring_system
output_dir = str(args.output)

##################################
###### DEFINE FUNCTIONS ##########
##################################

######### This works with the original dataset
def compute_sim_score(fastas_entries_1, fastas_entries_2, AA_pair_score_dict, scoring_system):
  #### remove ending gaps from both sequences
  if scoring_system == "NoEndGaps":
    first_gene = fastas_entries_1.seq.rstrip("-")
    second_gene = fastas_entries_2.seq.rstrip("-")
    #### get the actual aligned sequence
    aligned_sequence_len = min(len(first_gene), len(second_gene))
    #### select only the sequences in both genes that are actually aligned
    first_gene = first_gene[0:aligned_sequence_len]
    second_gene = second_gene[0:aligned_sequence_len]
  elif scoring_system == "NoStartEndGaps":
    first_gene = fastas_entries_1.seq.rstrip("-").lstrip("-")
    second_gene = fastas_entries_2.seq.rstrip("-").lstrip("-")
    #### get indexes of the aligned positions
    start_index = max((len(fastas_entries_1.seq) - len(fastas_entries_1.seq.lstrip("-"))), (len(fastas_entries_2.seq) - len(fastas_entries_2.seq.lstrip("-"))))
    end_index = start_index + min(len(fastas_entries_1.seq.rstrip("-").lstrip("-")), len(fastas_entries_2.seq.rstrip("-").lstrip("-")))
    #### select only the sequences in both genes that are actually aligned
    first_gene = fastas_entries_1.seq[start_index:end_index]
    second_gene = fastas_entries_2.seq[start_index:end_index]
  else:
    first_gene = fastas_entries_1.seq
    second_gene = fastas_entries_2.seq
  #### compute score as original
  final_score_first = 0
  final_score_second = 0
  ## This is now a bit redundat, but still
  for position in list(range(0, len(first_gene))):
    AA_pair = str(fastas_entries_1.seq[position]) + str(fastas_entries_2.seq[position])
    score = AA_pair_score_dict[AA_pair]
    final_score_first = final_score_first + score
  for position in list(range(0, len(second_gene))):
    AA_pair = str(fastas_entries_1.seq[position]) + str(fastas_entries_2.seq[position])
    score = AA_pair_score_dict[AA_pair]
    final_score_second = final_score_second + score
  first_gene_score = final_score_first/len(re.sub("-", "", str(first_gene)))
  second_gene_score = final_score_second/len(re.sub("-", "", str(second_gene)))
  return([first_gene_score, second_gene_score])


##################################
###### READ INPUTS ###############
##################################

#Orthogroups
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Generate geneID-species dict
geneID_species_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict()

#Read dataframe with scores for the AA pairs
score_df = pd.read_table(str(score_file), sep="\t", index_col=False, header=None, names=["First_AA", "Second_AA", "Score"])
## replace * with -
score_df["First_AA"] = [element if element != "*" else "-" for element in list(score_df["First_AA"])]
score_df["Second_AA"] = [element if element != "*" else "-" for element in list(score_df["Second_AA"])]
#Generate dictionary with key=AA_pair and value=Score
score_df["AA_pair"] = score_df["First_AA"]+score_df["Second_AA"]
AA_pair_score_dict = pd.Series(score_df.Score.values, index=score_df.AA_pair).to_dict()


##################################
###### MAIN ######################
##################################
#Select the orthogroup of interest
filtered_orthogroup_df = orthogroups_df.loc[orthogroups_df["OG_ID"]==my_orthogroup]

#Get alignments within the orthogroup
all_aln_fastas = list(list(SeqIO.parse(aln_file, "fasta")))

#Generate all possible pair of genes between genes belonging to the orthogroup
all_genes = list(filtered_orthogroup_df["GeneID"])
all_gene_pairs = list(itertools.combinations(all_genes, 2))
#Filter out pairs where genes come from the same species
filt_gene_pairs = [gene for gene in all_gene_pairs if geneID_species_dict[gene[0]] != geneID_species_dict[gene[1]]]
#Transform tuple to a list
filt_gene_pairs = [[gene[0], gene[1]] for gene in filt_gene_pairs]

#Compute sequence similarity based on blosum (alignments are already saved)
for gene_pair in filt_gene_pairs:
  geneID1 = gene_pair[0]
  geneID2 = gene_pair[1]
  species1 = geneID_species_dict[geneID1]
  species2 = geneID_species_dict[geneID2]

  ### there are 4 exceptions, but these are not in bilaterian conserved orthogroups, so I am just going to skip them
  #### almost all of these problems are within a single orthogroup
  if gene_pair==["GB55483", "LOC590007"] or gene_pair==['BLAG17000895', 'ENSMUSG00000051747'] or geneID1=="GB47977" or geneID2=="GB47977" or geneID1=="ENSG00000155657" or geneID2=="ENSG00000155657" or gene_pair==['CD03878', 'EBAG0002097'] or geneID1 == "ENSMUSG00000051747" or geneID2=="ENSMUSG00000051747" or gene_pair == ['CD03878', 'TC030701'] or gene_pair == ['EBAG0002097', 'TC030701']:
    continue #skip to next iteration
  else:
    #Isolate pre-computed alignments. The header is in the format QueryGene_TargetGene, where the sequence is the sequence of QueryGene in the alignment with TargetGene
    fastas_entries_1 = [element for element in all_aln_fastas if element.id == geneID1 + "_" + geneID2][0]
    fastas_entries_2 = [element for element in all_aln_fastas if element.id == geneID2 + "_" + geneID1][0]
  
    #Compute sim score
    sim_score = compute_sim_score(fastas_entries_1, fastas_entries_2, AA_pair_score_dict, scoring_system) #This is a list containing the sim_score of geneID1 and the sim_score of geneID2
  
    #Print the similarity score in a reciprocal way.
    output_sim_score_file = "%s/%s-%s-%s-sim_scores" % (output_dir, geneID1, geneID2, my_orthogroup)
    with open(output_sim_score_file, "w") as output_sim_score:
      output_sim_score.write("%s\t%s\t%s\t%s\t%s\t%f\n" % (my_orthogroup, species1, species2, geneID1, geneID2, sim_score[0]))
      output_sim_score.write("%s\t%s\t%s\t%s\t%s\t%f\n" % (my_orthogroup, species2, species1, geneID2, geneID1, sim_score[1]))
    output_sim_score.close()


##################################
######## CLEANING UP #############
##################################

##### Join scores for all pairs in the same file
cat_score_command = "cat %s/*-%s-sim_scores > %s/%s-all_sim_scores" % (output_dir, my_orthogroup, output_dir, my_orthogroup)
os.system(cat_score_command)
##### Remove single input files
single_sim_scores_rm_command = "rm %s/*-%s-sim_scores" % (output_dir, my_orthogroup)
os.system(single_sim_scores_rm_command)
