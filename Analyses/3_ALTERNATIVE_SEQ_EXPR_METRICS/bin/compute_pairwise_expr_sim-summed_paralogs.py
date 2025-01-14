#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import math
from scipy.stats import spearmanr
import itertools
import numpy as np

parser = argparse.ArgumentParser(description="Script to compute expression similarities between all pairs of genes from different species in the same gene orthogroup")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--expr_table_dir", "-e", required=True, metavar="expr_table_dir", help="Path to directory containing the tissue median/average expression tables by species")
parser.add_argument("--expr_file_suffix", "-ef", required=True, metavar="expr_file_suffix", help="Suffix for expression files to be uploaded")
parser.add_argument("--clade_species", "-c", required=True, metavar="clade_species", help="Comma separated list of all the species represented in the orthogroups. They will be used to read in the expression tables")
parser.add_argument("--tissues", "-t", required=True, metavar="tissues", help="Comma separated list of all tissues")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
expr_table_dir = args.expr_table_dir
expr_file_suffix = args.expr_file_suffix
clade_species_string = args.clade_species
tissue_string = args.tissues
output_file = args.output

##################################
###### READ INPUTS ###############
##################################
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
## Transform BmA to Bmo. I have been propagating this for too long...
orthogroups_df["Species"] = [species if species != "BmA" else "Bmo" for species in list(orthogroups_df["Species"])]

### Reduce the structure to OG_ID, Species
orthogroups_df = orthogroups_df[["OG_ID", "Species"]].drop_duplicates()

clade_species_list = clade_species_string.split(",")
all_tissues = tissue_string.split(",")

#Read the expression tables in. Join them in a unique table, adding columns of NAs when one tissue is missing.
joint_expr_table = pd.DataFrame()
for species in clade_species_list:
  species_expr_table = pd.read_table(expr_table_dir+"/" + species + expr_file_suffix, sep="\t", header=0, index_col=0)
  #rename index by adding the species
  species_expr_table = species_expr_table.rename(index = lambda x: x+"-"+species)
  #reorder entries
  species_expr_table = species_expr_table[all_tissues]
  #join to final table
  joint_expr_table = pd.concat([joint_expr_table, species_expr_table])


##################################
###### DEFINE FUNCTION ###########
##################################

##### Define logit function
def logit(p):
  logit_res = np.log(p/(1-p))
  return(logit_res)

##### Define function to compute pairwise expr_sim
def compute_pairwise_expr_sim(geneID1_expr, geneID2_expr):
  expr_sim = logit(1-sum([abs(a - b) for a,b in zip(geneID1_expr, geneID2_expr)])/len(geneID1_expr))
  return(expr_sim)

##################################
###### MAIN ######################
##################################
final_df = pd.DataFrame()

#group dataframe and cycle on groups
grouped_df = orthogroups_df.groupby("OG_ID")
for OG_ID, group in grouped_df:
  #Generate all possible pair of species between species belonging to the same orthogroup
  all_species = list(group["Species"])
  all_species_pairs = list(itertools.combinations(all_species, 2))
  #Filter out pairs with same species
  species_pairs_list = [[species[0], species[1]] for species in all_species_pairs]

  for pair in species_pairs_list:
    species1 = pair[0]
    species2 = pair[1]
    #check if the gene expression has been quantified (this part should be superflous)
    if OG_ID+"-"+species1 not in list(joint_expr_table.index.values) or OG_ID+"-"+species2 not in list(joint_expr_table.index.values):
      expr_sim = np.nan
    else:
      species1_expr = list(joint_expr_table.loc[OG_ID+"-"+species1,:])
      species2_expr = list(joint_expr_table.loc[OG_ID+"-"+species2,:]) 
      #Compute similarities
      expr_sim = compute_pairwise_expr_sim(species1_expr, species2_expr)
      #Join to final_df
    final_df = pd.concat([final_df, pd.DataFrame({"OG_ID" : [OG_ID, OG_ID], "Species1" : [species1, species2], "Species2" : [species2, species1], "GeneID1" : [OG_ID+"-"+species1, OG_ID+"-"+species2], "GeneID2" : [OG_ID+"-"+species2, OG_ID+"-"+species1], "Expr_sim" : [expr_sim, expr_sim]})])

#Save to output file
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
