#!/usr/bin/env python3

import argparse
import pandas as pd
from ete3 import PhyloTree
from ete3 import Tree
import os
from pathlib import Path

parser = argparse.ArgumentParser(description="Script to data gene duplications of a know gene tree by reconciling it with the relative species tree.")
parser.add_argument("--orthogroups_file", "-og", required=True, metavar="orthogroups_file", help="Broccoli output with: col1=OG_DI, col2=Species, col3=GeneID.")
parser.add_argument("--selected_orthogroups", "-so", required=True, metavar="selected_orthogroups", help="Comma-separated list of orthogroups for which to date gene duplications. It has to be provided as a single string.")
parser.add_argument("--input_dir", "-id", required=True, metavar="input_dir", help="Path to directory containing tree.")
parser.add_argument("--input_suffix", "-s", required=True, metavar="input_suffix", help="Suffix for the input gene trees (e.g., OG_ID{input_suffix}).")
parser.add_argument("--species_tree_file", "-is", required=True, metavar="species_tree_file", help="File containing the species tree in newick format.")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file.")

####################################
##### DEFINE ARGUMENTS #############
####################################
args = parser.parse_args()
orthogroups_file = args.orthogroups_file
selected_orthogroups = args.selected_orthogroups.split(",")
input_dir = args.input_dir
input_suffix = args.input_suffix
species_tree_file = args.species_tree_file
output = args.output


################################################################
################ FUNCTIONS #####################################
################################################################

def is_duplication(node):
  # Extract species from both child nodes
  left_species = {leaf.name.split('_')[0] for leaf in node.children[0].iter_leaves()}
  right_species = {leaf.name.split('_')[0] for leaf in node.children[1].iter_leaves()}
  # Check if there's any overlap in species
  return bool(left_species & right_species)


## lca stands for "most recent common ancestor"
def get_lca_in_species_tree(species_tree, species_set):
  return species_tree.get_common_ancestor(species_set)

def map_duplications_to_species_tree(gene_tree, species_tree, OG_ID):
  duplications = []
  for node in gene_tree.traverse():
    if not node.is_leaf() and is_duplication(node):
      # Get species sets from both child nodes of the duplication node
      left_species = {leaf.name.split('_')[0] for leaf in node.children[0].iter_leaves()}
      right_species = {leaf.name.split('_')[0] for leaf in node.children[1].iter_leaves()}
      # Union or Overlap of species from both clades under the duplication node
      union_species = left_species | right_species
      intersect_species = left_species & right_species
      # Find the MRCA of these species in the species tree
      if len(list(union_species)) > 1:
        lca_node_union = get_lca_in_species_tree(species_tree, union_species).name
        if lca_node_union == species_tree.get_tree_root().name:
          ##### Set the union to root only if there is reasonable ratio (i.e., less than 1:4) between the number of species on the two sides
          if (len(left_species) == 1 or len(right_species) == 1 or len(intersect_species) == 1) or (len(left_species)/len(intersect_species) >= 4 or len(left_species)/len(intersect_species) <= 0.25):
            lca_node_union = "Undetermined"
          else: 
            lca_node_union = "root"
      else:
        lca_node_union = list(union_species)[0]
      if len(list(intersect_species)) > 1:
        lca_node_intersect = get_lca_in_species_tree(species_tree, intersect_species).name
        if lca_node_intersect == species_tree.get_tree_root().name:
          lca_node_intersect = "root"
      else:
        lca_node_intersect = list(intersect_species)[0]
      # Record duplication mapping
      duplications.append({
          "orthogroup_ID": OG_ID,
          #### for some reason, the bootstrap gets read as the node name
          "bootstrap_support": node.name,
          "lca_union": lca_node_union,
          "lca_intersect": lca_node_intersect,
          "species_union": ",".join(list(left_species | right_species)),
          "species_intersect":  ",".join(list(left_species & right_species)),
          "species_left": ",".join(list(left_species)),
          "species_right": ",".join(list(right_species)),
          "genes_left": ",".join([leaf.name.split('_')[1] for leaf in node.children[0].iter_leaves()]),
          "genes_right": ",".join([leaf.name.split('_')[1] for leaf in node.children[1].iter_leaves()])
      })
  return pd.DataFrame(duplications)


###############################################################
############### MAIN ###########################################
################################################################
 
if __name__ == '__main__':

  ####################################
  ##### READ ARGUMENTS ###############
  ####################################
  
  ### Species tree
  species_tree = PhyloTree(species_tree_file, format=1)
  species_tree.resolve_polytomy(recursive=True) #Resolve polynomies, even if the inputs look totally fine to me
  
  ### GeneID-species dictionaries
  orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
  geneID_speciesID_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict() 
  
  ### Initialize dataframe
  all_duplications_df = pd.DataFrame()
  
  ####################################
  ##### DATE DUPLICATIONS ############
  ####################################
  
  for OG_ID in selected_orthogroups:
    ### Print OG_ID for log purposes
    print(OG_ID)
    ### Compose input name (this is not ideal)
    input_gene_tree = input_dir + OG_ID + "-" + input_suffix
    #################################################################
    ### For the moment, if file does not exist, move on. ############
    #################################################################
    if not os.path.exists(input_gene_tree) or Path(input_gene_tree).stat().st_size == 0:
      print("This file "+input_gene_tree+" is empty.")
      continue

    ### Read gene tree
    gene_tree = PhyloTree(input_gene_tree, format=1)
    ### Resolve polytomies and root
    gene_tree.resolve_polytomy(recursive=True) #Resolve polynomies, even if the inputs look totally fine to me 
    ### Root the tree (This part was copied from Broccoli step2 script)
    mid = gene_tree.get_midpoint_outgroup()
    try:
        gene_tree.set_outgroup(mid)
    except:
        pass
    #### Modify leaf names so that they include the speciesID
    for leaf in gene_tree.iter_leaves():
      leaf.name = geneID_speciesID_dict[leaf.name] + "_" + leaf.name
    #### Date duplications with a species overlap approach
    duplications = map_duplications_to_species_tree(gene_tree, species_tree, OG_ID)
    #### Join to final dataframe
    all_duplications_df = pd.concat([all_duplications_df, duplications])
  
  ####################################
  ##### SAVE TO FILE #################
  ####################################
  all_duplications_df.replace(r"^\s*$", "NA", regex=True, inplace=True) 
  all_duplications_df.to_csv(output, sep="\t", header=True, index=False, na_rep="NA")
