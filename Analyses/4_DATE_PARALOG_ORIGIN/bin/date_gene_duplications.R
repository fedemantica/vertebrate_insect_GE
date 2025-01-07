###################################################################################
##### Script to associate each event to the relative geneID and OGID ##############
###################################################################################

library(argparser)
library(tidyverse)

#########################################
########## DEFINE ARGUMENTS #############
#########################################

#Create a parser
p = arg_parser("Parse file with gene duplication info to return the most recent node where each gene was duplicated, with the relative bootstrap")

#Add command line arguments
p = add_argument(p, arg="--orthogroups_file", help="File containing orthogroups as in Broccoli output. Col1=OG_ID, Col2=Species, Col3=GeneID")
p = add_argument(p, arg="--tree_levels_file", help="File with Col1=node/leaf name, Col2=tree level. This file is useful to select the most recent duplication from which a given gene derives")
p = add_argument(p, arg="--duplications_info_file", help="Output of parse_gene_trees.py")
p = add_argument(p, arg="--lca_type", help="Type of approach to call the gene duplication. It can be lca_union or lca_intersect")
p = add_argument(p, arg="--output", help="Path to output file")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
orthogroups_file = argv$orthogroups_file
tree_levels_file = argv$tree_levels_file
duplications_info_file = argv$duplications_info_file
lca_type = argv$lca_type
output_file = argv$output

#########################################
########## MAIN #########################
#########################################

#####################
## READ FILES #######
#####################

#### Read in orthogroups file and get GeneID - Species correspondences
orthogroups_df = read.delim(orthogroups_file, header=FALSE, col.names=c("orthogroup_ID", "Species", "GeneID"))
species_geneID_df = orthogroups_df %>% dplyr::select(-orthogroup_ID)

#### Read in tree levels file
tree_levels_df = read.delim(tree_levels_file, header=TRUE)

#### Read in duplication input file
all_dated_duplications_df = read.delim(duplications_info_file, header=TRUE, sep="\t")

#####################
## PROCESS ##########
#####################

expanded_genes_df = all_dated_duplications_df %>%
  dplyr::select(orthogroup_ID, lca_type, bootstrap_support, genes_left) %>%
  separate_rows(genes_left, sep=",") %>%
  rename(GeneID = genes_left, lca = lca_type) %>%
  dplyr::bind_rows(all_dated_duplications_df %>% 
                     dplyr::select(orthogroup_ID, lca_type, bootstrap_support, genes_right) %>%
                     separate_rows(genes_right, sep=",") %>%
                     rename(GeneID = genes_right, lca = lca_type)) %>%
  ### Fix mistake (this should be done upstream in the snakemake)
  #But not harm is done if it runs after it has been fixed
  dplyr::mutate(lca = ifelse(lca == "Deuterostomia", "Deuterostoma", lca)) %>%
  dplyr::mutate(lca = ifelse(lca == "Protostomia", "Protostoma", lca)) %>%
  #### put NA if there is low support for the duplication branch. I am setting the minimum at 0.5
  dplyr::mutate(lca = ifelse(bootstrap_support <= 0.5 | is.na(bootstrap_support), "low_support", lca))
  

#### For each gene, pick the most recent ancestor where the gene appears within a duplication
selected_dup_expanded_genes_df = expanded_genes_df %>%
  left_join(tree_levels_df, join_by("lca"=="ancestor")) %>%
  group_by(GeneID) %>%
  #### If there are only NA, pick a random one. If there are other tree_level values, pick the lowest one (most recent duplication)
  mutate(all_na_tree_level = all(is.na(tree_level))) %>%
  filter((all_na_tree_level & row_number() == 1) | (!all_na_tree_level & tree_level == min(tree_level, na.rm = TRUE))) %>%
  ungroup() %>%
  distinct() %>%
  dplyr::select(-all_na_tree_level) %>%
  group_by(GeneID) %>%
  #### Maybe this part is not necessary anymore, I have to check
  slice_max(bootstrap_support, with_ties = TRUE) %>%
  ####### Add species from GeneID
  left_join(species_geneID_df, by="GeneID") %>%
  ungroup()

############################
## SAVE TO OUTPUT ##########
############################
write.table(selected_dup_expanded_genes_df, output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
