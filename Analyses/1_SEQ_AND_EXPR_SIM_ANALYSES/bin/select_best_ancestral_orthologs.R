##### Upload libraries
library(argparser)
library(tidyverse)

##### Define arguments
#Create a parser
p = arg_parser("Select Best-Ancestral paralogs within vertebrates and insects. Importantly, the outpupt orthology will not include genes for the outgroups")

#Add command line arguments
p = add_argument(p, arg="--gene_orthogroups", help="Input orthogroups (already filtered by expression). Col1=orthogroup, Col2=Species, Col3=Gene")
p = add_argument(p, arg="--seq_sim", help="All pairwise sequence similarities. Col1=orthogroup, Col2=Species1, Col3=Species2, Col4=GeneID1, Col5=GeneID2, Col6=Seq_Sim")
p = add_argument(p, arg="--vertebrata", help="Comma separated list of vertebrate species IDs")
p = add_argument(p, arg="--insecta", help="Comma separated list of insect species IDs")
p = add_argument(p, arg="--output", help="Output file with best-ancestral orthogroups")
p = add_argument(p, arg="--output_complete_BA", help="Output file with best-ancestral orthogroups and seq sim values")
p = add_argument(p, arg="--output_complete", help="Output file with complete orthogroups and seq sim values")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
all_orthogroups_file = argv$gene_orthogroups
all_pairwise_seq_sim_file = argv$seq_sim
Vertebrata = strsplit(argv$vertebrata, ",")[[1]]
Insecta = strsplit(argv$insecta, ",")[[1]]
output_file = argv$output
output_complete_BA_file = argv$output_complete_BA
output_complete_file = argv$output_complete

#########################################
########## MAIN #########################
#########################################

#######################
##### Read in files ###
#######################
all_orthogroups_df = read.delim(all_orthogroups_file, header=FALSE, col.names=c("OG_ID", "Species", "GeneID"))
all_pairwise_seq_sim_df = read.delim(all_pairwise_seq_sim_file, header=FALSE, col.names=c("OG_ID", "Species1", "Species2", "GeneID1", "GeneID2", "Seq_sim"))

#######################
### Average seq sim ###
#######################

#### Compute average by seq sim by gene considering both pairwise comparisons among species in the same clade
input_for_average_df = all_pairwise_seq_sim_df %>%
  filter(Species1 %in% c(Vertebrata, Insecta), Species2 %in% c(Vertebrata, Insecta)) %>%
  ### Add species with correct order
  dplyr::select("OG_ID", "GeneID1", "GeneID2", "Seq_sim") %>%
  left_join(all_orthogroups_df, by=c("OG_ID", "GeneID1"="GeneID")) %>%
  dplyr::rename(Species1 = Species) %>%
  left_join(all_orthogroups_df, by=c("OG_ID", "GeneID2"="GeneID")) %>%
  dplyr::rename(Species2 = Species) %>%
  ### Have everything in a long format
  dplyr::select(-GeneID2) %>%
  dplyr::rename(GeneID = GeneID1, Species=Species1, Target_species = Species2) %>%
  dplyr::bind_rows(all_pairwise_seq_sim_df %>%
                     filter(Species1 %in% c(Vertebrata, Insecta), Species2 %in% c(Vertebrata, Insecta)) %>%
                     ### Add species with correct order
                     dplyr::select("OG_ID", "GeneID1", "GeneID2", "Seq_sim") %>%
                     left_join(all_orthogroups_df, by=c("OG_ID", "GeneID1"="GeneID")) %>%
                     dplyr::rename(Species1 = Species) %>%
                     left_join(all_orthogroups_df, by=c("OG_ID", "GeneID2"="GeneID")) %>%
                     dplyr::rename(Species2 = Species) %>%
                     ### Have everything in a long format
                     dplyr::select(-GeneID1) %>%
                     dplyr::rename(GeneID = GeneID2, Species=Species2, Target_species = Species1))

#### Compute average seq_sim among all species within clade
average_seq_sim_df = input_for_average_df %>%
  mutate(Clade = ifelse(Species %in% Vertebrata, "Vertebrata", "Insecta"),
         Target_Clade = ifelse(Target_species %in% Vertebrata, "Vertebrata", "Insecta")) %>%
  filter(Clade==Target_Clade) %>%
  group_by(OG_ID, GeneID, Species, Clade) %>%
  summarize(Average_seq_sim = mean(Seq_sim))

#### Select the paralog with lower average seq sim per species and per clade
best_ancestral_complete_df = average_seq_sim_df %>%
  group_by(OG_ID, Species, Clade) %>%
  filter(Average_seq_sim == max(Average_seq_sim)) %>%
  ungroup() %>%
  #For the few cases with a tie, select one random paralog
  group_by(OG_ID, Species) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(OG_ID, Species, GeneID, Clade, Average_seq_sim)

best_ancestral_df = best_ancestral_complete_df %>%
  dplyr::select(OG_ID, Species, GeneID)

#######################
##### Save output #####
#######################
write.table(best_ancestral_df, output_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(best_ancestral_complete_df, output_complete_BA_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(average_seq_sim_df, output_complete_file, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
