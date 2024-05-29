##### Upload libraries
library(argparser)
library(tidyverse)

##### Define arguments
#Create a parser
p = arg_parser("Transfer GOs starting from representative ortholog of the query species in each orthogroup")

#Add command line arguments
p = add_argument(p, arg="--GO", help="GO annotation of query species. Col1=Gene_stable_ID, Col2=GO_domain, Col3=GO_term_accession, Col4=GO_term_name, Col5=GO_term_evidence_code")
p = add_argument(p, arg="--orthogroups", help="Input orthogroups (already filtered by expression and with one representative gene per species). Col1=orthogroup, Col2=Species, Col3=Gene")
p = add_argument(p, arg="--query_species", help="Identifier of the query species as in Col2 of the orthogroups file")
p = add_argument(p, arg="--output", help="Output file with transfered GO annotations")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
GO_file = argv$GO
orthogroups_file = argv$orthogroups
query_species = argv$query_species
output_file = argv$output

#########################################
########## MAIN #########################
#########################################

##### Format GOs
complete_GO_df = read.delim(GO_file, header=TRUE, col.names = c("GeneID", "Domain", "Accession", "Name", "Code")) %>%
  select(GeneID, Accession, Name) %>%
  mutate(Name = gsub(" ", "_", Name)) %>%
  filter(Accession != "")

##### Transfer GO
transferd_GO_annot_df = read.delim(orthogroups_file, col.names=c("OG_ID", "Species", "GeneID")) %>%
  filter(Species==query_species) %>%
  select(-Species) %>%
  left_join(complete_GO_df, by="GeneID") %>%
  select(-GeneID) %>%
  distinct()

##### Save to output
write.table(transferd_GO_annot_df, output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
