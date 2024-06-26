#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to retrieve coordinates of CDS exons for each gene")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Reference GTF (i.e. just one representative isoform for each gene)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

#######################
##### Main ############
#######################

#Read gtf in and add geneID column
my_gtf_df = pd.read_table(my_input, sep="\t", names=["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"])
intermediate_gene_list = [re.sub(";.*", "", element) for element in list(my_gtf_df["attributes"])]
my_gtf_df["geneID"] = pd.Series([re.sub(" ", "", re.sub("gene_id ", "", re.sub('"', '', element))) for element in intermediate_gene_list])

#### Filter only for the exons entry and subset columns
my_filtered_gtf = my_gtf_df.loc[my_gtf_df.feature=="CDS"]
final_df = my_filtered_gtf.loc[:,["geneID","chr","start","stop"]]

#### Save to output file
final_df.to_csv(my_output, sep="\t", index=False, header=False)
