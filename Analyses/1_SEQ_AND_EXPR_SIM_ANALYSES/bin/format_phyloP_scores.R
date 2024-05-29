#########################################
##### FORMAT PHYLOP SCORES #############
#########################################

##### Upload libraries
library(argparser)
library(rtracklayer)

####### Define arguments
#Create a parser
p = arg_parser("Format phyloP scores to tabular format (easier to parse)")

#Add command line arguments
p = add_argument(p, arg="--input", help="wigFix file containing the PhyloP score for each base of a chr (starting from a certain point)")
p = add_argument(p, arg="--output", help="Output file: tabular file containing the same information")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
input_file = argv$input
output_file = argv$output

#########################
###### MAIN #############
#########################

###### Read in file and transform to da
wig_data = import(input_file, format = "wig")
wig_df = as.data.frame(wig_data)

###### Save to output
gz_file <- gzfile("output_file", "w")
write.table(wig_df, gz_file, sep="\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
close(gz_file)