#########################################
##### AVERAGE PHYLOP SCORES #############
#########################################

##### Upload libraries
library(rtracklayer, lib.loc="/users/mirimia/fmantica/R/x86_64-pc-linux-gnu-library/4.2")
library(data.table)
#library(fuzzyjoin)
#library(GenomicRanges)

##### Define arguments
#Create a parser
p = arg_parser("Compute average phyloP score based on genomic coordinates of inputed genes")

#Add command line arguments
p = add_argument(p, arg="gene_coordinates", help="Bed file with coordinates of the genes (fourth column) for which to compute the average score")
p = add_argument(p, arg="wigFix", help="file containing the PhyloP score for each base of a chr (starting from a certain point)")
p = add_argument(p, arg="output", help="Output file")

##### Read in files
wig_file = "/no_backup/mirimia/fmantica/vertebrate_insects_GE/phyloP/chr1.phastCons100way.wigFix.top5000.gz"
gene_coordinates_file = "/users/mirimia/fmantica/projects/vertebrate_insect_GE/data/DB/gene_extra_info/Hs2_start_stop.tab"

#########################
###### MAIN #############
#########################

###### Read in files
wig_data = import(wig_file, format = "wig")
wig_df = as.data.frame(wig_data)
gene_coordinates_df = read.delim(gene_coordinates_file, header=FALSE, col.names=c("GeneID", "Chr", "Start", "Stop"))

##### Assign GeneID to all its coordinate scores
chr_gene_coordinates_df = gene_coordinates_df %>% filter(Chr == "chr1")

test_df = foverlaps(wig_df, chr_gene_coordinates_df, by.x=c("seqnames", "start", "end"), by.y=c("Chr", "Start", "Stop"))

all_genes_average_scores_df = data.frame()
#### Cycle on all the geneIDs (it is maybe faster this way)
for (geneID in as.vector(chr_gene_coordinates_df$GeneID)) {
  gene_start = chr_gene_coordinates_df %>% filter(GeneID==geneID) %>% .$Start
  gene_stop = chr_gene_coordinates_df %>% filter(GeneID==geneID) %>% .$Stop
  gene_scores_df = wig_df %>% filter(start >= gene_start, end <= gene_stop)
  ### Join to final dataframe
  all_genes_average_scores_df = dplyr::bind_rows(all_genes_average_scores_df, gene_scores_df)
}

test_df = fuzzy_left_join(chr_gene_coordinates_df, wig_df, 
                          by=c(Chr="seqnames", Start="start", Stop="end"), match_fun=list(`==`, `>=`, `<=`))
  
  
# Define your genomic ranges here
# granges <- GRanges(seqnames = "chr1",
#                    ranges = IRanges(start = c(11500, 13400),
#                                     end = c(12150, 13600)))



# Calculate average coverage in each specified range
average_values = sapply(granges, function(gr) {
  mean(coverage_data[[as.character(seqnames(gr))]][start(gr):end(gr)])
})