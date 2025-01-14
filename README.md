# Evolutionary diversification of ancestral genes across vertebrates and insects

This github repository contains the code use to perform the analyses and generate the figures associated with the manuscript titled **"Evolutionary diversification of ancestral genes across vertebrates and insects"** by Federica Mantica and Manuel Irimia (submitted for publication).

Table of contents
-------
* [Study Overview](#study-overview)  
* [Analyses](#analyses)  
  + [All pairwise similarities](#all-pairwise-similarities)
  + [Seq and expr sim analyses](#seq-and-expr-sim-analyses)
  + [Paralogous sets selection](#paralogous-sets-selection)
  + [Alternative expr metrics](#alternative-seq-expr-metrics)
  + [Date paralog origin](#data-paralog-origin)
* [Figures](#figures)  
  + [Main figures](#main-figures)   
  + [Alternative setups figures](#alternative-setups-figures)  

### Study Overview 

Vertebrates and insects diverged approximately 700 million years ago, and yet they retain a large core of conserved genes from their last common ancestor. These ancient genes likely present strong evolutionary constraints, which limit their overall sequence and expression divergence. However, evolutionary constraints can greatly vary across ancestral gene families. In few cases, sequence and expression changes resulting from looser constraints mediate functional evolution, as it has been demonstrated for many paralogous genes in closely related species. However, overall patterns of sequence and expression divergence and their potential functional outcome have never been explored in a genome-wide manner across large evolutionary distances in animals. In this study, we focus on approximately 7,000 highly conserved genes shared between vertebrates and insects, investigating global patterns of molecular diversification driven by changes in sequence and expression profiles. We identified molecular features generally linked to higher or lower diversification rates, together with gene groups with similar diversification profiles in both clades. However, we also discovered that some genes underwent differential diversification during vertebrate and insect evolution, potentially contributing to the emergence of unique phenotypes in each clade.

<br>

<figure>
  <img align="center" src="./docs/Fig1.png" width=700 height=400 />
  <figcaption><b>Framework overview and definition of sequence and expression similarities.</b> a. Left: phylogenetic tree including the scientific acronyms of the 8 vertebrate and 8 insect species considered in this study. Hsa: human, Mmu: mouse, Bta: cow, Mdo: opossum, Gga: chicken, Xtr: tropical clawed frog, Dre: zebrafish, Cmi: elephant shark. Dme: fruit fly, Eba: marmalade hoverfly, Aae: yellow fever mosquito, Bmo: domestic silk moth, Tca: red flour beetle, Ame: honey bee, Bge: cockroach, Cdi: mayfly. (see Methods for corresponding scientific names). Top right: example of one the 6,787 considered gene orthogroups. Bottom right: tissues represented in our bulk RNA-seq dataset. Evolutionary distances were derived from timetree (Kumar et al. 2022) (MYA: million years ago) and animal silhouettes were generated through Bing Chat by MIcrosoft (2023) https://www.bing.com/search. b, c. Scheme for the computation of the sequence (b) and expression (c) similarity measures. The depicted procedure was performed separately for vertebrates and insects, returning two values per orthogroup. d, e. Distributions of the sequence (d) and expression (e) similarity values for all gene orthogroups (n=6,788) within vertebrates (purple) and insects (orange).</figcaption>
<figure>

<br>

### Analyses  

This folder contains 4 subfolders, each containing a different snakemake pipeline and relative files.   

#### All pairwise similarities 
This folder contains snakemake pipeline (i.e., **Snakefile**) used to run the majority of computational analyses, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the bin subfolder. The Snakefile includes rules to:  

* Compute the **relative expression** across tissues for all the bilaterian-conserved, protein-coding genes of the 16 vertebrates and insect species included in our dataset.  
* Compute **all pairwise sequence similarities** between all possible pairs of genes (from different species) in an orthogroup.  
* Compute **all pairwise expression similarities** between all possible pairs of genes (from different species) in an orthogroup, starting from their relative expression across tissues.  

#### Seq and expr sim analyses
This folder contains snakemake pipeline (i.e., **Snakefile**) used to run the majority of computational analyses, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the bin subfolder. The Snakefile includes rules to:  

* Compute average protein **sequence similarities** for each gene orthogroups, starting from all the all species-pairwise similarities between its representative orthologs.
* Compute pairwise and average **expression similarities** for each gene orthogroups, starting from the relative expression across tissues of all its representative orthologs.  
* Define groups containing **500 most and least conserved gene orthogroups** in terms of sequence or expression similarities between either vertebrates or insects.  
* Generate inputs for **GO enrichments** of common lowly or highly diversified genes (see Fig.4 from Figures/All_Figures.Rmd).  
* NB: These rules are also use to generate the corresponding files upon all the alternative setup scenarios (see Analyses/3_ALTERNATIVE_SEQ_EXPR_METRICS and Figure/Alternative_setups_Figures.Rmd).  

#### Paralogous sets selection
This folder contains snakemake pipeline (i.e., **Snakefile**) used to run the majority of computational analyses, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the bin subfolder. The Snakefile includes rules to:  

* Compute **z-scores** of the average sequence or expression similarity per orthogroup and clade computed in Analyses/2_PARALOGOUS_SETS_SELECTION.  
* Select all extra paralogous sets (**BA-seq, BD-seq, BA-expr, BD-exr**) based on the most extreme similarity values across orthogroups: the highest for BA sets (Best-Ancestral) and the lowest for BD sets (Best-Divergent).  
* Define vertebrate housekeeping genes.  

#### Alternative seq and expr metrics  

This folder contains snakemake pipeline a (i.e., **Snakefile**) used to run the majority of computational analyses, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the bin subfolder. The Snakefile includes rules to:  

* Compute **all pairwise expression correlations** (spearman's rho) between all possible pairs of genes (from different species) in an orthogroup, starting from their relative expression across tissues.  
* Compute **all pairwise expression distances** (euclidean distance) between all possible pairs of genes (from different species) in an orthogroup, starting from their relative expression across tissues.  
* Compute **all pairwise sequence similarities** between all possible pairs of genes (from different species) in an orthogroup, directly using values from *blosum62*.    
* Compute **all pairwise sequence similarities** between all possible pairs of genes (from different species) in an orthogroup, directly using values from *blosum45*.  
* Compute **all pairwise sequence similarities** between all possible pairs of genes (from different species) in an orthogroup, as in the main approach (see Analyses/0_ALL_PAIRWISE_SIMILARITIES) but removing leading/trailing gaps from the alignments.  
* Compute **all pairwise expression similarities** between all possible pairs of genes (from different species) in an orthogroup, starting from their relative expression across tissues computed from **summed paralogs expression**.  
* Compute the average **phastCons** across the coding sequence of all human genes.  

#### Date paralog origin  
This folder contains snakemake pipeline a (i.e., **Snakefile**) used to run the majority of computational analyses, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the bin subfolder. The Snakefile includes rules to:  

* Generate **multiple alignments** of all the genes included in the considered orthogroups.  
* Build **gene trees** starting from these alignments.  
* Parse the gene trees to **date the origin** of all paralogs through a species-overal approach.  


### Figures

This folder contains two Rmarkdown files used to generate the majority of the figure panels associated with the publication, and perform very few of the less demanding computational analyses. 

#### Main figures    

The code to generate the following figures is available in **All_Figures.Rmd**.  

* **Figure 1**: Framework overview and definition of sequence and expression similarities.  
* **Figure 2**: Correlation of sequence and expression similarities between and within clades.  
* **Figure 3**: Characterization of fast and slow evolving genes.  
* **Figure 4**: Functional categories with common molecular diversification patterns between clades.  
* **Figure 5**: Functional categories with molecular diversification biases between clades.   
* **Supplementary Figure 1**: Schematics for computation of sequence/expression similarity and definition of extra 1:1 orthogroup sets.  
* **Supplementary Figure 2**: GO enrichments and GSEA performed with fruit fly based GO annotations.  
* **Supplementary Figure 3**: Deltas of sequence and expression similarities between vertebrates and insects.  
* **Supplementary Figure 4**: Comparison between alternative measures of sequence conservation.  
* **Supplementary Figure 5**: Comparison between sequence similarities computation methods.  
* **Supplementary Figure 9**: Comparison of sequence similarities computed upon different phylogenetic setups.  
* **Supplementary Figure 10**: Influence of paralog number and paralog origin on sequence similarities.  
* **Supplementary Figure 11**: Comparison between alternative measures of expression conservation.  

#### Alternative setups figures  

The code to generate the following figures is available in **Alternative_setups_Figures.Rmd**.  

* **Supplementary Figure 6**: Key results repeated by computing sequence similarities with the blosum62 computation method.  
* **Supplementary Figure 7**: Key results repeated by computing sequence similarities with the blosum45 computation method.  
* **Supplementary Figure 8**: Key results repeated by computing sequence similarities not considering leading and trailing gaps.  
* **Supplementary Figure 12**: Key results repeated by computing expression similarities from summed paralog expression.  
* **Supplementary Figures 13-16**: Key results repeated with extra paralogous sets (BA-seq, BD-seq, BA-expr, BD-expr).   


