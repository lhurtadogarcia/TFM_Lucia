################################################################
#                    Common prevalent DEGs                     #
#         Lucía Hurtado de Mendoza García de La Chica          #      
################################################################


## In this script, we'll create a Venn Diagram to visualize the number of common prevalent genes in all 4 studies.
## First, we'll set the working directory and load the package 'VennDiagram'

#setwd("~/TFM/your_path") 
library(VennDiagram)

## Now we'll fetch the names of prevalent genes from all studies: 

# DMD_hiseq4k_single

DMD_hiseq4k_single <- read.table('DMD_hiseq4k_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

## Now we fetch the prevalent genes: 

prev_DMD_hiseq_genes <- DMD_hiseq4k_single[DMD_hiseq4k_single$genes_tag == "PREVALENT_DEG", ]
dim(prev_DMD_hiseq_genes)

prev_DMD_hiseq_names <- row.names(prev_DMD_hiseq_genes) # prevalent genes' names
length(prev_DMD_hiseq_names)

# We'll print it out in a .txt file:

write(prev_DMD_hiseq_names, file = "prevalent_DMD_hiseq4k_single_genes.txt", sep = '\n')

# Nextseq500_paired

DMD_nextseq500_paired <- read.table('DMD_nextseq500_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

prev_DMD_nextseq_paired_genes <- DMD_nextseq500_paired[DMD_nextseq500_paired$genes_tag == "PREVALENT_DEG", ]
dim(prev_DMD_nextseq_paired_genes)


prev_DMD_nextseq_paired_names <- row.names(prev_DMD_nextseq_paired_genes) # names
length(prev_DMD_nextseq_paired_names)

# We'll print it out in a .txt file:

write(prev_DMD_nextseq_paired_names, file = "prevalent_nextseq500_paired_genes.txt", sep = '\n')

# Nextseq500_single

DMD_nextseq500_single <- read.table('DMD_nextseq500_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

prev_DMD_nextseq_single_genes <- DMD_nextseq500_single[DMD_nextseq500_single$genes_tag == "PREVALENT_DEG", ]
dim(prev_DMD_nextseq_single_genes)


prev_DMD_nextseq_single_names <- row.names(prev_DMD_nextseq_single_genes) # names
length(prev_DMD_nextseq_single_names)

# We'll print it out in a .txt file:

write(prev_DMD_nextseq_single_names, file = "prevalent_nextseq500_single_genes.txt", sep = '\n')

# ALS_hiseq4k_paired

ALS_hiseq4k_paired <- read.table('ALS_hiseq4k_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

prev_ALS_hiseq4k_paired_genes <- ALS_hiseq4k_paired[ALS_hiseq4k_paired$genes_tag == "PREVALENT_DEG", ]
dim(prev_ALS_hiseq4k_paired_genes)


prev_ALS_hiseq4k_paired_names <- row.names(prev_ALS_hiseq4k_paired_genes) # names
length(prev_ALS_hiseq4k_paired_names)

# We'll print it out in a .txt file:

write(prev_ALS_hiseq4k_paired_names, file = "prevalent_ALS_hiseq4k_paired_genes.txt", sep = '\n')

# Now we'll generate the Venn Diagram:

# First we prepare a palette of 4 colors with RColorBrewer:

library(RColorBrewer)

colors <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(prev_DMD_hiseq_names, prev_DMD_nextseq_names, prev_DMD_nextseq_single_names, prev_ALS_hiseq4k_paired_names),
  category.names = c("DMD_hiseq4k_single", "DMD_nextseq500_paired", "DMD_nextseq500_single", "ALS_hiseq4k_paired"),
  filename = 'common_prev_genes_venn_diagram.png',
  output = T,

  # Output features:
  imagetype = 'png',
  main = 'Common Prevalent Genes',

  # Circles:
  lwd = 2,
  lty = 'blank',
  fill = colors,

  # Numbers:
  cex = .6,
  fontface = 'bold',
  fontfamily = 'sans',

  # # Set names:
  cat.cex = 0.5,
  cat.fontface = 'bold',
  cat.default.pos = 'outer'

)

# By looking at the diagram, we can see that the 3 DMD studies have a gene in common.
# All 4 studies don't share any common prevalent genes. 

# We'll first take a look at the common gene between all 3 DMD studies: 
# First we'll compare DMD_hiseq4k_single and DMD_nextseq500_paired:

Reduce(intersect, list(prev_DMD_nextseq_paired_names, prev_DMD_hiseq_names, prev_DMD_nextseq_single_names))

# The gene ENSG00000137573 corresponds to SULF1. 

# We will now check whether this gene is up- or downregulated in all three datasets.

DMD_hiseq4k_single["ENSG00000137573", ] # upregulated
DMD_nextseq500_paired["ENSG00000137573", ] # upregulated
DMD_nextseq500_single["ENSG00000137573", ] # downregulated


## THE END ##




