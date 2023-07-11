################################################################
#                    Common Possible DEGs                      #
#         Lucía Hurtado de Mendoza García de La Chica          #      
################################################################


## In this script, we'll create a Venn Diagram to visualize the number of common possible genes in all 4 studies.
## First, we'll set the working directory and load the package 'VennDiagram'

#setwd("~/TFM/your_path")
library(VennDiagram)

## Now we'll fetch the names of possible genes from all studies: 

# DMD_hiseq4k_single

DMD_hiseq4k_single <- read.table('DMD_hiseq4k_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

## Now we fetch the possible genes: 

poss_DMD_hiseq_genes <- DMD_hiseq4k_single[DMD_hiseq4k_single$genes_tag == "POSSIBLE_DEG", ]
dim(poss_DMD_hiseq_genes)

poss_DMD_hiseq_names <- row.names(poss_DMD_hiseq_genes) # prevalent genes' names
length(poss_DMD_hiseq_names) 

# We'll print it out in a .txt file:

write(poss_DMD_hiseq_names, file = "possible_DMD_hiseq4k_single_genes.txt", sep = '\n')

# DMD_nextseq500_paired

DMD_nextseq500_paired <- read.table('DMD_nextseq500_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

poss_DMD_nextseq_paired_genes <- DMD_nextseq500_paired[DMD_nextseq500_paired$genes_tag == "POSSIBLE_DEG", ]
dim(poss_DMD_nextseq_paired_genes)


poss_DMD_nextseq_paired_names <- row.names(poss_DMD_nextseq_paired_genes) # names
length(poss_DMD_nextseq_paired_names)

# We'll print it out in a .txt file:

write(poss_DMD_nextseq_names, file = "possible_DMD_nextseq500_paired_genes.txt", sep = '\n')

# DMD_nextseq500_single

DMD_nextseq500_single <- read.table('DMD_nextseq500_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

poss_DMD_nextseq_single_genes <- DMD_nextseq500_single[DMD_nextseq500_single$genes_tag == "POSSIBLE_DEG", ]
dim(poss_DMD_nextseq_single_genes)


poss_DMD_nextseq_single_names <- row.names(poss_DMD_nextseq_single_genes) # names
length(poss_DMD_nextseq_single_names)

# We'll print it out in a .txt file:

write(poss_DMD_nextseq_single_names, file = "possible_DMD_nextseq500_single_genes.txt", sep = '\n')

# ALS_hiseq4k_paired

ALS_hiseq4k_paired <- read.table('ALS_hiseq4k_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

poss_ALS_hiseq4k_paired_genes <- ALS_hiseq4k_paired[ALS_hiseq4k_paired$genes_tag == "POSSIBLE_DEG", ]
dim(poss_ALS_hiseq4k_paired_genes)


poss_ALS_hiseq4k_paired_names <- row.names(poss_ALS_hiseq4k_paired_genes) # names
length(poss_ALS_hiseq4k_paired_names)

# We'll print it out in a .txt file:

write(poss_ALS_hiseq4k_paired_names, file = "possible_ALS_hiseq4k_paired_genes.txt", sep = '\n')

# Now we'll generate the Venn Diagram:

# First we prepare a palette of 4 colors with RColorBrewer:

library(RColorBrewer)

colors <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(poss_DMD_hiseq_names, poss_DMD_nextseq_names, poss_DMD_nextseq_single_names, poss_ALS_hiseq4k_paired_names),
  category.names = c("DMD_hiseq4k_single", "DMD_nextseq500_paired", "DMD_nextseq500_single", "ALS_hiseq4k_paired"),
  filename = 'common_possible_genes_venn_diagram.png',
  output = T,

  # Output features:
  imagetype = 'png',
  main = 'Common Possible Genes',

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

# By looking at the diagram, we can see that all studies have 1 gene in common. 


# Common possible genes
Reduce(intersect, list(poss_DMD_nextseq_paired_names, poss_DMD_hiseq_names, poss_DMD_nextseq_single_names, poss_ALS_hiseq4k_paired_names))

# The gene ENSG00000152527 corresponds to PLEKHH2.

# We will now fetch the respective cluster where  PLEKHH2 is in each study: 

DMD_hiseq4k_single["ENSG00000152527", ] # cluster 1
DMD_nextseq500_paired["ENSG00000152527", ] # cluster 13
DMD_nextseq500_single["ENSG00000152527", ] # cluster 3
ALS_hiseq4k_paired["ENSG00000152527", ] # cluster  7

# We'll get the clusters:

DMD_hiseq4k_ME1 <- DMD_hiseq4k_single[DMD_hiseq4k_single$Cluster_ID == '1', ]
dim(DMD_hiseq4k_ME1)

DMD_nextseq500_paired_ME13 <- DMD_nextseq500_paired[DMD_nextseq500_paired$Cluster_ID == '13', ]
dim(DMD_nextseq500_paired_ME13)

DMD_nextseq500_single_ME3 <- DMD_nextseq500_single[DMD_nextseq500_single$Cluster_ID == '3', ]
dim(DMD_nextseq500_single_ME3)

ALS_hiseq4k_ME7 <- ALS_hiseq4k_paired[ALS_hiseq4k_paired$Cluster_ID == '7', ]
dim(ALS_hiseq4k_ME7)

# We want to check whether TNNT3, MBNL1 and TNNT2 are present in these clusters alongside
# PLEKHH2. The relationship between these 4 genes was seen on the STRING database by creating a network
# of less than 20 interactors.
# TNNT3: ENSG00000130595
# MBNL1: ENSG00000152601
# TNNT2: ENSG00000118194

# We'll start with TNNT3:
DMD_hiseq4k_ME1["ENSG00000130595", ] # prevalent
DMD_nextseq500_paired_ME13["ENSG00000130595", ]
DMD_nextseq500_single_ME3["ENSG00000130595", ]
ALS_hiseq4k_ME7["ENSG00000130595", ]

# MBNL1
DMD_hiseq4k_ME1["ENSG00000152601", ] #not_deg
DMD_nextseq500_paired_ME13["ENSG00000152601", ]
DMD_nextseq500_single_ME3["ENSG00000152601", ]
ALS_hiseq4k_ME7["ENSG00000152601", ]

# TNNT2
DMD_hiseq4k_ME1["ENSG00000118194", ] # prevalent
DMD_nextseq500_paired_ME13["ENSG00000118194", ]
DMD_nextseq500_single_ME3["ENSG00000118194", ]
ALS_hiseq4k_ME7["ENSG00000118194", ]

# As seen above, TNNT3 is present in cluster 1 alongside PLEKHH2 and it's differentially expressed, 
# as it has been detected by all four packages (DESeq2, limma, edgeR and NOISeq). MBNL1 is also found 
# in cluster 1 but no packages detected it as a differentially expressed gene. Lastly, TNNT2 is also 
# categorized as prevalent in cluster 1 of the DMD_hiseq4k_single study, with a considerable mean_logFCs of 3.59.
# Neither of those genes were found in the rest of the respective clusters for each study. 

# We will now check if these two genes are present at all in the rest of the studies and, if present,
# the type of gene classification they belong to: 

DMD_nextseq500_paired["ENSG00000130595", ] # cluster 2, not_deg
DMD_nextseq500_paired["ENSG00000152601", ] # cluster 1, not_deg
DMD_nextseq500_paired["ENSG00000118194", ] # cluster 2, not_deg

DMD_nextseq500_single["ENSG00000130595", ] # filtered out
DMD_nextseq500_single["ENSG00000152601", ] # cluster 1, not_deg
DMD_nextseq500_single["ENSG00000118194", ] # cluster 5, not_deg

ALS_hiseq4k_paired["ENSG00000130595", ] # cluster 3, possible_deg
ALS_hiseq4k_paired["ENSG00000152601", ] # cluster 23, not_deg
ALS_hiseq4k_paired["ENSG00000118194", ] # cluster 1, possible_deg

## THE END ##








