################################################################
#                         Common CEGs                          #
#         Lucía Hurtado de Mendoza García de La Chica          #      
################################################################

# We will first load the datasets:

DMD_hiseq4k_single <- read.table('DMD_hiseq4k_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

DMD_nextseq500_paired <- read.table('DMD_nextseq500_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

DMD_nextseq500_single <- read.table('DMD_nextseq500_single/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

ALS_hiseq4k_paired <- read.table('ALS_hiseq4k_paired/hunter_results_table.txt', header = T, sep = '\t', dec = '.')

# We will now fetch the genes from each CEG, considering as such the closest cluster to the control/treat vectors. 

# DMD_hiseq4k_single dataset:
DMD_hiseq4k_ME2 <- DMD_hiseq4k_single[DMD_hiseq4k_single$Cluster_ID == '2', ]
dim(DMD_hiseq4k_ME2)
DMD_hiseq4k_ME2_names <- row.names(DMD_hiseq4k_ME2) # names
length(DMD_hiseq4k_ME2_names)

# DMD_nextseq500_paired:
DMD_nextseq500_paired_ME26 <- DMD_nextseq500_paired[DMD_nextseq500_paired$Cluster_ID == '26', ]
dim(DMD_nextseq500_paired_ME26)
DMD_nextseq500_paired_ME26_names <- row.names(DMD_nextseq500_paired_ME26) # names
length(DMD_nextseq500_paired_ME26_names)

# DMD_nextseq500_single:
DMD_nextseq500_single_ME1 <- DMD_nextseq500_single[DMD_nextseq500_single$Cluster_ID == '1', ]
dim(DMD_nextseq500_single_ME1)
DMD_nextseq500_single_ME1_names <- row.names(DMD_nextseq500_single_ME1) # names
length(DMD_nextseq500_single_ME1_names)

DMD_nextseq500_single_ME4 <- DMD_nextseq500_single[DMD_nextseq500_single$Cluster_ID == '4', ]
dim(DMD_nextseq500_single_ME4)
DMD_nextseq500_single_ME4_names <- row.names(DMD_nextseq500_single_ME4) # names
length(DMD_nextseq500_single_ME4_names)

# ALS_hiseq4k_paired:
ALS_hiseq4k_ME25 <- ALS_hiseq4k_paired[ALS_hiseq4k_paired$Cluster_ID == '25', ]
dim(ALS_hiseq4k_ME25)
ALS_hiseq4k_ME25_names <- row.names(ALS_hiseq4k_ME25) # names
length(ALS_hiseq4k_ME25_names)

# We will now visualize these CEGs in a Venn Diagram:

library(VennDiagram)
library(RColorBrewer)

colors <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(DMD_hiseq4k_ME2_names, DMD_nextseq500_paired_ME26_names, DMD_nextseq500_single_ME1_names, ALS_hiseq4k_ME25_names),
  category.names = c("DMD_hiseq4k_single", "DMD_nextseq500_paired", "DMD_nextseq500_single", "ALS_hiseq4k_paired"),
  filename = 'common_CEGs_venn_diagram.png',
  output = T,
  
  # Output features:
  imagetype = 'png',
  main = 'Common CEGs',
  
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

# Common genes in CEGs: 

common_CEGs <- Reduce(intersect, list(DMD_hiseq4k_ME2_names, DMD_nextseq500_paired_ME26_names, DMD_nextseq500_single_ME4_names, 
                                      ALS_hiseq4k_ME25_names ))
length(common_CEGs)