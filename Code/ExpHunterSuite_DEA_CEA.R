################################################################
#              DEA, CEA and Functional Analysis                #
#         Lucía Hurtado de Mendoza García de La Chica          #  
################################################################

# This script was created by taking the example pipeline from https://github.com/seoanezonjic/ExpHunterSuite
# and adjusting it to the newer version of ExpHunter Suite, as well as setting the parameters needed for the study.

library('ExpHunterSuite')

# First we need to set the design and count table to their respective variables:

target <- target_generation(from_file='control_vs_dmd_target.txt') # Read experiment design which describes the sample groups
raw_count_table <- read.table('final_counts.txt', header=TRUE, row.names=1, sep="\t") # Read the table counts with the number of reads per gene.

# Then we apply the function main_degenes_Hunter to obtain the differentially expressed genes (DEGs) and differentially
# coexpressed genes (CEGs), we write them to different files and create an html report. 

final_results <- main_degenes_Hunter( 
  target=target,
  raw=raw_count_table,
  modules='WDELN',
  p_val_cutoff = 0.05,
  lfc = 1, 
  WGCNA_mergecutHeight = 0.1,
  WGCNA_min_genes_cluster = 15,
  WGCNA_detectcutHeight = 0.995
  )

write.table(final_results[['raw_filter']], "filtered_count_data.txt", quote=FALSE, col.names=NA, sep="\t") # Raw table filtered by minreads parameter
write_df_list_as_tables(final_results[['all_data_normalized']], prefix = 'Normalized_counts_') #Normalized table by each expresion package used in modules argument.
write_df_list_as_tables(final_results[['all_counts_for_plotting']], prefix = 'allgenes_')
write.table(final_results[['DE_all_genes']], "hunter_results_table.txt", quote=FALSE, row.names=TRUE, sep="\t") # Table with all the expresion packages and the integrated results (gene_tag, combined p-value and combined log2FC)
write_expression_report(final_results) # Generate friendly html report with expresion data.

# For the funtional analysis, we use main_functional_hunter, adjust the parameters and get enrichment tables and an html report
# as output. 

func_results <- main_functional_hunter( #Perform enrichment analysis
  final_results,
  input_gene_id = "ENSEMBL",
  model_organism = 'Human', #Use specified organism database
  enrich_dbs = c("Reactome", "MF", "BP", "CC"), 
  enrich_methods = "ORA", #Perform enrichment  using only Overepresentation analysis (Not GSEA)
  pthreshold = 0.1,
  output_files = 'results_FEnr'
)

write_enrich_files(func_results) #Write enrichment tables
write_functional_report(final_results, func_results) # Generate friendly html report with functional data.


## THE END ##

