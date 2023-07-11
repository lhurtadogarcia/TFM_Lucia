################################################################
#              ExpHunter Suite format preparation              #
#         Lucía Hurtado de Mendoza García de La Chica          #  
################################################################

# This script exemplifies the preparation of the data for ExpHunter Suite once we have the count table. 
# We need another file that sets whether each sample belongs to the 'control' or 'treatment' group. 
# In order to do so, we create a table with each sample in the first column and the respective tag 'Ctrl' or 
#  'Treat' in the second column.
# This example in particular was made for the DMD_nextseq500_paired dataset, which consists of 18 samples:
# 9 from the control group and 9 from the treated group. 

# Samples from both conditions are labeled 'control' or 'treat' and a consecutive number. Thus, we'll create the first 
# column of the tables by doing so:

ctrl <- paste0("control", seq(1:9), sep = '') 

treated <- paste0("treat", seq(1:9), sep = '') 

# Then, we'll join all of the samples' labels:

labels <- c(ctrl, treated)

# Now we'll create the respective labels for each sample:

treat <- c(rep('Ctrl', 9), rep('Treat', 9))

# We'll create the table and name the columns as 'sample' and 'treat' for ExpHunter Suite:

table <- cbind(labels, treat)
colnames(table) <- c('sample', 'treat')


# Lastly, we'll export the table:

write.table (table,  
             file= "ctrl_vs_treat_target.txt", 
             sep = "\t", row.names=F, col.names = T, quote = F)
