# Adam Berman
# Princeton University
# 7 May 2018

# Change this path as appropriate
setwd("/Users/adamberman/metabolite-scorer")

library(futile.matrix)
library(reshape2)
library(KEGGREST)
library(plyr)
library(ROCR)
library(ggplot2)



#----------Read in raw data----------------------------------------------------------------
#----------Read in raw data----------------------------------------------------------------
#----------Read in raw data----------------------------------------------------------------

# Read raw expression data into a data frame, storing only the gene name, the expression value,
# and the indicator of whether a sample is from cancer tissue or healthy tissue
# "exp_array.BRCA-US.tsv", the raw expressional data, can be downloaded at the 
# following link: https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-US/exp_array.BRCA-US.tsv.gz
rnaseq_data = read.delim("exp_array.BRCA-US.tsv", colClasses = c("NULL","NULL","NULL","NULL","NULL","NULL","NULL",NA,NA,"NULL","NULL","NULL","NULL","NULL","NULL",NA,"NULL"))

# Read in metabolites and their associated genes
# The 'hmdb_metabolite_table.txt' file is created by pipelining the output of 
# make_hmdb_metabolites.py into a text file, where make_hmdb_metabolites_mod.py filters out the 20,440 unique endogenous metabolites  
# (as well as kegg_metabolite_table.txt itself) as can be found at: https://github.com/agberman/metabolite-scorer-v2.git
# Number of endogenous metabolites: 20,440
# Number of endogenous metabolites after merging synonyms: 1,120
# 19,522 rows / 1,120 metabolites = 17.4 genes associated with each metabolite on average
metabolites_with_associated_genes <- read.table('hmdb_metabolite_table.txt', header = TRUE, sep = '\t', quote = "", comment.char = "") 

# Read in genes and their gene lengths
# The amino_acid_gene_lengths file can be found in the Github repository at 
# the following link: https://github.com/agberman/metabolite-scorer-v2.git
gene_lengths <- read.table('amino_acid_gene_lengths', header = TRUE, sep = '\t', quote = "", comment.char = "")

# Read in raw TCGA mutation data
# TCGA_BRCA_mutect.txt can be downloaded in the Github repo above, or at the 
# following link: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/
all_tcga_data <- read.table("TCGA_BRCA_mutect.txt", sep="\t", quote = "", header=T)

# Get just the rows of the data that we want
partial_tcga_data <- data.frame(all_tcga_data$Hugo_Symbol, all_tcga_data$Gene, all_tcga_data$Tumor_Sample_Barcode, 
                                all_tcga_data$HGVSp, all_tcga_data$Protein_position, all_tcga_data$Variant_Classification)
colnames(partial_tcga_data) <- c("Hugo", "Ensemble", "Sample", "HGVSp", "Protein_position", "Variant_classification")
nrow(partial_tcga_data) # 120,987 mutations total

# Filter out all mutations that are not missense mutations
partial_tcga_data_filtered <- partial_tcga_data[partial_tcga_data$Variant_classification == "Missense_Mutation", ]
nrow(partial_tcga_data_filtered) # 62,315 missense mutations / 98,596 non-silent mutations










#----------Retrieve KEGG data----------------------------------------------------------------
#----------Retrieve KEGG data----------------------------------------------------------------
#----------Retrieve KEGG data----------------------------------------------------------------

# Read in raw KEGG data: all metabolic pathways found in humans ("hsa")
paths <- names(keggList("pathway", "hsa"))
length(paths) # 327 paths

# Initialize vectors to hold pathways and their constituent 
# metabolites (paths_with_metabolites) and metabolite names (names)
paths_with_metabolites <- vector(mode="list", length=length(paths))
names <- vector(mode="character", length=length(paths))

# Initialize lists to hold metabolites found in all pathways (met_pool)
# and metablites found only in cancer-linked pathways (cancer_pool)
met_pool <- character()
cancer_pool <- character()

# Iterate over all pathways
for (i in 1:length(paths)){
  print(paste("i:", i))
  # Parse out all metabolites in given pathway
  dat <- keggGet(paths[i])
  names[i] <- dat[[1]]["NAME"][[1]]
  metabolites_kegg <- names(dat[[1]]["COMPOUND"][[1]])
  
  # Add all metabolites to met_pool
  met_pool <- c(met_pool, metabolites_kegg)
  
  # Add metabolites to cancer_pool if pathway is cancer-linked (pathways 279-301 are cancer-linked: 23 of 327 paths)
  if ((i >= 279) && (i <= 301)) {
    cancer_pool <- c(cancer_pool, metabolites_kegg)
  }
  
  # Update paths_with_metabolites with new pathway
  paths_with_metabolites[[i]] <- metabolites_kegg
}

# Associate names with the stored pathways
pwp <- paths_with_metabolites
names(pwp) <- names

# Remove repeated metabolites in cancer_pool
cancer_pool <- unique(cancer_pool)
length(cancer_pool) # 166 unique metabolites in KEGG cancer-linked pathays
length(unique(met_pool)) # 3,454 unique metabolites in all KEGG pathways (5,905 with repeats)










#----------Compute mutational scores----------------------------------------------------------------
#----------Compute mutational scores----------------------------------------------------------------
#----------Compute mutational scores----------------------------------------------------------------

# Copy data
partial_data_filtered <- partial_tcga_data_filtered

# Add mutation count row (all ones because each row represents one mutation)
partial_data_filtered$Mutation_count <- rep(1,nrow(partial_data_filtered))

# Sum together total number of mutations to each gene
gene_scores <- aggregate(partial_data_filtered$Mutation_count, by=list(gene_symbol=partial_data_filtered$Hugo), FUN=sum)

# Rename gene_scores columns
colnames(gene_scores) <- c('gene_symbol', 'total_mutations')

# Match genes to gene lengths
gene_scores$gene_length <- gene_lengths$LENGTH[match(gene_scores$gene_symbol, gene_lengths$Hugo_Symbol)]
nrow(gene_scores)
# Note that there are 15,523 distinct genes mutated in the whole new TCGA data set
# There are an estimated 19,000-20,000 protein-coding genes in the human genome

# Remove genes with no gene length
gene_scores <- na.omit(gene_scores)
nrow(gene_scores)
# Note that there are 14,916 distinct genes mutated in the whole new TCGA data set that have gene lengths in my data set associated with them

# Compute mutational gene subscores by dividing total mutations by gene length
gene_scores$mutational_gene_subscore <- (gene_scores$total_mutations / gene_scores$gene_length)

# Sort mutational gene subscores
gene_scores_sorted <- gene_scores[with(gene_scores, order(-mutational_gene_subscore)), ]
head(gene_scores_sorted, 10) # TP53: 0.52162850, PIK3CA: 0.32022472, BRIP1: 0.08737864 all known to be associated to cancer

# Plot a histogram of mutational gene subscores
hist(gene_scores$mutational_gene_subscore, xlab = 'Mutational Gene Subscore', main = 'Distribution of Mutational Gene Subscores')

# Match mutational gene subscores to genes in metabolite list
metabolites <- metabolites_with_associated_genes
metabolites$mutational_gene_subscore <- gene_scores$mutational_gene_subscore[match(metabolites$Gene.Name, gene_scores$gene_symbol)]
nrow(metabolites) # 19,522 associated genes over all metabolites

# Remove rows where gene has no mutational gene subscore
metabolites <- na.omit(metabolites)
nrow(metabolites) # 16,184 associated genes over all metabolites with mutational gene subscores

# Sort metabolites by decreasing mutational gene subscore
#metabolites_sorted <- metabolites[with(metabolites, order(-mutational_gene_subscore)), ]

# Count number of genes associated with each metabolite
metabolites$genecount <- 1
metabolites_genecounts <- aggregate(genecount~Metabolite, metabolites, FUN=sum)
mean(metabolites_genecounts$genecount) # Average of 15.08 genes associated with each metabolite

# Find average of mutational gene subscores associated with each metabolite
metabolites_scored <- aggregate(mutational_gene_subscore~Metabolite, metabolites, FUN=mean)

# Re-append gene counts
metabolites_scored$genecount <- metabolites_genecounts$genecount

# CLEAN UP: Rename columns and return KEGG IDs to each metabolite after aggregation, then reorder columns
colnames(metabolites_scored) <- c('metabolite', 'mutational_score', 'gene_count')
metabolites_scored$KEGG_ID <- metabolites$KEGG[match(metabolites_scored$metabolite, metabolites$Metabolite)]
metabolites_scored <- metabolites_scored[,c(1,4,3,2)]

# Sort metabolites by descending mutational score
metabolites_scored_sorted <- metabolites_scored[with(metabolites_scored, order(-mutational_score)), ]
nrow(metabolites_scored_sorted) # 1,073 metabolites in total
head(metabolites_scored_sorted, 10)

# Remove all metabolites with fewer than three associated genes
metabolites_scored_sorted_threecut <- metabolites_scored_sorted[metabolites_scored_sorted$gene_count >= 3, ]
# 688 metabolites with three or more associated genes in TCGA-BRCA

# Normalize mutational scores to be between 0 and 1
metabolites_scored_sorted_threecut$normalized_mutational_score <- ((metabolites_scored_sorted_threecut$mutational_score - min(metabolites_scored_sorted_threecut$mutational_score)) / 
                                                                     (max(metabolites_scored_sorted_threecut$mutational_score) - min(metabolites_scored_sorted_threecut$mutational_score)))
metabolites_scored_sorted$normalized_mutational_score <- ((metabolites_scored_sorted$mutational_score - min(metabolites_scored_sorted$mutational_score)) / 
                                                            (max(metabolites_scored_sorted$mutational_score) - min(metabolites_scored_sorted$mutational_score)))
head(metabolites_scored_sorted_threecut, 10)


# Copy results into new data frames
mutational_results <- metabolites_scored_sorted_threecut # only metabolites with a minimum of three associated genes
mutational_results_full <- metabolites_scored_sorted # all metabolites

# Change column names
colnames(mutational_results) <- c("Metabolite", "KEGG_ID", "Gene_Count", "Mutation_Score", "Normalized_Score")
colnames(mutational_results_full) <- c("Metabolite", "KEGG_ID", "Gene_Count", "Mutation_Score", "Normalized_Score")

# Generate a histogram of the mutation scores
hist(mutational_results$Normalized_Score, xlab = 'Normalized Mutational Subscore', main = 'Distribution of Normalized Mutational Subscores')
hist(mutational_results_full$Normalized_Score, xlab = 'Normalized Mutational Subscore', main = 'Distribution of Normalized Mutational Subscores')

# Write out data for long-term storage
write.csv(mutational_results, file = "mutational_results.csv")
write.csv(mutational_results_full, file = "mutational_results_full.csv")










#----------Compute structural scores----------------------------------------------------------------
#----------Compute structural scores----------------------------------------------------------------
#----------Compute structural scores----------------------------------------------------------------

# Copy data
partial_data_filtered_struct <- partial_tcga_data_filtered

# Add column with mutation cleaned positions
partial_data_filtered_struct$Start_position <- NA
partial_data_filtered_struct$End_position <- NA

# Get start position and end position of mutations
for(i in 1:nrow(partial_data_filtered_struct)) 
{
  split <- strsplit(toString(partial_data_filtered_struct[i,5]), "[[:punct:]]")
  if (length(split[[1]]) > 0) {
    if (length(split[[1]]) == 2) {
      partial_data_filtered_struct[i,7] <- split[[1]][1]
      partial_data_filtered_struct[i,8] <- split[[1]][1]
    } else if (length(split[[1]]) == 3) {
      if (nchar(split[[1]][1]) > 0) {
        partial_data_filtered_struct[i,7] <- split[[1]][1]
        partial_data_filtered_struct[i,8] <- split[[1]][2]
      }
    }
  }
}

# Save and load data frame for future use
save(partial_data_filtered_struct,file="partial_data_filtered_struct.Rda")
load("partial_data_filtered_struct.Rda")
nrow(partial_data_filtered_struct) # 62,315 missense mutations with positional information

# Remove rows with no positional information
partial_data_filtered_struct <- partial_data_filtered_struct[!is.na(partial_data_filtered_struct$Start_position),]
nrow(partial_data_filtered_struct) # 62,288 missense mutations with positional information

# Simplify data to be written out
writeout_data <- data.frame(Ensemble = partial_data_filtered_struct$Ensemble, Start_position = partial_data_filtered_struct$Start_position, End_position = partial_data_filtered_struct$End_position)

# Write positional data to a CSV file
write.csv(writeout_data, file = "mutations_with_positions.csv")


#- - - - - - - - - - - - - - RUN genes_and_bindingweights.py (assuming that you have already run forsean_bindingweights_mod.py to generate genes_and_bindingweights.json) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Read in CSV containing genes and the number of mutations from the TCGA data that were considered positionally significant
# (binding weights at that position >= 0.1) according to Shilpa's positional Biolip data
genes_with_num_sig_mutations <- read.csv("genes_with_num_sig_mutations.csv", header = FALSE)
head(genes_with_num_sig_mutations)
colnames(genes_with_num_sig_mutations) <- c("Ensemble", "Number_of_positionally_significant_mutations")

# Basic specs of data
genes_with_num_sig_mutations[which.max(genes_with_num_sig_mutations$Number_of_positionally_significant_mutations),] # Titin gene has 115 positionally significant mutations
nrow(genes_with_num_sig_mutations) # 15,486 genes in all
hist(genes_with_num_sig_mutations$Number_of_positionally_significant_mutations)

# Check
genes_with_num_sig_mutations <- na.omit(genes_with_num_sig_mutations)
nrow(genes_with_num_sig_mutations) # still 15,486 genes in all

# Re-attach Hugo gene names
genes_with_num_sig_mutations$Hugo <- partial_data_filtered_struct$Hugo[match(genes_with_num_sig_mutations$Ensemble, 
                                                                             partial_data_filtered_struct$Ensemble)]

genes_with_num_sig_mutations <- genes_with_num_sig_mutations[c(3,2)]
head(genes_with_num_sig_mutations)

#---------------------------Above does not need to be rerun-------------------------------------

# Copy genes_with_num_sig_mutations
gene_scores_struct <- genes_with_num_sig_mutations

# Rename gene_scores_struct columns
colnames(gene_scores_struct) <- c('gene_symbol', 'total_sig_mutations')

# Match genes to gene lengths
gene_scores_struct$gene_length <- gene_lengths$LENGTH[match(gene_scores_struct$gene_symbol, gene_lengths$Hugo_Symbol)]
nrow(gene_scores_struct)
# Note that there are 15,486 distinct genes mutated in the whole new TCGA data set
# There are an estimated 19,000-20,000 protein-coding genes in the human genome

# Remove genes with no gene length
gene_scores_struct <- na.omit(gene_scores_struct)
nrow(gene_scores_struct) 
# Note that there are 14,891 distinct genes mutated in the whole new TCGA data set that have gene lengths associated with them

# Compute structural gene subscores by dividing total mutations by gene length
gene_scores_struct$structural_gene_subscore <- (gene_scores_struct$total_sig_mutations / gene_scores_struct$gene_length)

# Sort structural gene subscores
gene_scores_struct_sorted <- gene_scores_struct[with(gene_scores_struct, order(-structural_gene_subscore)), ]
head(gene_scores_struct_sorted, 10) # TP53, HIST1H2AE, AKT1, HIST1H2AK  / CDHI (Cadherin-1), GATA3 all known to be associated to breast cancer specifically

# Plot a histogram of structural gene subscores
hist(gene_scores_struct$structural_gene_subscore, xlab = 'Structual Gene Subscore', main = 'Distribution of Structural Gene Subscores')

# Copy fresh data
metabolites_struct <- metabolites_with_associated_genes

# Match structural gene subscores to genes in metabolite list
metabolites_struct$structural_gene_subscore <- gene_scores_struct$structural_gene_subscore[match(metabolites_struct$Gene.Name, 
                                                                                                 gene_scores_struct$gene_symbol)]
nrow(metabolites_struct) # 19,522 associated genes over all metabolites

# Remove rows where gene has no structural gene subscore
metabolites_struct <- na.omit(metabolites_struct)
nrow(metabolites_struct) # 16,184 associated genes over all metabolites with structural gene subscores

# Append count of genes
metabolites_struct$genecount <- 1
metabolites_struct_genecounts <- aggregate(genecount~Metabolite, metabolites_struct, FUN=sum)
mean(metabolites_struct_genecounts$genecount) # Average of 15.08 genes associated with each metabolite
head(metabolites_struct_genecounts)

# Find average of structural gene subscores associated with each metabolite
metabolites_struct_scored <- aggregate(structural_gene_subscore~Metabolite, metabolites_struct, FUN=mean)
head(metabolites_struct_scored)

# Append gene counts
metabolites_struct_scored$genecount <- metabolites_struct_genecounts$genecount
head(metabolites_struct_scored)

# CLEAN UP: Rename columns and return KEGG IDs to each metabolite after aggregation, then reorder columns
colnames(metabolites_struct_scored) <- c('metabolite', 'structural_score', 'gene_count')
metabolites_struct_scored$KEGG_ID <- metabolites_struct$KEGG[match(metabolites_struct_scored$metabolite, metabolites_struct$Metabolite)]
metabolites_struct_scored <- metabolites_struct_scored[,c(1,4,3,2)]
head(metabolites_struct_scored)

# Sort metabolites by descending mutational score
metabolites_struct_scored_sorted <- metabolites_struct_scored[with(metabolites_struct_scored, order(-structural_score)), ]
nrow(metabolites_struct_scored_sorted) # 1,073 metabolites in total
head(metabolites_struct_scored_sorted, 10) # Notably, Sialyl-Lewis X is known to be cancer-associated, and is female-specific

# Remove all metabolites with fewer than three associated genes
metabolites_struct_scored_sorted_threecut <- metabolites_struct_scored_sorted[metabolites_struct_scored_sorted$gene_count >= 3, ]
nrow(metabolites_struct_scored_sorted_threecut) # 688 metabolites with three or more associated genes
head(metabolites_struct_scored_sorted_threecut, 10) # Cortisone (suppresses immune system, reduces inflammation), Sepiapterin, Sumatriptan (seratonin associated), Superoxide, Clopidogrel, Zolpidem, Cortisol, L-Lactic acid (Warburg effect), N1-Acetylspermidine, Dimethylallylpyrophosphate

# Normalize structural scores to be between 0 and 1
metabolites_struct_scored_sorted_threecut$normalized_structural_score <- ((metabolites_struct_scored_sorted_threecut$structural_score - min(metabolites_struct_scored_sorted_threecut$structural_score)) / (max(metabolites_struct_scored_sorted_threecut$structural_score) - min(metabolites_struct_scored_sorted_threecut$structural_score)))
metabolites_struct_scored_sorted$normalized_structural_score <- ((metabolites_struct_scored_sorted$structural_score - min(metabolites_struct_scored_sorted$structural_score)) / (max(metabolites_struct_scored_sorted$structural_score) - min(metabolites_struct_scored_sorted$structural_score)))
head(metabolites_struct_scored_sorted_threecut, 10)
head(metabolites_struct_scored_sorted, 10)

# Copy results into new data frames
structural_results <- metabolites_struct_scored_sorted_threecut # only metabolites with a minimum of three associated genes
structural_results_full <- metabolites_struct_scored_sorted # all metabolites

# Change column names
colnames(structural_results) <- c("Metabolite", "KEGG_ID", "Gene_Count", "Structural_Score", "Normalized_Score")
colnames(structural_results_full) <- c("Metabolite", "KEGG_ID", "Gene_Count", "Structural_Score", "Normalized_Score")

# Generate a histogram of the structural scores
hist(structural_results$Normalized_Score, xlab = 'Normalized Structural Subscore', main = 'Distribution of Normalized Structural Subscores')
hist(structural_results_full$Normalized_Score, xlab = 'Normalized Structural Subscore', main = 'Distribution of Normalized Structural Subscores')

# Write out data for long-term storage
write.csv(structural_results, file = "structural_results.csv")
write.csv(structural_results_full, file = "structural_results_full.csv")










#----------Compute DiffMut score----------------------------------------------------------------
#----------Compute DiffMut score----------------------------------------------------------------
#----------Compute DiffMut score----------------------------------------------------------------

# Import DiffMut data
diffmut <- read.csv("TCGA_BRCA_mutect.txt-DiffMut.txt", sep = ' ')
diffmut_sorted <- diffmut[with(diffmut, order(-uEMDscore)), ]
head(diffmut_sorted, 15)

nrow(diffmut) # 19,460 genes in TCGA-BRCA
hist(diffmut$uEMDscore)
diffmut[diffmut$qVal != 1,] # 39

# Graph to show why we shreshold with uEMD scores instead of q-values
qplot(diffmut$uEMDscore, diffmut$qVal) + scale_x_log10() + labs(title = "DiffMut Results: q-Values vs. Log of uEMD Scores", x = "log10(uEMD Score)", y = "q-Value")

nrow(diffmut[diffmut$uEMDscore == 0,]) # 5428 of 19,460      / 14,032
nrow(diffmut[diffmut$uEMDscore >= 0.15,]) # 1,178 / 19,460 = 0.0605

# Copy data
metabolites_diffmut <- metabolites_with_associated_genes

# Match mutational gene subscores to genes in metabolite list
metabolites_diffmut$uEMDscore <- diffmut$uEMDscore[match(metabolites_diffmut$Gene.Name, diffmut$protNames)]

nrow(metabolites_diffmut) # 19,522
nrow(metabolites_diffmut[metabolites_diffmut$uEMDscore >= 0.15,]) # 1,586

# Remove rows where gene has no diffmut score
metabolites_diffmut <- na.omit(metabolites_diffmut)
nrow(metabolites_diffmut) # 19,042

# Assign 1 to genes with DiffMut uEMD >= 0.15, or 0 otherwise
nrow(metabolites_diffmut[metabolites_diffmut$uEMDscore >= 0.15,]) # 1,106
metabolites_diffmut$uEMDthresholded <- metabolites_diffmut$uEMDscore
metabolites_diffmut$uEMDthresholded[metabolites_diffmut$uEMDthresholded >= 0.15] <- 1
metabolites_diffmut$uEMDthresholded[metabolites_diffmut$uEMDthresholded < 0.15] <- 0
head(metabolites_diffmut[metabolites_diffmut$uEMDthresholded == 1,], 30)

# Count number of genes associated with each metabolite
metabolites_diffmut$gene_count <- 1

# Aggregate total number of genes and total number of genes with significant uEMD scores associated with each metabolite
metabolites_diffmut_scored <- aggregate(uEMDthresholded~Metabolite, metabolites_diffmut, FUN=sum)
head(metabolites_diffmut_scored)
metabolites_diffmut_genecount <- aggregate(gene_count~Metabolite, metabolites_diffmut, FUN=sum)
head(metabolites_diffmut_genecount)
nrow(metabolites_diffmut_scored) # 1,118 metabolites
nrow(metabolites_diffmut_genecount) # 1,118 metabolites

# Compute per-metabolite DiffMut scores
metabolites_diffmut_scored$gene_count <- metabolites_diffmut_genecount$gene_count
metabolites_diffmut_scored$diffmut_score <- (metabolites_diffmut_scored$uEMDthresholded / metabolites_diffmut_scored$gene_count)

# Sort by DiffMut score
metabolites_diffmut_scores_sorted <- metabolites_diffmut_scored[with(metabolites_diffmut_scored, order(-diffmut_score)), ]
head(metabolites_diffmut_scores_sorted, 10)

# Reassociate KEGG values and clean up
metabolites_diffmut_scores_sorted$KEGG_ID <- metabolites_diffmut$KEGG[match(metabolites_diffmut_scores_sorted$Metabolite, metabolites_diffmut$Metabolite)]
head(metabolites_diffmut_scores_sorted, 10)
colnames(metabolites_diffmut_scores_sorted) <- c('metabolite', 'num_sig_genes', 'gene_count', 'diffmut_score', 'kegg_id')
metabolites_diffmut_scores_sorted <- metabolites_diffmut_scores_sorted[,c(1,5,2,3,4)]
head(metabolites_diffmut_scores_sorted, 10)

# Remove all metabolites with fewer than three associated genes
metabolites_diffmut_scores_sorted_threecut <- metabolites_diffmut_scores_sorted[metabolites_diffmut_scores_sorted$gene_count >= 3, ]
head(metabolites_diffmut_scores_sorted_threecut, 30)

# Compute normalized DiffMut scores
metabolites_diffmut_scores_sorted$normalized_diffmut_score <- ((metabolites_diffmut_scores_sorted$diffmut_score - min(metabolites_diffmut_scores_sorted$diffmut_score)) / (max(metabolites_diffmut_scores_sorted$diffmut_score) - min(metabolites_diffmut_scores_sorted$diffmut_score)))
head(metabolites_diffmut_scores_sorted, 10)
metabolites_diffmut_scores_sorted_threecut$normalized_diffmut_score <- ((metabolites_diffmut_scores_sorted_threecut$diffmut_score - min(metabolites_diffmut_scores_sorted_threecut$diffmut_score)) / (max(metabolites_diffmut_scores_sorted_threecut$diffmut_score) - min(metabolites_diffmut_scores_sorted_threecut$diffmut_score)))
head(metabolites_diffmut_scores_sorted_threecut, 30)

colnames(metabolites_diffmut_scores_sorted) <- c("Metabolite", "KEGG_ID", "Num_Sig_Genes", "Gene_Count", "DiffMut_Score", "Normalized_Score")
colnames(metabolites_diffmut_scores_sorted_threecut) <- c("Metabolite", "KEGG_ID", "Num_Sig_Genes", "Gene_Count", "DiffMut_Score", "Normalized_Score")

# Final data structures
diffmut_results <- metabolites_diffmut_scores_sorted_threecut
head(diffmut_results, 10) # Carbamoyl phosphate, Prostaglandin I2, Cortisol all known to be cancer associated
diffmut_results_full <- metabolites_diffmut_scores_sorted
head(diffmut_results_full, 10)

# Generate a histogram of the structural scores
hist(diffmut_results$Normalized_Score, xlab = 'Normalized DiffMut Subscore', main = 'Distribution of Normalized DiffMut Subscores')
hist(diffmut_results_full$Normalized_Score, xlab = 'Normalized DiffMut Subscore', main = 'Distribution of Normalized DiffMut Subscores')

# Write out data for long-term storage
write.csv(diffmut_results, file = "diffmut_results.csv")
write.csv(diffmut_results_full, file = "diffmut_results_full.csv")










#----------Compute expressional score----------------------------------------------------------------
#----------Compute expressional score----------------------------------------------------------------
#----------Compute expressional score----------------------------------------------------------------

# Take expression data
df <- rnaseq_data

# Convert raw data into a matrix with one sample per row and one gene per column, where the value 
# at a row-column intersection is the expression value of the column's gene in the row's sample
matrix <- acast(df, raw_data_accession ~ gene_id, value.var="normalized_expression_value")

row_names <- rownames(matrix)
abbrev_row_names <- data.frame(strsplit(row_names, "-"))[4,]
abbrev_row_names <- t(abbrev_row_names)
abbrev_row_names <- as.character(as.vector(abbrev_row_names[,1]))


# 0 indicates a cancer sample (row), 1 indicates a healthy sample (row)
bin_row_names <- as.numeric(substring(abbrev_row_names, 1, 1))
matrix_copy <- matrix
matrix_copy <- cbind(matrix_copy, bin_row_names)

# 1. Split matrix into cancer and healthy submatrices
cancer_matrix <- matrix_copy[matrix_copy[,17815] == 0, c(1:17814)]
healthy_matrix <- matrix_copy[matrix_copy[,17815] == 1, c(1:17814)]


# 2. Run t-test between each corresponding pair of rows in the two submatrices

# Define a function to run your t.test, grab the p value, and put them in a data.frame
f <- function(x,y){
  test <- t.test(x, y)
  #out <- data.frame(pval = test$p.value)
  out <- test$p.value
  return(out)
}

# Iterate over columns in cancer_matrix and healthy_matrix via sapply
p_vals <- sapply(seq(ncol(cancer_matrix)), function(x) f(cancer_matrix[,x], healthy_matrix[,x]))

# Adjust p-values for multiple hypothesis testing
p_vals <- p.adjust(p_vals, method = "bonferroni")

# Binarize p-values according to the threshold of 0.01
binarized_p_vals <- as.integer(p_vals < 0.01) #######

# Create dataframe corresponding binarized p-values to genes
bpv <- data.frame(colnames(matrix), binarized_p_vals)
colnames(bpv) <- c("gene_symbol", "binarized_p_val")

# Histogram of p_vals
hist(p_vals, xlab = "Bonferroni-Corrected Per-Gene P-Values", main = "Distribution of Expressional Gene P-Values")
hist(binarized_p_vals)

# Sort p-values for examination
pv <- data.frame(colnames(matrix), p_vals)
colnames(pv) <- c("gene_symbol", "p_val")
pv_sorted <- pv[with(pv, order(-p_val)), ]
head(pv_sorted, 50)

# Import metabolite data
metabolites_exp <- metabolites_with_associated_genes
nrow(metabolites_exp)  # 19,522 total associated genes across all metabolites

# Append binarized p-value of each gene in metabolites to metabolites
metabolites_exp$binarized_p_val <- bpv$binarized_p_val[match(metabolites_exp$Gene.Name, bpv$gene_symbol)]

# Remove rows with no binarized p-value
metabolites_exp <- na.omit(metabolites_exp)
nrow(metabolites_exp)  # 17,422 associated genes across all metabolites that show up in RNAseq data

# Append count of genes per metabolite
metabolites_exp$genecount <- 1
metabolites_exp_genecounts <- aggregate(genecount~Metabolite, metabolites_exp, FUN=sum)
mean(metabolites_exp_genecounts$genecount) # Average of 15.7951 genes associated with each metabolite
head(metabolites_exp_genecounts)

# Aggregate total number of genes and total number of genes with significant uEMD scores associated with each metabolite
metabolites_exp_scored <- aggregate(binarized_p_val~Metabolite, metabolites_exp, FUN=sum)
head(metabolites_exp_scored)
head(metabolites_exp_genecounts)
nrow(metabolites_exp_scored) # 1,103 metabolites
nrow(metabolites_exp_genecounts) # 1,103 metabolites

# Compute per-metabolite expression scores
metabolites_exp_scored$gene_count <- metabolites_exp_genecounts$genecount
metabolites_exp_scored$expressional_score <- (metabolites_exp_scored$binarized_p_val / metabolites_exp_scored$gene_count)
head(metabolites_exp_scored)

# Sort by expressional score
metabolites_exp_scores_sorted <- metabolites_exp_scored[with(metabolites_exp_scored, order(-expressional_score)), ]
head(metabolites_exp_scores_sorted, 165) # There are 164 / 1,103 metabolites with perfect scores: this is why trimming at 3 associated genes is necessary

# Reassociate KEGG values and clean up
metabolites_exp_scores_sorted$KEGG_ID <- metabolites_exp$KEGG[match(metabolites_exp_scores_sorted$Metabolite, metabolites_exp$Metabolite)]
head(metabolites_exp_scores_sorted, 10)
colnames(metabolites_exp_scores_sorted) <- c('metabolite', 'num_sig_genes', 'gene_count', 'expressional_score', 'kegg_id')
metabolites_exp_scores_sorted <- metabolites_exp_scores_sorted[,c(1,5,2,3,4)]
head(metabolites_exp_scores_sorted, 10)

# Remove all metabolites with fewer than three associated genes
metabolites_exp_scores_sorted_threecut <- metabolites_exp_scores_sorted[metabolites_exp_scores_sorted$gene_count >= 3, ]
head(metabolites_exp_scores_sorted_threecut, 30) # 7-Dehydrocholesterol (6/6) already breast cancer associated, Lathosterol (5/5) also cholesterol-linked, 4-Hydroxyproline (10/11) important metabolic regulator already linked to cancer, N4-Acetylaminobutanal (7/8), 5-Methylcytosine (4/4), (S)-3-Hydroxyisobutyric acid (4/4), D-Fructose 2,6-bisphosphate (4/4), 
# 7-Dehydrocholesterol helps synthesize vitamin D, and is found in human breast milk (6 out of 6 associated metabolites significantly overexpressed)
# 4-Hydroxyproline (10 out of 11) plays a key role in collagen stability, and was completely absorbed by the synonym issue from last year

# Compute normalized expressional scores
metabolites_exp_scores_sorted$normalized_expressional_score <- ((metabolites_exp_scores_sorted$expressional_score - min(metabolites_exp_scores_sorted$expressional_score)) / (max(metabolites_exp_scores_sorted$expressional_score) - min(metabolites_exp_scores_sorted$expressional_score)))
head(metabolites_exp_scores_sorted, 10)
metabolites_exp_scores_sorted_threecut$normalized_expressional_score <- ((metabolites_exp_scores_sorted_threecut$expressional_score - min(metabolites_exp_scores_sorted_threecut$expressional_score)) / (max(metabolites_exp_scores_sorted_threecut$expressional_score) - min(metabolites_exp_scores_sorted_threecut$expressional_score)))
head(metabolites_exp_scores_sorted_threecut, 25)

# Rename columns
colnames(metabolites_exp_scores_sorted) <- c("Metabolite", "KEGG_ID", "Num_Sig_Genes", "Gene_Count", "Expression_Score", "Normalized_Score")
colnames(metabolites_exp_scores_sorted_threecut) <- c("Metabolite", "KEGG_ID", "Num_Sig_Genes", "Gene_Count", "Expression_Score", "Normalized_Score")

# Final data structures
expressional_results <- metabolites_exp_scores_sorted_threecut
head(expressional_results, 30)
expressional_results_full <- metabolites_exp_scores_sorted
head(expressional_results_full, 30)

# Generate a histogram of the structural scores
hist(expressional_results$Normalized_Score, xlab = 'Normalized Expressional Subscore', main = 'Distribution of Normalized Expressional Subscores')
hist(expressional_results_full$Normalized_Score, xlab = 'Normalized Expressional Subscore', main = 'Distribution of Normalized Expressional Subscores')

# Write out data for long-term storage
write.csv(expressional_results, file = "expressional_results.csv")
write.csv(expressional_results_full, file = "expressional_results_full.csv")










#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------total_score = (mutation_score*alpha) + (expression_score*(1-alpha))---------------------------------------------------

nrow(mutational_results) # 688 for TCGA-BRCA
nrow(structural_results) # 688 for TCGA-BRCA
nrow(diffmut_results) # 749 for TCGA-BRCA
nrow(expressional_results) # 727 for TCGA-BRCA

# Create data frames containing only scores from metabolites that are present in all four different score data frames
common_metabolites <- intersect(intersect(intersect(mutational_results$Metabolite, structural_results$Metabolite), expressional_results$Metabolite), diffmut_results$Metabolite)
length(common_metabolites) # 673 common metabolites among the four data frames

mutational_results_common <- mutational_results[mutational_results$Metabolite %in% common_metabolites, ]
structural_results_common <- structural_results[structural_results$Metabolite %in% common_metabolites, ]
diffmut_results_common <- diffmut_results[diffmut_results$Metabolite %in% common_metabolites, ]
expressional_results_common <- expressional_results[expressional_results$Metabolite %in% common_metabolites, ]

# Sort common data frames by metabolite name
mutational_results_common <- mutational_results_common[with(mutational_results_common, order(Metabolite)), ]
structural_results_common <- structural_results_common[with(structural_results_common, order(Metabolite)), ]
diffmut_results_common <- diffmut_results_common[with(diffmut_results_common, order(Metabolite)), ]
expressional_results_common <- expressional_results_common[with(expressional_results_common, order(Metabolite)), ]

# Change name of score to be "Type_Score" for common reference in function below
names(mutational_results_common)[names(mutational_results_common) == 'Mutation_Score'] <- 'Type_Score'
names(structural_results_common)[names(structural_results_common) == 'Structural_Score'] <- 'Type_Score'
names(diffmut_results_common)[names(diffmut_results_common) == 'DiffMut_Score'] <- 'Type_Score'
names(expressional_results_common)[names(expressional_results_common) == 'Expression_Score'] <- 'Type_Score'


# Function to compute the appropriate matrix between any two data sets
create_matrix <- function(data_1, data_2){
  
  # Filter out metabolites that don't have KEGG ID's in any KEGG human metabolism pathways
  data_1 <- data_1[data_1$KEGG_ID %in% met_pool, ]
  data_2 <- data_2[data_2$KEGG_ID %in% met_pool, ]
  
  # Remove last few rows of mutation and expression tables to make number of metabolites a multiple of 5 (e.g. 2,799 rows -> 2,795 rows)
  # This is for breaking up the data into evenly sized training and holdout sets (see below)
  modulo <- nrow(data_1) %% 5

  if (modulo != 0) {
    data_1 <- head(data_1, -modulo)
    data_2 <- head(data_2, -modulo)
  }
  
  # Re-normalize scores
  data_1$Normalized_Score <- ((data_1$Type_Score - min(data_1$Type_Score)) / (max(data_1$Type_Score) - min(data_1$Type_Score)))
  data_2$Normalized_Score <- ((data_2$Type_Score - min(data_2$Type_Score)) / (max(data_2$Type_Score) - min(data_2$Type_Score)))
  
  # Create sequence of all alphas from 0 to 1 by 0.01
  alphas <- seq(from = 0, to = 1, by = 0.01)
  
  # Create matrix for storing overall scores at all alpha values
  total_scores <- matrix(data = NA, nrow = nrow(data_1), ncol = length(alphas))
  rownames(total_scores) <- data_1$KEGG_ID
  colnames(total_scores) <- alphas
  
  # Calculate overall scores based on (mutation_score*alpha) + (expression_score*(1-alpha)) 
  # for all alphas 0-1 by 0.01 using normalized scores
  for(column in 1:ncol(total_scores)){
    alpha <- as.numeric(colnames(total_scores)[column])
    total_scores[, column] <- ((data_1$Normalized_Score * alpha) 
                               + (data_2$Normalized_Score * (1 - alpha)))
  }
  
  return(total_scores)
}

# Create matrices showing the combined scores of each of the six unique score pairings at all alpha values from 0 to 1 by 0.01 
# NOTE: the first score in each matrix's name is the one multiplied by alpha, and the second is the one multiplied by (1 - alpa)
total_scores_mutational_structural <- create_matrix(mutational_results_common, structural_results_common)
total_scores_mutational_expressional <- create_matrix(mutational_results_common, expressional_results_common)
total_scores_mutational_diffmut <- create_matrix(mutational_results_common, diffmut_results_common)
total_scores_structural_expressional <- create_matrix(structural_results_common, expressional_results_common)
total_scores_structural_diffmut <- create_matrix(structural_results_common, diffmut_results_common)
total_scores_expressional_diffmut <- create_matrix(expressional_results_common, diffmut_results_common)










#----------Perform machine learning to identify best alpha value----------------------------------------------
#----------Perform machine learning to identify best alpha value----------------------------------------------
#----------Perform machine learning to identify best alpha value----------------------------------------------

find_optimal_alpha_and_auc <- function(total_scores){
  
  # Create sequence of all alphas from 0 to 1 by 0.01
  alphas <- seq(from = 0, to = 1, by = 0.01)
  
  # Initialize list to store, for each metabolite, a 1 if the metabolite 
  # is present in a KEGG cancer-linked pathway, or a 0 if not
  in_cancer_pool <- c()
  
  # Iterate through all metabolites in the overall score matrix
  for (i in 1:length(rownames(total_scores))) {
    if (rownames(total_scores)[i] %in% cancer_pool) {
      # Mark a 1 in the list if the metabolite appears in a KEGG cancer-linked pathway
      in_cancer_pool = c(in_cancer_pool, 1)
    } else {
      # Mark a 0 in the list if the metabolite does not appear in a KEGG cancer-linked pathway
      in_cancer_pool = c(in_cancer_pool, 0)
    }
  }
  
  # Add the list of 1's and 0's for each metabolite to the overall score matrix
  total_scores <- cbind(total_scores, in_cancer_pool)
  
  
  # Randomly permute data row-wise
  total_scores_random <- total_scores[sample(nrow(total_scores)),]
  
  # Split random data into 5 training sets and 5 hold out sets
  set_size <- (nrow(total_scores_random) / 5)
  
  set1 <- total_scores_random[(set_size+1):nrow(total_scores_random), ]
  set1_holdout <- total_scores_random[1:set_size, ]
  
  set2 <- total_scores_random[1:set_size, ]
  set2 <- rbind(set2, total_scores_random[((2*set_size)+1):nrow(total_scores_random), ])
  set2_holdout <- total_scores_random[(set_size+1):(2*set_size), ]
  
  set3 <- total_scores_random[1:(2*set_size), ]
  set3 <- rbind(set3, total_scores_random[((3*set_size)+1):nrow(total_scores_random), ])
  set3_holdout <- total_scores_random[((2*set_size)+1):(3*set_size), ]
  
  set4 <- total_scores_random[1:(3*set_size), ]
  set4 <- rbind(set4, total_scores_random[((4*set_size)+1):nrow(total_scores_random), ])
  set4_holdout <- total_scores_random[((3*set_size)+1):(4*set_size), ]
  
  set5 <- total_scores_random[1:(4*set_size), ]
  set5_holdout <- total_scores_random[((4*set_size)+1):nrow(total_scores_random), ]
  
  sets <- list(set1, set2, set3, set4, set5)
  holdout_sets <- list(set1_holdout, set2_holdout, set3_holdout, set4_holdout, set5_holdout)
  
  # Run machine learning simulation on all training set-holdout set pairs, 
  # computing a maximum alpha and corresponding maximum holdout AUC for each
  alpha_maxes <- c()
  final_aucs <- c()
  holdout_perfs <- c()
  i = 0
  for (set in sets) {
    
    i = i + 1
    
    auc_vals <- matrix(data = NA, nrow = 1, ncol = length(alphas))
    rownames(auc_vals) <- c("AUC_Values")
    colnames(auc_vals) <- alphas
    
    for(column in 1:(ncol(total_scores_random)-1)) {
      prediction.obj <- prediction(set[, column], set[, "in_cancer_pool"])
      perf <- performance(prediction.obj, measure = 'auc')
      auc <- perf@y.values[[1]]
      auc_vals[1, column] <- auc
    }
    
    alpha_max <- colnames(auc_vals)[apply(auc_vals,1,which.max)]
    alpha_maxes <- c(alpha_maxes, as.numeric(alpha_max))
    
    holdout_pred <- prediction(holdout_sets[[i]][, alpha_max], holdout_sets[[i]][, "in_cancer_pool"])
    holdout_perf <- performance(holdout_pred, measure = 'auc')
    holdout_auc <- holdout_perf@y.values[[1]]
    
    # Plot ROC curve
    graph_perf <- performance(holdout_pred,"tpr","fpr")
    plot(graph_perf)
    
    # Store values
    holdout_perfs <- c(holdout_perfs, holdout_perf)
    final_aucs <- c(final_aucs, holdout_auc)
  }
  
  # Set final alphas to be those found in my writeup (due to the random nature 
  # of the set divisions, alpha values will vary each time the code is run)
  # final_alpha <- mean(c(1.0, 1.0, 0.13, 1.0, 1.0))
  
  # Compute average of alpha maxes (i.e. the optimal alpha)
  final_alpha <- mean(alpha_maxes)
  alpha_maxes
  final_alpha # I got 1, indicating a complete favoring of the mutation score (no weight at all is assigned to the expression score)
  
  # Compute the average of the final AUC values
  final_auc <- mean(final_aucs)
  final_aucs
  final_auc
  
  results <- c(final_alpha, final_auc, final_aucs)
  
  return(results)
}


# Compute final alpha and final AUC value for each scoring pairing 
# (alpha values closer to 1 favor first score, whereas alpha values closer to zero favors second score)
alpha_and_auc_mutational_structural <- find_optimal_alpha_and_auc(total_scores_mutational_structural)
alpha_and_auc_mutational_expressional <- find_optimal_alpha_and_auc(total_scores_mutational_expressional)
alpha_and_auc_mutational_diffmut <- find_optimal_alpha_and_auc(total_scores_mutational_diffmut)
alpha_and_auc_structural_expressional <- find_optimal_alpha_and_auc(total_scores_structural_expressional)
alpha_and_auc_structural_diffmut <- find_optimal_alpha_and_auc(total_scores_structural_diffmut)
alpha_and_auc_expressional_diffmut <- find_optimal_alpha_and_auc(total_scores_expressional_diffmut)

# Print results
alpha_and_auc_mutational_structural[[1]]
alpha_and_auc_mutational_structural[[2]]
alpha_and_auc_mutational_expressional[[1]]
alpha_and_auc_mutational_expressional[[2]]
alpha_and_auc_mutational_diffmut[[1]]
alpha_and_auc_mutational_diffmut[[2]]
alpha_and_auc_structural_expressional[[1]]
alpha_and_auc_structural_expressional[[2]]
alpha_and_auc_structural_diffmut[[1]]
alpha_and_auc_structural_diffmut[[2]]
alpha_and_auc_expressional_diffmut[[1]]
alpha_and_auc_expressional_diffmut[[2]]

# Get box-and-whiskers plots for final AUC values
boxplot_df <- data.frame("Structural_Expressional" = tail(alpha_and_auc_structural_expressional, 5), "Mutational_Structural" = tail(alpha_and_auc_mutational_structural, 5), "Mutational_Expressional" = tail(alpha_and_auc_mutational_expressional, 5), "Mutational_DiffMut" = tail(alpha_and_auc_mutational_diffmut, 5), "Expressional_DiffMut" = tail(alpha_and_auc_expressional_diffmut, 5), "Structural_DiffMut" = tail(alpha_and_auc_structural_diffmut, 5))
auc_boxplots <- ggplot(stack(boxplot_df), aes(x = factor(ind, levels = names(boxplot_df)), y = values)) + geom_boxplot() + xlab("Pairwise Subscore Aggregation") + ylab("Final AUC Value") + ggtitle("Pairwise Final AUC Value Boxplots")
auc_boxplots










#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------total_score = (mutation_score*alpha) + (expression_score*(1-alpha))---------------------------------------------------

# Set "best" variables (will be changed below if structural-DiffMut is not the best pairwise aggregation)
best_alpha_and_auc <- alpha_and_auc_structural_diffmut
best_results_1 <- structural_results
best_results_2 <- diffmut_results
score1_score2 <- "structural_diffmut"
colnames(best_results_1)[colnames(best_results_1)=="Structural_Score"] <- "Type_Score"
colnames(best_results_2)[colnames(best_results_2)=="DiffMut_Score"] <- "Type_Score"

# Identify best pairwise results and assign "best" variables appropriately
if (alpha_and_auc_mutational_structural[[2]] == max(alpha_and_auc_mutational_structural[[2]], alpha_and_auc_mutational_expressional[[2]], alpha_and_auc_mutational_diffmut[[2]], alpha_and_auc_structural_expressional[[2]], alpha_and_auc_structural_diffmut[[2]], alpha_and_auc_expressional_diffmut[[2]])) {
  best_alpha_and_auc <- alpha_and_auc_mutational_structural
  best_results_1 <- mutational_results
  best_results_2 <- structural_results
  score1_score2 <- "mutational_structural"
  colnames(best_results_1)[colnames(best_results_1)=="Mutation_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="Structural_Score"] <- "Type_Score"
} else if (alpha_and_auc_mutational_expressional[[2]] == max(alpha_and_auc_mutational_structural[[2]], alpha_and_auc_mutational_expressional[[2]], alpha_and_auc_mutational_diffmut[[2]], alpha_and_auc_structural_expressional[[2]], alpha_and_auc_structural_diffmut[[2]], alpha_and_auc_expressional_diffmut[[2]])) {
  best_alpha_and_auc <- alpha_and_auc_mutational_expressional
  best_results_1 <- mutational_results
  best_results_2 <- expressional_results
  score1_score2 <- "mutational_expressional"
  colnames(best_results_1)[colnames(best_results_1)=="Mutation_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="Expression_Score"] <- "Type_Score"
} else if (alpha_and_auc_mutational_diffmut[[2]] == max(alpha_and_auc_mutational_structural[[2]], alpha_and_auc_mutational_expressional[[2]], alpha_and_auc_mutational_diffmut[[2]], alpha_and_auc_structural_expressional[[2]], alpha_and_auc_structural_diffmut[[2]], alpha_and_auc_expressional_diffmut[[2]])) {
  best_alpha_and_auc <- alpha_and_auc_mutational_diffmut
  best_results_1 <- mutational_results
  best_results_2 <- diffmut_results
  score1_score2 <- "mutational_diffmut"
  colnames(best_results_1)[colnames(best_results_1)=="Mutation_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="DiffMut_Score"] <- "Type_Score"
} else if (alpha_and_auc_structural_expressional[[2]] == max(alpha_and_auc_mutational_structural[[2]], alpha_and_auc_mutational_expressional[[2]], alpha_and_auc_mutational_diffmut[[2]], alpha_and_auc_structural_expressional[[2]], alpha_and_auc_structural_diffmut[[2]], alpha_and_auc_expressional_diffmut[[2]])) {
  best_alpha_and_auc <- alpha_and_auc_structural_expressional
  best_results_1 <- structural_results
  best_results_2 <- expressional_results
  score1_score2 <- "structural_expressional"
  colnames(best_results_1)[colnames(best_results_1)=="Structural_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="Expression_Score"] <- "Type_Score"
} else if (alpha_and_auc_structural_diffmut[[2]] == max(alpha_and_auc_mutational_structural[[2]], alpha_and_auc_mutational_expressional[[2]], alpha_and_auc_mutational_diffmut[[2]], alpha_and_auc_structural_expressional[[2]], alpha_and_auc_structural_diffmut[[2]], alpha_and_auc_expressional_diffmut[[2]])) {
  best_alpha_and_auc <- alpha_and_auc_structural_diffmut
  best_results_1 <- structural_results
  best_results_2 <- diffmut_results
  score1_score2 <- "structural_diffmut"
  colnames(best_results_1)[colnames(best_results_1)=="Structural_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="DiffMut_Score"] <- "Type_Score"
} else {
  best_alpha_and_auc <- alpha_and_auc_expressional_diffmut
  best_results_1 <- expressional_results
  best_results_2 <- diffmut_results
  score1_score2 <- "expressional_diffmut"
  colnames(best_results_1)[colnames(best_results_1)=="Expression_Score"] <- "Type_Score"
  colnames(best_results_2)[colnames(best_results_2)=="DiffMut_Score"] <- "Type_Score"
}

# Print best pairing and corresponding optimal alpha and AUC values
score1_score2
best_alpha_and_auc[[1]]
best_alpha_and_auc[[2]]

# Optimal scores: Mutational and DiffMut alpha = 0.382; AUC = 0.75
final_alpha <- best_alpha_and_auc[[1]]

# Create data frames containing only scores from metabolites that are present in all four different score data frames
common_metabolites_score1_score2 <- intersect(best_results_1$Metabolite, best_results_2$Metabolite)
length(common_metabolites_score1_score2) 

score1_results_final <- best_results_1[best_results_1$Metabolite %in% common_metabolites_score1_score2, ]
score2_results_final <- best_results_2[best_results_2$Metabolite %in% common_metabolites_score1_score2, ]

# Alphabetize final data frames by metabolite name
score1_results_final <- score1_results_final[with(score1_results_final, order(Metabolite)), ]
score2_results_final <- score2_results_final[with(score2_results_final, order(Metabolite)), ]
all(score1_results_final$Metabolite == score2_results_final$Metabolite)

# Renormalize scores
score1_results_final$Normalized_Score <- ((score1_results_final$Type_Score - min(score1_results_final$Type_Score)) / (max(score1_results_final$Type_Score) - min(score1_results_final$Type_Score)))
score2_results_final$Normalized_Score <- ((score2_results_final$Type_Score - min(score2_results_final$Type_Score)) / (max(score2_results_final$Type_Score) - min(score2_results_final$Type_Score)))

# Raw (unweighted) sum of normalized mutation and expression scores
total_results <- data.frame(Metabolite = score1_results_final$Metabolite, 
                            KEGG_ID = score1_results_final$KEGG_ID, 
                            Normalized_Type_Score_1 = score1_results_final$Normalized_Score, 
                            Normalized_Type_Score_2 = score2_results_final$Normalized_Score, 
                            Total_Score = (((score1_results_final$Normalized_Score)*final_alpha) + ((score2_results_final$Normalized_Score)*(1 - final_alpha))))

# Compute normalized total score
total_results$Normalized_Total_Score <- ((total_results$Total_Score - min(total_results$Total_Score)) / (max(total_results$Total_Score) - min(total_results$Total_Score)))

# Sort by total scores
total_results_sorted <- total_results[with(total_results, order(-Normalized_Total_Score)), ]

# Print the 20 metabolites with highest normalized 
# overall scores at the optimal alpha value
head(total_results_sorted, 20)

# Generate a histogram of normalized total scores, computed using the 
# best alpha value as computed by ROC machine learning techniques
hist(total_results$Normalized_Total_Score, xlab = 'Normalized Overall Score', main = 'Distribution of Normalized Overall Scores Between Structural and DiffMut Subscores at Alpha=0.305', cex.main=0.7)

# Write overall scores to a CSV file
write.csv(total_results, file = "overall_scores.csv")
