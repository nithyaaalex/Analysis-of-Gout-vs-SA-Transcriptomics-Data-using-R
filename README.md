Used statistical analysis in R to answer the question: Can Septic Arthritis and Gout be differentiated based on their transcriptional levels in blood?
Grade: B1

Input Required for code: 
1. Expression Table with sample names as column names and gene IDs as row names
2. Sample Metadata
3. Annotations table
4. Differential Tables between two conditions (generated for one comparison, 2 were given)

What the code does:
1. Draw histograms to compare expression profiles of each sample type (Healthy, Gout, SA)
2. Performs t-tests to see if any of the variables in the metadata has a significant impact
3. Plots variation of variables from metadata
4. Makes a differential expression table from scratch for Septic Arthritis vs Gout and adjusts for large no of tests
5. Find the significant genes from each comparison
6. Plots the expression of any chosen gene across three samples (and by the metadata)
7. Makes a generalised linear model for the gene and plots it
8. PCA analysis of the data
9. Hypergeometric test to see if SA and gout have more genes in common than expected 
