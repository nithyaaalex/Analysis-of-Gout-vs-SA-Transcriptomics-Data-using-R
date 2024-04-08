#importing tables from folder
annotations = read.table("D:/Lab/Annotations.csv", header=TRUE, row.names=1, sep="\t")
de_gout_vs_hc = read.table("D:/Lab/DE_Gout_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
de_sa_vs_hc = read.table("D:/Lab/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
expression_table = read.table("D:/Lab/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")
sample_information = read.table("D:/Lab/Sample_Information.csv", header=TRUE, row.names=1, sep="\t")

#run this every time before starting
library(ggplot2)
library(ggfortify)

#making annotated tables for using later
expressions_annotated = merge(annotations, expression_table, by.x=0, by.y=0)
row.names(expressions_annotated) = expressions_annotated[,1]
names(expressions_annotated)[1] = "GENE_ID"

#experimenting with some graphs on the annotated expression table
ggp = ggplot(expressions_annotated, aes(x=log10(Healthy_1))) + geom_histogram()+
  labs(x="Log10 of sample Healthy_1", y="Expression Levels of Genes present", title="Histogram of sample Healthy_1 gene expression")
ggp

ggp = ggplot(expressions_annotated, aes(x=log10(Gout_1))) + geom_histogram()+
  labs(x="Log10 of sample Gout_1", y="Expression Levels of Genes present", title="Histogram of sample Gout_1 gene expression")
ggp

ggp = ggplot(expressions_annotated, aes(x=log10(SA_1))) + geom_histogram()+
  labs(x="Log10 of sample SA_1", y="Expression Levels of Genes present", title="Histogram of sample SA_1 gene expression")
ggp

#make subset tables of sa, gout and healthy 
samples_healthy = row.names(subset(sample_information, GROUP == "HEALTHY"))
samples_gout = row.names(subset(sample_information, GROUP == "GOUT"))
samples_sa = row.names(subset(sample_information, GROUP == "SA"))

expression_healthy = expressions_annotated[,samples_healthy]
expression_gout = expressions_annotated[,samples_gout]
expression_sa = expressions_annotated[,samples_sa]

#for figuring out the demographic of our data
summary(sample_information)
str(sample_information)
sample_information$GROUP=as.factor(sample_information$GROUP)
sample_information$AGE= as.numeric(sample_information$AGE)
sample_information$SEX = as.factor(sample_information$SEX)

samp_hc = sample_information[samples_healthy,]
samp_sa = sample_information[samples_sa,]
samp_gout = sample_information[samples_gout,]
summary(samp_hc)
summary(samp_sa)
summary(samp_gout)


#for calc pvalues between groups at each clinical measurement
sampleinfo_sa_vs_hc = sample_information[c(samples_sa,samples_healthy),]
sampleinfo_gout_vs_hc = sample_information[c(samples_gout,samples_healthy),]
sampleinfo_sa_vs_gout = sample_information[c(samples_sa,samples_gout),]
sampleinfo_sa_vs_hc$GROUP = as.factor(sampleinfo_sa_vs_hc$GROUP)
sampleinfo_gout_vs_hc$GROUP = as.factor(sampleinfo_gout_vs_hc$GROUP)
sampleinfo_sa_vs_gout$GROUP = as.factor(sampleinfo_sa_vs_gout$GROUP)

t.test(sampleinfo_sa_vs_hc$NEUTROPHILS~sampleinfo_sa_vs_hc$GROUP)
t.test(sampleinfo_sa_vs_hc$MONOCYTES~sampleinfo_sa_vs_hc$GROUP)
t.test(sampleinfo_sa_vs_hc$AGE~sampleinfo_sa_vs_hc$GROUP)

t.test(sampleinfo_gout_vs_hc$NEUTROPHILS~sampleinfo_gout_vs_hc$GROUP)
t.test(sampleinfo_gout_vs_hc$MONOCYTES~sampleinfo_gout_vs_hc$GROUP)
t.test(sampleinfo_gout_vs_hc$AGE~sampleinfo_gout_vs_hc$GROUP)

t.test(sampleinfo_sa_vs_gout$NEUTROPHILS~sampleinfo_sa_vs_gout$GROUP)
t.test(sampleinfo_sa_vs_gout$MONOCYTES~sampleinfo_sa_vs_gout$GROUP)
t.test(sampleinfo_sa_vs_gout$AGE~sampleinfo_sa_vs_gout$GROUP)

#didn't analyse sex above because it's already very skewed - any t.test wouldn't be fair to conduct 

#histograms of clinical information across all samples - for my reference only - not putting in report
ggp = ggplot(sample_information, aes(x=AGE)) + geom_histogram()+
  labs(x="Age", y = "Samples" , title="Histogram of variation of age across samples")
ggp
ggp = ggplot(sample_information, aes(x=NEUTROPHILS)) + geom_histogram()+
  labs(x="Neutrophils", y = "Samples" , title="Histogram of variation of neutrophils across samples")
ggp
ggp = ggplot(sample_information, aes(x=MONOCYTES)) + geom_histogram()+
  labs(x="Monocytes", y = "Samples" , title="Histogram of variation of monocytes across samples")
ggp

#boxplots of clinical information, with groups as factor - put age and neutrophils in report - monocytes nothing to note so exclude
ggp = ggplot(sample_information, aes(x=GROUP,y=AGE, fill = GROUP)) + geom_boxplot()+
  labs(x="Age", y = "Samples" , title="Boxplot of variation of age across each group")
ggp
ggp = ggplot(sample_information, aes(x=GROUP,y=NEUTROPHILS, fill = GROUP)) + geom_boxplot()+
  labs(x="Neutrophils", y = "Samples" , title="Boxplot of variation of neutrophils across each group")
ggp
ggp = ggplot(sample_information, aes(x=GROUP,y=MONOCYTES, fill = GROUP)) + geom_boxplot()+
  labs(x="Monocytes", y = "Samples" , title="Boxplot of variation of monocytes across each group")
ggp

#make de table of sa and gout using expression subset tables
de_sa_vs_gout =  as.data.frame(matrix(0, ncol = 3, nrow = nrow(expressions_annotated)))
names(de_sa_vs_gout) = c("Log2fold","P", "p.adj")
row.names(de_sa_vs_gout) = row.names(expressions_annotated)

for (row in 1:nrow(expressions_annotated))
{
  gene_data_sa = as.numeric(expression_sa[row,])
  gene_data_gout = as.numeric(expression_gout[row,])
  
  mean_sa = mean(gene_data_sa)
  mean_gout = mean(gene_data_gout)
  log2fold = log2(mean_sa) - log2(mean_gout)

  p = t.test(gene_data_gout,gene_data_sa)
  p = p$p.value

  de_sa_vs_gout[row,"Log2fold"] = log2fold
  de_sa_vs_gout[row,"P"] = p
  
}

#adding adjusted p to gout vs sa
p_values= de_sa_vs_gout[,"P"]
#trying different methods for p.adjust
adjusted_p_values = p.adjust(p_values, method = "bonferroni", n = length(p_values)) #too stringent - all values become 1
adjusted_p_values = p.adjust(p_values, method = "BH", n = length(p_values)) #both BH and fdr give same values > 0.53 
adjusted_p_values = p.adjust(p_values, method = "fdr") #using FDR FINALLY
adjusted_p_values
de_sa_vs_gout$p.adj = adjusted_p_values


#checking if the padjust function is working properly - it is, but still getting very high values of p.adj in sa vs gout table
#p_values= de_sa_vs_hc[,"p"]
#adjusted_p_values = p.adjust(p_values, method = "FDR", n = length(p_values))
#adjusted_p_values
#de_sa_vs_hc$p.adjTRIAL = adjusted_p_values
#de_sa_vs_hc = de_sa_vs_hc[,c("log2fold","p", "p.adj")]


#making significant tables for hc gout and hc sa
de_sa_vs_hc_sig = subset(de_sa_vs_hc,p<0.05)
idvector = row.names(de_sa_vs_hc_sig)
samplevector = names(expression_table)
expression_geneid_sa_hc_sig= expressions_annotated[idvector, samplevector]

de_gout_vs_hc_sig = subset(de_gout_vs_hc,p<0.05)
idvector = row.names(de_gout_vs_hc_sig)
samplevector = names(expression_table)
expression_geneid_gout_hc_sig= expressions_annotated[idvector, samplevector]

#making significant genes only table - for gout and sa
de_sa_vs_gout_sig = subset(de_sa_vs_gout,P<0.01)
nrow(de_sa_vs_gout_sig)
idvector = row.names(de_sa_vs_gout_sig)
samplevector = names(expression_table)
expression_geneid_gout_sa_sig= expressions_annotated[idvector, samplevector]

#list of genes to add in report

#genes from hc and gout, table - sampleinfo_gout_vs_hc
#ENSG00000108950 #FAM20A
#ENSG00000170439 #METTL7B #going with this for report 

#genes from hc and sa, table - sampleinfo_sa_vs_hc
#ENSG00000149256 #TENM4
#ENSG00000137869 #CYP19A1 #going with this for report

# genes i picked out to analyse from de_sa_vs_hc_sig, table - sampleinfo_sa_vs_gout
#ENSG00000225364 - good p value and good log2fold change - most negative #ATP6V0E1P1
#ENSG00000262061 -  good p value and good log2fold change - most positive #AC129507.1
#ENSG00000267541 - best p value #MTCO2P2

#for comparing expression of single gene across all data - code to reuse
gene_data = expression_table["ENSG00000267541",]
transposed= t(gene_data)
gene_data= as.data.frame(transposed)
names(gene_data) <- 'gene'
ggp = ggplot(gene_data, aes(x=(gene))) + geom_histogram(colour="darkgreen",fill="blue",size=1, alpha = 0.5) +labs(x = "ENSG00000267541", title = "Expression of MTCO2P2")
ggp #add log10 if needed

#adding clinical info
gene_data$GROUP = sample_information$GROUP
gene_data$SEX = sample_information$SEX
gene_data$AGE = sample_information$AGE
gene_data$NEUTROPHILS = sample_information$NEUTROPHILS
gene_data$MONOCYTES = sample_information$MONOCYTES

#gene across three groups compared
summary(gene_data)
ggp = ggplot(gene_data, aes(x=GROUP,y=gene, fill = GROUP)) +
  geom_boxplot() + 
  labs(x="GROUP", y="Expression of MTCO2P2")
ggp

#code for comparing across only two groups
# for gout and hcgene_data = gene_data[1:28,]
# for sa and hc gene_data = gene_data[c(1:14,29:42),]
gene_data = gene_data[c(15:42),]
ggp = ggplot(gene_data, aes(x=GROUP,y=gene, fill = GROUP)) +
  geom_boxplot() + 
  labs(x="GROUP", y="Expression of MTCO2P2")
ggp


#gene compared across sex
summary(gene_data)
ggp = ggplot(gene_data, aes(x=SEX,y=gene, fill = SEX)) + 
  geom_boxplot() + 
  labs(x="SEX", y="Expression of MTCO2P2")
ggp

#gene compared to age
ggp = ggplot(gene_data, aes(x=AGE,y=gene)) + 
  geom_point() + 
  labs(x="Age", y="Expression of MTCO2P2")
ggp

#gene compared to neutrophils
ggp = ggplot(gene_data, aes(x=NEUTROPHILS,y=gene)) + 
  geom_point() + 
  labs(x="Neutrophils", y="Expression of MTCO2P2")
ggp

#gene compared to monocytes
ggp = ggplot(gene_data, aes(x=MONOCYTES,y=gene)) + 
  geom_point() + 
  labs(x="Monocytes", y="Expression of MTCO2P2")
ggp

#for checking correlation #pearson for linear, spearman for ordinal add ,method = "pearson" in the brackets for that 
gene_data$GROUP=as.numeric(gene_data$GROUP)
gene_data$SEX=as.numeric(gene_data$SEX)
cor(gene_data$GROUP,gene_data$gene, method = "pearson")
cor(gene_data$AGE,gene_data$gene, method ="pearson")
cor(gene_data$SEX,gene_data$gene, method ="pearson")
cor(gene_data$NEUTROPHILS,gene_data$gene, method = "pearson")
cor(gene_data$MONOCYTES,gene_data$gene, method = "pearson")

#creating a GLM with gene data and the clinical information
model1 = lm(gene_data$gene~gene_data$GROUP+gene_data$SEX+gene_data$AGE+gene_data$NEUTROPHILS+gene_data$MONOCYTES) #remove group and do pls
anova(model1)
summary(model1)
hist(rstandard(model1)) #checking residuals assumption
autoplot(model1)

#not using below graph outputs - use these models only for like 
#plotting the model - use this scatter plot for covariates
ggp = ggplot(gene_data, aes(x=NEUTROPHILS, y=gene)) +
  geom_point() +
  labs(x="NEUTROPHILS", y="Expression of MTCO2P2") +
  geom_smooth(method = "lm", se = FALSE)
ggp

#plotting the model - use violin for factors
ggp = ggplot(gene_data, aes(x=GROUP, y=gene, fill=GROUP)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  stat_summary(fun=mean, colour="red") +
  labs(x="GROUP", y="expression")
ggp

for (x in 15:28) #setting the factor level for sa from 3 to zero to allow binomial comparison
{
  gene_data[x,"GROUP"] = 0
}
#now sa is 0 and gout is 1

#if autoplot results are not good/ skewed distribution, plotting binomial/poisson/inverse.gaussian whatever fits best
ggp = ggplot(gene_data, aes(x=gene, y=GROUP)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)
ggp

#binomial fit is good hence, making binomial generalized linear model
model2 = glm(gene_data$GROUP~gene_data$gene, family = binomial)
summary(model2)
ratio = 23.054/26 #residual deviance by residual degrees of freedom (for over or under dispersion check)

#PCA wrt grouping of disease
#reload all tables after doing this - healthy samples removed fully
expression_table = expression_table[,c(15:42)]
expression_geneid_gout_sa_sig = expression_geneid_gout_sa_sig[,c(15:42)]
expression_table.nm = as.matrix(sapply(expression_geneid_gout_sa_sig, as.numeric))
PCA = prcomp(t(expression_table.nm))
pca_coordinates = data.frame(PCA$x)
sample_information = sample_information[c(15:42),]
ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = sample_information$GROUP)) +
  geom_point()
ggp

#decent grouping according to disease found using pca 

#below tests only for sa vs gout
#t-test for comparing mean gene expression between SA and GOUT
t.test(gene_data$gene~gene_data$GROUP)
gene_data_SA = subset (gene_data, GROUP == "SA")
hist(gene_data_SA$gene)
gene_data_gout = subset (gene_data, GROUP == "GOUT")
hist(gene_data_gout$gene)
#not interesting DISTRIBUTION IN BOTH SA AND GOUT

#wilcox test to see if p value is better incase data is ordinal
wilcox.test(gene_data$gene~gene_data$GROUP, alternative = "two.sided")

#hypergeometric test between sa and gout common genes 
de_merged = merge(de_sa_vs_hc, de_gout_vs_hc, by=0)
tempvariable <- de_merged[,-1]
rownames(tempvariable) <- de_merged[,1]
de_merged <- tempvariable
rm(tempvariable)
names(de_merged) = c("log2fold_sa_vs_hc","p_sa_vs_hc","padj_sa_vs_hc","log2fold_gout_vs_hc","p_gout_vs_hc","padj_gout_vs_hc")

sa = nrow(subset(de_merged, padj_sa_vs_hc < 0.05))
gout = nrow(subset(de_merged, padj_gout_vs_hc < 0.05))
overlap = nrow(subset(de_merged, padj_sa_vs_hc < 0.05 & padj_gout_vs_hc < 0.05))

sa_only = sa - overlap
gout_only = gout - overlap
total = nrow(de_merged)

phyper(overlap-1, sa, total-sa, gout,lower.tail= FALSE)
#it is <<0.001 so it the overlap in genes is way more than expected

