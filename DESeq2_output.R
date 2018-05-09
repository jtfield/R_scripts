# Written by Jasper Toscani Field and Ethan Low
# 4/15/18
# takes two bam files from HISAT2 and runs differential expression using DESEQ2
# produces list of up and downregulated genes

library('readxl')
library('Rsubread')
library('DESeq2')
library(gplots)

setwd('/Users/Cristie/Desktop/School/Merced/Class/RNASeq Course- Spring 2018/Analysis')

#counts genomic features on each count file
features<-featureCounts(c('Galaxy62-[HISAT2_on_WT_new].bam', 'Galaxy63-[HISAT2_on_KO_new].bam'), 
                        nthreads = 8, annot.ext = "mm10.annot", useMetaFeatures=TRUE)

#you manually create excel file with sequence info, then this loads it in
columnnames <- read_excel('Book1.xlsx')

#reads feature counts and sequence info into DESeq object
deseq <- DESeqDataSetFromMatrix(features$counts, columnnames, ~condition)

#performs analysis on DESeq object
deAnalysis <- DESeq(deseq)

#takes results of analysis
res <- results(deAnalysis)

head(res)
plotMA(res, main="Picture")
summary(res)

n = 50
resOrdered <- res[order(res$padj),]

#Omit NA values
re_resOrdered <- na.omit(resOrdered)

#concatenate by row
topResults <- rbind( re_resOrdered[ re_resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     re_resOrdered[ re_resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
#create table of top 50 results for mean, log2FC and adj p-value
top <- topResults[c(1:50,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]

#create csv file from table see above
write.csv(top, file = "DESeq2_top_log_fold_genes.csv")
