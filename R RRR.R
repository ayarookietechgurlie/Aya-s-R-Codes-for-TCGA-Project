library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesameData)
library(sesame)

#list of projects
gdcprojects <- getGDCprojects()

#build first query for TCGA-BRCA Breast 
query_methyl1 <- GDCquery(project = 'TCGA-BRCA',
         data.category = 'DNA Methylation',
         platform = 'Illumina Human Methylation 450',
         access = 'open',
         data.type = 'Methylation Beta Value',
         barcode = c('TCGA-A7-A26E-01B-06D-A27B-05', 'TCGA-AR-A1AU-01A-11D-A12R-05'))

output_querymethyl1 <- getResults(query_methyl1)
GDCdownload(query_methyl1)

#TCGA-BRCA Breast Prepare 
dna.meth1 <- GDCprepare(query_methyl1, summarizedExperiment = TRUE)
assay(dna.meth1)

write.csv(assay(dna.meth1), file = "TCGA-BRCA1.csv")