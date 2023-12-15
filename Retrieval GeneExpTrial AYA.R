library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesameData)
library(sesame)
library(biomaRt)

#build gene expression query
query_exp1 <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling')
output_queryexp1 <- getResults(query_exp1)

#Retrieve Gene Expression Data Query
query_exp1 <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c("TCGA-A7-A26E-01", "TCGA-PL-A8LV-01", "TCGA-BH-A0BC-01", "TCGA-AR-A1AX-01",
                                   "TCGA-AC-A2FO-01", "TCGA-AQ-A0Y5-01", "TCGA-AC-A3EH-01", "TCGA-AC-A5EH-01",
                                   "TCGA-D8-A1XB-01", "TCGA-D8-A1JJ-01", "TCGA-E2-A10B-01", "TCGA-BH-A1FJ-01",
                                   "TCGA-OL-A66P-01", "TCGA-D8-A1XJ-01", "TCGA-AR-A24P-01", "TCGA-AO-A1KQ-01",
                                   "TCGA-AR-A24V-01", "TCGA-E9-A54X-01", "TCGA-BH-A0BA-01", "TCGA-EW-A1PH-01",
                                   "TCGA-BH-A209-01", "TCGA-E9-A2JS-01", "TCGA-E9-A1R3-01", "TCGA-E9-A1R6-01",
                                   "TCGA-BH-A1FU-01", "TCGA-D8-A1Y2-01", "TCGA-A2-A4RX-01", "TCGA-AO-A12B-01",
                                   "TCGA-XX-A899-01", "TCGA-C8-A26W-01", "TCGA-E9-A249-01", "TCGA-B6-A0RV-01",
                                   "TCGA-D8-A1JT-01", "TCGA-OL-A5RV-01", "TCGA-Z7-A8R5-01", "TCGA-S3-A6ZG-01",
                                   "TCGA-AR-A24Q-01", "TCGA-A2-A1G6-01", "TCGA-D8-A1XL-01", "TCGA-A2-A1G0-01",
                                   "TCGA-S3-AA12-01", "TCGA-BH-A0AZ-01", "TCGA-E9-A1RA-01", "TCGA-UL-AAZ6-01",
                                   "TCGA-BH-A1ES-01", "TCGA-EW-A1PA-01", "TCGA-EW-A1P4-01", "TCGA-D8-A1JH-01",
                                   "TCGA-B6-A409-01", "TCGA-AC-A7VB-01", "TCGA-AN-A0XN-01", "TCGA-GM-A2DN-01",
                                   "TCGA-AO-A12G-01", "TCGA-BH-A0B9-01", "TCGA-C8-A273-01", "TCGA-AR-A1AY-01",
                                   "TCGA-AR-A1AO-01", "TCGA-E9-A1RC-01", "TCGA-A2-A0EN-01", "TCGA-E9-A1RD-01",
                                   "TCGA-A8-A08O-01", "TCGA-AC-A2FE-01", "TCGA-B6-A0RN-01", "TCGA-BH-A0HK-01",
                                   "TCGA-D8-A1XZ-01", "TCGA-OL-A66K-01", "TCGA-AR-A24N-01", "TCGA-E2-A1IU-01",
                                   "TCGA-AR-A0TT-01", "TCGA-AO-A03L-01", "TCGA-E9-A54Y-01", "TCGA-A2-A0T2-01",
                                   "TCGA-AC-A3W6-01", "TCGA-A2-A3XV-01", "TCGA-BH-A0B3-01", "TCGA-5L-AAT1-01",
                                   "TCGA-E2-A15J-01", "TCGA-BH-A0H7-01", "TCGA-A2-A0YT-01", "TCGA-AN-A0XO-01",
                                   "TCGA-EW-A1P3-01", "TCGA-LL-A7SZ-01", "TCGA-B6-A0IK-01", "TCGA-AO-A03M-01",
                                   "TCGA-BH-A0W4-01", "TCGA-D8-A27R-01", "TCGA-B6-A0WZ-01", "TCGA-D8-A27H-01",
                                   "TCGA-AO-A03U-01", "TCGA-A7-A5ZV-01", "TCGA-BH-A1EY-01", "TCGA-AC-A4ZE-01",
                                   "TCGA-BH-A0H9-01", "TCGA-LD-A7W6-01", "TCGA-GM-A2D9-01", "TCGA-BH-A42T-01",
                                   "TCGA-A2-A0SX-01", "TCGA-A1-A0SI-01", "TCGA-A7-A0D9-01", "TCGA-BH-A8FZ-01"))
getResults(query_exp1)

#download data
GDCdownload(query_exp1)

#TCGA-BRCA Breast Prepare GENE EXP
gene.exp1 <- GDCprepare(query_exp1, summarizedExperiment = TRUE)
assay(gene.exp1)

#Extract more info for gene names
ensembl1 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
your_ensg_codes <- c("ENSG0000139618", "ENSG0000141510", "ENSG0000141510")
gene_names <- mapIds(ensembl, keys = your_ensg_codes, column = "ensembl_gene_id", 
                     keytype = "ensembl_gene_id", mart = ensembl)

# Replace ENSG codes with gene names on Y axis
your_y_axis <- c("ENSG00000139618", "ENSG00000141510", "ENSG00000141510")
gene_names <- mapIds(ensembl, keys = your_y_axis, column = "ensembl_gene_id", 
                     keytype = "ensembl_gene_id", mart = ensembl)
names(your_y_axis) <- gene_names$external_gene_name


#Retrieval
write.csv(assay(gene.exp1), file = "TCGA-BRCAGeneExp3.csv")


