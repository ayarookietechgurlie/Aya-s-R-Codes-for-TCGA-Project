query_met <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform = c("Illumina Human Methylation 450")
)
query_exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(
  substr(getResults(query_met, cols = "cases"), 1, 100),
  substr(getResults(query_exp, cols = "cases"), 1, 100)
)

# Only seelct the first 5 patients
query_metBRCA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = c("Illumina Human Methylation 450"),
  barcode = common.patients[1:100])

query_expBRCA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode = common.patients[1:100])

#Extract Query Codes
write.csv((query_expBRCA[[1]][[1]]), file = "TCGABRCAExp_Cases.csv")
write.csv((query_metBRCA[[1]][[1]]), file = "TCGABRCAMet_Cases.csv")