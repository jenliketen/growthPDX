library(TCGAbiolinks)
library(SummarizedExperiment)


# Download data from TCGA, for now we will download the 5 tissue types found in the PDXE training set
## The nomenclature in the two datasets is slightly different. We note them here:
### TCGA-BRCA, breast invasive carcinoma (name in training set: BRCA, breast carcinoma)
### TCGA-SKCM, skin cutaneous melanoma (name in training set: CM, cutaneous melanoma)
### TCGA-COAD, colon adenocarcinoma (name in training set: CRC, colorectal cancer)
### TCGA-READ, rectum adenocarcinoma, (name in training set: CRC, colorectal cancer)
#### TCGA divides colorectal cancer into colon and rectum adenocarcinomas; we will be downloading both
### TCGA-LUAD, lung squamous cell carcinoma (name in training set: NSCLC, non-small cell lung carcinoma)
### TCGA-PAAD, pancreatic adenocarcinoma (name in training set: PDAC, pancreatic ductal carcinoma)
get_geneX <- function(tissue) {
  query <- GDCquery(project=paste("TCGA", tissue, sep="-"), data.category="Gene expression",
                    data.type="Gene expression quantification",
                    experimental.strategy="RNA-Seq",
                    platform="Illumina HiSeq",
                    file.type="normalized_results",
                    legacy=TRUE)
  GDCdownload(query=query, method="api")
  rnaseqSE <- GDCprepare(query)
  expMatrix <- assay(rnaseqSE, "normalized_count")
  
  return(expMatrix)
}

tissueTypes <- c("BRCA", "SKCM", "COAD", "READ", "LUAD", "PAAD")

expMatrices <- list()
for (tissueType in tissueTypes) {
  rnaSeq_tissue <- get_geneX(tissueType)
  
  expMatrices[[paste0(tissueType)]] <- rnaSeq_tissue
}

saveRDS(expMatrices, file="~/Desktop/growth/data/results/tcga_exp_matrices.Rda")