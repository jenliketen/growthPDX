library(TCGAbiolinks)
library(SummarizedExperiment)


# Download data from TCGA, for now we will download the 5 tissue types found in the PDXE training set
## TCGA-BRCA, breast invasive carcinoma (name in training set: BRCA, breast carcinoma)
query_BRCA <- GDCquery(project="TCGA-BRCA", data.category="Gene expression",
                  data.type="Gene expression quantification",
                  experimental.strategy="RNA-Seq",
                  platform="Illumina HiSeq",
                  file.type="normalized_results",
                  legacy=TRUE)
GDCdownload(query=query_BRCA, method="api")
BRCARnaseqSE <- GDCprepare(query_BRCA)
BRCAMatrix <- assay(BRCARnaseqSE, "normalized_count")

## TCGA-SKCM, skin cutaneous melanoma (name in training set: CM, cutaneous melanoma)
query_SKCM <- GDCquery(project="TCGA-SKCM", data.category="Gene expression",
                       data.type="Gene expression quantification",
                       experimental.strategy="RNA-Seq",
                       platform="Illumina HiSeq",
                       file.type="normalized_results",
                       legacy=TRUE)
GDCdownload(query=query_SKCM, method="api")
SKCMRnaseqSE <- GDCprepare(query_SKCM)
SKCMMatrix <- assay(SKCMRnaseqSE, "normalized_count")

## TCGA-COAD, colon adenocarcinoma (name in training set: CRC, colorectal cancer)
### TCGA divides colorectal cancer into colon and rectum adenocarcinomas; we will be downloading both
query_COAD <- GDCquery(project="TCGA-COAD", data.category="Gene expression",
                       data.type="Gene expression quantification",
                       experimental.strategy="RNA-Seq",
                       platform="Illumina HiSeq",
                       file.type="normalized_results",
                       legacy=TRUE)
GDCdownload(query=query_COAD, method="api")
COADRnaseqSE <- GDCprepare(query_COAD)
COADMatrix <- assay(COADRnaseqSE, "normalized_count")

## TCGA-READ, rectum adenocarcinoma, (name in training set: CRC, colorectal cancer)
### TCGA divides colorectal cancer into colon and rectum adenocarcinomas; we will be downloading both
query_READ <- GDCquery(project="TCGA-READ", data.category="Gene expression",
                       data.type="Gene expression quantification",
                       experimental.strategy="RNA-Seq",
                       platform="Illumina HiSeq",
                       file.type="normalized_results",
                       legacy=TRUE)
GDCdownload(query=query_READ, method="api")
READRnaseqSE <- GDCprepare(query_READ)
READMatrix <- assay(READRnaseqSE, "normalized_count")

## TCGA-LUAD, lung squamous cell carcinoma (name in training set: NSCLC, non-small cell lung carcinoma)
query_LUAD <- GDCquery(project="TCGA-LUAD", data.category="Gene expression",
                       data.type="Gene expression quantification",
                       experimental.strategy="RNA-Seq",
                       platform="Illumina HiSeq",
                       file.type="normalized_results",
                       legacy=TRUE)
GDCdownload(query=query_LUAD, method="api")
LUADRnaseqSE <- GDCprepare(query_LUAD)
LUADMatrix <- assay(LUADRnaseqSE, "normalized_count")

## TCGA-PAAD, pancreatic adenocarcinoma (name in training set: PDAC, pancreatic ductal carcinoma)
query_PAAD <- GDCquery(project="TCGA-PAAD", data.category="Gene expression",
                       data.type="Gene expression quantification",
                       experimental.strategy="RNA-Seq",
                       platform="Illumina HiSeq",
                       file.type="normalized_results",
                       legacy=TRUE)
GDCdownload(query=query_PAAD, method="api")
PAADRnaseqSE <- GDCprepare(query_PAAD)
PAADMatrix <- assay(PAADRnaseqSE, "normalized_count")