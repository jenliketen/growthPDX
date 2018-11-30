setwd("~/Desktop/growth/data/results")

library(Biobase)


# Load the PDXE dataset
pdxun <- readRDS("~/Desktop/growth_arvind/data/pdxe_untreated.Rda")


# Compute the Pearson correlation for gene expression vs. each of the 3 growth features
## All across the board; not tissue-specific
corGrowthAll <- function(feature, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  return(cor(x=t(expMatrix), y=eSet[, feature], method="pearson"))
}

cor.survival.allTissues <- corGrowthAll("time.last_published")
cor.slope.allTissues <- corGrowthAll("slope")
cor.doubleTime.allTissues <- corGrowthAll("timeToDouble_published")

cor.growth.allTissues <- data.frame(cor.survival.allTissues, cor.slope.allTissues, cor.doubleTime.allTissues)
colnames(cor.growth.allTissues) <- paste(c("Survival", "Slope", "DoublingTime"), "AllTissues", sep=".")

## Tissue-specific
corGrowthTissue <- function(feature, tissue, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  expTissue <- expMatrix[, pdxun@phenoData$tumor.type==tissue]
  featureTissue <- eSet[pdxun@phenoData$tumor.type==tissue, ]
  return(cor(x=t(expTissue), y=featureTissue[, feature], method="pearson"))
}

tmpMatrix <- NULL
for (tissue in unique(pdxun$tumor.type)) {
  cor.survival.tissue <- corGrowthTissue("time.last_published", tissue)
  cor.slope.tissue <- corGrowthTissue("slope", tissue)
  cor.doubleTime.tissue <- corGrowthTissue("timeToDouble_published", tissue)
  tmpMatrix <- cbind(tmpMatrix, cor.survival.tissue, cor.slope.tissue, cor.doubleTime.tissue)
}

cor.growth.tissue <- as.data.frame(tmpMatrix)

growthTissueNames <- c()
for (name in unique(pdxun$tumor.type)) {
  for (parameter in c("Survival", "Slope", "DoublingTime")) {
    x <- paste(parameter, name, sep=".")
    growthTissueNames <- c(growthTissueNames, x)
  }
}

colnames(cor.growth.tissue) <- growthTissueNames

## All tissues as well as tissue-specific
cor.growth <- cbind(cor.growth.allTissues, cor.growth.tissue)

saveRDS(cor.growth, file="cor_geneGrowth.Rda")


# Diagnostic plots for inspection
## All across the board; not tissue-specific
corPlotAll <- function(gene, feature, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  plot(expMatrix[gene, ], eSet[, feature],
       ylab=paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep=""),
       xlab="Expression",
       main=gene)
}

## Tissue-specific
corPlotTissue <- function(gene, feature, tissue, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  expTissue <- expMatrix[, pdxun@phenoData$tumor.type==tissue]
  featureTissue <- eSet[pdxun@phenoData$tumor.type==tissue, ]
  plot(expTissue[gene, ], featureTissue[, feature],
       ylab=paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep=""),
       xlab="Expression",
       main=gene)
}