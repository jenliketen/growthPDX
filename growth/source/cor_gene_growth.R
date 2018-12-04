setwd("~/Desktop/growth/data/results")

library(Biobase)


# Load the PDXE dataset
pdxun <- readRDS("~/Desktop/growth/data/pdxe_untreated.Rda")


# Compute the Pearson correlation for gene expression vs. each of the 3 growth features:
# doubling time (timeToDouble_published), survival (time.last_published), and slope
## All across the board; not tissue-specific
corGrowthAll <- function(feature, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  return(cor(x=t(expMatrix), y=eSet[, feature], method="pearson"))
}

cor.doublingTime.allTissues <- corGrowthAll("timeToDouble_published")
cor.survival.allTissues <- corGrowthAll("time.last_published")
cor.slope.allTissues <- corGrowthAll("slope")

cor.growth.allTissues <- data.frame(cor.doublingTime.allTissues, cor.survival.allTissues, cor.slope.allTissues)
colnames(cor.growth.allTissues) <- paste(c("DoublingTime", "Survival", "Slope"), "AllTissues", sep=".")

## Tissue-specific
corGrowthTissue <- function(feature, tissue, expMatrix=exprs(pdxun), eSet=pData(pdxun)) {
  expTissue <- expMatrix[, pdxun$tumor.type==tissue]
  featureTissue <- eSet[eSet$tumor.type==tissue, ]
  return(cor(x=t(expTissue), y=featureTissue[, feature], method="pearson"))
}

tmpMatrix <- NULL
for (tissue in unique(pdxun$tumor.type)) {
  cor.doublingTime.tissue <- corGrowthTissue("timeToDouble_published", tissue)
  cor.survival.tissue <- corGrowthTissue("time.last_published", tissue)
  cor.slope.tissue <- corGrowthTissue("slope", tissue)
  tmpMatrix <- cbind(tmpMatrix, cor.doublingTime.tissue, cor.survival.tissue, cor.slope.tissue)
}

cor.growth.tissue <- as.data.frame(tmpMatrix)

growthTissueNames <- c()
for (name in unique(pdxun$tumor.type)) {
  for (parameter in c("DoublingTime", "Survival", "Slope")) {
    x <- paste(parameter, name, sep=".")
    growthTissueNames <- c(growthTissueNames, x)
  }
}

colnames(cor.growth.tissue) <- growthTissueNames

## All tissues as well as tissue-specific
cor.growth <- cbind(cor.growth.allTissues, cor.growth.tissue)

saveRDS(cor.growth, file="cor_gene_growth.Rda")


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
  expTissue <- expMatrix[, pdxun$tumor.type==tissue]
  featureTissue <- eSet[eSet$tumor.type==tissue, ]
  plot(expTissue[gene, ], featureTissue[, feature],
       ylab=paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep=""),
       xlab="Expression",
       main=gene)
}