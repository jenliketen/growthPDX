setwd("~/Desktop/growth/data/results")

library(BBmisc)
library(ggplot2)
library(edgeR)
library(limma)


# Load the PDXE dataset
pdxun <- readRDS("~/Desktop/growth/data/pdxe_untreated.Rda")


# Group the growth data into high and low classes, according to median
growthClass <- data.frame(id = rownames(pData(pdxun)),
                          DoublingTime=pData(pdxun)$timeToDouble_published,
                          Survival=pData(pdxun)$time.last_published,
                          Slope=pData(pdxun)$slope, stringsAsFactors=FALSE)
growthClass$dtClass <-  ifelse(growthClass$DoublingTime<median(growthClass$DoublingTime), "Low", "High")
growthClass$surClass <- ifelse(growthClass$Survival<median(growthClass$Survival), "Low", "High")
growthClass$slClass <- ifelse(growthClass$Slope<median(growthClass$Slope), "Low", "High")

## Visualize the classes
growthClass.boxplot <- function(growthClassData, varClass,varFeature) {
  colClass <- c("#BC3C29F2", "#0072B5F2")
  plot <- ggplot(data=growthClassData, aes_string(x=varClass, y=varFeature, fill=varClass))
  plot <- plot + geom_boxplot() + scale_fill_manual(values=colClass)
  plot <- plot + theme_bw() + theme(panel.border=element_blank(),
                              panel.grid.major=element_blank(),
                              panel.grid.minor=element_blank(),
                              axis.line=element_line(color="black"),
                              plot.title=element_text(hjust=0.5))
  plot <- plot + geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5)
  plot <- plot + theme(legend.position="none")
  return(plot)
}

pdf("growth_plot_classes.pdf",
    width=12, height=8)
growthClass.boxplot(growthClassData=growthClass, varClass="dtClass", varFeature="DoublingTime")
growthClass.boxplot(growthClassData=growthClass, varClass="surClass", varFeature="Survival")
growthClass.boxplot(growthClassData=growthClass, varClass="slClass", varFeature="Slope")
dev.off()


##### FROM NOW ON FOCUS ON DOUBLING TIME ONLY #####

# Check for heteroscedasticity
cpm.mat <- log(cpm(exprs(pdxun)))
mean.vec <- apply(cpm.mat, 1, mean)
sd.vec <- apply(cpm.mat, 1, sd)
doublingTimeSpread <- plot(mean.vec, sd.vec, pch=".",
                           ylab="sd", xlab="Average logCPM", main="PDXE")
## As expected, we see from the plot that the expression values are heteroscedastic


# Subset the highest and lowest 25 doubling time values
## 1x1x1 PDX experimental design; number of replicates is 1!!!
## We do not know how much error there is with the median doubling time, so we take values at
## the extremes to fall on the safe side
growthFeature <- "DoublingTime"
growthClass <- sortByCol(x=growthClass, col=growthFeature, asc=FALSE)
n=25
extremes <- growthClass[c(1:n, (nrow(growthClass)-n+1):nrow(growthClass)), ] ## For now we will take 25 from each end
eSet.subsetByN <- pData(pdxun)[extremes$id, ]
table(eSet.subsetByN$tumor.type, extremes$dtClass) ## Make sure all tissue types are included in our selection

pdf("dtClasses_plot.pdf",
    width=12, height=8)
growthClass.boxplot(growthClassData=extremes, varClass="dtClass", varFeature="DoublingTime")
dev.off()


# Preprocessing for the expression matrix; generate Z-scores
## The function assumes that the column means are zero
centerAndScale <- function(trData, tsData=NULL) {
  preProcValues <- caret::preProcess(trData, method=c("center", "scale")) 
  trS <- predict(preProcValues, trData)
  
  if (is.null(tsData)) {
    return(trS)
    }
  
  tsS <- predict(preProcValues, tsData) ## The tsData will be the TCGA data we feed in for ML later
  return(list(tr_data=trS, ts_data=tsS, scale_fact=preProcValues))
}

expS <- exprs(pdxun)[, extremes$id]
expS <- expS[-which(apply(expS, 1, var)==0), ] ## Remove genes with zero variance or function complains
expS <- t(centerAndScale(t(expS))) ## Tranpose twice because we want our row means to be zero
dim(expS)


# Limma baseline model for t-test
## Ideally we would use DESeq2 or edgeR to correct for heteroscedasticity, but alas, we do not
## have raw count data
doDiffExp <- function(mat, flab) {
  design <- model.matrix(~flab)
  fitTrtMean <- limma::lmFit(mat, design)
  efit <- limma::eBayes(fitTrtMean)
  diffLst <- topTable(efit, coef=2, number=dim(mat)[1], sort.by="logFC")
  diffLst$abs.logFC <- abs(diffLst$logFC)
  diffLst <- BBmisc::sortByCol(diffLst, c("abs.logFC"), asc=FALSE)
  
  return(diffLst)
}

doublingTime.doDiffExp <- doDiffExp(expS, extremes$dtClass)

saveRDS(doublingTime.doDiffExp, file="diff_gene_exp_doublingTime.Rda")