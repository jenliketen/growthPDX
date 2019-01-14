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
growthClass <- growthClass[-which(growthClass$Survival < 10), ] ## Include only samples with greater than 10 days survival
                                                                ## Less than 10 days indicates mouse died of other reasons
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


################## 1x1x1 PDX experimental design; number of replicates is 1!!! ###################
### We do not know how much error there is with the median doubling time, so we take values at ###
############################## the extremes to fall on the safe side #############################

# First, we try subsetting the highest and lowest N doubling time values
## Subset N samples
samples.subset <- function(growthFeature, n) {
  growthClass <- sortByCol(x=growthClass, col=growthFeature, asc=FALSE)
  subsettedSamples <- growthClass[c(1:n, (nrow(growthClass)-n+1):nrow(growthClass)), ]
  phenoDataSubsetted <- pData(pdxun)[subsettedSamples$id, ]
  tissuesIncluded <- table(phenoDataSubsetted$tumor.type, subsettedSamples$dtClass) ## Make sure all tissue types are included in our selection
  if (any(tissuesIncluded==0)) {
    warning("Not all tissue types are included. Would you like to change the value of N?")
  }
  return(subsettedSamples)
  }

## Subset expression data for the selected N samples
expression.subset <- function(n) {
  growthClass <- sortByCol(x=growthClass, col="DoublingTime", asc=FALSE)
  subsettedSamples <- growthClass[c(1:n, (nrow(growthClass)-n+1):nrow(growthClass)), ]
  expressionSubsetted <- exprs(pdxun)[, subsettedSamples$id]
  return(expressionSubsetted)
}

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

## Prepare sample and expression data
sampleDataAllN <- list()
expDataAllN <- list()
for (i in c(18, seq(20, 50, by=5))) {
  samples <- samples.subset("DoublingTime", n=i)
  expS <- expression.subset(n=i)
  expS <- expS[which(apply(expS, 1, var)!=0), ] ## Remove genes with zero variance or function complains
  expS <- t(centerAndScale(t(expS))) ## Tranpose twice because we want our row means to be zero
  sampleDataAllN[[paste0("n=", i)]] <- samples
  expDataAllN[[paste0("n=", i)]] <- expS
}

## Plots for inspection
pdf("dtClasses_plot.pdf",
    width=12, height=8)
for (i in 1:length(sampleDataAllN)) {
  print(growthClass.boxplot(growthClassData=sampleDataAllN[[i]], varClass="dtClass", varFeature="DoublingTime"))
}
dev.off()

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

doublingTime.doDiffExp.N <- list()
for (i in 1:length(sampleDataAllN)) {
  doublingTime.doDiffExp.N[[i]] <- doDiffExp(expDataAllN[[i]], sampleDataAllN[[i]]$dtClass)
}

howManySamples <- c()
for (i in c(18, seq(20, 50, by=5))) {
  howManySamples <- c(howManySamples, paste0("n=", i))
}

names(doublingTime.doDiffExp.N) <- howManySamples

saveRDS(doublingTime.doDiffExp.N, file="diff_gene_exp_doublingTime_N.Rda")


# We can also try taking the high values above the 3rd quartile and the low ones below the 1st quartile
## Prepare sample and expression data
above_3_quart <- growthClass[which(growthClass$DoublingTime > quantile(growthClass$DoublingTime, 0.75)), ]
below_1_quart <- growthClass[which(growthClass$DoublingTime < quantile(growthClass$DoublingTime, 0.25)), ]
samplesByQuartiles <- rbind(above_3_quart, below_1_quart)
dim(samplesByQuartiles) ## In this case N=43, similar to the case we did for N=45

expressionByQuartiles <- expression.subset(n=nrow(above_3_quart))
expressionByQuartiles <- expressionByQuartiles[-which(apply(expressionByQuartiles, 1, var)==0), ] ## Remove genes with zero variance or function complains
expressionByQuartiles <- t(centerAndScale(t(expressionByQuartiles))) ## Tranpose twice because we want our row means to be zero

doublingTime.doDiffExp.quartiles <- doDiffExp(expressionByQuartiles, samplesByQuartiles$dtClass)

saveRDS(doublingTime.doDiffExp.quartiles, file="diff_gene_exp_doublingTime_quartiles.Rda")