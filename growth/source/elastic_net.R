setwd("~/Desktop/growth/data/results")

library(Biobase)
library(caret)
library(doMC)
library(plotROC)


# Prepare matrix of predictors (gene expresison) and vector of response (doubling time)
pdxun <- readRDS("~/Desktop/growth/data/pdxe_untreated.Rda")
dgea <- readRDS("~/Desktop/growth/data/results/diff_gene_exp_doublingTime_N.Rda")
dgea <- dgea$`n=50`
dgea <- dgea[which(dgea$P.Value < 0.05), ]

## Gene expression matrix
expMatrix <- exprs(pdxun)
expMatrix <- expMatrix[which(apply(expMatrix, 1, var) != 0), ]
expMatrix <- expMatrix[rownames(dgea), ]
expMatrix <- t(expMatrix)

## Doubling time vector
growthClass <- data.frame(id = rownames(pData(pdxun)),
                          DoublingTime=pData(pdxun)$timeToDouble_published,
                          Survival=pData(pdxun)$time.last_published, stringsAsFactors=FALSE)
growthClass <- growthClass[-which(growthClass$Survival < 10), ] ## Include only samples with greater than 10 days survival
                                                                ## Less than 10 days indicates mouse died of other reasons

rownames(growthClass) <- growthClass$id

expMatrix <- expMatrix[growthClass$id, ] ## Subset the gene expression matrix accordingly

growthClass$dtClass <- ifelse(growthClass$DoublingTime > median(growthClass$DoublingTime), "high", "low")

dt <- factor(growthClass$dtClass); names(dt) <- growthClass$id

## Preprocess the expression matrix
preProcValues <- preProcess(expMatrix, method=c("center", "scale"))
expMatrix <- predict(preProcValues, expMatrix)


# Build the elastic net model using caret workflow
predictAccuracy <- function(model, data, response) {
  predictions <- predict(model, data.matrix(data), type="prob")
  predictions$class <- predict(model, data.matrix(data))

  accuracy <- caret::confusionMatrix(predictions$class, response)
  print(accuracy)

  return(predictions)
}

trainModel <- function(x.train, y.train, x.test=NULL, y.test=NULL, NCPU=3) {
  if (NCPU > 1) {registerDoMC(cores=NCPU)}
  
  trnCtrl <- trainControl(method="repeatedcv", number=4,
                         allowParallel=TRUE, savePredictions=TRUE)

  
  lambda.grid <- 10^seq(2, -10, length=10)
  alpha.grid <- seq(0, 0.5, length = 6)
    
  srchGrid = expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)
    
  my_train <- train(x=x.train, y=y.train,
                      method="glmnet", tuneGrid=srchGrid, trControl=trnCtrl,
                      standardize=FALSE, maxit=1000000)

  self.predict <- predictAccuracy(my_train, x.train, y.train)
  
  if(!is.null(x.test) && !is.null(y.test)) {
    test.predict <- predictAccuracy(my_train, x.test, y.test)
  }
  
  return(list(fit=my_train, training_predictions=self.predict, test_predictions=test.predict))
}


# Nested cross validation within the PDXE dataset (external test sets: TCGA and PDMR)
allMod <- list()
pdf("roc_curves.pdf",
    width=12, height=8)

set.seed(42)

testIndex <- createDataPartition(dt, times=4, p=0.1, list=TRUE)
for(ts in testIndex) {
  genes.train=expMatrix[-ts, ]; genes.test= expMatrix[ts, ]
  
  dt.train=dt[rownames(genes.train)]; dt.test=dt[rownames(genes.test)]
  
  elastic_net <- trainModel(genes.train, dt.train, genes.test, dt.test, NCPU=3)
  
  allMod[[length(allMod)+1]] <- elastic_net
  
  predictions <- elastic_net$test_predictions
  predictions$org.class <- ifelse(dt.test=="high", 1, 0)
  roc <- ggplot(predictions, aes(m=low, d=org.class)) + geom_roc() + style_roc() + coord_equal()
  auc <- calc_auc(roc)
  
  print(roc)
  print(auc)
}
dev.off()