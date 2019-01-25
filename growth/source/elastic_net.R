library(Biobase)
library(caret)


pdxun <- readRDS("~/Desktop/growth/data/pdxe_untreated.Rda")
dgea <- readRDS("~/Desktop/growth/data/results/diff_gene_exp_doublingTime_N.Rda")
dgea <- dgea$`n=50`
dgea <- dgea[which(dgea$P.Value < 0.05), ]

expMatrix <- exprs(pdxun)
expMatrix <- expMatrix[which(apply(expMatrix, 1, var) != 0), ]
expMatrix <- expMatrix[rownames(dgea), ]
expMatrix <- t(expMatrix)


growthClass <- data.frame(id = rownames(pData(pdxun)),
                          DoublingTime=pData(pdxun)$timeToDouble_published,
                          Survival=pData(pdxun)$time.last_published, stringsAsFactors=FALSE)
growthClass <- growthClass[-which(growthClass$Survival < 10), ] ## Include only samples with greater than 10 days survival
                                                                ## Less than 10 days indicates mouse died of other reasons

rownames(growthClass) <- growthClass$id

expMatrix <- expMatrix[growthClass$id, ] ## Subset the gene expression matrix accordingly

growthClass$dtClass <-  ifelse(growthClass$DoublingTime > median(growthClass$DoublingTime), "high", "low")

dtOutcome <- factor(growthClass$dtClass)


# Z-transform expression matrix
preProcValues <- preProcess(expMatrix, method=c("center", "scale"))
expMatrix <- predict(preProcValues, expMatrix)


# Cross validation and building the model
lambda.grid <- 10^seq(2, -10, length=10)
alpha.grid <- seq(0, 0.5, length=6)
searchGrid <- expand.grid(alpha=alpha.grid, lambda=lambda.grid)

trCtrl <- trainControl(method="repeatedcv", number=4, p=0.8, repeats=5, allowParallel=TRUE, savePredictions = TRUE)

set.seed(42)

model <- train(x=expMatrix[, 1:100], y=dtOutcome,
                  method="glmnet", tuneGrid=searchGrid, trControl=trCtrl,
                  standardize=FALSE, maxit=1000000)
plot.train(model)

# Self-predict
testIndex <- createDataPartition(response, times=4, p=0.1, list=TRUE)
response <- growthClass$dtClass; names(response) <- growthClass$id

allMod <- list()
for(ts in testIndex)
{
  trainMatrix=expMatrix[-ts, ]; testMatrix= expMatrix[ts, ]
  
  trainResponse=response[rownames(trainMatrix)]; testResponse=response[rownames(testMatrix)]
  
  model <- train(x=trainMatrix, y=trainResponse,
                 method="glmnet", tuneGrid=searchGrid, trControl=trCtrl,
                 standardize=FALSE, maxit=1000000)
  
  self_predict <- predict(model, newdata=data.matrix(testMatrix), type="prob")
  self_predict$class <- predict(model, newdata=data.matrix(testMatrix))
  accuracy <- confusionMatrix(self_predict$class, testResponse)
  
  
  allMod[[length(allMod)+1]] <- accuracy
}

self_predict <- predict(model, newdata=data.matrix(expMatrix[, 1:100]), type="prob")
self_predict$class <- predict(model, newdata=data.matrix(expMatrix[, 1:100]))
accuracy <- confusionMatrix(self_predict$class, dtOutcome) ## sum(self_predict$class==dtOutcome)/length(dtOutcome)

# top positive and negative genes
