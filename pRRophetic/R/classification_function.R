#' Dichotimize a training expression set and fit a logistic ridge regression model which is applied to the test expression matirx.
#' 
#' Dichotimize a training expression set and fit a logistic ridge regression model which is applied to the test expression matirx.
#' This function will return a set of probabilities.
#' 
#' @import genefilter
#' @import ridge
#' 
#' @param trainingExprData - Gene expression matrix for samples for which we the phenotype is already known.
#' @param trainingPtype The known phenotype, a vector in the same order as the columns of "trainingExprData" or with the same names as colnames of "trainingExprData".
#' @param testExprData Gene expression matrix for samples on which we wish to predict a phenotype. Gene names as rows, samples names as columns.
#' @param batchCorrect The type of batch correction to be used. Options are "eb", "none", .....
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#' @param numGenesSelected Specifies how genes are selected for "variableSelectionMethod". Options are "tTests", "pearson" and "spearman".
#' @param numSens The number of sensitive cell lines to be fit in the logistic regression model.
#' @param numRes The number of resistant cell lines fit in the logistic regression model.
#' @param minNumSamples The minimum number of test samples, print an error if the number of columns of "testExprData" is below this threshold. A large number of test samples may be necessary to correct for batch effects.
#' 
#' @export
classifySamples <- function(trainingExprData, trainingPtype, testExprData, batchCorrect="eb", minNumSamples=10, selection=-1, printOutput=TRUE, numGenesSelected=1000, numSens=15, numRes=55)
{
  sensInd <- order(trainingPtype)[1:numSens]
  resInd <- order(trainingPtype)[(length(trainingPtype)-numRes):length(trainingPtype)]

  homDataErlot <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)

  tTests <- rowttests(data.matrix(cbind(homDataErlot$train[, sensInd], homDataErlot$train[, resInd])), as.factor(c(rep("sens", length(sensInd)), rep("res", length(resInd)))))
  topCorVarGenes <- rownames(tTests[order(tTests[, "p.value"]),])[1:numGenesSelected] # this makes more sense than spearman of pearson....

  trainExpr <- t(homDataErlot$train[ topCorVarGenes, c(sensInd, resInd)])
  trainPtyle <- as.numeric(as.factor(c(rep("sens", length(sensInd)), rep("res", length(resInd))))) - 1

  trainDat <- data.frame(trainPtyle, trainExpr)
  cat("Fitting model, may take some time....\n")
  ridgeLogisticModel_all <- logisticRidge(trainPtyle~., data=trainDat)

  preDataLRR <- data.frame(t(homDataErlot$test[topCorVarGenes, ]))
  predsLRR <- predict(ridgeLogisticModel_all, preDataLRR)

  return(predsLRR)
}

## (under construction....)
## use X fold (default 10 fold) cross validataion to estimate the accuracy of our classifier.... 
#testClassification <- function(...., fold=10)