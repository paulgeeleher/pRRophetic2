% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/compute_phenotype_function.R
\name{predictionAccuracyByCv}
\alias{predictionAccuracyByCv}
\title{Cross validation on training dataset}
\usage{
predictionAccuracyByCv(trainingExprData, trainingPtype, testExprData = -1,
  cvFold = -1, powerTransformPhenotype = TRUE, batchCorrect = "eb",
  removeLowVaryingGenes = 0.2, minNumSamples = 10, selection = 1)
}
\arguments{
\item{trainingExprData}{The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"}

\item{trainingPtype}{The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".}

\item{testExprData}{The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".}

\item{cvFold}{Specify the "fold" requried for cross validation. "-1" will do leave one out cross validation (LOOCV)}

\item{powerTransformPhenotype}{Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.}

\item{batchCorrect}{How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.}

\item{removeLowVaryingGenes}{What proportion of low varying genes should be removed? 20 precent be default}

\item{minNumSamples}{How many training and test samples are requried. Print an error if below this threshold}

\item{selection}{How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.}

\item{printOutput}{Set to FALSE to supress output}
}
\value{
An object of class "pRRopheticCv", which is a list with two members, "cvPtype" and "realPtype", which correspond to the cross valiation predicted phenotype and the  user provided measured phenotype respectively.
}
\description{
This function does cross validation on a training set to estimate prediction accuracy on a training set.
If the actual test set is provided, the two datasets can be subsetted and homogenized before the
cross validation analysis is preformed. This may improve the estimate of prediction accuracy.
}
\keyword{phenotype}
\keyword{predict}

