#' Take an expression matrix and if duplicate genes are measured, summarize them by their means
#' 
#' This function accepts two expression matrices, with gene ids as rownames() and 
#' sample ids as colnames(). It will find all duplicate genes and summarize their
#' expression by their mean.
#'
#' @param exprMat a gene expression matrix with gene names as row ids and sample names as column ids.
#'
#' @return a gene expression matrix that does not contain duplicate gene ids
#'
#' @keywords summarize duplicate gene ids by their mean.
#'
#' @export
summarizeGenesByMean <- function(exprMat) {
  
  geneIds <- rownames(exprMat)
  t <- table(geneIds) # how many times is each gene name duplicated
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]
  
  # create a *new* gene expression matrix with everything in the correct order....
  # start by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), , drop = FALSE]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]
  
  # add all the duplicated genes to the bottom of "exprMatUnique", summarizing as you go
  for(numDups in allNumDups) {
    geneList <- names(which(t == numDups))
    
    for(i in 1:length(geneList)) {
      # Get the duplicated rows for this gene
      duplicated_rows <- exprMat[which(geneIds == geneList[i]), , drop = FALSE]
      
      # Calculate mean across duplicated rows
      if (nrow(duplicated_rows) == 1) {
        mean_expr <- duplicated_rows[1, ]
      } else {
        mean_expr <- colMeans(duplicated_rows, na.rm = TRUE)
      }
      
      # Add to the unique matrix
      exprMatUnique <- rbind(exprMatUnique, mean_expr)
      gnamesUnique <- c(gnamesUnique, geneList[i])
    }
  }
  
  # Fixed class check for R 4.0+ compatibility
  # Check if exprMatUnique is a vector (which happens when there's only one gene)
  if (is.vector(exprMatUnique) || (!is.matrix(exprMatUnique) && !is.data.frame(exprMatUnique))) {
    # Convert vector to matrix
    exprMatUnique <- matrix(exprMatUnique, ncol = length(exprMatUnique))
  }
  
  # Ensure it's a proper matrix
  exprMatUnique <- as.matrix(exprMatUnique)
  rownames(exprMatUnique) <- gnamesUnique
  
  return(exprMatUnique)
}

