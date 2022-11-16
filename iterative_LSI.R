
#####################################
# Cluster scRNA using iterative LSI
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(Matrix)
})

################
# Iterative LSI
################

getVarGenes <- function(mat, nvar = 2000, blacklist = NULL){
  # Get the top nvar variable genes present in mat (a gene x sample/cell matrix)
  # If blacklist is present, do not return any genes in the blacklist
  if(!is.null(blacklist)){
    ncount <- nrow(mat)
    mat <- mat[!rownames(mat) %in% blacklist,]
    message(sprintf("Removed %s genes overlapping blacklist prior to selecting variable genes...", ncount - nrow(mat)))
  }
  if(is(mat, "sparseMatrix")){
    varGenes <- rownames(mat)[head(order(sparseMatrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }else{
    varGenes <- rownames(mat)[head(order(matrixStats::rowVars(mat), decreasing = TRUE), nvar)]
  }
  return(varGenes)
}


runLSI <- function(mat, nComponents, binarize = FALSE){
  # TF-IDF LSI adapted from Jeff Granja, who adapted from flyATAC (i.e. Cusanovich et al. 2018)
  #
  # Calculate the Term Frequency - Inverse Document Frequency (TF-IDF) for a feature x cell counts
  # matrix, then calculate the Singular Value Decomposition of that matrix, which is then used as
  # input for Surat's SNN clustering
  if(binarize){
    message("Binarizing matrix...")
    # The 'x' slot of the dgCMatrix class contains the non-zero elements of the matrix
    mat@x[mat@x > 0] <- 1
  }
  #Calculate RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)

  # Calculate TF-IDF
  message(sprintf("Calculating TF-IDF with %s features (terms) and %s cells (documents)...", nrow(mat), ncol(mat)))
  start <- Sys.time()
  scaleTo <- 10^4
  tf <- t(t(mat) / colSm)
  idf <- as(ncol(mat) / rowSm, "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  # Log transform TF-IDF
  tfidf <- sparseLogX(tfidf, logtype="log", scale=TRUE, scaleFactor=scaleTo)
  
  # Clean up
  rm(tf)
  rm(idf)
  invisible(gc())
  
  # Calculate SVD for LSI
  message("Calculating SVD for LSI...")
  svd <- irlba::irlba(tfidf, nv=nComponents, nu=nComponents)
  svdDiag <- Matrix::diag(x=svd$d)
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC", seq_len(ncol(matSVD)))
  
  # Return matSVD and svd
  message(sprintf("LSI complete: %s minutes", round(difftime(Sys.time(), start, units="mins"), 3)))
  if(is.null(rownames(mat))){
    rownames(mat) <- 1:nrow(mat)
  }
  return(
    list(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm, 
        svd = svd, 
        binarize = binarize
        )
    )
}



