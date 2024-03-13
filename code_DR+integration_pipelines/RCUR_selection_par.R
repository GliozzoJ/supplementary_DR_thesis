library(doParallel)
library(cytominer)
#' Feature selection using clustering
#' 
#' @description This function performs unsupervised feature selection exploiting the RCUR algorithm
#' In practice, it considers subsets of features and applies RCUR
#' to select the features that allow to maintain most of the variance iin the dataset. 
#' After all subset
#' of features are processed, the procedure is repeated maxit times on the selected
#' features at the previous iteration. Thus, at each iteration the set of selected
#' features is refined. 
#'
#' @param df dataframe. Dataframe (examples x features).
#' @param maxD integer. Number of features estimated using the blocking ID method.
#' @param maxit integer. Maximum number of iterations (def. 10).
#'
#' @return Matrix (samples x features) with the selected set of variables.
#' @export
RCUR_selection_par <- function(df, maxD =  round((min(dim(df))-1)/2), maxit = 3, 
                                            patience = 1000, ncores = my_detectCores()){
  X = df
  cat("RCUR selection\n")
  cat("maxD = ", maxD, "\n")
 
  mat_votes = rep(0, times = ncol(X))
  names(mat_votes) = colnames(X)
  
  cl <- makePSOCKcluster(min(ncores, detectCores()-1))
  registerDoParallel(cl)
  
  # if I have missing values I impute
  Ximp = imputa(X)
  k = min(maxD, min(50,round(min(dim(Ximp)/5)-1)))   # the smaller the k, the faster the algorithm
  selected = foreach (numIts = 1:patience, .combine='c', .packages = c("rsvd")) %dopar% {
    r <- rcur(Ximp, k = min(maxD, k), rand = TRUE)
    return(colnames(r$C))
  }
  
  selected = selected[!is.na(selected)]
  tt = sort(table(selected), decreasing = TRUE)
  tt = tt[tt>0]
  selected = names(tt[1:min(length(tt),maxD)])
  
  stopCluster(cl)
  unregister_dopar()
  
  #return non imputed columns
  X = X[, match(selected, colnames(X))]

  cat('final dimension = ', dim(X), '\n')

  return(X)
}