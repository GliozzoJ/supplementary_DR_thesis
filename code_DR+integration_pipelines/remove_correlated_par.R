library(doParallel)



#' Parallel determination of high correlated variables
#' 
#' @description This function considers subsets of features and 
#' removes correlated features pair-wise from each subset. If two variables have 
#' a high correlation, the function looks at the mean absolute correlation of 
#' each variable and removes the variable with the largest mean absolute 
#' correlation. The computation on subsets is done in parallel.
#' Once a passage through the entire set of features is performed, the remaining 
#' features are shuffled and the process is repeated again (till maxiter is 
#' reached).
#'
#' @param df dataframe. Dataframe examples x features.
#' @param dim_split integer. Number of features considered at each iteration 
#' (def. 1000).
#' @param maxiter integer. Maximum number of iterations (def. 10).
#' @param method string. Method used to compute correlation. Options are
#' c("pearson", "kendall", "spearman"). Default "kendall".
#' @param cutoff numeric.A numeric value for the pair-wise absolute correlation 
#' cutoff (def. 0.6).
#'
#' @return A matrix (samples x features) having only the selected subset of 
#' variables.
#' @export
remove_correlated_par <- function(mat, dim_split = 1000, maxiter = 10, 
                                  method = "kendall", cutoff = 0.6, ncores = my_detectCores()){
  
  X = t(mat)
  niter = 0
  if (!is.finite(dim_split)) dim_split = nrow(X)
  cat("dim(X) before starting remove correlation: ", dim(X), "\n")
  while(niter < maxiter ){
    cat("niter = ", niter, '\n')
    filtered_X = NULL
    
    cl <- makeCluster(min(ncores, detectCores()-1))
    registerDoParallel(cl)
    
    filtered_X = foreach(nR = seq(1, nrow(X), by=dim_split), 
                         .combine='rbind', .packages = c("caret", "ggplot2")) %dopar% {
                           
         subX = X[nR:min(nrow(X), (nR+dim_split-1)), ]
         if (is.null(dim(subX))) return(subX)
         cc = cor(t(subX), use = "pairwise.complete.obs", method = method)
         if (any(is.na(cc))) return(subX)
         select_corr = caret::findCorrelation(cc, cutoff = cutoff, exact = FALSE)
         #print(length(select_corr))
         
         if (length(select_corr)>0){ 
           subX = subX[-select_corr, ]
         }
         
         return(subX)
           
     }
    
    stopCluster(cl)
    
    print(names(filtered_X))
    #    if (length(unique(filtered_X[,1]))>1) cat('PROBLEMA!')
    no_removed = nrow(X)-nrow(filtered_X)
    cat('Removed = ', no_removed, '\n')
    cat('dim filetered_X =', dim(filtered_X), '\n')

    X = filtered_X[sample(nrow(filtered_X)), ]
    
    niter = niter + 1
    cat("nrow(X) =", dim(X), '\n')
    if (no_removed ==0) break;
  }
  cat('final dimension = ', dim(t(X)), '\n')
  
  
  
  return(t(X))
}

