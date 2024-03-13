library(doParallel)
library(cytominer)
#' Feature selection using clustering
#' 
#' @description This function performs unsupervised feature selection exploiting the svd entropy algorithm
#' In practice, it considers subsets of features and applies svd entropy
#' to select the features that most contribute to the dataset entropy. After all subset
#' of features are processed, the procedure is repeated maxit times on the selected
#' features at the previous iteration. Thus, at each iteration the set of selected
#' features is refined. 
#'
#' @param df dataframe. Dataframe (examples x features).
#' @param dim_split integer. Number of features considered at each iteration 
#' (def. 1000).
#' @param maxD integer. Number of features estimated using the blocking ID method.
#' @param M integer. minimum split dimension for applying the svd entropy selection 
#' @param maxit integer. Maximum number of iterations (def. 10).
#'
#' @return Matrix (samples x features) with the selected set of variables.
#' @export
hierarchical_svd_entropy_par <- function(df, dim_split = 1000, 
                                         M = (min(dim(df))-1),
                                         maxD =  round((min(dim(df))-1)/2), maxit = 3, ncores = my_detectCores()){
  X = df
  cat("CLUSTERING FEATURES\n")
  cat("maxD = ", maxD, "\n")
  
  nc_per_split =  (maxD/ncol(X))^(1/maxit)
  
  print(nc_per_split)
  
  iter = 0
  while((ncol(X)>maxD) & (iter<maxit)){
    iter = iter+1
    cat("ITERATION=", iter, "\n")
    cat("dim_split =", dim_split, "nc_per_split=", nc_per_split , "\n")
    
    # cl <- makePSOCKcluster(ncores)
    # registerDoParallel(cl)
    selected_X = NULL
    
  #  selected_X = foreach(nC = seq(1, ncol(X), by=dim_split), .combine='cbind', .packages = c("cytominer")) %dopar% {
    for (nC in seq(1, ncol(X), by=dim_split)){
      if ((min(ncol(X), (nC+dim_split-1))-nC)>M){
        subX = X[, nC:min(ncol(X), (nC+dim_split-1))]
        
        
          #cat('dimensions of subX = ', dim(subX), '\n')
          k = round(ncol(subX)*nc_per_split)
          if ((k > 1) & (k < ncol(subX))){
            
            cat("svd_entropy selection\n")
            tbl <- tibble::as_tibble(as.data.frame(subX))
            entropy <- cytominer::svd_entropy(colnames(subX), tbl,  detectCores()-1)
            entropy <- entropy[order(entropy$svd_entropy, decreasing = T), ]
            
            fs_names <- entropy$variable[1:k]
            sel_sub <- subX[, colnames(subX) %in% fs_names]
            
          }else{
            #cat('no clusters = ', k, 'too little; selecting all elements\n')
            sel_sub = subX
          }
          
        
      }else{
        sel_sub = subX
      }
      if (nC == 1) selected_X = sel_sub
      else selected_X = cbind(selected_X, sel_sub)
    }
    
    X = selected_X[,sample(ncol(selected_X))]
    selected_X = NULL
    
    cat("dim di X", dim(X), "\n")
    #stopCluster(cl)
    
  }


  cat('final dimension = ', dim(X), '\n')
   
  
  return(X)
}