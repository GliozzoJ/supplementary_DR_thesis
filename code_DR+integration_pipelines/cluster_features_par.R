library(doParallel)
library(genieclust)


#' Feature selection using clustering
#' 
#' @description This function performs feature selection exploiting the clustering
#' algorithm Genie. In practice, it considers subsets of features and applies Genie
#' to select the medoids of the considered clusters as features. After all subset
#' of features are processed, the procedure is repeated maxit times on the selected
#' features at the previous iteration. Thus, at each iteration the set of selected
#' features is refined. 
#'
#' @param df dataframe. Dataframe (examples x features).
#' @param dim_split integer. Number of features considered at each iteration 
#' (def. 1000).
#' @param maxD integer. Number of features estimated using the blocking ID method.
#' @param M integer. Smoothing factor for Genie algorithm. Note that for each
#' considered subset of features M can be set to min(M, round(nrow(subX)/2)).
#' @param maxit integer. Maximum number of iterations (def. 10).
#' @param method string. Method used to compute correlation. Options are
#' c("pearson", "kendall", "spearman"). Default "kendall".
#'
#' @return Matrix (samples x features) with the selected set of variables.
#' @export
cluster_features_par <- function(df, dim_split = 1000, 
                             maxD = 250,  M = 10, maxit = 3, method = "pearson", ncores = my_detectCores()){
  X = t(df)
  cat("CLUSTERING FEATURES\n")
  cat("maxD = ", maxD, "\n")

  nc_per_split =  (maxD/nrow(X))^(1/maxit)
  
  print(nc_per_split)
  
  iter = 0
  while((nrow(X)>maxD) & (iter<maxit)){
    iter = iter+1
    cat("ITERATION=", iter, "\n")
    cat("dim_split =", dim_split, "nc_per_split=", nc_per_split , "\n")
    
    cl <- makePSOCKcluster(min(ncores, detectCores()-1))
    registerDoParallel(cl)
    clustered_X = NULL
    
    clustered_X = foreach(nR = seq(1, nrow(X), by=dim_split), .combine='rbind', .packages = c("genieclust")) %dopar% {
      if ((min(nrow(X), (nR+dim_split-1))-nR)>0){
        subX = X[nR:min(nrow(X), (nR+dim_split-1)), ]
        
        if (nrow(subX)>(M*2)){
          #cat('dimensions of subX = ', dim(subX), '\n')
          k = round(nrow(subX)*nc_per_split)
          if ((k > 1) & (k < nrow(subX))){
            cc = cor(t(subX), use = "pairwise.complete.obs", method = method)
            dist_mat = 1- abs(cc)
            
            h <- gclust(dist_mat, M = min(M, round(nrow(subX)/2)))
            res = cutree(h, k)
            
            cluster_centroids = NULL
            rnames = NULL
            for (nc in min(res):max(res)){
              idx_cluster = which(res == nc)
              #new_feats[nc] = paste(features_in_cluster[idx_cluster], collapse = ' + ')
              if (length(idx_cluster)>1){
                centroidval = colMeans(subX[idx_cluster, ], na.rm = TRUE)
                xxx = rbind(centroidval,subX[idx_cluster, ])
                #prendi solo la prima riga per avere tutte e sole le correlazioni con il centroide
                xcor = cor(t(xxx), use = "pairwise.complete.obs")[1, ]
                #elimina la correlazione del centroide con se stesso
                idx_best_feat = which.max(abs(xcor[-1]))
                idx_cluster = idx_cluster[idx_best_feat]
                rnames = c(rnames, idx_cluster)
              }else{
                # TENGO O BUTTO?? per ora tengo
                #if (keep_isolated){
                rnames = c(rnames, idx_cluster)
                #}
              }
            }
            print(rnames)
            cluster_centroids = subX[rnames, ]
            rownames(cluster_centroids) = names(rnames)
          }else{
            #cat('no clusters = ', k, 'too little; selecting all elements\n')
            cluster_centroids = subX
          }
          
        }else{
          cluster_centroids =  subX
        }
      }else{
        cluster_centroids = NULL
      }
      return(cluster_centroids)
    }
    
    X = clustered_X[sample(nrow(clustered_X)), ]
    clustered_X = NULL
    
    cat("dim di X", dim(X), "\n")
    stopCluster(cl)
    unregister_dopar()
    
  }
  # cat('Too many features, (' , ncol(projected_data), ') I will cluster them!\n')
  # h <- gclust(emst_mlpack(X), distance = 'cosine')
  # #h <- gclust(X, distance = 'cosine')
  # 
  
  # 
  # cat('with clustering I have ', ncol(projected_data),  ' clustered features\n') 
   
   X = t(X)
   cat('final dimension = ', dim(X), '\n')
   
  
   return(X)
}