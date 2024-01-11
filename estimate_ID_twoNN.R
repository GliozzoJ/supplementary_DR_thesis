library(doParallel)

#' Computation of an unbiased two-nn ID estimate
#'
#' @description Computation of an unbiased two-nn ID estimate by averaging all the 
#' two-nn ID estimates computed on "maxit" under-sampled versions of the dataset, 
#' where the under-sampling randomly selects (with repetition) a fixed 
#' percentage, "perc_points", of the dataset points.
#' 
#' The two-nn ID estimator is fully described in the paper:
#' 
#' Facco, Elena, et al. "Estimating the intrinsic dimension of datasets by a 
#' minimal neighborhood information." Scientific reports 7.1 (2017): 12140.
#' 
#' @param mat_data matrix. Matrix (samples x features).
#' @param dist_fun_twoNN matrix. Distance matrix for twonn.
#' @param maxit numeric. Maximum number of under-sampled versions of the 
#' dataset (def. 11).
#' @param perc_points numeric. Percentage of data points to select in each
#' undersampled version of the dataset (def. 0.9).
#' @param verbose numeric. If verbose=1, print informative messages (def. 0).
#' @param ncores numeric. Number of cores used, set as 
#' ncores = min(ncores, detectCores()-1).
#'
#' @return List with two components:
#' - "id": estimated ID.
#' - "sd_id": standard deviation of the estimated ID.
#' @export
estimate_ID_twoNN <- function(mat_data = NULL , dist_fun_twoNN = 'canberra', 
                              maxit = 11, 
                              perc_points = 0.9, verbose = 0, ncores = 8){
  #mat_data e' una matrice con features sulle colonne
  ID_twoNN = NULL
  sd_twoNN = NULL
  #if (perc_points > 0.9999) maxit = 1
  cl <- makePSOCKcluster(min(ncores, detectCores()-1))
  registerDoParallel(cl)
  
  ID_twoNN = foreach(nit = 1:maxit, .combine='c', .packages = c("intRinsic"), .export=c("compute_distance") ) %dopar% {
  #for(nit in 1:maxit){
    cat('.')
    if (perc_points < 0.9999){
      sub_mat_data = mat_data[sample(nrow(mat_data), size = round(nrow(mat_data)*perc_points), replace = FALSE), ]  
    }else{
      sub_mat_data = mat_data
    }
    
    dm = compute_distance(sub_mat_data, dist_fun_twoNN)
    
    if (verbose == 1 & nit ==1) cat('distance mat with shape rxc = ', nrow(dm), ' x ', ncol(dm), '\n')
    id_twonn = twonn(dist_mat = dm, method = "mle")
    
    id_twonn_est = id_twonn$est[2]

    ID_twoNN = c(ID_twoNN, id_twonn_est)
  }  
  
  stopCluster(cl)
  
  sd_ID = sd(ID_twoNN)
  ID_est = mean(ID_twoNN)
 
  if (verbose == 1){ 
    cat('\nestimated ID on ', maxit, ' iterations on the ', round(perc_points*100), '% of points = \n', 
      ID_est, '[', -sd_ID,', ', sd_ID, '] \n')
  }
  return(list(id = ID_est, sd_id = sd_ID))
  
}
