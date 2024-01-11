library(doParallel)

source("estimate_ID_twoNN.R")
source("compute_distance.R")

#' Blocking-id estimate
#'
#' @param mat matrix. Matrix (samples x features) representing a single
#' data view in a multi-modal/multi-omic dataset.
#' @param ID_orig numeric. Unbiased original estimate of the intrinsic
#' dimensionality for the considered view 
#' (i.e. $\hat{d}_{twonn}\left(\mathbf{X}\right)$)
#' @param factor numeric. Factor used for the computation of LO (see below).
#' @param max_blocking_runs numeric. Maximum number of blocks of increasing 
#' size to evaluate (def. 51).
#' @param L0 numeric. Dimension of the smallest block. Default value is
#' round(ID_orig * factor)-
#' @param ID_estimator_fun string. Function to estimate ID, default is 
#' "estimate_ID_twonn".
#' @param ntry numeric. Number of resampling for each block (def. 31).
#' @param args_ID list. List of parameters to use by 'estimate_ID_twonn'.
#' In our experiments we used:
#' args_ID= list(dist_fun_twoNN = 'canberra', perc_points = 0.9, maxit = 11, 
#' ncores = min(ncores, detectCores()-1)).
#' @param task character. Name of the considered predictive task (def. NULL).
#' @param str_desc character. Descriptive name for the evaluated view 
#' (def. 'data').
#' @param ncores numeric. Save number of cores used, which is 
#' ncores = min(ncores, detectCores()-1).
#' @param no_els_to_monitor numeric. Minimum number of steps/blocks to evaluate.
#' @param nDigits numeric. Number of digits to save id estimates results (def. 2).
#' @param verbose numeric. If verbose=1, print informative messages (def. 0).
#'
#' @return Matrix containing the estimated IDs and their variance for each block. 
#' @export
blocking_ID <- function(mat = NULL, ID_orig = NULL, # ID_orig MUST BE PROVIDED
                        factor = 3,
                        max_blocking_runs = 51,
                        L0 = round(ID_orig*factor), 
                        ID_estimator_fun = 'estimate_ID_twonn', 
                        ntry = 31, args_ID = NULL, 
                        task = NULL, str_desc ='data', 
                        ncores = 8,
                        no_els_to_monitor = 10, 
                        nDigits = 2, 
                        verbose = 0) {
  
  
  # N = number of blocks of increasing size I'm going to evaluate
  # maxDim is the maximum blocking Dimension; that is, to speed up we will 
  # analyze blocks with dimension <= L0, ..., maxDim
  N = max_blocking_runs
  maxDim = min(ncol(mat), max(ID_orig*75, nrow(mat)*75), L0*N) # max block dim can't be greater than number of columns! 
  step_inc = L0
  
  # if there are less than 2*no_els_to_monitor steps of increment, then recompute step_inc to have them 
  if (round((maxDim-L0)/step_inc)<(no_els_to_monitor*2)) step_inc = max(round(ID_orig/2), floor((maxDim-L0)/(no_els_to_monitor*2))) # increment should at least be > ID_orig/2
  cat('Blocking method on', str_desc, 'in', task, '\n')
  cat('[L0, + N x step, maxDim] = [', L0, ', +', step_inc, ', ' , maxDim,']\n')
  
  
  # for each block, id_est is a matrix containing, the mean of the id_estimates 
  # over that block, and the standard error (= sd/sqrt(ntry)) 
  colnames_id_mat =   c('no_vars_in_block', 'ntry', 'block_id', 'var_block_id', 'sd_block_id',  
                        'mean_id', 'within_var', 'within_sd', 
                        'between_var', 'between_sd',
                        'pooled_var', 'pooled_sd', 'total_var', 'total_sd')
  
  id_est = matrix(0, nrow = 0, ncol = length(colnames_id_mat))
  colnames(id_est) = colnames_id_mat
  
  
  list_block_dims = unique(c(seq(L0, maxDim, by = step_inc), ncol(mat)))
  stop_at = NA
  stop_at_Lj = NA
  start_count_last = FALSE
  stop_after = 10
  
  # for each block
  for(Lj in list_block_dims) {
      #cat('Block ', i, '\n')
      block_ids = NULL
      # estimate all the ntry ids
      args_ID[['ncores']] = (min(ncores, detectCores()-1))
      for (i in 1:ntry){  
        # randomly sample Lj columns 
          elem = sample(ncol(mat), Lj, replace = FALSE)
        # create B_j(i), that is the i-th block with size Lj   
          Bj_i = mat[, elem]
          args_ID[['mat_data']] <- Bj_i
        # compute id on B_j(i)
          id = do.call(ID_estimator_fun, args_ID)
          id_res = ifelse(is.finite(id[['id']]), id[['id']], -1)
        #return(id_res)
          block_ids = c(block_ids, id_res)
      }
      
      block_ids = block_ids[block_ids>0]
      id_est_Bj = mean(block_ids, na.rm = TRUE)
      var_id_estBj = ifelse(length(block_ids)>1, var(block_ids, na.rm = TRUE),0)
      sd_id_estBj = ifelse(length(block_ids)>1, sd(block_ids, na.rm = TRUE),0)
  
      
      mean_est = mean(c(id_est[, 'block_id'], id_est_Bj), na.rm = TRUE)
      within_var = mean(c(id_est[, 'within_var'], var_id_estBj), na.rm = TRUE)
      within_sd = mean(c(id_est[, 'within_sd'], sd_id_estBj), na.rm = TRUE)
      
      between_var = ifelse(nrow(id_est)>0, var(c(id_est[, 'block_id'], id_est_Bj), na.rm = TRUE),0)
      between_sd = ifelse(nrow(id_est)>0, sd(c(id_est[, 'block_id'], id_est_Bj), na.rm = TRUE),0)
      
      pooled_var = mean(between_var, within_var)
      pooled_sd = mean(between_sd, within_sd)
      
      tot_var = between_var + within_var # law of total variance
      tot_sd = sqrt(tot_var) 
      
      id_est = rbind(id_est, cbind('no_vars_in_block' = Lj, 'ntry' = ntry, 'block-id' = round(id_est_Bj, nDigits), 
                                              'var_block_id' = round(var_id_estBj, nDigits), 'sd_block_id' = round(sd_id_estBj, nDigits), 
                                              'mean_id' = round(mean_est, nDigits), 
                                              'within_var' = round(within_var, nDigits), 'within_sd' = round(within_sd, nDigits), 
                                              'between_var' = round(between_var, nDigits),   'between_sd' = round(between_sd, nDigits), 
                                              'pooled_var' = round(pooled_var, nDigits),  'pooled_sd' = round(pooled_sd, nDigits), 
                                              'total_var' = round(tot_var, nDigits), 'total_sd' = round(tot_sd, nDigits)))
      
      
  
      if (verbose ==1) cat(paste(colnames(id_est), collapse = '-'), '\n', paste(id_est[nrow(id_est), ], collapse='\t'), '\n')
    
      if ((nrow(id_est)>no_els_to_monitor) & (!start_count_last)){
          # compute the std of the last 5 elements. When the std falls below 0.1 
          # stop computation and take the last block id as the last dataset id
          last_els = id_est[(nrow(id_est)-no_els_to_monitor):nrow(id_est),'block_id']
          std_last_els = sd(last_els, na.rm = TRUE)
          lag = 3
          sum_deriv = sum(diff(last_els, lag = lag)/lag)
          if (verbose == 1) cat('std_last_els = ', std_last_els, ' - sum_deriv = ' , sum_deriv, '\n')
          if ((std_last_els < 0.25) | (sum_deriv < 0.05)){ 
              if (verbose == 1)  cat('convergence reached after ', nrow(id_est), ' steps, that is when the block has dimension ', id_est[nrow(id_est), 'no_vars_in_block'], '\n')
              stop_at = nrow(id_est)
              stop_at_Lj = id_est[nrow(id_est), 'no_vars_in_block']
              
              # start counting last ten iterations and then stop computing
              start_count_last = TRUE 
          }
      }
    
    # if you are counting the last iterations then decrease this counter
    if (start_count_last){ 
      stop_after = stop_after-1
      if (verbose == 1) cat(stop_after, ' iterations left\n')
    }
    if (stop_after==0) break
    
  } #end for each block
  

  return(id_est)
}


