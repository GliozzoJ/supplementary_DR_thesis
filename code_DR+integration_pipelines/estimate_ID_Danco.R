library(doParallel)
estimate_ID_danco <- function(mat_data, 
                              maxit =  11, maxD = 100,
                              list_k = c(4,8,16,24,32,64), perc_points = 0.95, niter = 1, verbose = 0, ncores = my_detectCores()){
  #mat_data ? una matrice con features sulle colonne e sono min-max normalized 
  
  # per ogni k itero niter volte e aggrego le niter iterazioni con fagg
  if ((perc_points >0.999) & (niter ==1) ) niter = maxit
  niter = 1
  
  ids = NULL
  sd_ids = NULL
  npoints = round(nrow(mat_data)*perc_points)
  for (k in list_k){
    idsk = NULL
    
    for (extit in 1:maxit){
        if (perc_points<0.999){ 
        mm = mat_data[sample(nrow(mat_data), size = npoints, replace = FALSE ), ]
      }else{mm = mat_data}
      if (verbose ==1) cat('ID estimation with num neigh =  ', k, '\n')
      idsk_on_mm = NULL
      
      cl <- makePSOCKcluster(min(ncores, detectCores()-1))
      registerDoParallel(cl)
      
      idsk_on_mm = foreach(it = 1:niter, .combine='c', .packages = c("intrinsicDimension") ) %dopar% {
      #for (it in 1:niter){
        ID_danco <- tryCatch({
            ID_danco = globalIntrinsicDimension(.data = mm, .method = 'dancoDimEst', k = k, D = maxD)
          },
          error = function(cond){
            message("Error:")
            message(cond)
            ID_danco = list('dim.est' = 0)
          }
        )    
        return(ID_danco$dim.est)
      }
      
      idsk = c(idsk,mean(idsk_on_mm))
      stopCluster(cl)
    }
    ids = c(ids,  mean(idsk))
    sd_ids = c(sd_ids, sd(idsk)) # standard deviation
    if (verbose ==1) cat('k = ', k,  '*** estimate= ', ids[length(ids)], '\n')
  }
  

      
  #cat('danco ID = ', ID_danco$dim.est, '\n')
  #cat('Mind MLi and MLk IDs= ', ID_MLi$dim.est, ' - ', ID_MLk$dim.est, '\n')
  ID_danco = mean(ids)
  
  if (verbose ==1) cat('ID estimates with danco = ', ids, '\n --> estimated ID = ', ID_danco, '\n')
  return(list(id = ids, sd_id = sd_ids))
}