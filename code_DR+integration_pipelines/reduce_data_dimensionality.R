reduce_data_dimensionality <- function(mat_data = NULL, str_desc = 'data', df_pt = NULL, data_path = NULL, dim_red = dim_red, 
                         minD = 100, best_red = NULL,
                         dist_fun_twoNN = 'euclidean', ID_col = NULL, color_col = NULL,
                         target_id = NULL, ID_estimator_fun = 'estimate_ID_danco', args_ID = NULL, verbose = 1){
  
  library(rsvd)  
  if ((!is.null(dim_red)) & (ncol(mat_data)>minD)){
      
      if (is.null(best_red)){
        if (verbose == 1) cat('dimensionality reduction - choosing between ', names(dim_red), '\n ')
      }else{
        if (verbose == 1) cat('dimensionality reduction with ', best_red, '- evaluationg also ', 
            names(dim_red)[!(names(dim_red) %in% best_red)], '\n ')
      }
      red_data = list()
      IDs = NULL

      for (algo in names(dim_red)){
        if (verbose ==1) cat('*****', algo, '***** \n')
        
        red_data[[algo]] = do.call(dim_red[[algo]], list(x= mat_data, k = round(minD)))
        if (is.null(row.names(red_data[[algo]]))){row.names(red_data[[algo]]) = row.names(mat_data)}
        
        args_ID[['mat_data']] = red_data[[algo]]
        id_est = do.call(ID_estimator_fun, args_ID)
        
        IDs = c(IDs, mean(id_est[['id']]))
      
        # if (!is.null(data_path)){
        #    
        #      plot_mat = red_data[[algo]]
        #      algostr = algo
        #    
        #   do_scatter(mat_data = plot_mat, df_pt = df_pt, ID_col_in_pt = ID_col[['pt']], outcome_col = color_col, 
        #              str_desc = str_desc, task ='red', save_path = data_path, str_par = algo, 
        #              str_add_title =
        #                paste('DR with ', algostr, ' - Danco ID = ', round(id_est[['id']],2), ' (+/-', 
        #                      round(id_est[['sd_id']],2),')' , sep =''))
        #     
        # }      
      }
      
      names(IDs) = names(red_data)
      
      if (verbose ==1) {
        cat('estimated IDs = \n')
        print(IDs)
      }
      
      if (is.null(best_red)){
        if(is.null(target_id)) {
          dists = abs(IDs-ID_estimate)
        }
        else {
          dists = abs(IDs-target_id)
        }
          
        gt_0 = which(dists>0)
        if (length(gt_0)>0){
          dists = dists[gt_0]
          IDs = IDs[names(dists)]
          best_red = names(dists)[which(dists == min(dists))]
        }else{
          dists = abs(dists)
          best_red = names(IDs)[which(IDs == max(IDs))]
        }
      }
      
      projected_data = red_data[[best_red]]
      ID_estimate = IDs[best_red]
      cat(str_desc, ':', best_red, 'chosen for dimensionality reduction (ID = ',ID_estimate,'\n')
   }else{
     projected_data = mat_data
     best_red = 'identity'
   }
     
   return(list(mat_data = projected_data, best_red = best_red)) 
}