filter_cols_by_p <- function(mat_data = NULL, df_p_vals = NULL, fracVar = 0.5, p_thresh = NA){
    # mat_data è la matrice con i dati (features in colonna) 
  
    if (!(is.infinite(fracVar))){
      df_p_vals  = df_p_vals[rownames(df_p_vals) %in% colnames(mat_data), ]
      
      if (fracVar < 1){
        novars = round(ncol(mat_data)*fracVar)
      }else{    
        novars = fracVar 
      }
      
      
      if (is.finite(p_thresh))
        df_p_vals = df_p_vals[df_p_vals[['padj']]< p_thresh, ]
      
      all_ps = df_p_vals[['padj']] 
      names(all_ps) = row.names(df_p_vals)
      
      #all_ps = all_ps[all_ps<thr_on_p]
      all_ps = sort(all_ps)
      
      
      idx_match = match(names(all_ps), colnames(mat_data))
      cat('no unmatched pvalues = ', sum(is.na(idx_match)), '\n')
      cat('no matched pvalues = ', sum(is.finite(idx_match)), '\n')
      
      idx_match = idx_match[is.finite(idx_match)]
      
      # resort columns of miRna based on p-values
      mat_data = mat_data[, idx_match[1:min(novars, length(idx_match))]]
      all_ps = all_ps[1:min(novars, length(idx_match))]
      cat("after filtering by p-values the new dimension is: ", dim(mat_data), '\n',
          '(min,max) p-value = (', min(all_ps), ', ', max(all_ps), ')\n')
      
    
    }else{ cat('no selection based on p\n')}
    
    return(mat_data)
}
