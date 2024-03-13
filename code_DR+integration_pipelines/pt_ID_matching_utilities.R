match_pt_ID <-function(mat = NULL, df_pt = df_pt, ID_col = 'ID',  other_cols_from_df = NULL){
  
  # prende mat e sostituisce all'ID in id_col_in_mat, l'ID in df_pt
  
  df_data = data.frame(mat)
  df_data[[ID_col]] = row.names(df_data)
  
  
  pt_id_index_in_pheno = match(df_pt[[ID_col]], df_data[[ID_col]])
  
  cat('matching row.names(data) to data_ID in df_pt - unmatched names =: \n')
  
  print(df_pt[[ID_col]][is.na(pt_id_index_in_pheno)])
  
  
  
  df_data_merge = merge(df_data, df_pt[, names(df_pt) %in% c(ID_col, other_cols_from_df)], 
                        by = ID_col, all.x = FALSE, all.y = FALSE, sort = TRUE)
  rownames(df_data_merge) = df_data_merge[[ID_col]]

  return(df_data_merge)
}


