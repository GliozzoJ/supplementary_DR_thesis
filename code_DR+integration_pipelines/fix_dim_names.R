fix_dim_names <- function(mat_data = NULL, df_pt = NULL,  str_desc = 'data', ID_col = NULL){
  
  #mat_data è già la matrice ribaltata con features sulle colonne!
  
  source('./pt_ID_matching_utilities.R')
  
  
  
  # normalizza per poter filtrare per sd
  #qui df_miRna diventa una matrice normalizzata con min-max normalization e 
  # le features sono sulle colonne!!
  
  # if (!is.null(normalize_f)){  
  #   cat('norm with\n', normalize_f)
  #   mat_data = do.call(normalize_f, list(df_data = df_data))
  #   
  # }
  
  if (str_desc == 'SNP'){
    
    tronca = function(x) { 
      res = str_locate(x, "_")
      x_cut = str_sub(x, start = (res[1]+1)) 
      return(x_cut)
    }

    rnames = gsub('_', '',tronca(gsub('X', '',row.names(mat_data))) )
    row.names(mat_data) = rnames
  }

  cat('fixing patient ids...\n')
  rnames = gsub('X', '', gsub('-', '_', gsub('\\.','_',row.names(mat_data))))
  row.names(mat_data) = rnames
  
  cat('matching patient ids in data type to patient ID in patient file...\n')  
  names_in_data = df_pt[[ID_col[[str_desc]]]]
  names_pts = row.names(df_pt)
  idx_existent = (names_in_data != "NO") & (!is.na(names_in_data))
  
  names_in_data = names_in_data[idx_existent]
  names_pts = names_pts[idx_existent]
  
  expected_no_elements = length(names_pts)
  
  pt_id_index_in_pheno = match(row.names(mat_data), names_in_data)

  # così elimini quelli che non vengono dallo stesso centro!
  mat_data = mat_data[!(is.na(pt_id_index_in_pheno)), ]
  
  row.names(mat_data) = names_pts[pt_id_index_in_pheno[!is.na(pt_id_index_in_pheno)]]
  if (nrow(mat_data)<expected_no_elements){ cat("WARNING!\nWARNING!\nWARNING!\n patients removed because names were unmatching!!!\n")}
  
  return(mat_data)
}
