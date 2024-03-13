create_mofa_embedding <- function(df_outcome = NULL, df_pt_use = NULL, 
                                  embedded_data_types = NULL, 
                                  use_types = names(embedded_data_types), 
                                  embed_pt = FALSE,
                                  no_mofa_factors = 5, # default
                                  sparse_weights = TRUE, # default
                                  sparse_factors = FALSE, # default
                                  seed = 3, 
                                  python_path = python_path, 
                                  save_path_task = './', res_dir = './', force_recreation = FALSE){

  list_data_types = list()
  #reticulate::use_python(python = python_path)
  all_pts = row.names(df_outcome)
  if ((length(use_types)==1) & (use_types[1] == 'pt')){
    df_use = merge(df_outcome, df_pt_use,  by = 'row.names')
    row.names(df_use) = df_use[["Row.names"]]
    df_use = df_use[, !(names(df_use) %in% "Row.names")]
    return(df_use)    
  }else{
    
    for (ut in use_types){
  
      if (ut !='pt') {
        if (!(is.data.frame(embedded_data_types[[ut]]))) {
          load(embedded_data_types[[ut]])
        }else mat_data = embedded_data_types[[ut]]
        
        mat_data_all = matrix(NA, nrow = length(all_pts), ncol=ncol(mat_data))
        row.names(mat_data_all) = all_pts
        colnames(mat_data_all) = colnames(mat_data)
        idx_names = match(row.names(mat_data), all_pts)
        mat_data_all[idx_names[!is.na(idx_names)], ] = mat_data[!(is.na(idx_names)), ]
        
        mat_data_all = t(mat_data_all)
        
        list_data_types[[ut]] = mat_data_all # FOR MOFA, samples are stored in columns and features in rows of the matrix
        print(dim(list_data_types[[ut]]))
      }else if (embed_pt){
        df_pt_use = df_pt_use[match(all_pts, row.names(df_pt_use)), ]
        df_pt_use_mat = do.call(rbind, lapply(df_pt_use, as.numeric))
        colnames(df_pt_use_mat) = row.names(df_pt_use)
        row.names(df_pt_use_mat) = names(df_pt_use)
        
        if (any(row.names(df_pt_use_mat) %in% 
                c('years_to_birth', 'patient.age_at_initial_pathologic_diagnosis'))){
            int_cols = df_pt_use_mat[ row.names(df_pt_use_mat) %in% 
                                    c('years_to_birth', 'patient.age_at_initial_pathologic_diagnosis'), ]
            list_data_types[['pt_int']] = int_cols
        }
        
        if (any(!(row.names(df_pt_use_mat) %in% 
                  c('years_to_birth', 'patient.age_at_initial_pathologic_diagnosis')))){
            factor_cols = df_pt_use_mat[ !(row.names(df_pt_use_mat) %in% 
                                    c('years_to_birth', 'patient.age_at_initial_pathologic_diagnosis')), ]
            list_data_types[['pt_factor']] = factor_cols
        }
        # mat_data_all = matrix(NA, nrow = length(all_pts), ncol=3)
        # row.names(mat_data_all) = all_pts
        # idx_names = match(row.names(df_pt_use), all_pts)
        # mat_data_all[idx_names[!is.na(idx_names)], ] = df_pt_use[!is.na(idx_names), c('race', 'ethnicity', 'gender')]
        # 

      }
  
    }
  
    names_aligned_maybe = do.call(rbind, lapply(list_data_types, colnames))
  
    data_str =   paste(use_types, collapse = '_')
  
    mofa_outfile = file.path(save_path_task, paste('MOFA', ifelse(embed_pt, 'PT', ''), '_', data_str, '.Rda', sep = '') )
  
  
#     if (file.exists(mofa_outfile)){
#         load(mofa_outfile)
#         Z = mofaobj@expectations$Z # latent factor matrix
#         Z = Z$group1
#         W = mofaobj@expectations$W # individual weights
# 		    
#         if (any(!is.finite(Z))){ 
# 			    force_recreation = TRUE 
# 		    }else{
# 			    factor_mean = colMeans(Z, na.rm = TRUE)
# 			    factor_std = apply(Z,2, sd, na.rm = TRUE)
# 			    if (sum(factor_std> 10^(-20)) < 2) force_recreation = TRUE
# 		    }      
#     } 
    
    if ((force_recreation) | (!file.exists(mofa_outfile))){
  
      mofaobj = create_mofa_from_matrix(list_data_types)
  
      data_opts <- get_default_data_options( mofaobj  )
      #data_opts$scale_views = TRUE
  
      model_opts <- get_default_model_options( mofaobj  )
      model_opts$num_factors <- no_mofa_factors
      model_opts$spikeslab_factors = sparse_factors
      model_opts$spikeslab_weights = sparse_weights
      model_opts$ard_factors = sparse_factors
      model_opts$ard_weights = sparse_weights
  
      train_opts <- get_default_training_options( mofaobj  )
      train_opts$seed <- seed
     # train_opts$convergence_mode <- "slow"
  
      mofaobj <- prepare_mofa( mofaobj,
                                  training_options = train_opts,
                                  model_options = model_opts,
                                  data_options = data_opts )
      mofaobj <- run_mofa( mofaobj,
                              outfile = file.path(res_dir, 'MOFA.hdf5'), save_data = FALSE, use_basilisk = TRUE)
  
  
      Z = mofaobj@expectations$Z # latent factor matrix
      Z = Z$group1
      W = mofaobj@expectations$W # individual weights
  
      save( mofaobj , file = mofa_outfile )
      
      # sometimes it returns all Inf or all NA values 
      if (any(!is.finite(Z))) return(NULL)
      # sometimes it returns a Z matrix with all values set to 0 
      factor_mean = colMeans(Z, na.rm = TRUE)
      factor_std = apply(Z,2, sd, na.rm = TRUE)
      if (sum(factor_std> 10^(-20)) < 2) return(NULL)

      
    }else{
        cat('loading MOFA outfile: ', mofa_outfile, '\n')
        load(mofa_outfile)
        Z = mofaobj@expectations$Z # latent factor matrix
        Z = Z$group1
        W = mofaobj@expectations$W # individual weights
        
        # sometimes it returns all Inf or all NA values 
        if (any(!is.finite(Z))) return(NULL)
        # sometimes it returns a Z matrix with all values set to 0 
        factor_mean = colMeans(Z, na.rm = TRUE)
        factor_std = apply(Z,2, sd, na.rm = TRUE)
        if (sum(factor_std> 10^(-20)) < 2) return(NULL)
      
    }
    
    df_emb = data.frame(Z)
    if ('pt' %in% use_types){
      # se voglio viste complete comunque ho giá risolto togliendo a priori i pz a cui mancava una vista all.x / all.y = TRUE
      df_emb = merge(df_emb, df_pt_use,  by = 'row.names', all.x = TRUE, all.y = TRUE)
      row.names(df_emb) = df_emb[["Row.names"]]
      df_emb = df_emb[, !(names(df_emb) %in% "Row.names")]
    }
    # se voglio viste complete comunque ho giá risolto togliendo a priori i pz a cui mancava una vista all.x / all.y = TRUE
    df_emb =  merge(df_emb, df_outcome,  by = 'row.names', all.x = TRUE, all.y = TRUE) 
    row.names(df_emb) = df_emb[["Row.names"]]
    df_emb = df_emb[, !(names(df_emb) %in% "Row.names")]
  
    return(df_emb)
  }
}