library(parallel)
library(impute)

my_detectCores <- function(detect = FALSE){
  return(min(11, detectCores()))
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

select_points_by_name <- function(x, name_list = NULL){
  return(x[match(name_list, rownames(x)), ]) 
} 

open_data_fun <- function(str_desc = NULL, fn = NULL, normalize_fun = standard_scaling , feat_select = NULL, centers = c('CHUT', 'OSR')){
  df_data_all = data.frame()
  for (cc in centers){
    df_data = read.csv(file = fn[grepl(cc,fn, ignore.case = TRUE)], 
                     header = TRUE, sep = '\t')
  
    if (!is.null(var_with_rownames[[str_desc]])){
      rownames(df_data) = df_data[[var_with_rownames[[str_desc]]]]
      df_data = df_data[, !(names(df_data) %in% c(var_with_rownames[[str_desc]],drop_columns[[str_desc]]))]
    }
    if (!(is.null(feat_select))){
      df_data = df_data[feat_select, ]
    }
    if (!is.null(normalize_fun)){
      mat_data = as.matrix(df_data)
      #normalize_fun assumes variables are on columns!
      mat_data = do.call(normalize_fun, list(mat = t(mat_data)))
      df_data = data.frame(t(mat_data))
    }
    if (sum(dim(df_data_all))==0){ df_data_all = df_data 
    }else{ df_data_all = cbind(df_data_all, df_data[match(row.names(df_data), row.names(df_data_all)), ] )}
    
  }
  
  return(df_data_all)
}

my_quantile_normalization <- function(df_data){ 
  mat_data = t(as.matrix(df_data))
  mat_data = imputa(mat_data)
  mat_data = normalize.quantiles(mat_data, copy = TRUE, keep.names = TRUE) 
  return(mat_data)
}



create_pt_df <- function(df_pt = NULL, ID_patient_col = 'ID_Patient', 
                         neg_label = 'EDA', pos_label = 'NEDA', 
                         outcome_col = 'Activity'){
  colNames = names(df_pt)
  df_pt = df_pt[-which(duplicated(df_pt[[ID_patient_col]])),]
  
  cols_ID = colNames[grepl('ID', colNames)]
  for (cc_ID in cols_ID){
    df_pt[[cc_ID]] =  gsub('X', '',gsub('-', '_', gsub('\\.','_',df_pt[[cc_ID]])))
  }
  
  cat('creating new column = ', outcome_col, '\n')
  
  
  cat(names(df_pt), '\n')
  rownames(df_pt) =  df_pt[[ID_patient_col ]]
  df_pt = df_pt[, !(names(df_pt) %in% ID_patient_col)]
  
#  df_pt$SEX = factor(df_pt$SEX)
  df_pt$Sex = factor(df_pt$Sex)
  df_pt$Therapy_pre = factor(df_pt$Therapy_pre=='Yes')
  
  df_pt$AgeAtOnset = NULL #as.numeric(df_pt$AgeAtOnset)
  df_pt$Age_at_t0 = as.numeric(df_pt$Age_at_t0)
  df_pt$DiseaseDuration_at_T0 = as.numeric(df_pt$DiseaseDuration_at_T0)
  df_pt$EDSS_t0 = as.numeric(df_pt$EDSS_t0)
  df_pt$VitaminD = as.numeric(df_pt$VitaminD)
#  df_pt$Therapy_Pre = as.factor(df_pt$Therapy_Pre)
  
  df_pt[[outcome_col]] = factor(df_pt[[outcome_col]], levels = c(neg_label, pos_label))
  
  df_pt[, !(names(df_pt) %in% c(outcome_col, cols_ID))] = imputa(df_pt[, !(names(df_pt) %in% c(outcome_col, cols_ID))])
  
  
  return(df_pt)
}


do_scatter <- function(mat_data = NULL, df_pt = NULL, ID_col_in_pt = NULL, outcome_col = NULL, str_desc = 'data', task ='', 
                       save_path = './', str_par = '', str_add_title = NULL){
  cat('scatterplot using of ', str_desc, ' on ', task, '\n')
  
  mat_data = mat_data[complete.cases(mat_data), ]
  idx = match(rownames(mat_data), rownames(df_pt))
  image_fn = file.path(save_path, paste(str_desc, '_', task, str_par, '_', 
                                        nrow(mat_data), 'x', ncol(mat_data), '.png', sep =''))
  cat('saving scatter in path = ', image_fn, '\n')
  png(filename = image_fn)
  
  if (!is.null(outcome_col)){
    df = data.frame(x = mat_data[, 1], y =  mat_data[, 2], activity = df_pt[[outcome_col]][idx])
    p = ggplot(df,aes(x=x,y=y,col=activity))+geom_point()
  }else{
    df = data.frame(x = mat_data[, 1], y =  mat_data[, 2])
    p = ggplot(df,aes(x=x,y=y))+geom_point()
  }
  p = p+ggtitle(paste(str_desc, str_add_title) )
  print(p)
  
  dev.off()
}



imputa <- function(mat_data, fun = ifelse(ncol(mat_data)<500, 'missRanger', 'knn'), 
                   k = ifelse(fun== 'missRanger', 1, 5), 
                   my_f = as.formula('. ~ .'), num.trees = 101){
  # fun may be knn or missRanger
  
  #mat_data has features on columns
  miss_indicator = apply(apply(mat_data, 2, is.na),2, as.numeric)
  
  if (any(miss_indicator>0)){
    no_pt_with_missing = sum(rowMeans(miss_indicator)>0)
    no_vars_with_missing = sum(colMeans(miss_indicator)>0)
    
    cat('no. cases with missing values = ', no_pt_with_missing, 
        '(', round(100*(no_pt_with_missing/nrow(mat_data)),3), '%) - average (min,max) missingness =', 
        mean(rowMeans(miss_indicator)), ' (',
        min(rowMeans(miss_indicator)),',',max(rowMeans(miss_indicator)),')\n')
    
    cat('no. variables with missing values = ', no_vars_with_missing, 
        '(', round(100*(no_vars_with_missing/ncol(mat_data)),3), '%) - average (min,max) missingness =', 
        mean(colMeans(miss_indicator)), ' (',
        min(colMeans(miss_indicator)),',',max(colMeans(miss_indicator)),')\n')
   
    
    if (fun == 'knn'){
      cat('imputation with knn (k = ', k, ')\n')
      # mat_imp has genes on columns - impute knn wants it on rows
      imputed_list = impute.knn(t(mat_data), k = k) 
      imputed_data = t(imputed_list$data)
      # appico pmm di missranger per pescare valori che già esistevano!
      # nel caso degli SNP almeno faccio rounding
      imputed_data[(is.na(mat_data))] = pmm(xtrain = mat_data[!(is.na(mat_data))], imputed_data[(is.na(mat_data))], ytrain = imputed_data[!(is.na(mat_data))])
    }else{
        if (ncol(mat_data) > 500){ 
          max.depth = 9
          num.trees = 51
        }else{ max.depth = max(round(nrow(df)/4), 31) }
        cat('imputing with ranger settings: max.depth =', max.depth, ' - num.trees = ', num.trees)
        cat('imputation with missRanger (k = ', k, ')\n')
        if (!is.data.frame(mat_data)) df_data = data.frame(mat_data)
        else df_data = mat_data
        names(df_data) = str_replace_all(pattern = '-', replacement = '_', names(df_data))
        imputed_data =  missRanger(df_data, my_f, pmm.k = k, num.trees = num.trees, max.depth = max.depth, num.threads = 31, maxiter = 5)
        if (!is.data.frame(mat_data)) imputed_data = as.matrix(imputed_data)
    
    }
  }else{
    cat('no missing values\n')
    imputed_data = mat_data
  }
  
  return(imputed_data)
  
}

robust_scaler_var <-function(x, probs = c(0.25, 0.75)){ 
  quants = quantile(x, probs = probs, na.rm = TRUE)
  sd_r = max(quants)-min(quants)
  xx = (x-median(x, na.rm = TRUE))/sd_r
  if (any(!is.finite(xx))){
    cat('qualcosa non va\n')
  }
  return(xx)
}


robust_scaler = function(mat = NULL, probs = c(0.25, 0.75)){
  # assumes variables are on columns!
  # mat_norm = apply(mat, 2, 
  #                  )
  mat_norm = mat
  for (cc in 1:ncol(mat_norm)){
    mat_norm[,cc] = robust_scaler_var(mat_norm[,cc], probs = probs)
  } 
  
  return(mat_norm)
}


standard_scaler = function(mat = NULL){
  # assumes variables are on columns!
  mat_norm = apply(mat, 2, function(x){ xx = (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)})
  return(mat_norm)
}


min_max_scaler = function(mat = NULL){
  mat_norm = apply(mat, 2, function(x){ xx = (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
  return(mat_norm)
}





cv <- function(x, ...) {
  # computes the coefficient of variation = sd/mean
  return(sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
}


# 
# del_cols_too_many_NA <- function(df, thr_del = 0.15){
#   # removes columns that have more than the thr_del % of missing values
#   
#   sumNA = apply(apply(apply(df,2,is.na),2,as.numeric),2,sum)
#   print(sumNA[1:10])
#   del_cols = sumNA[sumNA > round(nrow(df)*thr_del)]
#   if (length(del_cols)>0){
#     cat('Removing ', length(del_cols), ' columns with more that the ', thr_del, '% of missing data \n')
#   }else{
#     print(del_cols)
#   }
#   cat('initial number of columns = ', ncol(df), '\n')
#   df = df[, !(names(df) %in% names(del_cols))]
#   cat('number of columns after removing columns with, ', thr_del*100, '% of NA vals= ', ncol(df), '\n')
#   return(df)
#   
# }

del_cols_too_many_NA <- function(mat_data, thr_del = 0.1){
  # removes COLUMNS that have more than the thr_del % of missing values
  
  sumNA = apply(apply(mat_data,2,is.na),2,as.numeric)
  
  del_cols = colMeans(sumNA) > thr_del
  
  if (sum(del_cols)>0){
    cat('Removing ', sum(del_cols), ' features (on columns) with more that the ', thr_del, '% of missing data \n')
  
    cat('initial number of rows = ', ncol(mat_data), '\n')
    mat_data = mat_data[, !(del_cols)]
    cat('number of COLUMNS after removing those with, ', thr_del*100, '% of NA vals = ', ncol(mat_data), '\n')
  }else{
    cat('no features to delete; matrix still has ', nrow(mat_data), ' features\n')
    
  }
  return(mat_data)
  
}


del_cols_low_sd <- function(mat, thr_sd = 0.05, str_source = 'data', save_path = '.\\results\\'){
  # mat has features on columns and it is min-max normalized!
  
  print(dim(mat))
  
  sds_mat = apply(mat,2,sd, na.rm = TRUE)
  png(paste( save_path, 'hist_sds_', str_source, '.png', sep =''))
  hist(sds_mat,breaks = seq(min(sds_mat), max(sds_mat), length.out =20))
  dev.off()
  
  cat('mean - [min, max] sd = ', mean(sds_mat), 
      ' [', min(sds_mat), ' - ' , max(sds_mat),']\n')
  
  del_cols = sds_mat < thr_sd
  if (sum(del_cols)>0){
    cat('Removing ', sum(del_cols), ' features (on columns) with sd < ', thr_sd, '\n')
    
    cat('initial number of  = columns', ncol(mat), '\n')
    mat = mat[, !(del_cols)]
    cat('number of features (on columns) after removing rows with low variability = ', ncol(mat), '\n')
  }else{
    cat('no features to delete; dataframe still has ', ncol(mat), ' features\n')
  }
  return(mat)
}


intelligent_random_p <- function(df,  nsplit = 0.5, nshuffle = 3, ID_col = "Patient_ID", outcome_col = "Activity", 
                          boolean_cols = NULL, integer_cols = NULL, double_cols = NULL, factor_cols = NULL, 
                         verbose =0, use_only_not_NA = FALSE){
  # restituisce un dataframe in cui ogni colonna colntiene il p_value per la corrispondente colonna di df
  # nsplit ? la PROPORZIONE del numero di campioni della MAJORITY class che voglio considerare per ogni random sampling; 
  # per ogni sottocampionamento della MAJORITY class, la MINORITY class 
  # viene allo stesso modo sottocampionata per avere lo split bilanciati
  # se nsplit = NA, allora la proporzione ? calcolata in modo da avere lo stesso numero di campioni della minority class
  # nshuffle = numero di resampling della minority e della majority
  # se nsplit = 1 allora anche viene settato anche nshuffle = 1; e viene calcolato lo standard p usando tutti gli elementi
  
  if (use_only_not_NA){
    df = df[complete.cases(df), ]
  }
  
  df[[outcome_col]] = factor(df[[outcome_col]])

  class_cardinalities = table(df[[outcome_col]])
  minority_class = names(which.min(class_cardinalities))[1]
  majority_class = names(which.max(class_cardinalities))[1]
  
  idx_min = which(df[[outcome_col]] == minority_class)
  idx_max = which(df[[outcome_col]] == majority_class)
  
  df_min = df[idx_min, ]
  df_max = df[idx_max, ]
  # calculate the number of splits
  if (is.na(nsplit)){
    nsplit =  max(floor(nrow(df_min) /(nrow(df_max))),1)
  }
  
  dim_split = floor(nrow(df_max)*nsplit)
  
  if (nsplit==1){ 
    nshuffle = 1
    dim_split = nrow(df_max)
  }
  
  
  cat('prop of split = ', nsplit, ' (dim max = ', dim_split, ')\n')
  
  
  if (is.null(factor_cols) &  is.null(boolean_cols) & is.null(integer_cols) & is.null(double_cols) ){
    cat("no column type provided; use all columns as double columns", "\n")
    double_cols = names(df)[!(names(df) %in% c(ID_col, outcome_col))]
  }
  
  p_values_shuffled = data.frame(matrix(NA, nrow= nshuffle , ncol = length(c(boolean_cols, integer_cols, double_cols))))
  names(p_values_shuffled) = c(boolean_cols, integer_cols, double_cols)
  
  for (nits in 1:nshuffle){
      cat("*******************************", "\n", "start shuffle = ", nits , "\n")

      if (nsplit ==  1){
        df_split = df
      }else{
        df_split = rbind(df_max[sample(class_cardinalities[majority_class][1], dim_split), ], 
                         df_min[sample(class_cardinalities[minority_class][1], min(dim_split, class_cardinalities[minority_class][1])), ])
      } 
    
      if(!(is.null(boolean_cols))){
        for (cc in boolean_cols){
          p_values_shuffled[nits, cc] = cor(as.integer(df_split[[cc]]), as.integer(df_split[[outcome_col]]==minority_class),  method = "pearson", use = "complete.obs")
          #  cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no boolean cols \n')}
      
      if(!(is.null(factor_cols))){
        for (cc in factor_cols){
          res_chi = chisq.test(df_split[[cc]], df_split[[outcome_col]], correct=TRUE)
          p_values_shuffled[nits, cc] = res_chi$p.value
          #cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no factor cols \n')}
      
      if (!(is.null(c(integer_cols, double_cols)))){
        for (cc in c(integer_cols, double_cols)){
          pk = kruskal.test(df_split[[cc]], df_split[[outcome_col]], na.omit = TRUE)
          p_values_shuffled[nits, cc] = pk$p.value
          #  cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no integer or double cols \n')}
      
    

  }
  
  return(colMeans(p_values_shuffled))
  
}




filter_by_p <- function(df_p = NULL, mat = NULL, fracVar = Inf){
  # p ? il df con i pvalue
  # mat ? la matrice di cui prendere le colonne

  # INNANZI TUTTO CONSIDERA SOLO I miRna con pvalue (di ferdinando) sotto una certa soglia
 # df_p = read.csv(file = 'dati_elena_Discovery_rinominati/association_statistics/miRNA_p.txt', 
  #                      header = TRUE, sep = '\t')
  if (!(is.infinite(fracVar))){
    all_ps = df_p[['padj']] 
    names(all_ps) = row.names(df_p)
    all_ps = sort(all_ps)
    
    
    #match returns a vector of the positions of (first) matches of its first argument in its second.
    # x = c(1,2,3,21) 
    # y = c(8,9,7,10,5,4,1,2 ,6,3)
    # match(x,y)
    # 7  8  10  NA
  
    idx_match = match(names(all_ps), colnames(mat))
    cat('no unmatched pvalues = ', sum(is.na(idx_match)), '\n')
    cat('no matched pvalues = ', sum(is.finite(idx_match)), '\n')
    
    idx_match = idx_match[is.finite(idx_match)]
  
    if (fracVar < 1){
      novars = round(ncol(mat)*fracVar)
    }else{    novars = fracVar }
  
    # resort columns of miRna based on p-values
    mat = mat[ , idx_match[1:min(novars, length(idx_match))]]
  }
  return(mat)
}



intelligent_p <- function(df, nshuffle = 3, ID_col = "Patient_ID", outcome_col = "Activity", majority_class = "EDA", minority_class = "NEDA",
                          boolean_cols = NULL, integer_cols = NULL, double_cols = NULL, factor_cols = NULL, 
                          nsplit = NA, verbose =0){
  # restituisce un dataframe in cui ogni colonna colntiene il p_value per la corrispondente colonna di df
  
  
  idx_pos = which(df[[outcome_col]] == minority_class)
  idx_neg = which(df[[outcome_col]] == majority_class)
  
  df_pos = df[idx_pos, ]
  df_neg = df[idx_neg, ]
  # calculate the number of splits
  if (is.na(nsplit)){
    nsplit =  max(floor(nrow(df_neg) /(nrow(df_pos))),1)
  }
  
  dim_split = ceiling(nrow(df_neg)/nsplit)
  
  if (nsplit==1){ 
    nshuffle = 1
    dim_split = nrow(df_neg)
  }
  
  
  cat('number of splits = ', nsplit, ' (dim = ', dim_split, ')\n')
  
  
  if (is.null(factor_cols) &  is.null(boolean_cols) & is.null(integer_cols) & is.null(double_cols) ){
    cat("no column type provided; use all columns as double columns", "\n")
    double_cols = names(df)[!(names(df) %in% c(ID_col, outcome_col))]
  }
  
  p_values_shuffled = data.frame(matrix(NA, nrow= nshuffle , ncol = length(c(boolean_cols, integer_cols, double_cols))))
  names(p_values_shuffled) = c(boolean_cols, integer_cols, double_cols)
  
  for (nits in 1:nshuffle){
    #        cat("*******************************", "\n", "shuffle = ", nits , "\n")
    p_values = data.frame(matrix(NA, nrow= nsplit , ncol = ncol(p_values_shuffled)))
    names(p_values) = names(p_values_shuffled)
    df_neg = df_neg[sample(1:nrow(df_neg)), ]
    for (ns in 1:nsplit){
      #            cat("-----------------", "\n", "split = ", ns , "\n")
      if (nsplit > 1){ 
        sub_neg = df_neg[(((ns-1)*dim_split)+1):min(dim_split*ns, nrow(df_neg)), ]
        df_split = rbind(cbind(df_pos, "outcome_num" = 1), cbind(sub_neg, "outcome_num" = 0))
      }else{
        df_split = rbind(cbind(df_pos, "outcome_num" = 1), cbind(df_neg, "outcome_num" = 0))
      }
      if(!(is.null(boolean_cols))){
        for (cc in boolean_cols){
          p_values[[cc]][ns] = cor(as.integer(df_split[[cc]]), df_split$outcome_num,  method = "pearson", use = "complete.obs")
          #  cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no boolean cols \n')}
      
      if(!(is.null(factor_cols))){
        outcome = factor(df_split$outcome_num)
        for (cc in factor_cols){
          res_chi = chisq.test(df_split[[cc]], outcome, correct=TRUE)
          p_values[[cc]][ns] = res_chi$p.value
          #cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no factor cols \n')}
      
      if (!(is.null(c(integer_cols, double_cols)))){
        for (cc in c(integer_cols, double_cols)){
          
          pk = kruskal.test(df_split[[cc]], df_split$outcome_num, na.omit = TRUE)
          p_values[[cc]][ns] = pk$p.value
          #  cat("p_value of ", cc, " = ", p_values[[cc]][ns],"\n")
        } 
      }else{ cat('no integer or double cols \n')}
      
    }
    #        cat("computing means for shuffle = ", nits, "\n")
    xx = colMeans(p_values, na.rm = TRUE)
    if (verbose ==1){
      print(xx)
    }
    cat("nits = ", nits, 'done! \n')
    
    p_values_shuffled[nits, names(xx)[which(!is.na(xx))]] = xx[which(!is.na(xx))]
  }
  
  return(colMeans(p_values_shuffled))
}
