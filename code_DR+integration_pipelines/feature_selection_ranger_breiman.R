my_cor_fun = function(v, outcome, method = "pearson", thrconf = '0.05'){ 
  res = cor.test(v, outcome, method = method)
  cor_v = res$estimate
  cor_p = res$p.value
  cor_v = ifelse(cor_p < thrconf, cor_v, 0)
  return(cor_v)
}

#' Feature selection on training set data
#' 
#' This function performs feature selection on the training using a multiple
#' holdout procedure. During each holdout, the training of the feature selection
#' model is done using a balanced number of positives and negatives (by
#' downsampling of the majority class). Random Forest used the feature importance
#' as mean descrease in accuracy and the mean values across holdouts is used 
#' as indication of features importance. 
#'
#' @param df_train dataframe. Dataframe (samples X features) with the training
#' set data. idvar and outcome are columns of this dataset.
#' @param feature_selection_method string. Method chosen to perform feature
#' selection. Options are "RF_importance" (def.), "boruta", "glmnet".
#' @param minNumVar integer. Minimum number of features to select.
#' @param maxNumVar integer. Maximum number of features to select.
#' @param perc_samples double. Percentage of samples in the training set to 
#' create balanced "n_internal_iters" bootstraps.
#' @param outcome string. Name of the column in df_train used as target.
#' @param idvar string. Name of column in df_train containing samples ID. 
#' NULL if there is no variable with IDs
#' @param n_internal_iters integer. Number of stratified internal multiple 
#' holdouts used for feature selection.
#' @param pos_label integer. It represents the positive label (e.g. 1)
#' @param thrconf double. Threshold to select most informative features 
#' (feat >= thrconf).
#'
#' @return Vector containing the names of the selected features.
#' @export
feature_selection_ranger_breiman <- function(df_train, RF_method = 'RF', 
                                             feature_selection_method = 'RF_importance', 
                                             balanced_bootstraps = TRUE, 
                                             class_weights = 'class_priors',
                              minNumVar = 1, maxNumVar = Inf, perc_samples = 0.8, 
                              outcome = 'label', idvar = NULL, n_internal_iters = 51, 
                              pos_label = levels(df_train[[outcome]])[2], 
                              thrconf = 0.05, default_feature_selection = 'RF_importance', 
                              ncores = my_detectCores(), 
                              timeout_per_int_iter = 5,
                              verbose = 0){
  
    # do not care about the ID of cases
    df_train = df_train[, !(names(df_train) %in% idvar)]
    names_vars = names(df_train[, !(names(df_train) %in% c(idvar, outcome))])
    
    if (verbose == 1){ cat('starting feature selection with pos_label = ', pos_label, '\n')}
    
    if (feature_selection_method=='NULL') return(names_vars)
    
    #fn_log = file.path(getwd(), 'log_fs.log')
    #cat('feature selection saving logs in file = ', fn_log, '\n')
    
    cl <- makeCluster(min(ncores, detectCores()-1))#, outfile = fn_log)
    registerDoParallel(cl)
    
    internal_folds <- caret::createDataPartition(as.factor(df_train[[outcome]]), 
                                               times = n_internal_iters, 
                                               p = perc_samples, list = TRUE)
    # if (any(unlist(lapply(df_train,class))=="factor")){ 
    #   cat('Data contains some factors: using ', default_feature_selection, '\n')
    #   feature_selection_method = default_feature_selection
    # } 
    
    all_non_important = 0
   
    if (feature_selection_method == 'RF_importance'){
    #cat('for feature selection using ', RF_method, '\n')
        df_imp_feat_select = foreach(n_internal_it = 1:n_internal_iters, 
                                  .combine = 'cbind', .packages = c('ranger', 'randomForest', 'R.utils')) %dopar%{
      
          cat('n_internal_it = ', n_internal_it, '\n')                                      
          sub_train_idx = internal_folds[[n_internal_it]]
          df_t = df_train[sub_train_idx, ]
          x = df_t[, !(names(df_t) %in% c(outcome, idvar))] 
          y = df_t[[outcome]]      
          tab_class = table(y)      
          
          cat('tab_class = ', tab_class, '\n')                                      
          
          if (class_weights == 'class_priors'){ 
              # in this way, the less represented class will lower priors
              classwt = tab_class/sum(tab_class)
          }else if (class_weights == 'inv_class_priors'){ 
              # in this way, the less represented class will have higher priors
              classwt = sum(tab_class)/tab_class
          } else if (class_weights == 'none'){
              classwt = c(1,1)
          }
         
          if ((RF_method == 'extratrees') | (RF_method == 'ranger')){
              cat('train of ', RF_method, 'models\n')
              if (balanced_bootstraps) sample.fraction = (min(tab_class))* 1./ tab_class
              else sample.fraction = 0.632
              cat('sample fraction = ', sample.fraction, '\n')
              rf = withTimeout({ranger(x = x, y = y, 
                               sample.fraction = sample.fraction, importance = 'permutation',
                               splitrule = ifelse(RF_method == 'extratrees', 'extratrees', 'gini'),
                               num.random.splits = ifelse(splitrule == 'gini',1,5),
                               class.weights = classwt,
                               probability = TRUE)}, timeout = timeout_per_int_iter, onTimeout = "warning")
              if (is.null(rf)) return(NULL)
              else return(rf$variable.importance)
              
          }else if (RF_method == 'RF'){
              #cat('Breiman RF chosen \n')    
              cat('train of RF models\n')  
              if (balanced_bootstraps){
                  nsamp = floor(min(tab_class)*perc_samples)
                  sampsize = c(nsamp,nsamp)
                  names(sampsize) = names(tab_class)
              }else{
                  sampsize = 0.632*nrow(x)
              }
              cat('sampsize = ', sampsize, '\n')
              rf = withTimeout({randomForest(x = x, y = y, 
                                     sampsize = sampsize, 
                                     norm.votes=TRUE,
                                     strata = y, 
                                     classwt = classwt, 
                                     importance=TRUE, 
                                     probability = TRUE)}, timeout = timeout_per_int_iter, onTimeout = "warning")
              if (is.null(rf)) return(NULL)
              else return(rf$importance[, "MeanDecreaseAccuracy"])
              
          }
          
                                  }
        
        if (is.null(df_imp_feat_select)){ 
          cat("foreach runs always took too long to run selection - existing with NULL\n")
          return(NULL)
        }else if (!is.vector(df_imp_feat_select))  df_imp_feat_select = apply(df_imp_feat_select,1, mean)
        
        if (!any(df_imp_feat_select>0)){
          cat('variables were all non-important. using all of them\n')
          all_non_important = 1 
          selected_vars = names(df_imp_feat_select)
        }else{ 
          df_imp_feat_select = sort(df_imp_feat_select,decreasing = TRUE)
          min_num_vars = names(df_imp_feat_select)[1:minNumVar]
          df_imp_feat_select = df_imp_feat_select[df_imp_feat_select>0]
          if (length(df_imp_feat_select)< minNumVar){
              selected_vars = min_num_vars
          }else{
              cat('selecting features based on importance\n')
              
              df_imp_feat_select = df_imp_feat_select/sum(abs(df_imp_feat_select))
              df_imp_feat_select = cumsum(df_imp_feat_select)
              idx = which(df_imp_feat_select < 0.95)
              
              df_imp_feat_select = df_imp_feat_select[1: max(idx, minNumVar)]
              selected_vars = names(df_imp_feat_select)
          }
        }
    }else if (feature_selection_method=='boruta'){    
      
        df_imp_feat_select = foreach(n_internal_it = 1:n_internal_iters, 
                                     .combine = 'cbind', .packages = c('Boruta', 'caret')) %dopar%{
          
          sub_train_idx = internal_folds[[n_internal_it]]
          df_t = df_train[sub_train_idx, ]
        
          # balancing the training set via downsampling
            balanced_df = caret::downSample(df_t, df_t[[outcome]], list = FALSE)
            x = balanced_df[, !(names(balanced_df) %in% c(outcome, 'Class'))]
            y = balanced_df[['Class']]
          
            boruta_decision = withTimeout({Boruta(x = x, y = y, pValue = (thrconf*2), maxRuns = n_internal_iters)}, 
                                            timeout = timeout_per_int_iter, onTimeout = "warning")
            if (is.null(boruta_decision)) return(NULL)
            selected_vars = getSelectedAttributes(boruta_decision, withTentative = TRUE)
          return(selected_vars)
        }
        if (is.null(df_imp_feat_select)){ 
          cat("foreach runs always took too long to run selection - existing with NULL\n")
          return(NULL)
        }else if (!is.vector(df_imp_feat_select)) df_imp_feat_select = apply(df_imp_feat_select,1, mean)
        df_imp_feat_select = sort(df_imp_feat_select,decreasing = TRUE)
        df_imp_feat_select = df_imp_feat_select[1:min(maxNumVar, length(df_imp_feat_select))]
        
        if (!any(df_imp_feat_select>=(thrconf^2))){
          cat('variables were all non-important. using all of them')
          all_non_important = 1 
          selected_vars = names(df_imp_feat_select)
        }else{
          df_imp_feat_select = df_imp_feat_select[df_imp_feat_select>=(thrconf^2)]
          if(length(df_imp_feat_select)<minNumVar){
            selected_vars = names_vars
          }else {
            selected_vars = names(df_imp_feat_select)
          }
        }
    
    }else if(feature_selection_method=='glmnet'){
      
        df_imp_feat_select = foreach (n_internal_it = 1:n_internal_iters, .combine = 'cbind', 
                                      .packages = c('glmnet', 'caret')) %dopar%{
          
          sub_train_idx = internal_folds[[n_internal_it]]
          df_t = df_train[sub_train_idx, ]
          
          # balancing the training set via downsampling
          balanced_df = caret::downSample(df_t, df_t[[outcome]], list = FALSE)
          
          x = balanced_df[, !(names(balanced_df) %in% c(outcome, 'Class'))]
          if (any(unlist(lapply(x,class))=="factor")){ x = as.matrix(data.frame(lapply(x, as.numeric)))
          }else{ x = as.matrix(x)}      
          
          y = as.numeric(balanced_df[['Class']] == pos_label) 
          x_norm = apply(x, 2, function(xx) (xx - mean(xx, na.rm = TRUE))/(sd(xx, na.rm = TRUE)))
          cvfit = cv.glmnet(x_norm,y, alpha = 0.5, family = 'binomial')
          vals = coef(cvfit, s = "lambda.min")
          
          return(vals[,1])
        }
        
        if (is.null(df_imp_feat_select)){ 
          cat("foreach runs always took too long to run selection - existing with NULL\n")
          return(NULL)
        }else if (!is.vector(df_imp_feat_select)) {
          df_imp_feat_select = df_imp_feat_select[!grepl('Intercept', row.names(df_imp_feat_select)), ]
          df_imp_feat_select = apply(abs(df_imp_feat_select),1, mean)
        }else{
          df_imp_feat_select = df_imp_feat_select[!grepl('Intercept', names(df_imp_feat_select)) ]
        }
        df_imp_feat_select = sort(df_imp_feat_select,decreasing = TRUE)
        df_imp_feat_select = df_imp_feat_select[1:min(maxNumVar, length(df_imp_feat_select))]
        
        if (!any(df_imp_feat_select>=(thrconf^2))){
          cat('variables were all non-important. using all of them')
          all_non_important = 1 
          selected_vars = names(df_imp_feat_select)
        }else{
          df_imp_feat_select = df_imp_feat_select[df_imp_feat_select>=(thrconf^2)]
          if(length(df_imp_feat_select)<minNumVar){
            selected_vars = names_vars
          }else {
            selected_vars = names(df_imp_feat_select)
          }
        }
    
    }else  if(grepl('corr-', feature_selection_method)){
        method_corr = strsplit(feature_selection_method,'corr-')[[1]][2]
        
        
        df_imp_feat_select = foreach(n_internal_it = 1:n_internal_iters, 
                                     .combine = 'cbind', .export = 'my_cor_fun', .packages = 'caret') %dopar%{
          sub_train_idx = internal_folds[[n_internal_it]]
          df_t = df_train[sub_train_idx, ]
          
          # balancing the training set via downsampling
          balanced_df = caret::downSample(df_t, df_t[[outcome]], list = FALSE)
          
          x = balanced_df[, !(names(balanced_df) %in% c(outcome, 'Class'))]
          if (any(unlist(lapply(x,class))=="factor")) x = as.matrix(data.frame(lapply(x, as.numeric)))
          else x = as.matrix(x)
          
          y = as.numeric(balanced_df[['Class']] == pos_label) 
          
          cor_vals = apply(x,2, my_cor_fun, outcome = y, method = method_corr )
          
          return(cor_vals)
        }
        
        if (is.null(df_imp_feat_select)){ 
          cat("foreach runs always took too long to run selection - existing with NULL\n")
          return(NULL)
        }else{ 
          if (!is.vector(df_imp_feat_select)) {  df_imp_feat_select = apply(abs(df_imp_feat_select),1, mean)}
        }
        df_imp_feat_select = sort(df_imp_feat_select,decreasing = TRUE)
        df_imp_feat_select = df_imp_feat_select[1:min(maxNumVar, length(df_imp_feat_select))]
        
        if (!any(df_imp_feat_select>=0.1)){
          cat('variables were all non-important. using all of them\n')
          all_non_important = 1 
          selected_vars = names(df_imp_feat_select)
        }else{
          df_imp_feat_select = df_imp_feat_select[df_imp_feat_select>=0.1]
          if(length(df_imp_feat_select)<minNumVar){
            selected_vars = names_vars
          }else {
            selected_vars = names(df_imp_feat_select)
          }
        }
    
    }
    
    
    
    stopCluster(cl)
    unregister_dopar()
    
    return(list('selected_vars' = selected_vars, 'all_non_important' = all_non_important)) 
}
