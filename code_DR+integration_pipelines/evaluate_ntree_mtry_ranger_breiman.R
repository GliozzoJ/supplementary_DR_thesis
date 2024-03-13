evaluate_ntree_mtry_ranger_breiman <- function(df_train = NULL, 
                                               RF_method = 'RF',
                                               outcome = 'Activity', 
                                               idvar = NULL, 
                                               use_vars = names(df_train)[!(names(df_train) %in% c(outcome, idvar))],
                                               pos_label = levels(df_train[[outcome]])[1], 
                                               neg_label = levels(df_train[[outcome]])[2],
                                               n_internal_iters = 101, 
                                               ntree_mtry_candidate = NULL, 
                                               opti_meas = 'aucpr', 
                                               class_weights = 'none',
                                               balanced_bootstraps = TRUE,
                                               perc_samples = 0.9, 
                                               ratio_train = 0.9, 
                                               verbose = 0, 
                                               timeout_per_int_iter = 5,
                                               default_ntree = 501, ncores = my_detectCores()){
    
  # uses n_internal_iters holdouts to compute, for each parameter combination, the performance measure opti_meas
  # as the average across the optimeas computed over all the holdouts 
  # returns the first combination maximizing the opti_meas value
    
    # do not care about the ID of cases
    df_train = df_train[, (names(df_train) %in% c(use_vars, outcome))]
    for (opti_m in opti_meas) ntree_mtry_candidate[[opti_m]] = NA
    folds <- caret::createDataPartition(df_train[[outcome]], times = n_internal_iters, 
                                      p = ratio_train, list = TRUE);
    cat('parameter tuning with ', RF_method, '\n')
    if (verbose == 1){ cat('pos_label = ', pos_label, ' - neg_label = ', neg_label,'\n using variables: ', use_vars, '\n')}
    
    for (nr in 1:nrow(ntree_mtry_candidate)){
        if ('ntree' %in% names(ntree_mtry_candidate)) 
          nt = ntree_mtry_candidate[['ntree']][nr]
        else nt = default_ntree
        
        if ('mtry' %in% names(ntree_mtry_candidate)){ 
          mt = ntree_mtry_candidate[['mtry']][nr]
        }else{ 
          if (!is.factor(df_train[[outcome]])) mt = max(floor((ncol(df_train)-1)/3), 1) 
          else mt = floor(sqrt(ncol(df_train)-1))
        }
        if (verbose ==1){ cat('comb', nr,' =  ', unlist(ntree_mtry_candidate[nr, ]), ':') }
        cl <- makePSOCKcluster(min(ncores, detectCores()-1))
        registerDoParallel(cl)
        
        res_df = foreach(n_internal_it = 1:n_internal_iters, .combine = 'rbind', 
                         .packages = c('ranger','randomForest', 'ROCR', 'PerfMeas', 'R.utils', 'cutpointr'), 
                          .export = 'compute_perf_eval') %dopar%{
                           
        #for (n_internal_it in 1:n_internal_iters){
            
              idx = folds[[n_internal_it]]
              
              df_t = df_train[idx, ]
              df_test = df_train[-idx, ]
              # there is no more idvar here!
              x_train = df_t[, !(names(df_t) %in% c(outcome))] 
              y_train = df_t[[outcome]]
              tab_class = table(y_train)
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
                  if (balanced_bootstraps) sample.fraction = (min(tab_class) * perc_samples)* 1./ tab_class
                  else sample.fraction =  0.632
                  
                  rf = withTimeout({ranger(x = x_train, y = y_train, 
                         sample.fraction = sample.fraction, importance = 'permutation',
                         splitrule = ifelse(RF_method == 'extratrees', 'extratrees', 'gini'),
                         num.random.splits = ifelse(splitrule == 'gini',1,5),
                         class.weights = classwt,
                         num.trees =  nt, mtry = mt, 
                         probability = TRUE)}, timeout = timeout_per_int_iter, onTimeout = "warning")
                  if (is.null(rf)) return(NULL)
                  
                  test_probs = stats::predict(rf, data = df_test[,  !(names(df_test) %in% outcome)], type = 'response')
                  test_probs = test_probs$predictions
                  
              }else if (RF_method == 'RF'){
                  if (balanced_bootstraps){ 
                      nsamp = floor(min(tab_class)*perc_samples)
                      sampsize = c(nsamp,nsamp)
                      names(sampsize) = names(tab_class)
                  }else{ sampsize =  ceiling(0.632*nrow(x_train)) }
                  
                  rf_best =  withTimeout({randomForest(x = x_train, y = y_train, 
                                         xtest = df_test[, names(df_test) %in% names(x_train)], 
                                         sampsize = sampsize, norm.votes=TRUE,
                                         strata = y_train, classwt = classwt,
                                         importance=TRUE, 
                                         ntree = nt, mtry = mt)}, timeout = timeout_per_int_iter, onTimeout = "warning")
                  if (is.null(rf_best)) return(NULL)
                  test_probs = rf_best$test$votes
              }
              
              
              df_perf = compute_perf_eval(prediction_probs = test_probs, labels = df_test[[outcome]], 
                                          neg_label = neg_label, pos_label = pos_label, measures = opti_meas)
              
              return(df_perf)
          }
    
        stopCluster(cl)
        unregister_dopar()
        if (is.null(res_df)) return(NULL)
        
        if (length(opti_meas)==1) ntree_mtry_candidate[nr, opti_meas ] = mean(res_df[[opti_meas]])
        else{
            vecMeans = colMeans(res_df[, match(names(res_df), opti_meas)])
            ntree_mtry_candidate[nr, opti_meas ] = vecMeans
        }
        
        
        if (verbose ==1){ cat(' (', unlist(opti_meas), ' =) ', unlist(ntree_mtry_candidate[nr, opti_meas ]), '\n')}
    }
    
    #idx_best = which.max(ntree_mtry_candidate[[opti_meas]])
    
    return(ntree_mtry_candidate)
}