analisi_RF_myOpti_ranger_breiman <- function(df, exp_desc_str = '',
                                             recompute_RF = FALSE,
                                        RF_method = 'RF', # if RF_method == 'extratrees' use ranger and use as split.criterion the extratrees
                                        maxNumVar = Inf, minNumVar = 3,
                                        pos_label = levels(df[[outcome]])[2], 
                                        neg_label = levels(df[[outcome]])[1], 
                                        outcome = 'Activity', idvar = 'Patient_ID',
                                        str_desc ='data', 
                                        save_path = './', 
                                        perc_samples = 0.8, 
                                        ratio_train = 0.9,
                                        # almeno posso compararmi a quando ne avevo 101!
                                        max_n_external_iters = 101, 
                                        n_internal_iters = 10, 
                                        n_external_iters = 10, 
                                        feature_selection_method = 'RF_importance', 
                                        # if feature_selection_method == 'none' no feature selection is applied
                                        seed = Inf, thr_prob = 0.5, 
                                        tuneRF = TRUE, 
                                        class_weights = 'none',
                                        balanced_bootstraps = TRUE,
                                        ntree = ifelse(ncol(df)>100, 1001,ifelse(ncol(df)<20, 101,501)), 
                                        mtry = sqrt(ncol(df)-1),
                                        measures = c("auc", "aucpr"), 
                                        opti_meas = c('aucpr', 'auc'), 
                                        weights_opti_meas = c(0.5,0.5), 
                                        feat_selection_patience = round(n_internal_iters/2), 
                                        just_collect_results = FALSE,
                                        nDigits = 5, verbose = 0,
                                        timeout_per_int_iter = 5,
                                        timeout_per_ext_iter = 5,
                                        ncores = my_detectCores()){

    opti_meas_str = ifelse(length(opti_meas)==1, opti_meas, paste(opti_meas, collapse ='_'))
    base_save_str = paste0(exp_desc_str, paste(str_desc, class_weights, balanced_bootstraps, sep = '_'))
    if (!(RF_method == 'RF'))  base_save_str = paste0(base_save_str, '_', RF_method)    
    if (!(feature_selection_method == 'RF_importance')) base_save_str =  paste0(base_save_str, '_', feature_selection_method)
    if (!(opti_meas_str == 'aucpr')) base_save_str =  paste0(base_save_str, '_', opti_meas_str)
    fn_res_list = file.path(save_path, paste0('listRF_', base_save_str, '.Rda', sep=''))
    cat("result file will be saved/load from ", fn_res_list, '\n')
    

    if ((!recompute_RF) & (file.exists(fn_res_list) & (file.info(fn_res_list)$size>0))){
        cat('loading classification results\n')
        load(fn_res_list)
        if (!exists("res_list")){ 
            res_list = NULL
            save(res_list, file = fn_res_list)
            return(res_list)
        }
        if (is.null(res_list)){ 
          return(res_list)
        }
        
        cat(names(res_list), '\n-----\n')
        saved_measures = measures[measures %in% names(res_list$all_holdouts)]
        if ((sum(measures %in% names(res_list$all_holdouts)) < length(measures)) | (n_external_iters != res_list$mean_over_holdouts$n_external_iters)) {
            cat('some performance measures are missing; I\'ll try to reestimate them!\n')
            res_mega_perf = res_list$all_holdouts
            res_perf = res_list$mean_over_holdouts
            all_predictions = res_list$all_predictions_for_explanation
            var_importances = res_list$var_importances
            
            if (ratio_train != res_perf$ratio_train) {
                warning('ratio_train is different than saved results!', 
                            'ratio_train = ', ratio_train, 
                            'saved ratio_train = ', res_perf$ratio_train, '\n')
            }
            if (n_external_iters != res_perf$n_external_iters) {
                warning('no_holdouts is different than saved results!', 
                        'no_holdouts = ', n_external_iters, 
                        'saved external_iters= ', res_perf$n_external_iters, '\n')
            }
            
            expected_no_test_per_holdout = nrow(all_predictions)/n_external_iters
            all_perf = data.frame()
            all_predictions_new = all_predictions
            if (!any(names(all_predictions) %in% 'holdout_no')){
                all_predictions_new[['holdout_no']] = NA
            }
            for (niter in 1:n_external_iters){
                if (any(names(all_predictions) %in% 'holdout_no')){
                    test_probs = all_predictions[all_predictions$holdout_no == niter,]
                }else{
                    start = ((niter-1)*expected_no_test_per_holdout)+1
                    stop = start+expected_no_test_per_holdout-1
                    test_probs = all_predictions[start:stop,]
                    all_predictions_new[['holdout_no']][start:stop] = niter 
                }
                
                test_probs[[pos_label]] = test_probs$pos_pred
                test_probs[[neg_label]] = 1-test_probs$pos_pred
                
                df_perf = compute_perf_eval(prediction_probs = test_probs, 
                                            labels = factor(test_probs[['labels']], levels = c(neg_label, pos_label)), 
                                            neg_label = neg_label, 
                                            pos_label = pos_label,
                                            measures = measures)
                all_perf = rbind(all_perf, cbind(holdout_no = niter, df_perf))
                if (verbose ==1){ 
                    print(df_perf)
                    already_saved_measures = names(res_list$all_holdouts)[names(res_list$all_holdouts) %in% measures]
                    cat(niter, ']', sum(df_perf[already_saved_measures] == 
                                res_list$all_holdouts[niter, already_saved_measures]) == 
                            length(already_saved_measures), ' matching measures: \n')
                    print(df_perf[already_saved_measures])  
                    print(res_list$all_holdouts[niter, already_saved_measures])
                        
                }
                
            }
            mean_meas = apply(all_perf, 2, mean, na.rm = TRUE)
            mean_meas_df = data.frame(t(mean_meas))
            names(mean_meas_df) = paste(names(mean_meas_df), '_avg', sep ='')
            
            sd_meas = apply(all_perf, 2, sd, na.rm = TRUE)
            sd_meas_df = data.frame(t(sd_meas))
            names(sd_meas_df) = paste(names(sd_meas_df), '_std', sep ='')
            
            #### SISTEMA IL SALVATAGGIO!!!
            
            #SISTEMA IL SALVATAGGIO!!!
            names(all_perf)[names(all_perf) %in% saved_measures] = 
                paste0(names(all_perf)[names(all_perf) %in% saved_measures], '_copy') 
            res_mega_perf_new = cbind(res_mega_perf, all_perf)
            
            names(mean_meas_df)[names(mean_meas_df) %in% paste0(saved_measures, '_avg')] = 
                paste0(names(mean_meas_df)[names(mean_meas_df) %in% paste0(saved_measures, '_avg')], '_copy')             
            names(sd_meas_df)[names(sd_meas_df) %in% paste0(saved_measures, '_std')] = 
                paste0(names(sd_meas_df)[names(sd_meas_df) %in% paste0(saved_measures, '_std')], '_copy')             
            
            res_perf_new = cbind(res_perf, cbind(mean_meas_df, sd_meas_df))
            
            if ((sum(abs(res_mega_perf_new$aucpr_copy - res_mega_perf_new$aucpr)) > 0) | 
                (abs(res_perf_new$aucpr_avg_copy - res_perf_new$aucpr_avg) > 0)){
                cat('merda qualcosa non va!\n')
                stop()
            }
            
            res_list = list('mean_over_holdouts' = res_perf_new, 
                            'all_holdouts' = res_mega_perf_new, 
                            'all_predictions_for_explanation' =  all_predictions_new, 
                            'var_importances' = var_importances) 
            
            
            
            cat('Resaving RF results in file = ', fn_res_list, '\n')
            save(res_list, file = fn_res_list)
            file.remove(fn_res_list)
        }
        
    }else{
        if (just_collect_results) return(NULL)
        # se seed ? Inf non viene settato e gli holdout cambiano sempre!
        #feature_selection may be = 'RF_importance', 'boruta', 'glmnet'
        
        if (class_weights == 'class_priors'){ 
            # in this way, the less represented class will lower priors
            cat('using class priors for weighing\n')
        }else if (class_weights == 'inv_class_priors'){ 
            # in this way, the less represented class will have higher priors
            cat('using inverse class priors for weighing\n')
        } else if (class_weights == 'none'){
            cat('unweighted classification\n') 
        }
        
        if (balanced_bootstraps){
            cat('balanced bootstraps\n')
        }else{ 
            cat('unbalanced bootstraps\n')
        }
        
        all_combs = NULL
        nneg = sum(df[[outcome]]!= pos_label)
        npos = sum(df[[outcome]]== pos_label)
        
        cat('\n\nstarting (not optimized) analysis with matrix ', str_desc, ' - dim = ', dim(df),  '\n')
        cat('pos label =', pos_label, ' (npos =', npos,' ) -  neg label = ', neg_label, ' (nneg = ', nneg, ')\n')
        
        
        names_vars = names(df)[!(names(df) %in% c(idvar, outcome))]
        df_imp <- data.frame(matrix(0, nrow = length(names_vars), ncol = n_external_iters))
        row.names(df_imp) <- names_vars
        
        all_perf = data.frame()  
        if (is.finite(seed)) set.seed(seed) # to reproduce the same holdouts
        folds <- caret::createDataPartition(df[[outcome]], times = max_n_external_iters, 
                                          p = ratio_train, list = TRUE);
        
        mega_df_for_explanations = data.frame()
        prop_vars_selected = NULL
        no_all_non_important = 0
        for (niter in 1:min(n_external_iters, max_n_external_iters)) {
            if (verbose ==1) cat('external holdout no = ', niter, '\n')
            
            train_idx = folds[[niter]]
            df_test = df[-train_idx, ]
            df_train = df[train_idx, ]
            
            nneg_train = sum(df_train[[outcome]]!= pos_label)
            npos_train = sum(df_train[[outcome]]== pos_label)
            
            if (verbose ==1) cat(niter, '] nneg_train =', nneg_train,'- npos_train = ', npos_train, '\n') 
            
            if (!(feature_selection_method == 'none')){
              
              if (ncol(df_train)>minNumVar){ 
                        res_feat_sel = feature_selection_ranger_breiman(df_train, 
                                                          RF_method = RF_method, 
                                                          feature_selection_method = feature_selection_method,
                                                          minNumVar = 2, 
                                                          maxNumVar = nrow(df)-1, 
                                                          perc_samples = perc_samples, 
                                                          outcome = outcome, 
                                                          idvar = idvar, 
                                                          n_internal_iters = n_internal_iters,
                                                          pos_label = pos_label, 
                                                          balanced_bootstraps = balanced_bootstraps,
                                                          class_weights = class_weights,
                                                          ncores = min(ncores, detectCores()-1), 
                                                          timeout_per_int_iter = timeout_per_int_iter,
                                                          verbose = verbose)
                        if (is.null(res_feat_sel)) return(NULL)
              }else{ 
                res_feat_sel= list('selected_vars' = names_vars, 'all_non_important' = 0 )
              }
              selected_vars = res_feat_sel[['selected_vars']]
              
            }else{
              res_feat_sel= list('selected_vars' = names_vars, 'all_non_important' = 0 )
              selected_vars = names_vars
              
            }
            no_all_non_important = no_all_non_important + res_feat_sel[['all_non_important']]
            if (niter >= min(round(n_external_iters/3), feat_selection_patience)){
                 if (no_all_non_important == niter){
                     return(NULL)
                 }
            } 
            nvar = length(selected_vars)
            prop_vars_selected = c(prop_vars_selected, round(nvar/(ncol(df_train)-1), 2))
            if (verbose==1){
              cat(nvar , ' features selected \n')
            }
            
            # use only selected vars
            df_train = df_train[, names(df_train) %in% c(selected_vars, outcome)]
            #df_test = df_test[, names(df_test) %in% c(selected_vars, outcome)]
            if (tuneRF){
                if (nvar > 2){ 
                    if (nvar > 10){
                        ntree_candidate = c(501, 751, 1001)
                        mtry_candidate = unique(c(floor(nvar/10), floor(nvar/5), floor(sqrt(nvar))))
                        mtry_candidate = mtry_candidate[mtry_candidate>0 & mtry_candidate<nvar]                        
                    }else{
                        ntree_candidate = c(101,501)
                        mtry_candidate = floor(sqrt(nvar))
                    }
                    ntree_mtry_candidate = cbind(expand.grid(ntree_candidate,mtry_candidate))
                    names(ntree_mtry_candidate) = c("ntree", "mtry")
                    
                    #set.seed(seed)
                    
                    evals_ntree_mtry = evaluate_ntree_mtry_ranger_breiman(
                        df_train = df_train, 
                        RF_method = RF_method,
                        outcome = outcome, 
                        neg_label = neg_label, 
                        pos_label = pos_label,
                        use_vars = selected_vars,
                        n_internal_iters = n_internal_iters, 
                        idvar = idvar,
                        ntree_mtry_candidate = ntree_mtry_candidate, 
                        opti_meas = opti_meas, 
                        perc_samples = perc_samples,
                        ratio_train = ratio_train, 
                        class_weights = class_weights,
                        balanced_bootstraps = balanced_bootstraps,
                        verbose = verbose, 
                        timeout_per_int_iter = timeout_per_int_iter,
                        ncores = min(ncores, detectCores()-1)) 
                    if (is.null(evals_ntree_mtry)) return(NULL)
                    
                    if (length(opti_meas)<=1){
                        fscore = evals_ntree_mtry[[opti_meas]]
                    }else{
                        # columnwise multiplication of columns in opti_meas and weights_opti_meas
                        int_score = t(t(evals_ntree_mtry[, names(evals_ntree_mtry) %in% opti_meas]) * 
                                          weights_opti_meas[match(opti_meas, names(evals_ntree_mtry[, names(evals_ntree_mtry) %in% opti_meas]))])
                        fscore = rowSums(int_score)
                        
                    }
                    evals_ntree_mtry[['final_score']] = fscore
                    idx_best = which.max(evals_ntree_mtry[['final_score']])
                    best_comb = evals_ntree_mtry[idx_best, ]                              
                    
                }else{
                    best_comb = c(101,1,rep(NA, length(opti_meas)), NA)
                    names(best_comb) = c('ntree', 'mtry',opti_meas, 'final_score')
                    best_comb = data.frame(t(best_comb))
                }
            }else{
                if (is.na(mtry)){
                    cat('setting default parameter for mtry\n')
                    mtry = ifelse((ncol(df_train)-1)<3,1,sqrt(ncol(df_train)-1))
                }
                if (is.na(ntree)){
                    cat('setting default parameter for ntree\n')
                    ntree = 500
                }
                cat('using user defined values for RF hyper-params: ntree = ', ntree, ' - mtry = ', mtry, '\n')
                best_comb = c(ntree,mtry,rep(NA, length(opti_meas)), NA)
                names(best_comb) = c('ntree', 'mtry',opti_meas, 'final_score')
                best_comb = data.frame(t(best_comb))
            }
            
            all_combs = rbind(all_combs, best_comb)
            if (verbose ==1) {
                print(best_comb) 
                cat( ' selected as best parameters\n')
            }
            
            x_train = df_train[, !(names(df_train) %in% outcome)] 
            y_train = df_train[[outcome]]
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
            
            if (verbose == 1) cat('class weights = ', classwt, '\n')
            
            if ((RF_method == 'extratrees') | (RF_method == 'ranger')){
                cat('train extratrees models\n')
                
                if (balanced_bootstraps) sample.fraction = (min(tab_class))* 1./ tab_class
                else sample.fraction = 0.632
                cat('sample fraction = ', sample.fraction, '\n')
                rf_best = withTimeout({ranger(x = x_train, y = y_train, 
                                 sample.fraction = sample.fraction, importance = 'permutation',
                                 splitrule = ifelse(RF_method == 'extratrees', 'extratrees', 'gini'),
                                 num.random.splits = ifelse(splitrule == 'gini',1,5),
                                 class.weights = classwt,
                                 num.trees =  best_comb[,'ntree'], mtry = best_comb[,'mtry'], 
                                 probability = TRUE)}, timeout = timeout_per_ext_iter, onTimeout = "warning")
                if (is.null(rf_best)) return(NULL)
                imp = rf_best$variable.importance
                #setting the correct response type for prediction on test set
                test_probs = stats::predict(rf_best, data = df_test, type = 'response')
                test_probs = test_probs$predictions
                rm(rf_best)
            }else if (RF_method == 'RF'){
                #cat('Breiman RF chosen \n')    
                cat('train of RF models\n')  
                if (balanced_bootstraps){
                    nsamp = floor(min(tab_class)*perc_samples)
                    sampsize = c(nsamp,nsamp)
                    names(sampsize) = names(tab_class)
                }else{
                    sampsize = 0.632*nrow(x_train)
                }
                rf_best =  withTimeout({randomForest(x = x_train, y = y_train, xtest = df_test[, (names(df_test) %in% colnames(x_train))] , 
                                       sampsize = sampsize, norm.votes=TRUE,
                                       strata = y_train, 
                                       classwt = classwt, 
                                       importance=TRUE, ntree = best_comb[,'ntree'], mtry = best_comb[,'mtry'],
                                       probability = TRUE)}, timeout = timeout_per_ext_iter, onTimeout = "warning")
                if (is.null(rf_best)) return(NULL)
                imp = rf_best$importance[, "MeanDecreaseAccuracy"]
                #setting the correct prediction on test set
                test_probs = rf_best$test$votes
            }
            if (verbose==1) cat('training fnished\n')
            
            df_imp[names(imp), niter] = imp
            
            
            if (niter >= min(round(n_external_iters/3), feat_selection_patience)){
              mean_imp = apply(as.matrix(df_imp[,1:niter]),1,mean)
              if (!(any(mean_imp>0))) return(NULL)
            } 
            
            if (verbose==1) cat('evaluation of test data\n')
            
            
            mega_df_for_explanations = rbind(mega_df_for_explanations, data.frame(patient_id = row.names(test_probs), holdout_no = niter, 
                                                                                  pos_pred = test_probs[, colnames(test_probs) %in% pos_label],
                                                                                  labels = df_test[[outcome]]))
            df_perf = compute_perf_eval(prediction_probs = test_probs, labels = df_test[[outcome]], 
                                        neg_label = neg_label, pos_label = pos_label,
                                        measures = measures)
            all_perf = rbind(all_perf, cbind(holdout_no = niter, df_perf))
            if (verbose ==1) print(df_perf)
        }
        
        mean_meas = apply(all_perf, 2, mean, na.rm = TRUE)
        mean_meas_df = data.frame(t(mean_meas))
        names(mean_meas_df) = paste(names(mean_meas_df), '_avg', sep ='')
        
        sd_meas = apply(all_perf, 2, sd, na.rm = TRUE)
        sd_meas_df = data.frame(t(sd_meas))
        names(sd_meas_df) = paste(names(sd_meas_df), '_std', sep ='')
        
        # str_perfs_df = apply(rbind(mean_meas, sd_meas), 2, function(x){paste(round(x[1],3), '+/-', round(x[2],3))})
        # str_perfs_df = data.frame(t(str_perfs_df))
        # names(str_perfs_df) = paste(names(str_perfs_df), '_str', sep ='')
        df_comb_save = data.frame(all_combs)
        nn = names(df_comb_save)
        nn[match(opti_meas, nn)] = paste(opti_meas, 'train', sep ='_')
        names(df_comb_save) = nn
        opti_meas_str = paste(opti_meas, collapse='_')
        weight_meas_str = paste(weights_opti_meas, collapse='_')
        
        res_mega_perf = cbind(data.frame(opti_meas = opti_meas_str, weight_meas = weight_meas_str, 
                                         nrow_data = nrow(df), ncol_data = ncol(df), 
                                        prop_vars_selected = prop_vars_selected, 
                                        perc_samples = perc_samples, ratio_train = ratio_train, 
                                        feature_selection_method = feature_selection_method, 
                                        balanced_bootstraps = balanced_bootstraps,
                                        class_weights = class_weights,
                                        n_external_iters = niter, 
                                        n_internal_iters= n_internal_iters), 
                                        df_comb_save, all_perf)
        
        res_perf = cbind(data.frame(opti_meas = opti_meas_str, 
                                weight_meas = weight_meas_str, 
                                nrow_data = nrow(df), ncol_data = ncol(df),
                                prop_vars_selected = mean(prop_vars_selected), 
                                perc_samples = perc_samples, 
                                ratio_train = ratio_train, 
                                feature_selection_method = feature_selection_method, 
                                n_external_iters = niter, n_internal_iters= n_internal_iters, 
                                maxvoted_no_tree = names(which.max(table(all_combs[, 'ntree']))), 
                                maxvoted_mtry = names(which.max(table(all_combs[, 'mtry'])))),
                                mean_meas_df, 
                                sd_meas_df, 
                                path_save = fn_res_list)
        
        
        if (verbose ==1){
            cat('****** perf on ', str_desc,' *******\n')
            print(res_perf)
        }
        
        df_imp_save = data.frame('name' = row.names(df_imp), 'permutation_importance' = apply(df_imp,1,mean), 'permutation_importance_sd' = apply(df_imp,1,sd)) 
        # fn_importances = file.path(save_path, paste(base_save_str,'_importances', sep=''))
        # write_xlsx(x=df_imp_save, path = paste(fn_importances, '.xlsx', sep =''))
        
        
        #save(mega_df_for_explanations, file = file.path(save_path, paste(base_save_str, '_all_predictions.Rda', sep='')))
        #mega_df_for_explanations_csv = cbind('pts_ID'=row.names(mega_df_for_explanations), mega_df_for_explanations)
        #write_xlsx(x = mega_df_for_explanations, path = file.path(save_path, paste(base_save_str,'_all_predictions.xlsx', sep='')))
        
        
        res_list = list('mean_over_holdouts' = res_perf, 
                        'all_holdouts' = res_mega_perf, 
                        'all_predictions_for_explanation' = mega_df_for_explanations, 
                        'var_importances' = df_imp_save) 
        
        save(res_list, file = fn_res_list)
    }
    
    if (!any(names(res_list$mean_over_holdouts) %in% 'path_save')) {
        res_list$mean_over_holdouts = cbind(res_list$mean_over_holdouts, path_save = fn_res_list)
        save(res_list, file = fn_res_list)
    }
    if (verbose == 1) print(res_list[['mean_over_holdouts']])
  
  return(res_list)
}