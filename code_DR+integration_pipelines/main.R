gc()
rm(list = ls())

# #### settings for reticulate (interface to Python)
# python_path = '/home/jessica/miniconda3/bin' #insert python path
# setwd('.') #set the working directory
# Sys.setenv(RETICULATE_PYTHON = python_path)
# library(reticulate)
# 
# use_python(python = python_path)


# MKL : https://academic.oup.com/bioinformatics/article/34/6/1009/4565592 
#Unsupervised multiple kernel learning for heterogeneous data integration

# #### settings to set the conda environment
setwd(".") # set working directory

library(reticulate)
use_condaenv(condaenv="r_4.3") #set path to conda environment
# ####

# Load all required packages (reticulate was loaded before)
library(Boruta)
library(caret)
library(cutpointr)
library(dplyr)
library(doParallel)
library(entropy) # for unsup feat selection via entropy estimation
library(foreach)
library(genieclust)
library(ggplot2)
library(glmnet)
library(kernlab)
library(igraph)
library(impute)
library(intrinsicDimension)
library(intRinsic)
library(limma)
library(loe)
library(missRanger)
library(mixKernel) # for MKL di Mariette and Vollaneaux 
library(SNFtool)
library(MOFA2)
library(parallel)
library(PerfMeas)
library(plyr)
library(preprocessCore)
library(randomForest)
library(ranger)
library(Rdimtools)
library(readr)
library(readxl)
library(RMKL) # per supervised multiple kernel learning (install from source file, it is archived)
library(ROCR)
library(R.utils)
library(rsvd)
library(StatMatch)
library(stringr)
library(tsne)
library(umap)
library(tableone)
library(writexl)

source('./analisi_RF_myOpti_ranger_breiman.R')
source('./blocking_ID.R')
source('./cluster_features_par.R')
source('./compute_distance.R')
source('./compute_perf_eval.R')
source('./create_mofa_embedding.R')
source('./estimate_ID_Danco.R')
source('./estimate_ID_twoNN.R')
source('./evaluate_ntree_mtry_ranger_breiman.R')
source('./feature_selection_ranger_breiman.R')
source('./filter_cols_by_p.R')
source('./fix_dim_names.R')
source('./hierarchical_svd_entropy_par.R')
source('./RCUR_selection_par.R')
source('./reduce_data_dimensionality.R')
source('./remove_correlated_par.R')
source('./utilities.R')

#reticulate::source_python('embedding.py')

# defined path to data
dati_rinominati <- file.path('./data')


ext_no_neighs <<- 5
ext_seed <<- 3




dim_red <- list(
    "umap" = function(x=NULL, k=NULL, method = 'euclidean', n_neighbors = ext_no_neighs){
        set.seed(3)
        cat('umap with ', method, ' distance\n')
        if ((sum(is.na(x)))>0){ 
            x_dist    = compute_distance(x, method = method)
            proj_data = umap(x_dist, n_components = k,   
                             n_neighbors = n_neighbors, input="dist", method = 'umap-learn')
        }else{
            proj_data = umap(x, n_components = k,    
                             n_neighbors = n_neighbors, method = 'umap-learn')
     }
        return(as.matrix(proj_data$layout))},
    "tsne" = function(x=NULL, k=NULL){
        set.seed(3)
        cat('tsne\n')
        initial_dims = ifelse(ncol(x)<(nrow(x)-1), min(k*10, round(ncol(x)*2/3)), (nrow(x)-1))
        return(tsne(imputa(x), k = k, initial_dims = initial_dims))
     },
    "rpca" = function(x=NULL, k=NULL){
        cat('RPCA\n')
        s <- rpca(imputa(x),k = k,
                            center = TRUE,
                            scale = TRUE,
                            retx = TRUE,
                            rand = TRUE)
        return(s$x)
     },
    "rcur" = function(x=NULL, k=NULL){
        r <- rcur(imputa(x),k = k,
                            rand = TRUE)
        x = x[ ,  r$C.idx]
        cat(row.names(x))
        cat('returning mat with dimension:', dim(x), '\n')
        return(x)
    },
    "cur" = function(x=NULL, k=NULL){
        r <- rcur(imputa(x),k = k,
                            rand = FALSE)
        x = x[ ,  r$C.idx]
        cat(row.names(x))
        cat('returning mat with dimension:', dim(x), '\n')
        return(x)
    },
    "laplacianEigenmaps_dimRedtools" = function(x = NULL, k = NULL, type = c("knn", ext_no_neighs)) {
        # type tells how to connect the graph. type proportion means that the 10% of edges are kept
        # type = c("knn", n_neighbors) composes knn neighborhoods and guarantess connectivity
        # usa gli eigenvectors, quindi k deve per forza essere il minimo tra k, ncol(x)-1, nrow(x)-1
        cat('LaplacianEigenmaps\n')
        set.seed(3)
        xx = do.lapeig(imputa(x), ndim = min(k, nrow(x)-1), type = type)
        return(xx$Y)
    },
    "mofa" = function(x = NULL, k = 15, seed = ext_seed, no_mofa_factors = 15, 
                                        sparse_weights = TRUE, sparse_factors = FALSE){
        k = min(k, no_mofa_factors)
        
        mofaobj = create_mofa_from_matrix(list(t(x)))
        
        data_opts <- get_default_data_options( mofaobj )
        #data_opts$scale_views = TRUE
        
        model_opts <- get_default_model_options( mofaobj )
        model_opts$num_factors <- k
        # model_opts$spikeslab_factors = sparse_factors
        # model_opts$spikeslab_weights = sparse_weights
        # model_opts$ard_factors = sparse_factors
        # model_opts$ard_weights = sparse_weights
        
        train_opts <- get_default_training_options( mofaobj    )
        train_opts$seed <- seed
        train_opts$convergence_mode <- "slow"
        
        mofaobj <- prepare_mofa( mofaobj,
                                 training_options = train_opts,
                                 model_options = model_opts,
                                 data_options = data_opts )
        mofaobj <- run_mofa( mofaobj, outfile = file.path('./', 'MOFA.hdf5'), save_data = FALSE, use_basilisk = TRUE)
        
        
        Z = mofaobj@expectations$Z # latent factor matrix
        Z = Z$group1
        return(Z)
    }

)


codes = unique(unlist(lapply(lapply(list.files(file.path(dati_rinominati)), str_split_1, pattern ='_'), function(x) x[1])))
#codes = c('BLCA2', 'BRCA1', 'BRCA2', 'KIRC1', 'LUAD2', 'LUSC3', 'OV1', 'PRAD1', 'SKCM1')




# To compare ID and NoID
# main_best_models(str_add_save = 'compNew_speedy_', 
#      cohort_list = codes,
#      use_id_list = list(noID_1 = list(use_id = FALSE, frac_samples = 1, minNumViews = 4), 
#                         noID_2 = list(use_id = FALSE, frac_samples = 2, minNumViews = 4), 
#                         ID = list(use_id = TRUE, frac_samples = NA, minNumViews = 3)),
#      just_collect_results = FALSE,
#      n_internal_iters = 1, n_external_iters = 15, 
#      perc_samples = 1,
#      ntree = NA, mtry = NA, # use default parameters for RF
#      timeout_per_int_iter = 10, 
#      feature_selection_method_list = 'none',
#      unsup_feature_selection_list = c('feature_clustering', 'par_rcur', 'NULL'), 
#      best_red_list =  c('rcur', 'rpca','laplacianEigenmaps_dimRedtools', 'tsne', 'umap', 'NULL'), 
#      regenerate_DR_data = list(miRna = FALSE, RnaSeq = FALSE, RPPA = FALSE, methy = FALSE), # to regenerate some data
#      data_integration_list = c('concatenation', 'SNF_nonorm','mofa', 'mkl', 'mofaPT',  'SNF_nonormPT', 'mklPT'), # se Lenovo 'here we ask to also integrate pt data
#      regenerate_integration_list = list('concatenation'=FALSE, 'SNF_nonorm'=FALSE, 'SNF_nonormPT'=FALSE,  'mofa'=FALSE, 
#                                         'mofaPT'=FALSE, 'mkl'=FALSE,  'mklPT'=FALSE),
#      recompute_RF = FALSE, 
#      only_embed_dataset = FALSE, only_integrate_views = FALSE)


# For mofa
# only ID
# main_best_models(str_add_save = 'compNew_speedy_', 
#      use_id_list = list(ID = list(use_id = TRUE, frac_samples = 3, minNumViews = 3)),
#      just_collect_results = FALSE,
#      n_internal_iters = 1, n_external_iters = 15, 
#      perc_samples = 1,
#      ntree = NA, mtry = NA, # use default parameters for RF
#      feature_selection_method_list = 'none',
#      unsup_feature_selection_list = c('feature_clustering', 'par_rcur', 'NULL'), 
#      best_red_list =  c('tsne','rpca','laplacianEigenmaps_dimRedtools', 'umap', 'NULL'), 
#      regenerate_DR_data = list(miRna = FALSE, RnaSeq = FALSE, RPPA = FALSE, methy = FALSE), # to regenerate some data
#      data_integration_list = c('concatenation', 'SNF_nonorm','mofa', 'mkl', 'mofaPT',  'SNF_nonormPT', 'mklPT'), # se Lenovo 'here we ask to also integrate pt data
#      regenerate_integration_list = list('concatenation'=FALSE, 'SNF_nonorm'=FALSE, 'SNF_nonormPT'=FALSE,  'mofa'=FALSE, 
#                                         'mofaPT'=FALSE, 'mkl'=FALSE,  'mklPT'=FALSE),
#      recompute_RF = FALSE, 
#      timeout_per_int_iter = 3,
#      timeout_per_ext_iter = 3,
#      only_embed_dataset = FALSE, only_integrate_views = FALSE)



# Main function for DR+data fusion pipelines
main_best_models <- function(str_add_save = 'BestModels_', data_types = c('miRna', 'RnaSeq', 'RPPA', 'methy'),  
                                     maxVars = 30000,  
                                     with_pt = TRUE,
                                     cohort_list = codes,
                                     normalize_fun = 'standard_scaler',
                                     use_id_list = list(noID_1 = list(use_id = FALSE, frac_samples = 1, minNumViews = 4), 
                                                        noID_2 = list(use_id = FALSE, frac_samples = 2, minNumViews = 4), 
                                                        ID = list(use_id = TRUE, frac_samples = NA, minNumViews = 3)),
                                     # each element of the list above is a list with three elements: use_id, frac_samples, minNunViews
                                      #  - frac_samples is related to use_id = FALSE 
                                      #  when the corresponding value in use_id_list is TRUE --> the ID is used and frac_samples is NA
                                      #  OTHERWISE 
                                      #  when the corresponding value in use_id_list is FALSE --> no ID is used to define the dimensionality of the reduced space
                                      #  In this case the lower-dim (minD) is defined as:
                                      #  if (frac_samples >= 1 ) minD = (min(N,D))-1 
                                      #  else minD = round(min(N,D)/frac_samples)
                                      #  - minNumViews is set to all views (4) when use_ID = FALSE, else we consider other combinations
                                     factor_SSS = 1, 
                                      # factor_SSS is the coefficient that I use to decide 
                                      # whether a view is to big and may benefit from the intermediate reduction. 
                                      # If D/N > factor_SSS then I use unsup FS , otherwise I don't
                                     unsup_feature_selection_list = c('feature_clustering', 'par_rcur', 'NULL'), # can't use simple RCUR here because it will select at most min(N,D) variables 
                                     best_red_list = c('rpca', 'laplacianEigenmaps_dimRedtools', 'tsne', 'umap', 'NULL'), 
                                     rcur_patience = 1000,
                                     no_neighs = ext_no_neighs,
                                     data_integration_list = c('concatenation', 'SNF_nonorm', 'mofa', 'mkl', 'mofaPT',  'SNF_nonormPT', 'mklPT'), # 'here we ask to also integrate pt data
                                     #default MOFA parameters
                                     no_mofa_factors = 15, sparse_weights = TRUE, sparse_factors = FALSE, #mofa default
                                     RF_method_list = c('RF'), #, 'extratrees'), 
                                     k_snf = 20, # default SNF parameters: k = number of neighbors
                                     sigma_snf = 0.5, # default SNF parameters: sigma = sigma value in affinity matrix
                                     t_snf = 11,  # default SNF parameters: t = number of iterations
                                     partial_views = TRUE, 
                                     feature_selection_method_list = c('RF_importance'),
                                     #feature_selection_method may also be c("corr-pearson", "corr-kendall", "corr-spearman", 'none'),
                                     perc_samples = 1, 
                                     ratio_train = 0.9,
                                     seed = ext_seed,
                                     outcome_col = 'os',
                                     neg_label = '0',
                                     pos_label = '1',
                                     measures =    c("aucpr", "auc", 'sens', 'spec', 'ppv', 'npv', 'acc', 'f'), 
                                     opti_meas =  'aucpr',
                                     weights_opti_meas = 1,
                                     n_internal_iters = 21,
                                     n_external_iters = 15,
                                     tuneRF = FALSE,
                                     mtry = NA, 
                                     ntree = NA,
                                     class_weights =  'none',
                                     balanced_bootstraps = TRUE,
                                     thr_conf = 0.05,
                                     vis_norm_data = FALSE, 
                                     ncores = max(min(my_detectCores(), detectCores()-5),1), 
                                     timeout_per_int_iter = 15,
                                     timeout_per_ext_iter = 15,
                                      just_collect_results = FALSE, # to only collect already run results
                                      only_compute_blocking = FALSE, # just compute blocking IDs
                                      only_embed_dataset = FALSE, # just arrive until embedding is done and do not perform integration and classification 
                                      only_integrate_views = FALSE, # just arrive until integration and do not perform classification 
                                      regenerate_blocking = list(miRna = FALSE, RnaSeq = FALSE, RPPA = FALSE, methy = FALSE),
                                      regenerate_DR_data = list(miRna = FALSE, RnaSeq = FALSE, RPPA = FALSE, methy = FALSE),
                                      regenerate_integration_list = list('concatenation'=FALSE, 'SNF_nonorm'=FALSE, 'SNF_nonormPT'=FALSE,  'mofa'=FALSE, 
                                                                         'mofaPT'=FALSE, 'mkl'=FALSE,  'mklPT'=FALSE),
                                      recompute_RF = FALSE){ 

        
    
        str_add_save = paste0(ifelse(str_add_save  != '', paste0(str_add_save ,'_'),''), 
                              'useID_', paste(names(use_id_list), collapse = '-')) 
    
        
        clin_vars = c('years_to_birth','gender', 'race', 'ethnicity', 'patient.age_at_initial_pathologic_diagnosis')
       # outcome_col <- 'pfi'
        center_col <- 'COHORT'
        corr_filter = TRUE
        
        dist_fun_twoNN = 'canberra'
        perc_points = 0.9
        maxit = 11
        factor_blocking = 3
        max_blocking_runs = list(miRna = 100, RnaSeq = 100, RPPA = 100, methy = 100)
        ntry = list(miRna = 21, RnaSeq = 21, RPPA = 21, methy = 21)
        no_els_to_monitor = 10
        
        ID_estimator_fun = 'estimate_ID_twoNN'
        args_list = list('estimate_ID_danco' = list(list_k = c(4,6,8,12,16), maxD = 100, perc_points = perc_points, maxit = maxit, 
                                                    ncores = min(ncores, detectCores()-1)) ,    
                                         'estimate_ID_twoNN' = list(dist_fun_twoNN = 'canberra', perc_points = perc_points, maxit = maxit, 
                                                                    ncores = min(ncores, detectCores()-1)))
        
        #serve per limitare il numero di variabili tolte dalla p_value_select
        limitVar = list(miRna = 350, RnaSeq = 3000, RPPA = 3000, methy = 3000)
        
        
        #normalize_data_fun = list(miRna = min_max_norm, RnaSeq =    min_max_norm, RPPA = NULL, methy =    min_max_norm)
        normalize_data_fun = list(miRna = my_quantile_normalization, 
                                                            RnaSeq =    my_quantile_normalization, 
                                                            RPPA = my_quantile_normalization, 
                                                            methy =    my_quantile_normalization
        )
        
        
        dim_split_corr = list(miRna = 505, RnaSeq = 1000, RPPA = 1000, methy = 1000)
        cutoff = list(miRna = 0.8, RnaSeq = 0.8, RPPA = 0.8, methy = 0.8)
        maxit_corr = list(miRna = 11, RnaSeq = 11, RPPA = 11, methy = 11)
        method_corr = list(miRna = 'spearman', RnaSeq = 'spearman', RPPA = 'spearman', methy = 'spearman')
        maxDim = list(miRna = NA, RnaSeq = NA, RPPA = NA, methy = NA)
        minD = list(miRna = NA, RnaSeq = NA, RPPA = NA, methy = NA)
        dim_split_feat_clustering = list(miRna = 75, RnaSeq = 75, RPPA = 75, methy = 75)
        maxit_unsup_feature_selection = list(miRna = 5, RnaSeq = 5, RPPA = 5, methy = 5)
        
        in_all_df = NULL
        # fracVar = 0.9
        
        #    embedded_data_types_fn = file.path('embedded_data_types.Rda')
        #    explain_data_types_fn = file.path('explain_data_types.Rda')
        all_centers_res = data.frame()
        all_centers_mega_res = data.frame()
        
        uncompleted = matrix(0, nrow = length(cohort_list)*length(use_id_list), ncol = length(data_integration_list))
        not_converged = matrix(0, nrow = length(cohort_list)*length(use_id_list), ncol = length(data_integration_list))
        colnames(uncompleted) = data_integration_list
        colnames(not_converged) = data_integration_list
        rownames(uncompleted) = unlist(lapply(cohort_list,function(x) paste0(x, '_', names(use_id_list))))
        rownames(not_converged) = unlist(lapply(cohort_list,function(x) paste0(x, '_', names(use_id_list))))
        
        
        for (centers in cohort_list){
            # if we have more centers, the data of the second, third, ... centers will be aligned to the data of the first center
            str_centers = paste(centers, sep = '', collapse = '_')    
            cat('\n******+++++++++++++++++++++++++++++++++*****\n')
            
            
            cat('starting to process ', str_centers, '\n')
            
            
            if (length(centers)==1){ norm_quantiles = FALSE }else{ norm_quantiles = TRUE }
            
            #res_dir = file.path('current_exp', 'data_exp', paste('June_res', str_centers,    sep = '_'))
            res_dir = file.path('./results', paste('res', str_centers, sep = '_'))
            if (!dir.exists(res_dir)) dir.create(res_dir)
            
            # riikka_dir = file.path(res_dir, 'riikka')
            # if (!dir.exists(riikka_dir)) dir.create(riikka_dir)
            
            
            ### files xlsx for saving RF results
            selected_file_mean_results_RF = file.path(res_dir, paste(str_add_save, 'mean_results_selected.xlsx', sep =''))
            selected_file_mega_results_RF = file.path(res_dir, paste(str_add_save, 'mega_results_selected.xlsx', sep =''))
            
            selected_mean_res_df = data.frame()
            selected_mega_res_df = data.frame()
            
            wkl_weights = data.frame()
           
            
            ###################
            
            df_pt = data.frame()
            new_names = NULL
            for (center_now in centers){
                if (file.exists(file.path(dati_rinominati, paste0(center_now, '_', outcome_col,'.rds')))){
                    df_c = data.frame(readRDS(file = file.path(dati_rinominati, paste0(center_now, '_', outcome_col,'.rds'))))
                }else{
                    warning(paste0('NESSUN OUTCOME TROVATO PER ', center_now, ' - SKIPPING CENTER\n'))
                    next
                }
                names(df_c) = outcome_col
                df_c[[center_col]] = center_now
                df_clin = data.frame(readRDS(file = file.path(dati_rinominati, paste0(center_now, '_clinics.rds'))))
                df_clin = df_clin[, names(df_clin) %in% clin_vars]
                df_c = merge(df_c, df_clin, by = 'row.names')
                row.names(df_c) = df_c$Row.names
                # plyr:: rbind.fill perde i rownames!
                new_names = c(new_names, row.names(df_c))
                df_pt = plyr::rbind.fill(df_pt, df_c[, !(names(df_c) %in% 'Row.names')])
                row.names(df_pt) = new_names
                df_c = NULL
                
            } 
            
            df_pt[[outcome_col]] = factor(df_pt[[outcome_col]], levels = c(neg_label, pos_label))
            if (any(grepl('gender', names(df_pt)))) df_pt$gender = factor(df_pt$gender, levels = c('male', 'female'))
            if (any(grepl('race', names(df_pt)))) df_pt$race = factor(df_pt$race, levels = c('white', 'black or african american', 'asian'))
            if (any(grepl('ethnicity', names(df_pt)))) df_pt$ethnicity = factor(df_pt$ethnicity, levels = c('not hispanic or latino','hispanic or latino'))
            
            for (cc in names(df_pt)){
                cat(cc, ' - ')
                if (length(unique(df_pt[[cc]]))==1) df_pt = df_pt[, !(names(df_pt) %in% cc)]
            } 
            
            rm(cc)
            
            df_outcome = data.frame(df_pt[[outcome_col]], row.names = row.names(df_pt))
            # dataframe che poi andra' a comporre il dataframe dei dati combinati da passare alla RF
            names(df_outcome) = outcome_col
            df_pt_use = df_pt[, !(names(df_pt) %in% c(outcome_col, 'COHORT'))]
            vars_pt = names(df_pt_use)[!names(df_pt_use) %in% outcome_col]
            
            
            for (exp_id_no in 1:length(use_id_list)){
                exp_id = names(use_id_list)[exp_id_no] # I'll need it later to record the uncompleted or unconverged experiments
                
                use_id = use_id_list[[exp_id]][['use_id']]
                frac_samples = use_id_list[[exp_id]][['frac_samples']]
                minNumViews = use_id_list[[exp_id]][['minNumViews']]
                for (unsup_feature_selection in unsup_feature_selection_list){
                    ID_string_unsupFS = ifelse(use_id, '_ID', paste0('_noID_', frac_samples)) # when I use the ID the unsup feature selection process is equal for all frac_samples    
                    
                    unsup_feature_selection = paste0(unsup_feature_selection, ID_string_unsupFS)
                    
                    for (best_red in best_red_list){    
                        ID_string_red = paste0('_', ifelse(use_id, '', 'no'), 'ID_', frac_samples)
                      
                        if (!use_id){ 
                            # if I'm not using the ID the last dimensionality will surely be lower than N; therefore, no need to use par_rcur
                            if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') =='par_rcur') next
                      
                            # when not using the ID to choose the dimensionality only 1-shot DR is chosen to avoid empirically choosing also the intermediate step  
                            if ((gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') != 'NULL') & (best_red != 'NULL')) next
                            
                        }
                      
                      # if both are null I skip the experiment bacause it is too demanding
                        if (( gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'NULL') & (best_red == 'NULL')) next
                      
                      # the combo:
                      # FS = 'rcur' (or 'cur') + FE = 'NULL' 
                      # is equal to the inverse combo:
                      # FS = 'NULL' + FE = 'rcur' (or cur); if both of them are scheduled I'm jumping the first of them 
                        if ((grepl('cur', gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '')) & (best_red == 'NULL')) & 
                            any(best_red_list == gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '')) & 
                                 any(gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'NULL')) next
                      
                        
                       
                        if (!(best_red == 'NULL')){
                            best_red = paste0(best_red, ID_string_red) 
                        }
                        
                        mat_data_all_names = data.frame('pt' = row.names(df_outcome), row.names =  row.names(df_outcome))
                        
                        
                        # per problemi di psazio non posso tenerle in memoria ma nelle liste seguenti salvo solo la posizione dei df
                        explain_data_types <- list()
                        embedded_data_types <- list()
                        
                        embedded_data_types[['pt']] = df_pt_use
                        explain_data_types[['pt']] = df_pt_use
                        
                        if (verbose ==1){ 
                          cat('******\n')
                          cat('start with setting =  - unsup feature selection = ', unsup_feature_selection, 
                                ' - best_red = ', best_red, ' - use_id = ', use_id, '\n')
                        }
                        
                        # dati gi? standardizzati!!
                        norm_omics = readRDS(file = file.path(dati_rinominati, paste0(centers, '_norm_omics.rds')))
                        
                        
                        for (str_desc in data_types){
                            
                            cat('\n\n*****************---------********************\n')
                            
                            name_data = names(norm_omics)[grepl(paste0('_', str_desc), names(norm_omics), ignore.case = TRUE)]
                            #str_desc = names(norm_omics)[grepl(str_desc, names(norm_omics), ignore.case = TRUE)]
                                                        
                            cat('processing ', str_desc, ' data type\n')
                           # unsup_feature_selection = unsup_feature_selection
                            best_red_now = best_red
                            
                            
                            save_ID_list = list()
                            save_ID_list[['ID_est']] = NULL
                            save_ID_list[['sd_ID_est']] = NULL
                            save_ID_list[['blocking_estimates']] = list()
                            
                            cat('-------------------------------------------------------\n')
                            task = 'clean_data'
                            
                            cat(task, '\n')
                            # ognuno segue il suo percorso - 
                            my_path_records = list()
                            my_path_records[[task]] = res_dir
                            
                            file_data_fn = file.path(my_path_records[[task]], paste(str_desc, '_', task,'.Rda' , sep ='' ))
                            
                            #################### OPENING AND CLEANING DATA
                            if (file.exists(file_data_fn)){ 
                              load(file_data_fn)
                              #save(mat_data, file = file_data_fn)
                            }
                            else mat_data = norm_omics[[name_data]]
                            # se ? troppo grosso riduco con entropia (dovrei poi ricomputare la blocking id)
                            if (ncol(mat_data)>maxVars){
                                cat('data is hyper-diper high dimensional!\n')
                                file_data_fn_entropy_reduction = file.path(my_path_records[[task]], 
                                                                           paste(str_desc, '_entropy_filt_after_opening.Rda' , sep ='' ))
                                if (file.exists(file_data_fn_entropy_reduction)){ 
                                    cat('loading lower dimensional data \n')
                                    load(file_data_fn_entropy_reduction)
                                }else{
                                    cat('entropy filtering to reduce dimensionality to ', maxVars, ' variables \n')
                                    mat_data = imputa(mat_data)
                                    entp_feat = apply(mat_data,2, entropy, method="Laplace")
                                    entp_feat = sort(entp_feat, decreasing = TRUE)
                                    mat_data = mat_data[ , match(names(entp_feat)[1:maxVars], colnames(mat_data))]
                                    save(mat_data, file = file_data_fn_entropy_reduction)
                                }
                            }
                            
                            
                            
                            if ((norm_quantiles) & (length(centers)>1)){
                                file_data_fn = file.path(my_path_records[[task]], paste(str_desc, '_', task, '_quantile_norm.Rda' , sep ='' ))
                                if (file.exists(file_data_fn)){
                                    load(file = file_data_fn)
                                    #    load(file = ID_fn)
                                    #    blocking_mat = blocking_estimates[[task]]
                                    cat('LOADED ', str_desc, ' mat_data after ', task, ' (quantile_normalized) = ', dim(mat_data), ' \n')
                                }else{
                                    
                                    mat_data = imputa(mat_data)
                                    #apply quantile norm to target center using the first center in the list of centers as target
                                    target_center = centers[1]
                                    idx_in_pt_df = match(row.names(mat_data), row.names(df_pt))
                                    center_data = df_pt[idx_in_pt_df, center_col]
                                    target_length = NULL
                                    trasf_data = function(x) return(x)
                                    norm_quant_data = TRUE
                                    if (str_desc == "RPPA"){ 
                                        norm_quant_data = FALSE 
                                    }
                                    
                                    cat('normalize centers\n')
                                    if (norm_quant_data){
                                        for (cc in 1:ncol(mat_data)){
                                            if (vis_norm_data){ 
                                                df_s1 = data.frame(data = mat_data[, cc], aes = paste(center_data, '_before'))
                                                h1 = ggplot(df_s1, aes(x = data, color=aes, fill = aes)) + geom_density(alpha = 0.2)+
                                                    ggtitle('before normalization')
                                                print(h1)
                                            }
                                            target_data = trasf_data(normalize.quantiles.determine.target(
                                                as.matrix(mat_data[(center_data == target_center),cc]),
                                                target.length = target_length))
                                            mat_data[!(center_data == target_center),cc] = trasf_data(
                                                normalize.quantiles.use.target(as.matrix(mat_data[!(center_data == target_center),cc]), target_data))
                                            if (vis_norm_data){ 
                                                df_s1 = rbind(df_s1,data.frame(data = mat_data[, cc], aes = (center_data)))
                                                h2 = ggplot(df_s1, aes(x = data, color=aes, fill = aes)) + geom_density(alpha = 0.1)+
                                                    ggtitle('after normalization')
                                                print(h2)
                                            }
                                        }
                                      rm(cc)
                                    }
                                    save(mat_data, file = file_data_fn)
                                    cat('SAVED ', str_desc, ' mat_data after ', task, '    (quantile normalized) = ', dim(mat_data), ' \n')    
                                }
                            }
                            
                            
                            #################### END OPENING AND CLEANING DATA
                            
                            ID_fn = file.path(my_path_records[[task]], 
                                                                paste0('ID_estimates', ifelse(use_id, '_ID', '_noID'), '_' , factor_blocking,'_', str_desc, '.Rda'))
                            
                            if (use_id){
                                
                                if (file.exists(ID_fn) & !(regenerate_blocking[[str_desc]])){
                                    load(ID_fn)
                                  
                                  ############## POI DA TOGLIERE!!!!
                                  # save_ID_list = save_ID_list[!(names(save_ID_list) %in% c('maxDim','minD'))]
                                  # save(save_ID_list, file = ID_fn)
                                  # 
                                  
                                    cat('********************\nRETRIEVING SAVED BLOCK ESTIMATES \n ')
                                    if (length(save_ID_list[['blocking_estimates']])>0){ 
                                        
                                        blocking_mat = save_ID_list[['blocking_estimates']][[task]]
                                        # # check whether a plateau has been found or, instead, there is no redundancy in the features
                                        # converged = FALSE
                                        # dim_convergence = NA
                                        # dim_lt_1std = NA
                                        
                                        converged = save_ID_list$blocking_estimates$convergence
                                        dim_convergence = save_ID_list$blocking_estimates$dim_convergence
                                        dim_lt_1std = save_ID_list$blocking_estimates$dim_lt_1std
                                        
                                        
                                        # if (nrow(blocking_mat)>no_els_to_monitor){  # should already be!
                                        #   # compute the std of the last 5 elements. When the std falls belo 0.1 stop computation and take the last block id as the last dataset id
                                        #   
                                        #   for (j in (no_els_to_monitor+1):nrow(blocking_mat)){
                                        #     last_els = blocking_mat[(j-no_els_to_monitor):j,'block_id']
                                        #     std_last_els = sd(last_els, na.rm = TRUE)
                                        #     lag = 3
                                        #     sum_deriv = sum(diff(last_els, lag = lag)/lag) # per la derivate non faccio il modulo perchè mi va bene se scende e poi risale
                                        #     if (verbose == 1) cat('std_last_els = ', std_last_els, ' - sum_deriv = ' , sum_deriv, '\n')
                                        #     if ((std_last_els < 0.25) | (sum_deriv < 0.05)){ 
                                        #       if (verbose == 1)  cat('convergence reached at ', j, ' steps, that is when the block has dimension ', blocking_mat[j, 'no_vars_in_block'], '\n')
                                        #       converged = TRUE
                                        #       dim_convergence = blocking_mat[j, 'no_vars_in_block']
                                        #       break
                                        #     }
                                        #   }
                                        # }
                                        
                                        
                                        
                                        if (converged){
                                          ####################### CHOOSING INTERMEDIATE DIMENSION
                                          if (verbose ==1) cat('BLOCKING ANALYSIS FOR ', str_desc,  '-', paste(centers, collapse = '+'),'CONVERGED - setting dim_convergence as maxDim\n')
                                          target = blocking_mat[nrow(blocking_mat), 'mean_id']
                                          target_sd = blocking_mat[nrow(blocking_mat), 'total_sd']
                                          lower_bound = target_sd/2
                                          idx_block = which((blocking_mat[, 'mean_id']> (target-lower_bound )) & 
                                                              (blocking_mat[, 'mean_id'] < (target+lower_bound )))[1]
                                          dim_lt_1std =  blocking_mat[idx_block, 'no_vars_in_block']
                                          # intermediate dimension is the max between the intermediate dimension estimated 
                                          # with blocking-distribution analysis and 
                                          # e the number of cases = N
                                          cat('dim_at_lt_1std\tdim_convergence\tmin(dim)\n')
                                          cat(dim_lt_1std, '\t', dim_convergence, '\t', min(dim(mat_data))-1, '\n')
                                          
                                          maxDim[[str_desc]] = max(dim_convergence, nrow(mat_data)-1)  # anyhow, I do not want it to be less than the sample size
                                          
                                        }else{
                                          dim_convergence = NA
                                          dim_lt_1std = NA
                                          
                                          cat('BLOCKING ANALYSIS FOR ', str_desc,  '-', paste(centers, collapse = '+'), ' did not converge; no FS applied\n' )  
                                          
                                          ID_orig = save_ID_list$ID_est[1]
                                          sd_ID_orig = save_ID_list$sd_ID_est[1]
                                          
                                          target = round(ID_orig, digits = 2)
                                          target_sd = round(sd_ID_orig, digits = 2)
                                          maxDim[[str_desc]] = ncol(mat_data) # do not apply intermediate reduction
                                        }
                                        
                                        # save_ID_list$blocking_estimates$convergence = converged
                                        # save_ID_list$blocking_estimates$dim_convergence = dim_convergence
                                        # save_ID_list$blocking_estimates$dim_lt_1std = dim_lt_1std
                                        # 
                                        # save(save_ID_list, file = ID_fn)
                                        
                                        
                                        ####### plot blocking_ID
                                        # centers_str = paste(centers, sep ='_')
                                        # png(filename = paste0('C:/DATI/Anacleto/similarities_PNet/missSNF/data_TCGA/AARisultatiAnalizzati_per_articolo/images/', 
                                        #                       centers_str, '_', str_desc, '.png'))
                                        # 
                                        # df_blocking = data.frame(blocking_mat[, match(c('no_vars_in_block', 'block_id', 'mean_id', 'sd_block_id' , 'total_sd'), colnames(blocking_mat))]  )
                                        # df_blocking = rbind(rep(0,ncol(df_blocking)), df_blocking)
                                        # upper = df_blocking[nrow(df_blocking), 'mean_id'] - df_blocking[nrow(df_blocking), 'total_sd'] 
                                        # lower = df_blocking[nrow(df_blocking), 'mean_id'] + df_blocking[nrow(df_blocking), 'total_sd']
                                        # 
                                        # names(df_blocking) = c('Lj', paste0('block_id (std) = ', round(df_blocking[nrow(df_blocking), 'block_id'],0),' (', round(df_blocking[nrow(df_blocking), 'sd_block_id'],1), ')'), 
                                        #                        paste0('cumulative mean (pooled std) = ', round(df_blocking[nrow(df_blocking), 'mean_id'],0),' (', round(df_blocking[nrow(df_blocking), 'total_sd'],1), ')'),  
                                        #                        'sd_block_id', 'total_sd')
                                        # 
                                        # 
                                        # df_block_id = cbind('block-IDs (std.)', df_blocking[, c(1,2,4)])
                                        # df_cumulative = cbind('cumulative mean (total std.)', df_blocking[, c(1,3,5)])
                                        # names(df_block_id) = c('estimate', 'Lj', 'ID', 'std(ID)')
                                        # names(df_cumulative) = c('estimate', 'Lj', 'ID', 'std(ID)')
                                        # 
                                        # df_ggplot = rbind( df_block_id, df_cumulative )
                                        # df_ggplot$upper = upper
                                        # df_ggplot$lower = lower
                                        # 
                                        # 
                                        # colors_to_use = NULL
                                        # out2 = ggplot(df_ggplot, aes(x=Lj, y=ID, col=estimate)) + geom_line() +geom_point()
                                        # out2 = out2 +  geom_errorbar(aes(ymin= .data[['ID']]-.data[['std(ID)']], 
                                        #                                  ymax= .data[['ID']]+.data[['std(ID)']], colour = estimate),
                                        #                              width=5,position=position_dodge(0.05))
                                        # 
                                        # # out2 <- out2 + geom_ribbon(aes(ymin= upper, 
                                        # #                                ymax= lower), alpha = 0.1, colour = NA)
                                        # out2 = out2 + labs(title = 
                                        #                      paste0('ID estimation: block analysis for ', str_desc, '\n(', centers_str,') after preprocessing'), 
                                        #                    x = 'number of randomly selected features', 
                                        #                    y = 'block-ID with cumulative mean and total std. error', color = 'Legend')
                                        # 
                                        # 
                                        # if (converged){
                                        #   vdf = data.frame(xintercept_name = c('dim at convergence', 'dim at 1 std away from ID estimate'), 
                                        #              xintercept = c(dim_convergence,dim_lt_1std), stringsAsFactors = FALSE)
                                        #   out2 = out2+geom_vline(mapping = aes(xintercept = xintercept, colour = xintercept_name), 
                                        #                          data = vdf, show.legend = TRUE)
                                        #           
                                        #   # out2 = out2+geom_vline(xintercept = dim_convergence, show.legend = TRUE, colour="red")
                                        #   # out2 = out2+geom_vline(xintercept = dim_lt_1std, show.legend = TRUE, colour = "grey")
                                        #   # out2 = out2 + annotate("text", dim_convergence, 1, hjust = -.01, 
                                        #   #                        label = 'dim at convergence') + annotate("text", dim_lt_1std, dim_convergence, 1, hjust = -.01, 
                                        #   #                                                                 label = 'dim at lt 1 std')
                                        # }
                                        # 
                                        # print(out2)
                                        # dev.off()
                                        
                                        # ############### end plotting
                                        
                                    }else{ 
                                        ID_orig = save_ID_list$ID_est[1]
                                        sd_ID_orig = save_ID_list$sd_ID_est[1]
                                        
                                        target = round(ID_orig, digits = 2)
                                        target_sd = round(sd_ID_orig, digits = 2)
                                        maxDim[[str_desc]] = ncol(mat_data) # do not apply intermediate reduction
                                    }
    
                                    
                                  
                                }else{
                                  
                                    cat('saving ID estimates in path ', ID_fn, '\n')
                                  
                                    ids_twonn = estimate_ID_twoNN(mat_data, 
                                                                  dist_fun_twoNN, 
                                                                  maxit = maxit, perc_points = perc_points)
                                    
                                    
                                    cat('estimated ids on ', str_desc, ' = ', ids_twonn[['id']], '\n')
                                    
                                    save_ID_list[['df_estimates']] = rbind(data.frame(k = 2, IDs = ids_twonn[['id']],
                                                                                        low = ids_twonn[['id']]-ids_twonn[['sd_id']], 
                                                                                        up = ids_twonn[['id']]+ids_twonn[['sd_id']], 
                                                                                        ID_method = ID_estimator_fun))
                                    save_ID_list[['ID_est']][task] =    ids_twonn[['id']]
                                    save_ID_list[['sd_ID_est']][task] = ids_twonn[['sd_id']]
                                    
                                    # data_norm = standardNormalization(mat_data)
                                    # dist_mat = (dist2(as.matrix(mat_data),as.matrix(mat_data)))^(1/2)
                                    # h = hist(dist_mat[upper.tri(dist_mat, diag=FALSE)])
                                    # ent = entropy(h$density)
                                    
                                    ID_orig = ids_twonn[['id']]
                                    sd_ID_orig = ids_twonn[['sd_id']]
                                    if (round(ncol(mat_data) / nrow(mat_data),1) > factor_SSS){
                                    # if we are in the case of small-sample-size                                
     
                                        #BLOCKING ID
                                        blocking_mat = blocking_ID(mat = mat_data, 
                                                                 ID_estimator_fun = ID_estimator_fun, 
                                                                 args_ID = args_list[[ID_estimator_fun]], 
                                                                 save_path = my_path_records[[task]], 
                                                                 task = task, 
                                                                 str_desc = str_desc, 
                                                                 ID_orig = ids_twonn[['id']],
                                                                 factor = factor_blocking,
                                                                 max_blocking_runs = max_blocking_runs[[str_desc]], 
                                                                 ntry = ntry[[str_desc]],
                                                                 no_els_to_monitor = no_els_to_monitor,
                                                                 verbose = 0, 
                                                                 ncores = min(ncores, detectCores()-1))
                                        save_ID_list[['blocking_estimates']][[task]] = blocking_mat
                                        save(save_ID_list, file = ID_fn)
                                        
                                        ################## END blocking ID
                                    
                                        # check whether a plateau has been found or, instead, there is no redundancy in the features
                                        converged = FALSE
                                        dim_convergence = NA
                                        dim_lt_1std = NA
                                        
                                        if (nrow(blocking_mat)>no_els_to_monitor){  # should already be!
                                          # compute the std of the last 5 elements. When the std falls belo 0.1 stop computation and take the last block id as the last dataset id
                                          
                                          for (j in (no_els_to_monitor+1):nrow(blocking_mat)){
                                            last_els = blocking_mat[(j-no_els_to_monitor):j,'block_id']
                                            std_last_els = sd(last_els, na.rm = TRUE)
                                            lag = 3
                                            sum_deriv = sum(diff(last_els, lag = lag)/lag) # per la derivate non faccio il modulo perchè mi va bene se scende e poi risale
                                            if (verbose == 1) cat('std_last_els = ', std_last_els, ' - sum_deriv = ' , sum_deriv, '\n')
                                            if ((std_last_els < 0.25) | (sum_deriv < 0.05)){ 
                                              if (verbose == 1)  cat('convergence reached at ', j, ' steps, that is when the block has dimension ', blocking_mat[j, 'no_vars_in_block'], '\n')
                                              converged = TRUE
                                              dim_convergence = blocking_mat[j, 'no_vars_in_block']
                                              break
                                            }
                                          }
                                        }
                                        
                                        
                                        if (converged){
                                          ####################### CHOOSING INTERMEDIATE DIMENSION
                                          target = blocking_mat[nrow(blocking_mat), 'mean_id']
                                          target_sd = blocking_mat[nrow(blocking_mat), 'total_sd']
                                          lower_bound = target_sd/2
                                          idx_block = which((blocking_mat[, 'mean_id']> (target-lower_bound )) & 
                                                              (blocking_mat[, 'mean_id'] < (target+lower_bound )))[1]
                                          dim_lt_1std =  blocking_mat[idx_block, 'no_vars_in_block']
                                          # intermediate dimension is the max between the intermediate dimension estimated 
                                          # with blocking-distribution analysis and 
                                          # e the number of cases = N
                                          cat('dim_at_lt_1std\tdim_convergence\tmin(dim)\n')
                                          cat(dim_lt_1std, '\t', dim_convergence, '\t', min(dim(mat_data))-1, '\n')
                                          
                                          maxDim[[str_desc]] = max(dim_convergence, nrow(mat_data)-1)  # anyhow, I do not want it to be less than the sample size
                                          
                                          
                                        }else{
                                          dim_convergence = NA
                                          dim_lt_1std = NA
                                          
                                          cat('BLOCKING ANALYSIS FOR ', str_desc,  '-', paste(centers, collapse = '+'), ' did not converge; no FS applied' )  
                                          
                                          ID_orig = save_ID_list$ID_est[1]
                                          sd_ID_orig = save_ID_list$sd_ID_est[1]
                                          
                                          target = round(ID_orig, digits = 2)
                                          target_sd = round(sd_ID_orig, digits = 2)
                                          maxDim[[str_desc]] = ncol(mat_data) # do not apply intermediate reduction
                                        }
                                        
                                        save_ID_list$blocking_estimates$convergence = converged
                                        save_ID_list$blocking_estimates$dim_convergence = dim_convergence
                                        save_ID_list$blocking_estimates$dim_lt_1std = dim_lt_1std
                                        
                                        save(save_ID_list, file = ID_fn)
                                        
                                        
                                        ####### plot blocking_ID
                                        centers_str = paste(centers, sep ='_')
                                        if (!dir.exists('./results/images/')) dir.create('./results/images/')
                                        png(filename = paste0('./results/images/', 
                                                              centers_str, '_', str_desc, '.png'))
                                        
                                        df_blocking = data.frame(blocking_mat[, match(c('no_vars_in_block', 'block_id', 'mean_id', 'sd_block_id' , 'total_sd'), colnames(blocking_mat))]  )
                                        df_blocking = rbind(rep(0,ncol(df_blocking)), df_blocking)
                                        upper = df_blocking[nrow(df_blocking), 'mean_id'] - df_blocking[nrow(df_blocking), 'total_sd'] 
                                        lower = df_blocking[nrow(df_blocking), 'mean_id'] + df_blocking[nrow(df_blocking), 'total_sd']
                                        
                                        names(df_blocking) = c('Lj', paste0('block_id (std) = ', round(df_blocking[nrow(df_blocking), 'block_id'],0),' (', round(df_blocking[nrow(df_blocking), 'sd_block_id'],1), ')'), 
                                                               paste0('cumulative mean (pooled std) = ', round(df_blocking[nrow(df_blocking), 'mean_id'],0),' (', round(df_blocking[nrow(df_blocking), 'total_sd'],1), ')'),  
                                                               'sd_block_id', 'total_sd')
                                        
                                        
                                        df_block_id = cbind('block-IDs (std.)', df_blocking[, c(1,2,4)])
                                        df_cumulative = cbind('cumulative mean (total std.)', df_blocking[, c(1,3,5)])
                                        names(df_block_id) = c('estimate', 'Lj', 'ID', 'std(ID)')
                                        names(df_cumulative) = c('estimate', 'Lj', 'ID', 'std(ID)')
                                        
                                        df_ggplot = rbind( df_block_id, df_cumulative )
                                        df_ggplot$upper = upper
                                        df_ggplot$lower = lower
                                        
                                        colors_to_use = NULL
                                        out2 = ggplot(df_ggplot, aes(x=Lj, y=ID, col=estimate)) + geom_line() +geom_point()
                                        out2 = out2 +  geom_errorbar(aes(ymin= .data[['ID']]-.data[['std(ID)']], 
                                                                         ymax= .data[['ID']]+.data[['std(ID)']], colour = estimate),
                                                                     width=5,position=position_dodge(0.05))
                                        
                                        # out2 <- out2 + geom_ribbon(aes(ymin= upper, 
                                        #                                ymax= lower), alpha = 0.1, colour = NA)
                                        out2 = out2 + labs(title = 
                                                             paste0('ID estimation: block analysis for ', str_desc, '\n(', centers_str,') after preprocessing'), 
                                                           x = 'number of randomly selected features', 
                                                           y = 'block-ID with cumulative mean and total std. error', color = 'Legend')
                                        
                                        
                                        if (converged){
                                          vdf = data.frame(xintercept_name = c('dim at convergence', 'dim at 1 std away from ID estimate'), 
                                                           xintercept = c(dim_convergence,dim_lt_1std), stringsAsFactors = FALSE)
                                          out2 = out2+geom_vline(mapping = aes(xintercept = xintercept, colour = xintercept_name), 
                                                                 data = vdf, show.legend = TRUE)
                                          
                                          # out2 = out2+geom_vline(xintercept = dim_convergence, show.legend = TRUE, colour="red")
                                          # out2 = out2+geom_vline(xintercept = dim_lt_1std, show.legend = TRUE, colour = "grey")
                                          # out2 = out2 + annotate("text", dim_convergence, 1, hjust = -.01, 
                                          #                        label = 'dim at convergence') + annotate("text", dim_lt_1std, dim_convergence, 1, hjust = -.01, 
                                          #                                                                 label = 'dim at lt 1 std')
                                        }
                                        
                                        print(out2)
                                        dev.off()
                                        
                                        # ############### end plotting
                                        
                                  }else{
                                      # there is no suspect of small-sample-size
                                      # cat(str_desc, ': number of columns ( D = ', ncol(mat_data),') is comparable to number of samples ( N = ', 
                                      #         nrow(mat_data), ' --> D/N = )', round(ncol(mat_data)/nrow(mat_data), 1) , ' - no need for intermediate step\n')
                                      # 
                                      warning(str_desc, ': number of columns ( D = ', ncol(mat_data),') is comparable to number of samples ( N = ', 
                                                 nrow(mat_data), ' --> D/N = )', round(ncol(mat_data)/nrow(mat_data), 1) , 
                                                ' - no need for intermediate step\n', immediate. = TRUE)
                                      
                                      target = round(ID_orig, digits = 2)
                                      target_sd = round(sd_ID_orig, digits = 2)
                                      maxDim[[str_desc]] = ncol(mat_data) # do not apply intermediate reduction
                                  }
                                } # end if (file.exists(ID_fn))
                                
                                # If there is no feature extraction I keep only the most representative features
                                if ((gsub(best_red_now, pattern = ID_string_red, replacement='') =='NULL') | 
                                    (grepl('cur', best_red_now))) minD[[str_desc]] = maxDim[[str_desc]]
                                else minD[[str_desc]] = ceiling(frac_samples*target +(target_sd*3))
                                
                            }else{ # NOT USING THE ID
                                cat('NOT USING ID\n')
                                blocking_mat = NULL
                                target = ifelse(frac_samples==1, min(dim(mat_data))-1, round(min(dim(mat_data))/frac_samples))
                                target_sd = 0        
                                maxDim[[str_desc]] = target # either apply FS or FE; anyhow go down to target
                                minD[[str_desc]] = target
                            } # end if (use_id)
                            
                            
                            
                            if (only_compute_blocking) next
                            
                            if (!(gsub(pattern = ID_string_unsupFS, replacement = '', unsup_feature_selection) == 'NULL')){
                                old_path = my_path_records[[task]]
                                str_add_fs = ''
                                maxDim_now = maxDim[[str_desc]]
                                if (use_id){
                                    cat('FOR ', str_desc, 'I should select ', maxDim[[str_desc]], ' features by FS!!\n')
                                }else{
                                    cat('no ID reduction maxDim == ncol(mat_data)?\n', maxDim_now == ncol(mat_data))
                                }

                                task = unsup_feature_selection
                                my_path_records[[task]] = file.path(old_path, task)
                                if (!file.exists(my_path_records[[task]])){ dir.create(my_path_records[[task]])}
                                
                                
                                file_data_fn = file.path(my_path_records[[task]], paste(str_desc, str_add_fs, '.Rda' , sep ='' ))
                                cat('-------------------------------------------------------\n')
                                
                                cat(task, '\n')
                                
                                if (file.exists(file_data_fn) & (!regenerate_DR_data[[str_desc]])){
                                    load(file = file_data_fn)
                                    cat('loading ', task,    ' data with dim', dim(mat_data), '\n')
                                }else{
                                    cat(task,    'until dimension is greater than ', maxDim_now,    '\n')
                                    if (maxDim_now == ncol(mat_data)){
                                            warning(str_desc, ": intermediate dimenson - ", maxDim_now, " == ", ncol(mat_data), 
                                                " - mat dimension ! mat_data remains the same!\n", immediate. = TRUE)
                                           # cat(str_desc, ": intermediate dimenson - ", maxDim_now, " == ", ncol(mat_data), 
                                          #    " - mat dimension ! mat_data remains the same!\n")
                                    }else if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'feature_clustering'){
                                            cat("feature clustering\n")
                                            mat_data = cluster_features_par(mat_data, dim_split = 1000, maxD = maxDim_now, 
                                                                            maxit = ifelse(use_id, 3, 1),
                                                                            method = method_corr[[str_desc]], 
                                                                            ncores = min(ncores, detectCores()-1))
                                    }else if(gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'svd_entropy'){
                                            cat("parallel svd selection\n")
                                            mat_data = hierarchical_svd_entropy_par(mat_data, 
                                                                                    dim_split = 10000, 
                                                                                    maxD = maxDim_now,
                                                                                    maxit = maxit_unsup_feature_selection[[str_desc]],
                                                                                    ncores = min(ncores, detectCores()-1))
                                    }else if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'par_rcur'){
                                            cat("parallel RCUR selection\n")
                                            mat_data = RCUR_selection_par(mat_data, patience = rcur_patience, maxD = maxDim_now, 
                                                                          maxit = maxit_unsup_feature_selection[[str_desc]])
                                    }else if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'rcur'){
                                        cat("RCUR selection\n")
                                        r <- rcur(mat_data, k = maxDim_now, rand = TRUE)
                                        mat_data = mat_data[ , r$C.idx]
                                        cat(dim(mat_data))
                                    }else if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'cur'){
                                        cat("CUR selection\n")
                                        r <- rcur(mat_data, k = maxDim_now, rand = FALSE)
                                        mat_data = mat_data[ , r$C.idx]
                                        cat(dim(mat_data))
                                    }else if (gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'entropy'){
                                        cat('entropy filtering \n')
                                        mat_data = imputa(mat_data)
                                        entp_feat = apply(mat_data,2, entropy, method="Laplace")
                                        entp_feat = sort(entp_feat, decreasing = TRUE)
                                        mat_data = mat_data[ , match(names(entp_feat)[1:maxDim_now], colnames(mat_data))]
                                    }
                                      
                                    
                                    cat('mat_data of ', str_desc, ' after ', task,
                                            ' = ', dim(mat_data), '\n')
                                    cat('saving in path ', file_data_fn, '\n')
                                    save(mat_data , file = file_data_fn)
                                    
                                }
                                explain_data_types[[str_desc]] = file_data_fn
                                
                            }
                            rm(maxDim_now)
                            # se anche non faccio unsup feat selection prendo il file dopo il clean data!
                            explain_data_types[[str_desc]] = file_data_fn
                            
                            ##################### DIM REDUCTION
                            
                            #imputazione e salvataggio 
                            file_data_fn = str_replace(file_data_fn, pattern = '.Rda', replacement = '_after_imputation.Rda')
                            if (file.exists(file_data_fn) & (!regenerate_DR_data[[str_desc]])){
                                cat('loading  imputed data  of ', str_desc, ' after feature selection\n')
                                load(file = file_data_fn)
                            }else{
                                cat('eventual imputation of ', str_desc, ' data\n')
                                mat_data = imputa(mat_data)
                                save(mat_data, file = file_data_fn)
                            }
        
                            
                            if (!(gsub(best_red_now, pattern = ID_string_red, replacement = '') == 'NULL')){
                              #use_python(python = python_path)
                                minD_now = minD[[str_desc]]
                                if (!use_id){ 
                                    if (verbose==1) cat('no use of ID for defining the dimension of the lower dimensional space: minD = ', minD_now, '\n')
                                }else{
                                    if (verbose==1) cat('use of ID for defining the dimension of the lower dimensional space: minD = ', minD_now, '\n')
                                }
                                if ((str_replace(best_red_now, pattern = ID_string_red, replacement = '') == 'tsne') 
                                                                              & (ncol(mat_data)>10000)){ 
                                    # # tsne crashes is the supploied matrix has a number of columns > (N-1)
                                    cat('D > N not proper in ', best_red_now,' reduction - using cur-tsne instead\n')
                                    best_red_now = paste0('cur-', best_red_now)
                                }
                              
                              
                                if (str_replace(best_red_now, pattern = ID_string_red, replacement = '') == 'umap'){
                                    # # umap wants a lower dimensional space  < (N-2)
                                    minD_now = min(minD_now, min(dim(mat_data))-2)
                                }
                                cat('minD for ', best_red_now,' = ', minD_now, '\n')
                                
                                red_algos = strsplit(str_replace(best_red_now, pattern = ID_string_red, replacement = ''), '-')[[1]]
                                if (length(red_algos) > 1)  minD_now = c(min(nrow(mat_data)-1, ncol(mat_data)-1), minD_now)
                                
                                
                                old_path = my_path_records[[task]]
                                task = best_red_now
                                my_path_records[[task]] = file.path(old_path, task)
                                if (!file.exists(my_path_records[[task]])){ dir.create(my_path_records[[task]])}
                                
                                file_data_fn = file.path(my_path_records[[task]], paste(str_desc, '.Rda' , sep ='' ))
                                if (file.exists(file_data_fn) & (!regenerate_DR_data[[str_desc]]) ){
                                    cat('loading data reduced with ', best_red_now, '\n')
                                    load(file = file_data_fn)
                                }else{
                                    cat(task,    'until dimension is greater than ', minD[[str_desc]],    '\n')
                                    
                                    for(idx_red_algo in 1:length(red_algos)){
                                        red_algo = red_algos[idx_red_algo]
                                        intermediate_fn= file.path(my_path_records[[task]], paste(str_desc, '_', idx_red_algo,'_', red_algo,'.Rda' , sep ='' ))
                                        if (file.exists(intermediate_fn)){
                                            load(intermediate_fn)
                                        }else{
                                            
                                            red_data  = do.call(dim_red[[red_algo]], list(x= mat_data, k = round(minD_now[idx_red_algo])))
                                            if (is.null(row.names(red_data))){row.names(red_data) = row.names(mat_data)}
                                            mat_data = red_data
                                            rm(red_data)
                                            
                                            if (!grepl('cur', red_algo)){ 
                                                colnames(mat_data) = 
                                                paste(str_desc, red_algo, as.character(1:ncol(mat_data)), sep = '_')
                                            }else{ 
                                                colnames(mat_data) = paste(str_desc, red_algo, colnames(mat_data), sep = '_')
                                            }
                                            if ((length(red_algos)>1) & (idx_red_algo<length(red_algos))){
                                                intermediate_fn= file.path(my_path_records[[task]], 
                                                                           paste(str_desc, '_', idx_red_algo,'_', red_algo,'.Rda' , sep ='' ))
                                                save(mat_data , red_algo, file = intermediate_fn)
                                            }
                                        }
                                    }
                                    rm(minD_now)
                                    cat('mat_data of ', str_desc, ' after ', task,
                                            ' = ', dim(mat_data), '\n')
                                    cat('saving in path ', file_data_fn, '\n')
                                    save(mat_data , file = file_data_fn)
                                   
                                    args_ll = args_list[[ID_estimator_fun]]
                                    args_ll[['mat_data']] = mat_data
                                    ids = do.call(ID_estimator_fun, args_ll)
                                    
                                    save_ID_list[['ID_est']][task] = mean(ids[['id']])
                                    save_ID_list[['sd_ID_est']][task] = mean(ids[['sd_id']])
                                    cat('estimated ids after ', task, 'on ', str_desc, ' = ', save_ID_list[['ID_est']][task], '\n')
                                    
                                    save(save_ID_list, file = ID_fn)
                               }
                            }else{
                                cat('no dimensionality reduction applied\n')
                            }                
        
                            embedded_data_types[[str_desc]] = file_data_fn
                            ############################ END DIM REDUCTION
                            
                            df_names = data.frame(row.names(mat_data), row.names = row.names(mat_data))
                            names(df_names) = str_desc 
                            mat_data_all_names = merge(mat_data_all_names, 
                                                       df_names, by = "row.names", all.x = TRUE, all.y = TRUE)
                            rm(df_names)
                            row.names(mat_data_all_names) = mat_data_all_names[["Row.names"]]
                            mat_data_all_names = mat_data_all_names[, !(names(mat_data_all_names) %in% "Row.names")]
                            
                        }#end for (str_desc in data_types)
        
                        

                        
                        if ((only_embed_dataset) | (only_compute_blocking)) next
                        
                        
                                        
                        cat('\n\n******\n******\n******\n')
                        cat('setting =    - feature selection = ', unsup_feature_selection, 
                                ' - best_red = ', best_red, '\n')
                    
                    #    cat('MEGA_DF dim:', dim(mega_df), ' - data-types    =', names(embedded_data_types), '\n')
                        n <- length(names(embedded_data_types))
                        #design_mat = rbind(diag(n), tril(matrix(1, nrow = n, ncol= n))[2:n,])
                        #colnames(design_mat) = names(embedded_data_types)
                    
                        if (n>1){
                            n <- length(names(embedded_data_types))
                            l <- rep(list(0:1), n)
                            
                            names(l) = names(embedded_data_types)
                            design_mat = as.matrix(expand.grid(l))
                            design_mat = design_mat[rowSums(design_mat)>0, ]
                            colnames(design_mat) = names(embedded_data_types)
                            # tolgo la riga che corrisponde solo al paziente
                            #design_mat = design_mat[-which((rowSums(design_mat)==1) & design_mat[,which(colnames(design_mat)=='pt')]),]
                            
                            # se vuoi avere almeno tre viste decommenta questa!
                            if (minNumViews>1){
                                cat('At least ', minNumViews, ' are requested! \n' )   
                                design_mat = design_mat[-which((rowSums(design_mat[, !(colnames(design_mat) %in% 'pt')])<(minNumViews))),]
                                if (!is.na(with_pt)){
                                    #se ? NA provi entrambe le conbinazioni
                                    if (!with_pt) design_mat = design_mat[which(design_mat[, (colnames(design_mat) %in% 'pt')] == 0 ), ]
                                    if (with_pt) design_mat = design_mat[which(design_mat[, (colnames(design_mat) %in% 'pt')] == 1 ), ]
                                } 
                                if (!is.matrix(design_mat)) design_mat = t(as.matrix(design_mat))
                            }
                            # togli i pz per cui hai solo i dati del paziente!
                         # present_omics = rowSums(!is.na(mat_data_all_names[, !(names(mat_data_all_names) %in% 'pt')]))
                            if (length(data_types)==1) 
                                present_omics = 
                                    as.numeric(!is.na(mat_data_all_names[, names(mat_data_all_names) %in% data_types]))
                            else 
                                present_omics = 
                                    apply(!is.na(mat_data_all_names[, names(mat_data_all_names) %in% data_types]),1,sum)
                        
                            mat_data_all_names = 
                                mat_data_all_names[!(rownames(mat_data_all_names) %in% names(present_omics[present_omics==0])), ]
                            
                            if (partial_views){ 
                                all_pts = row.names(mat_data_all_names)
                            }else{ 
                                missing_omics = rowSums(is.na(mat_data_all_names[, !(names(mat_data_all_names) %in% 'pt')]))
                                all_pts = names(present_omics[missing_omics==0])
                            } 
                            
                        }else{
                            design_mat = matrix(1, nrow = 1, ncol =1)
                            colnames(design_mat)=names(embedded_data_types)
                            all_pts = row.names(df_outcome)
                            if (names(embedded_data_types)=='pt'){ 
                                task = 'pt_open'
                                my_path_records = list('pt_open' = res_dir)
                            }
                        }
                        
                        for (j in 1:nrow(design_mat)){ 
                            use_types = colnames(design_mat)[design_mat[j,]==1]
                            
                            data_str = paste(use_types, collapse = '_')
                            #df_outcome[['rn']] = row.names(df_outcome) check that they match
                            
                            if (data_str == 'pt'){ 
                                unsup_feature_selection = 'none'
                                best_red = 'none'
                                data_integration = 'none'
                            }
                            
                            cat('*********************************************\n')
                            cat('center = ', centers, ' - tasks = ', unsup_feature_selection, ' + ', best_red, 
                                ' - data misc (j)= ',    data_str, '(', j ,') \n')    
                            
                            for (data_integration in data_integration_list){
                                cat('LINEA 931\n')
                                cat( ' - integration = ', data_integration, '\n' )
                                
                                # se non uso feature extraction per dimensionality reduction significa che ho le features
                                # originali che sono gi? opportunamente normalizzate. Ho gi? eseguito (con SNF) - inutile rifarlo
                                if ((str_sub(data_integration, start= -2) == 'PT') & (!any(use_types=='pt'))){ 
                                    #voglio un metodo che integri anche i dati dei pz ma non li uso in questa combinazione di dati
                                    next
                                }
                                if ((gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '') == 'NULL') & 
                                    (gsub(best_red, pattern = ID_string_red, replacement = '') == 'NULL') & 
                                    (data_integration == 'concatenation')){
                                    cat('no feature selection and no dim-red - simple RF will explode! skip them \n')
                                    next
                                }
                                if ((length(use_types[!(use_types %in% 'pt')])==1) & (data_integration== 'mofa')){ 
                                    cat(data_str, ' unique not-pt data type. Mofa might fail\n')
                                    
                                }
                                # if ((length(use_types[!(use_types %in% 'pt')])==1) & (data_integration== 'mofa')){ 
                                #     cat(data_str, ' unique not-pt data type. Do not integrate\n')
                                #     next
                                # }
                                
                                df_use = data.frame(df_outcome[match(all_pts, row.names(df_outcome)), ])
                                row.names(df_use) = all_pts
                                names(df_use) = names(df_outcome)
                                
                                df_use_explain = df_use
                                
                                data_integration_str = paste0(data_integration, '_')   
                                #se hai dovuto rigenerare qualche vista cambierà tutto, quindi devi rigenerare anche le integrazioni
                                regenerate_integrated_views = any(unlist(regenerate_DR_data)) | regenerate_integration_list[[data_integration]]
                                
                                if (str_sub(data_integration, end = 3)  == 'mkl') {
                                    fn_save_mkl_integrated = file.path(my_path_records[[task]], 
                                                                       paste0(data_integration, '_', outcome_col, '_', data_str,'.Rda' ))
                                    if (file.exists(fn_save_mkl_integrated) & (!regenerate_integrated_views) ){ 
                                         load(fn_save_mkl_integrated)
                                    }else{
                                        W = list()
                                        W_int = NULL
                                        flagStop = FALSE
                                        for (ut in use_types){
                                            
                                            if (ut !='pt'){
                                                if (!(is.data.frame(embedded_data_types[[ut]]))) {
                                                    load(embedded_data_types[[ut]])
                                                }else mat_data = embedded_data_types[[ut]]
                                                
                                                mat_data_all = matrix(NA, nrow = length(all_pts), ncol=ncol(mat_data))
                                                row.names(mat_data_all) = all_pts
                                                idx_names = match(row.names(mat_data), all_pts)
                                                mat_data_all[idx_names[!is.na(idx_names)], ] = mat_data[!is.na(idx_names), ]
                                                # now compute kernel with radial basis function and sigma =  median (1 /(x_i - x_j)^2 ) as 
                                                # used in https://academic.oup.com/bioinformatics/article/34/6/1009/4565592#393770555
                                                
                                                ## Calculate the pair-wise distance using the function from SNF;
                                                ## If the data is continuous, we recommend to use the function "dist2" as follows
                                                dist_mat = (dist2(as.matrix(mat_data_all),as.matrix(mat_data_all)))
                                                ## next, build kernel
                                                sigma = median(1/dist_mat)
                                                W_mat = tryCatch(expr = mixKernel::compute.kernel(mat_data_all, kernel.func= "gaussian.radial.basis",
                                                                                  sigma = sigma), error = function(e) NULL)
    
                                                ## next, we fuse all the graphs
                                                ## then the overall matrix can be computed by similarity network fusion(SNF):
                                                if (is.null(W_mat)){ 
                                                    cat("something went wrong - matrix is not positive semi-definite\n")
                                                    flagStop = TRUE
                                                    break
                                                }else{
                                                    W[[ut]] = W_mat
                                                }                                        
                                                print(dim(W[[ut]]$kernel))
                                                
                                            }else{
                                                # if there are pt data, merge it with df_use (now only contains outcome!)
                                                if (str_sub(data_integration, start= -2) != 'PT'){
                                                    df_use = merge(df_use, df_pt_use, by = 'row.names', all.x = TRUE, all.y = TRUE)
                                                    row.names(df_use) = df_use[['Row.names']]
                                                    df_use = df_use[, !(names(df_use) %in% 'Row.names')]
                                                }else if (str_sub(data_integration, start= -2) == 'PT'){
                                                    index_df_pt = match(all_pts, row.names(df_pt_use))
                                                    dist_mat = gower.dist(df_pt_use[index_df_pt, ],df_pt_use[index_df_pt, ])
                                                    rm(index_df_pt)
                                                    ## next, build kernel
                                                    sigma = median(1/dist_mat)
                                                    W_mat = tryCatch(expr = mixKernel::compute.kernel(dist_mat, kernel.func= "gaussian.radial.basis",
                                                                                                       sigma = sigma), error = function(e) NULL)
                                                    # 
                                                    # ## next, we fuse all the graphs
                                                    # ## then the overall matrix can be computed by similarity network fusion(SNF):
                                                    if (is.null(W_mat)){ 
                                                        cat("something went wrong - matrix is not positive semi-definite\n")
                                                        flagStop = TRUE
                                                        
                                                        break
                                                        
                                                    }else{
                                                        W[[ut]] = W_mat
                                                    }                                        
                                                    print(dim(W[[ut]]$kernel))
                                                }
                                            }
                                            
                                        }
                                        # integrate by running MKL
                                        if (!flagStop){ 
                                            W[['method']] = "full-UMKL"
                                            if (length(W)>1) W_int = tryCatch(expr = do.call('combine.kernels', args = W),
                                                                              error = function(e) NULL)
                                            else W_int = W[[1]]
                                        }
                                        
                                        if (flagStop | is.null(W_int)){
                                            W = NULL
                                            W_int = NULL
                                        }
                                        save(W, W_int, df_use, file = fn_save_mkl_integrated)
                                        rm(fn_save_mkl_integrated)
                                    }                                
                                    
                                    if (is.null(W_int)){
                                        cat(data_integration_str, ' FAILED!\n')
                                        W_int = NULL
                                    }else{
                                    
                                        cat(use_types, ' sparse UMKL weights:\n')
                                        cat(W_int$weights, '\n')
                                        weights_to_be_saved = W_int$weights
                                        names(weights_to_be_saved) = names(W)[names(W) != 'method']
                                        
                                        # CONTROLLA MA DOVRESTI AVERE GIA' RISOLTO IL FATTO DEI DATI PARZIALI??
                                        df_use = merge(df_use, W_int$kernel, by = 'row.names', all.x = TRUE, all.y = TRUE)
                                        row.names(df_use) = df_use[['Row.names']]
                                        df_use = df_use[, !(names(df_use) %in% 'Row.names')]
                                        
                                        cat(dim(df_use))
                                        # aggancia i pazienti se ci sono!
                                        rm(W, W_int)
                                    }
                                }else if (grepl('mofa', data_integration)){
    
                                    cat(data_str, ' INTEGRATION WITH MOFA \n')
                                    embed_pt = grepl('PT', data_integration)
                                    df_use = create_mofa_embedding(df_outcome = df_use, df_pt_use = df_pt_use, 
                                                                   embedded_data_types, use_types,
                                                                   embed_pt = embed_pt,
                                                                no_mofa_factors = no_mofa_factors, 
                                                                save_path_task = my_path_records[[task]], 
                                                                seed = seed, res_dir = res_dir, 
                                                                sparse_weights = sparse_weights, 
                                                                sparse_factors = sparse_factors,
                                                                python_path = python_path,
                                                                force_recreation = regenerate_integrated_views)

                                    if (is.null(df_use)){
                                        cat('ALL FACTORS RESULTED with no STD! - ARE NOT SUFFICIENT\n')
                                        insufficient_dim_data = ncol(df_use)
                                        
                                        # rimetti in df_use solo la colonna degli outcome!
                                        df_use = data.frame(df_outcome[match(all_pts, row.names(df_outcome)), ])
                                        row.names(df_use) = all_pts
                                        names(df_use) = names(df_outcome)
                                        
                                        cat(c(data_str, str_centers, unsup_feature_selection, best_red, data_integration_str, insufficient_dim_data), '\n')
                                        
                                    }else{
                                        if (any(use_types %in% 'pt') & (data_integration=='mofa')) min_num_var = 1+ncol(df_pt_use) # 1 per l'outcome!
                                        else min_num_var = 1   # 1 per l'outcome!
                                        if (ncol(df_use) < (min_num_var+2) ){ #min(3,(min_num_var)*2)){
                                            cat('SKIPPED INTEGRATION! NUMBER OF FACTORS ARE NOT SUFFICIENT\n')
                                            insufficient_dim_data = ncol(df_use)
                                            
                                            # rimetti in df_use solo la colonna degli outcome!
                                            df_use = data.frame(df_outcome[match(all_pts, row.names(df_outcome)), ])
                                            row.names(df_use) = all_pts
                                            names(df_use) = names(df_outcome)
                                            
                                            cat(c(data_str, str_centers, unsup_feature_selection, best_red, data_integration_str, insufficient_dim_data), '\n')
                                        }
                                    }
                                }else if (data_integration == 'concatenation'){
                                    cat('INTEGRATION WITH CONCATENATION \n')
                                    data_integration_str = 'concat_'
                                    fn_save_concatenated = file.path(my_path_records[[task]], 
                                                                       paste0(data_integration_str, '_', outcome_col, '_', data_str,'.Rda' ))
                                    if (file.exists(fn_save_concatenated) & (!regenerate_integrated_views)){
                                        cat('reloading concatenated data \n')
                                        load(fn_save_concatenated)
                                    }else{
                                        for (ut in use_types){
                                            if (!(is.data.frame(embedded_data_types[[ut]]))) {
                                                if (!is.null(embedded_data_types[[ut]])){ 
                                                    load(embedded_data_types[[ut]])
                                                    if (!is.data.frame(mat_data)) df_mat = data.frame(mat_data)
                                                }else{ df_mat = NULL }
                                                
                                                # carica i dati per explanation
                                                if (!is.null(explain_data_types[[ut]])){ 
                                                    load(explain_data_types[[ut]])
                                                    if (!is.data.frame(mat_data)) df_mat_explain = data.frame(mat_data)
                                                }else{ df_mat_explain = NULL}
                                            }else{ 
                                                df_mat = embedded_data_types[[ut]]
                                                df_mat_explain = explain_data_types[[ut]]
                                            }
                                            # se voglio viste complete comunque ho gia' risolto togliendo a priori i pz a cui mancava una vista MA
                                            # nella nuova vista da integrare ci potrebbero essere piu' pazienti 
                                            # perch? ho filtrato df_outcome in modo da togliere i pazienti che hanno solo la clinica! 
                                            # se metto all.x = TRUE creo un NA nell'outcome variable
                                            if(!is.null(df_mat)){
                                                df_use = merge(x = df_mat, y = df_use, by = 'row.names', all = partial_views)
                                                row.names(df_use)    =    df_use$Row.names
                                                df_use = df_use[, !(names(df_use) %in% "Row.names")]
                                            }
                                            
                                            if(!is.null(df_mat_explain )){
                                                df_use_explain = merge(x = df_mat_explain, y = df_use_explain, by = 'row.names', all = partial_views)
                                                row.names(df_use_explain) = df_use_explain$Row.names
                                                df_use_explain = df_use_explain[, !(names(df_use_explain) %in% "Row.names")]
                                            }
                                        }
                                        save(df_use_explain, df_use, file = fn_save_concatenated)
                                    }
                                }else if (str_sub(data_integration, end = 3) =='SNF'){
                                    
                                    fn_save_SNF_integrated = file.path(my_path_records[[task]], 
                                                                       paste0(data_integration, '_', outcome_col, '_', data_str,'.Rda' ))
                                    if (file.exists(fn_save_SNF_integrated) & (!regenerate_integrated_views)){ 
                                        load(fn_save_SNF_integrated)
                                    }else{
                                        W = list()
                                        for (ut in use_types){
                                            
                                            if (ut !='pt'){
                                                if (!(is.data.frame(embedded_data_types[[ut]]))) {
                                                    load(embedded_data_types[[ut]])
                                                }else mat_data = embedded_data_types[[ut]]
                                                
                                                mat_data_all = matrix(NA, nrow = length(all_pts), ncol=ncol(mat_data))
                                                row.names(mat_data_all) = all_pts
                                                idx_names = match(row.names(mat_data), all_pts)
                                                mat_data_all[idx_names[!is.na(idx_names)], ] = mat_data[!is.na(idx_names), ]
                                                # now compute affinity matrix
                                                
                                                ## Calculate distance matrices
                                                ## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
                                                ## If the data are all continuous values, we recommend the users to perform
                                                ## standard normalization before using SNF,
                                                ## though it is optional depending on the data the users want to use.
                                                if (data_integration == 'SNF'){
                                                    # as suggested by authors apply data normalization
                                                    data_norm = standardNormalization(mat_data_all)
                                                }else{
                                                    data_norm = mat_data_all
                                                }
                                                ## Calculate the pair-wise distance;
                                                ## If the data is continuous, we recommend to use the function "dist2" as follows
                                                dist_mat = (dist2(as.matrix(data_norm),as.matrix(data_norm)))^(1/2)
                                                ## next, construct similarity graphs
                                                W_mat = affinityMatrix(dist_mat, K = k_snf, sigma = sigma_snf)
                                                ## next, we fuse all the graphs
                                                ## then the overall matrix can be computed by similarity network fusion(SNF):
                                                W[[ut]] = W_mat
                                                
                                                print(dim(W[[ut]]))
                                            }else{
                                                
                                                if (str_sub(data_integration, start= -2) != 'PT'){
                                                # if there are pt data, merge it with df_use (now only contains outcome!)
                                                df_use = merge(df_use, df_pt_use, by = 'row.names', all.x = TRUE, all.y = TRUE)
                                                row.names(df_use) = df_use[['Row.names']]
                                                df_use = df_use[, !(names(df_use) %in% 'Row.names')]
                                                }else  if (str_sub(data_integration, start= -2) == 'PT'){
                                                    index_df_pt = match(all_pts, row.names(df_pt_use))
                                                    dist_mat = gower.dist(df_pt_use[index_df_pt, ],df_pt_use[index_df_pt, ])
                                                    rm(index_df_pt)
                                                    row.names(dist_mat) = all_pts
                                                    colnames(dist_mat) = all_pts
                                                    ## next, construct similarity graphs
                                                    W_mat = affinityMatrix(dist_mat, K = k_snf, sigma = sigma_snf)
                                                    ## next, we fuse all the graphs
                                                    ## then the overall matrix can be computed by similarity network fusion(SNF):
                                                    W[[ut]] = W_mat                    
                                                    print(dim(W[[ut]]))
                                                }
                                            }
        
                                        }
                                    # integrate by running SNF
                                        if (length(W)) W_int <- SNF(W, K=k_snf, t=t_snf)
                                        else W_int = W[[1]]
                                        save(W, W_int, df_use, file = fn_save_SNF_integrated)
                                        rm(fn_save_SNF_integrated)
                                    }   
                                    
                                    if (!is.null(W_int)){
                                        # partial data problem is solved already
                                        df_use = merge(df_use, W_int, by = 'row.names', all.x = TRUE, all.y = TRUE)
                                        row.names(df_use) = df_use[['Row.names']]
                                        df_use = df_use[, !(names(df_use) %in% 'Row.names')]
                                        cat(dim(df_use))
                                        rm(W, W_int)
                                    }else{
                                        cat(data_integration_str, ' failed\n')
                                    }
                                }

                                
                                #### WHEN I'm here if df_use contains only the outcome column, this means something went wrong with the data integration!
                                
                                # dal momento che ho tenuto tutti i pazienti (perch? nelle merge avevo all.x = TRUE e all.y = TRUE) ora 
                                # butto quelli che hanno outcome = NA
                                if (ncol(df_use)>1){
                                    if (any(is.na(df_use[[outcome_col]]))){
                                        idx = which(!is.na(df_use[[outcome_col]]))
                                        df_use = df_use[-idx, ]
                                    }
                                    
                                    cat('-----*****------****** \n')
                                    cat('-----*****------****** \n')
                                    cat('-----*****------****** \n')
                                    
                                    cat("outcome used = ", outcome_col, '\n')
                                    print(table(df_use[[outcome_col]]))
                                    
                                    
                                    cat(' data misc (j)= ',    data_str, '(', j , ') - center = ', centers, 
                                        ' -  tasks = ', unsup_feature_selection, ' + ', best_red , ' - integration = ', 
                                        data_integration_str, '\n')
                                    if (!just_collect_results){
                                      
                                      if (any(is.na(df_use))){
                                        fn_imputed = 
                                          file.path(my_path_records[[task]], paste(outcome_col, '_', data_str, '_', 
                                                                                   data_integration_str, '_',
                                                                                   ifelse(gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = '')=='NULL','', unsup_feature_selection), 
                                                                                   '_',
                                                                                   ifelse(gsub(pattern = ID_string_red, replacement = '', best_red)=='NULL','', best_red),
                                                                                   '.Rda', sep = ''))
                                        if ((!file.exists(fn_imputed)) | (regenerate_integrated_views)){
                                          cat('Found some missing views... Imputing them with missRanger\n')
                                          df_use = tryCatch(imputa(df_use, fun = 'missRanger', k = 1, 
                                                          my_f = as.formula(paste('. ~ . - ', outcome_col, sep = ''))), 
                                                    error = function(e) NULL)
                                          
                                          save(df_use, file = fn_imputed)
                                        }else{
                                          cat('Found some missing views... loading data imputed with missRanger\n')
                                          load(fn_imputed)
                                        }                                
                                      }
                                      if (!is.null(df_use)){                                      
                                          if (!is.factor(df_use[[outcome_col]]))
                                            df_use[[outcome_col]] = factor(df_use[[outcome_col]], levels = c(neg_label, pos_label))
                                          if (any(grepl('gender', names(df_use)))) 
                                            if (!is.factor(df_use$gender)) df_use$gender = factor(df_use$gender, levels = c('male', 'female'))
                                          if (any(grepl('race', names(df_use)))) 
                                            if (!is.factor(df_use$race)) df_use$race = factor(df_use$race, levels = c('white', 'black or african american', 'asian'))
                                          if (any(grepl('ethnicity', names(df_use)))) 
                                            if (!is.factor(df_use$ethnicity)) df_use$ethnicity = factor(df_use$ethnicity, levels = c('not hispanic or latino','hispanic or latino'))
                                      }else{
                                          cat('imputation did not suceed due to memory limits or convergence problems\n')
                                          df_use = data.frame(df_outcome[match(all_pts, row.names(df_outcome)), ])
                                          row.names(df_use) = all_pts
                                          names(df_use) = names(df_outcome)
                                      }
                                    }
                                    
                                    
                                }else{
                                    cat('integration method ', data_integration_str, ' did not converge!!\n')
                                    cat('FILLING OUTPUT results for data misc (j)= ',    data_str, '(', j , ') - center = ', centers, 
                                            ' - tasks = ', unsup_feature_selection, ' + ', best_red, ' - integration = ', 
                                            data_integration_str, '\n')
                                    
                                }                                
                                
                                
                                if (only_integrate_views) next
                                
                                for (RF_method in RF_method_list){
                                    str_RF_method = RF_method
                                    for (feature_selection_method in feature_selection_method_list){
                                        cat('Feature selection method = ', feature_selection_method, '\n',
                                            'RF_method = ', RF_method, '\n')
                                        exp_desc_str = paste0('BEST_models_90_10_' , ifelse(outcome_col == 'os', paste0(data_integration_str, 'os_'), data_integration_str ))
                                        
                                        # if I have only one column, it is because the integration failed!
                                        if (ncol(df_use)>1){
                                            
                                            if (recompute_RF) warning('********** RECOMPUTE RF SET!!! ******************\n')
                                            
                                            
                                             timeout = n_external_iters*timeout_per_int_iter*n_internal_iters
                                           
                                            # if you needed to recompute the integrated views, then you also need to recompute the RFs!
                                            results_RF_list = withTimeout({
                                               analisi_RF_myOpti_ranger_breiman(recompute_RF = recompute_RF | 
                                                                                    any(unlist(regenerate_DR_data)) | 
                                                                                    regenerate_integrated_views,
                                                                                                 df = df_use, 
                                                                                                 exp_desc_str = exp_desc_str,
                                                                                                 RF_method = RF_method, 
                                                                                                 verbose = 0,
                                                                                                 idvar = NULL, 
                                                                                                 outcome = outcome_col, 
                                                                                                 neg_label = neg_label,
                                                                                                 pos_label = pos_label,
                                                                                                 save_path = my_path_records[[task]], 
                                                                                                 str_desc = data_str,
                                                                                                 n_internal_iters = n_internal_iters, 
                                                                                                 n_external_iters = n_external_iters, 
                                                                                                 feature_selection_method = feature_selection_method,
                                                                                                 seed = seed,
                                                                                                 perc_samples = perc_samples, 
                                                                                                 ratio_train = ifelse(str_centers == 'PRAD1', 0.8,ratio_train),
                                                                                                 measures = measures, 
                                                                                                 tuneRF = tuneRF,
                                                                                                 class_weights = class_weights,
                                                                                                 balanced_bootstraps = balanced_bootstraps,
                                                                                                 opti_meas = opti_meas, 
                                                                                                 weights_opti_meas = weights_opti_meas, 
                                                                                                 ncores = min(ncores, detectCores()-1), 
                                                                                                 ntree = ntree, 
                                                                                                 mtry = mtry,
                                                                                                 timeout_per_int_iter = timeout_per_int_iter,
                                                                                                 timeout_per_ext_iter = timeout_per_ext_iter,
                                                                                                 just_collect_results = just_collect_results)  
                                            }, timeout = timeout, onTimeout = "warning")
                                            
                                        }else{
                                          save_path = my_path_records[[task]]
                                          
                                          opti_meas_str = ifelse(length(opti_meas)==1, opti_meas, paste(opti_meas, collapse ='_'))
                                          base_save_str = paste0(exp_desc_str, paste(data_str, class_weights, balanced_bootstraps, sep = '_'))
                                          if (!(RF_method == 'RF'))  base_save_str = paste0(base_save_str, '_', RF_method)    
                                          if (!(feature_selection_method == 'RF_importance')) base_save_str =  paste0(base_save_str, '_', feature_selection_method)
                                          if (!(opti_meas_str == 'aucpr')) base_save_str =  paste0(base_save_str, '_', opti_meas_str)
                                          fn_res_list = file.path(save_path, paste0('listRF_', base_save_str, '.Rda', sep=''))
                                          cat("something went wrong in the integration - saving NULL file in ", fn_res_list, '\n')
                                          
                                          results_RF_list = NULL
                                          res_list = NULL
                                          save(res_list, file = fn_res_list)
                                        }
                                        
                                        prop_pos = sum(df_use[[outcome_col]]==1)/nrow(df_use)
                                        prop_neg = sum(df_use[[outcome_col]]==0)/nrow(df_use)
                                        
                                        if (is.null(results_RF_list)) { 

                                            save_path = my_path_records[[task]]
  
                                            opti_meas_str = ifelse(length(opti_meas)==1, opti_meas, paste(opti_meas, collapse ='_'))
                                            base_save_str = paste0(exp_desc_str, paste(data_str, class_weights, balanced_bootstraps, sep = '_'))
                                            if (!(RF_method == 'RF'))  base_save_str = paste0(base_save_str, '_', RF_method)
                                            if (!(feature_selection_method == 'RF_importance')) base_save_str =  paste0(base_save_str, '_', feature_selection_method)
                                            if (!(opti_meas_str == 'aucpr')) base_save_str =  paste0(base_save_str, '_', opti_meas_str)
                                            fn_res_list = file.path(save_path, paste0('listRF_', base_save_str, '.Rda', sep=''))
                                            cat("saving NULL file in ", fn_res_list, '\n')
                                            if (!file.exists(fn_res_list)){  
                                                cat('RF took too long to complete\n')
                                                uncompleted[paste0(str_centers, '_', exp_id) , data_integration] = uncompleted[paste0(str_centers, '_', exp_id), data_integration]+1
                                            }else{
                                                cat('RF took too long to complete\n')
                                                not_converged[paste0(str_centers, '_', exp_id), data_integration] = not_converged[paste0(str_centers, '_', exp_id), data_integration]+1
                                              
                                            }
                                          #                                          res_list = NULL
#                                          save(res_list, file = fn_res_list)
                                          
                                            results_RF = data.frame(data = data_str, center = str_centers, 
                                                                    unsup_feature_selection = gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = ''),
                                                                    unsup_feature_extraction = gsub(pattern = ID_string_red, replacement = '', best_red),
                                                             outcome = outcome_col, 
                                                             use_id = use_id, frac_samples = frac_samples,
                                                             RF_method = str_RF_method, data_integration = data_integration_str, 
                                                             opti_meas = opti_meas, weight_meas = weights_opti_meas, 
                                                             nrow_data = nrow(df_use), ncol_data = ncol(df_use),
                                                             prop_vars_selected = 0, 
                                                             perc_samples = perc_samples, ratio_train = ratio_train, 
                                                             feature_selection_method = feature_selection_method, 
                                                             n_external_iters = n_external_iters, n_internal_iters = n_internal_iters,
                                                             maxvoted_no_tree = 0, maxvoted_mtry = 0,
                                                             holdout_no_avg = 0,
                                                             aucpr_avg = min(prop_pos, prop_neg), auc_avg = 0.5, 
                                                             delta_aucpr = prop_pos, pos_rate = prop_pos)
                                            
                                            mega_results_RF = cbind(data = data_str, center = str_centers, 
                                                                    unsup_feature_selection = gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = ''),
                                                                    unsup_feature_extraction = gsub(pattern = ID_string_red, replacement = '', best_red),
                                                                    outcome = outcome_col,  
                                                                    use_id = use_id, frac_samples = frac_samples,
                                                                    RF_method = str_RF_method, data_integration = data_integration_str, 
                                                                    opti_meas = opti_meas, weight_meas = weights_opti_meas, 
                                                                    nrow_data = nrow(df_use), ncol_data = ncol(df_use),
                                                                    prop_vars_selected = 0, 
                                                                    perc_samples = perc_samples, ratio_train = ratio_train, 
                                                                    feature_selection_method = feature_selection_method, 
                                                                    balanced_bootstraps = balanced_bootstraps, 
                                                                    class_weights = class_weights,
                                                                    n_external_iters = n_external_iters, n_internal_iters = n_internal_iters,
                                                                    ntree = 0, mtry = 0,
                                                                    data.frame(holdout_no = 1:n_external_iters), 
                                                                    aucpr = min(prop_pos, prop_neg), 
                                                                    auc = 0.5, delta_aucpr = prop_pos, pos_rate = prop_pos)
                                            
                                            
                                        }else{ #(!is.null(results_RF_list)) { 
                                            # se il valore ritornato ? null ? perch? il numero di volte in cui le var erano non significative
                                            # ha raggiunto la soglia
                                            
                                            results_RF = results_RF_list[['mean_over_holdouts']]
                                            mega_results_RF = results_RF_list[['all_holdouts']]
                                            # DA USARE PER EXPLAINABILITY
                                            all_predictions = results_RF_list[['all_predictions_for_explanation']] 
                                            
                                            # salva mean results
                                            if (data_integration == 'mkl'){
                                                results_RF = cbind(data = data_str, center = str_centers, 
                                                                   unsup_feature_selection = gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = ''),
                                                                   unsup_feature_extraction = gsub(pattern = ID_string_red, replacement = '', best_red),
                                                               outcome = outcome_col,  
                                                               use_id = use_id, frac_samples = frac_samples,
                                                               RF_method = str_RF_method, data_integration = data_integration_str,
                                                               results_RF,  data.frame(t(weights_to_be_saved)), 
                                                               delta_aucpr = results_RF$aucpr_avg - prop_pos, pos_rate = prop_pos)
                                            }else{
                                                results_RF = cbind(data = data_str, center = str_centers, 
                                                                   unsup_feature_selection = gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = ''),
                                                                   unsup_feature_extraction = gsub(pattern = ID_string_red, replacement = '', best_red),
                                                                   outcome = outcome_col,  
                                                                   use_id = use_id, frac_samples = frac_samples,
                                                                   RF_method = str_RF_method, data_integration = data_integration_str,
                                                                   results_RF, 
                                                                   delta_aucpr = results_RF$aucpr_avg - prop_pos, pos_rate = prop_pos)
                                            }
                                            
                                            # salva mega results                                    
                                            mega_results_RF = cbind(data = data_str, center = str_centers, 
                                                                    unsup_feature_selection = gsub(unsup_feature_selection, pattern = ID_string_unsupFS, replacement = ''),
                                                                    unsup_feature_extraction = gsub(pattern = ID_string_red, replacement = '', best_red),
                                                                    outcome = outcome_col,
                                                                    use_id = use_id, frac_samples = frac_samples,
                                                                    RF_method = str_RF_method, data_integration = data_integration_str, 
                                                                    mega_results_RF,  
                                                                    delta_aucpr = mega_results_RF$aucpr - prop_pos, pos_rate = prop_pos)
                                            
                                        }
                                        
                                        selected_mean_res_df = plyr::rbind.fill(selected_mean_res_df, results_RF)
                                        selected_mega_res_df = plyr::rbind.fill(selected_mega_res_df, mega_results_RF)
                                        cat('finishing and saving new RF results in the selected results\n')  
                                        cat('path mean:',  selected_file_mean_results_RF, '\n')
                                        writexl::write_xlsx(selected_mean_res_df, path = selected_file_mean_results_RF)
                                        cat('path mega:',  selected_file_mega_results_RF, '\n')
                                        writexl::write_xlsx(selected_mega_res_df, path = selected_file_mega_results_RF)
                                            
                                            # df_use_explain_with_pred = df_use_explain
                                            # df_use_explain_with_pred[['patient_id']] = row.names(df_use_explain_with_pred )
                                            # 
                                            # df_use_explain_with_pred = merge(all_predictions, df_use_explain_with_pred,  by = 'patient_id', all = FALSE)
                                            # columns_to_test =  names(df_use_explain_with_pred[, !(names(df_use_explain_with_pred) %in% 
                                            #                                                           c('patient_id', 'pos_pred', 'labels'))] )
                                            # pvals = rep(0, length(columns_to_test))
                                            # names(pvals) = columns_to_test
                                            # for (cc in columns_to_test){
                                            #     res.pvalue = kruskal.test(df_use_explain_with_pred[[cc]],df_use_explain_with_pred$pos_pred)
                                            #     pvals[cc] = c(res.pvalue$p.value)
                                            #     
                                            # }
                                            # rm(cc)
                                            # 
                                            # pvals_df = data.frame(pvals)
                                            # pvals_df[['var_names']] = names(pvals)
                                            # pvals_df[['select']] = pvals_df$pvals < thr_conf
                                            # save(pvals_df,  file = file.path(my_path_records[[task]], 
                                            #                                  paste0('explanations_', data_str, '_', 
                                            #                                         n_external_iters, '_', 
                                            #                                         n_internal_iters, '_', 
                                            #                                         RF_method, '_', 
                                            #                                         feature_selection_method, 
                                            #                                         '.Rda' ) ))
                                            # 
                                            # writexl::write_xlsx(x = pvals_df, path = file.path(my_path_records[[task]], 
                                            #                                                    paste0('explanations_', data_str, '_', 
                                            #                                                           n_external_iters, '_', 
                                            #                                                           n_internal_iters, '_', 
                                            #                                                           RF_method, '_', 
                                            #                                                           feature_selection_method, 
                                            #                                                           '.xlsx' ) ))
                                            # 
                                        #                                       
                                        # }
                                        
                                    }#for (feature_selection_method in feature_selection_method_list)        
                                }# for (RF_method in RF_method_list)
                            }# for (data_integration in data_integration_list)
                        }#for (j in 1:nrow(design_mat))
                    } # best_red
                }# unsup_feature_selection
                rm(use_id, frac_samples, minNumViews)
            } #for (exp_id in use_id_list){
                
            all_centers_res = rbind.fill(all_centers_res, selected_mean_res_df)
            all_centers_mega_res = rbind.fill(all_centers_mega_res, selected_mega_res_df)
            
            # writexl::write_xlsx(all_centers_res, path = paste0(paste(centers, collapse ='_'), '_mean_res.xlsx'))
            # writexl::write_xlsx(all_centers_mega_res, path = paste0(paste(centers, collapse ='_'), '_mega_res.xlsx'))
            
        }#centers in centers_list
        cat('FINITO!!\n') 
        
        
        # writexl::write_xlsx(selected_mean_res_df, path = selected_file_mean_results_RF)
        # writexl::write_xlsx(selected_mega_res_df, path = selected_file_mega_results_RF)
        # 
        if(!str_add_save==''){
            writexl::write_xlsx(all_centers_res, path =  paste0(str_add_save, '_all_mean_results.xlsx'))
            writexl::write_xlsx(all_centers_mega_res, path = paste0(str_add_save, '_all_mega_results.xlsx'))
            write.csv(uncompleted, file =  paste0(str_add_save, '_uncompleted.csv'))
            write.csv(not_converged, file =  paste0(str_add_save, '_not_converged.csv'))
        }else{
            writexl::write_xlsx(all_centers_res, path =  'all_mean_results.xlsx')
            writexl::write_xlsx(all_centers_mega_res, path = 'all_mega_results.xlsx')
            write.csv(uncompleted, file =  'uncompleted.csv')
            write.csv(not_converged, file =  'not_converged.csv')
        }
        
#        writexl::write_xlsx(all_centers_mega_res, path = 'all_mega_results.xlsx')


       
    
        return(list('mean_res' = all_centers_res, 'mega_res' = all_centers_mega_res))
}

    