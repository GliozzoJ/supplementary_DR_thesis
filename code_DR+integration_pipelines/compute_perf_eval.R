# compute_perf_eval <- function(prediction_probs = NULL, labels = NULL, cutoff = 0.5, # probability for positive class
#                               pos_label = '1', neg_label = '0', measures = c("aucpr", "auc") ){
#     #prediction_probs is a matrix where each row is a case and each column contains probabilities for a class
#     # labels are true labels - they should be factors 
#     
#     class_probs = prediction_probs[, colnames(prediction_probs) %in% pos_label]
#     
#     no_aucs = sum(measures %in% c('auc', 'aucpr'))
#     no_nonaucs = length(measures)-no_aucs # if some measures are not aucs I must find the threshold for the class prob and derive a classification
#     
#     
#     if (no_nonaucs>0){
#         classifications = ifelse(class_probs>cutoff, 1, 0 )
#         #pred_class = ROCR::prediction(classifications, labels,  label.ordering = c(neg_label,pos_label))
#     }
#     
#     perfs = NULL
#     numeric_labels_1 = as.numeric(labels==pos_label)
#     numeric_labels_0 = as.numeric(labels==neg_label)
#     
#     for (measure in measures){
#         #pred_ROCR = ROCR::prediction(class_probs, labels,  label.ordering = c(neg_label,pos_label))
#         pre_list=c('1'=list(precision.at.all.recall.levels(prediction_probs[, colnames(prediction_probs) %in% pos_label], numeric_labels_1)), 
#                    '0' =list(precision.at.all.recall.levels(prediction_probs[, colnames(prediction_probs) %in% neg_label], numeric_labels_0)))
#         
#         
#         if (measure == 'auc'){ 
#             #print(performance(pred, 'auc'))
#             cat(measure, '\n')
#             auc_perfMeas = PerfMeas::AUC.single(class_probs, numeric_labels_1)
#             perfs = c(perfs, auc_perfMeas)
#         }else if  (measure == 'aucpr'){ 
#             pre_list=c('1'=list(precision.at.all.recall.levels(prediction_probs[, colnames(prediction_probs) %in% pos_label], numeric_labels_1)), 
#                        '0' =list(precision.at.all.recall.levels(prediction_probs[, colnames(prediction_probs) %in% neg_label], numeric_labels_0)))
#             aucpr_perfMeas = PerfMeas::AUPRC(pre_list, comp.precision=TRUE)
#             perfs = c(perfs, aucpr_perfMeas[1])
#         }else{
#             # in this case you must use the thresholded classification and numeric labels
#             metrics = PerfMeas::F.measure.single(classifications, numeric_labels_1)
#             meas = ifelse(measure== 'prec', metrics[1], 
#                     ifelse((measure == 'rec') | (measure == 'sens'), metrics[2], 
#                         ifelse(measure == 'spec', metrics[3], 
#                             ifelse(measure == 'f', metrics[4],
#                                 ifelse(measure == 'acc', metrics[5],NA)))))
#             perfs = c(perfs, meas)
#             
#         }
#     }
#     names(perfs) = measures
#     res_df = data.frame(t(perfs))
#     
#     return(res_df)
# }
# 


compute_perf_eval <- function(prediction_probs = NULL, labels = NULL, cutoff = 0.5,
                              pos_label = '1', neg_label = '0', measures = c('aucpr', 'auc', 'acc', 'sens', 'spec', 'f', 'npv', 'ppv', 'mcc')){
  #prediction_probs is a matrix where each row is a case and each column contains probabilities for a class
  # labels are true labels - they should be factors

  class_probs = prediction_probs[, colnames(prediction_probs) %in% pos_label]

  no_aucs = sum(measures %in% c('auc', 'aucpr'))
  no_nonaucs = length(measures)-no_aucs # if some measures are not aucs I must find the threshold for the class prob and derive a classification


  
  perfs = NULL
  # le ordino per comodità mia
  idx = order(class_probs)
  cprobs = class_probs[idx]
  labs = as.numeric(labels[idx] == pos_label)
  pred = ROCR::prediction(cprobs, labs)
  cutpoints <- cutpointr(data.frame(x = cprobs, class = labs), pos_class = 1, neg_class = 0, 
                         x = 'x', class = 'class', direction = '>=', maximize_metric = 'youden' )
  opt_cutpoint = cutpoints$optimal_cutpoint
  # ho usato direzione >= quindi tutti i valori sotto soglia saranno la classe negativa
  classifications = ifelse(cprobs<cutpoints$optimal_cutpoint, 0, 1 )
  
  P = sum(labs)
  N = sum(labs==0)
  tp = sum((classifications==1) & (labs==1))
  tn = sum((classifications==0) & (labs==0))
  fp = sum((classifications==1) & (labs==0))
  fn = sum((classifications==0) & (labs==1))
  
  meas = c(sens = tp/P,
            spec = tn/N,
           prec = tp/(tp+fp),
           acc = (tp+tn)/(P+N),
           f = (2*tp)/(2*tp+fn+fp),
           ppv = tp/(tp+fp),
           npv = tn/(tn+fn),
           mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fn)*(tn+fp)))
  
  
  for (measure in measures){
    if ((measure == 'auc') | (measure == 'aucpr')){
      perfs = c(perfs, performance(pred, measure = measure)@y.values[[1]])
    }else{
        # meas = NA
        # if (measure== 'prec') meas = prec
        # if ((measure == 'rec') | (measure == 'sens')) meas = sens
        # if (measure == 'spec') meas = spec
        # if (measure == 'f') meas = f
        # if (measure == 'acc') meas = acc
        perfs = c(perfs, meas[measure])
    }
  }
  names(perfs) = measures
  res_df = data.frame(t(perfs))

  return(res_df)


}