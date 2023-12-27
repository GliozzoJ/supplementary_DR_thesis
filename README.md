### Supplementary files - Chapter 3 "Feature selection guided by intrinsic dimensionality on omics data"

This repository contains supplementary **tables** (.xlsx excel files), **figures** (.png) and 
**code** (.R) relative to the third chapter "Feature selection guided by intrinsic dimensionality on omics data" of the PhD thesis "Patient Similarity Networks-based methods for multimodal data integration and clinical outcome prediction".
Note that files directly cited in the thesis as "Supplementary" are named with the prefix
**S#** (S followed by the file number, i.e. S11-S15) and we report also additional result tables for
completeness. All the figures in this repository are explained in the appendices of the thesis.

Supplementary files are organized in the following folders:

- **1_NO_DR_at_all** contains the results of experiments that does not exploit the Dimensionality Reduction (DR) framework proposed in our work. We report results with (`NODR_tuning_mean_results.xlsx`, `noDR_tuning.png`) and without (`NODR_noTUNING_mean_results.xlsx`, `noDR_NO_tuning.png`) Random Forests optimization (supervised feature selection and hyper-parameter tuning). File `NODR_tuning_mean_resultsWITH_PT.xlsx` contains the results of tuned Random Forests 
with the integration of demographic descriptors.

- **2_Comp_HD1_HD2** reports the results comparing DR pipelines using the heuristics 
HD<sub>1</sub> and HD<sub>2</sub> to set the dimension of the reduced space:
    * `Comp_HD1_HD2_AUC.png`, `Comp_HD1_HD2_AUCPR.png` show extracts of the win-tie-loss
    tables performed to make pairwise comparisons between the evaluated pipelines for 
    AUC and AUCPR metrics, respectively. Tables `S1_comp_HD1_HD2_wtl_AUC.xlsx` and 
    `S2_comp_HD1_HD2_wtl_AUCPR.xlsx` reports the corresponding complete win-tie-loss tables.
    * Table `S3_comp_HD1_HD2_best_models.xlsx` collects the three DR+data-fusion pipelines 
    that obtained the best AUC and/or the best AUCPR.

- **3_Comp_HD2_useID** contains the comparisons between DR pipelines that use of
heuristic HD<sub>2</sub> and block-analysis to the set dimension of the reduced
space of each data view:
    * Tables `S5_comp_HD2_useID_wtl_AUC.xlsx` and `S6_comp_HD2_useID_wtl_AUCPR.xlsx`
    report the complete win-tie-loss results for AUC and AUCPR, respectively. 
    * Figures `Comp_HD2_ID1_AUC.png` and `Comp_HD2_ID1_AUCPR.png` depict extract
    of the previous two tables.
    * `S4_comp_HD2_useID_best_models.xlsx` collects the three DR+data-fusion pipelines 
    that obtained the best AUC and/or the best AUCPR.

- **4_Comp_ID1_4omics_noPT** contains the results of the proposed DR pipeline guided by
id-estimation for the integration of 4 -omics, without demographic data: 
    * Tables `S7_comp_useID_4omic_noPT_wtl_AUC.xlsx` and `S8_comp_useID_4omic_noPT_wtl_AUCPR.xlsx` 
    report the complete win-tie-loss tables among DR+data-fusion pipelines for
    AUC and AUCPR metrics, while 
    figures `Comp_SOLOID_4omics_noPT_AUC_detailed.png` 
    and `Comp_SOLOID_4omics_noPT_AUCPR_detailed.png` show the corresponding 
    extracts of the best 25 pipelines. 
    * Table `S9_comp_useID_4omics_noPT_best_models.xlsx` collects the three 
    DR+data-fusion pipelines that obtained the best AUC and/or the best AUCPR.

- **5_Comp_ID1_gt2omics_noPT** contains the results of the proposed DR pipeline guided by
id-estimation for the integration of at least 2 -omics, without demographic data:
     * Tables `S10_comp_gt2omics_NOPT_wtl_AUC.xlsx` and `S11_comp_gt2omics_NOPT_wtl_AUCPR.xlsx` 
    report the complete win-tie-loss tables among DR+data-fusion pipelines for
    AUC and AUCPR metrics, while figures `Comp_ID1_gt2omics_noPT_auc_detailed.png` 
    and `Comp_ID1_gt2omics_noPT_aucpr_detailed.png` show the corresponding 
    extracts of the best 25 pipelines. 
    * Table `S12_comp_gt2omics_NOPT_best_models.xlsx` collects the three 
    DR+data-fusion pipelines that obtained the best AUC and/or the best AUCPR.

- **6_Comp_ID1_4omics_withPT_selectedFS-FE** contains the results of the proposed DR 
pipeline guided by id-estimation for the integration of 4 -omics and demographic data:
     * Tables `S13_comp_4omics_withPT_AUC.xlsx` and `S14_comp_4omics_withPT_AUCPR.xlsx` 
    report the complete win-tie-loss tables among DR+data-fusion pipelines for
    AUC and AUCPR metrics, while figures `6_Comp_ID1_4omics_withPT_auc_detailed.png` 
    and `6_Comp_ID1_4omics_withPT_aucpr_detailed.png` show the corresponding 
    extracts of the best 25 pipelines. 
    * Table `S15_comp_4omics_withPT_best_models.xlsx` collects the three 
    DR+data-fusion pipelines that obtained the best AUC and/or the best AUCPR.

- **7_optimized_RF** contains the results of optimized RF by supervised feature
selection through RF importance and hyper-parameter tuning through internal
stratified holdout validation:
    * Table `Best_results_pair_auc_aucpr_new.xlsx` 
    reports the results for the best models, while tables 
    `BestModels__useID_ID_1_all_mean_results.xlsx` and 
    `BestModels__useID_ID_1_all_mega_results.xlsx` report the complete results.

- **supp_blocking_images** contains the images for each cancer of the
block-analysis performed using the two-nn estimator.

- **blocking_ID.R** contains the code for the proposed ID-estimation by 
block-analysis.
