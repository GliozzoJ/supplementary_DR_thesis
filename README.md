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

- Folder **code_DR+integration_pipelines** contains the code to run the proposed
DR + data fusion pipelines. Instructions to run the code are available below.

1. Install the following R packages:

```
# Packages installed from CRAN
install.packages(c("reticulate", "Boruta", "caret", "cutpointr", "dplyr", 
"doParallel", "entropy", "foreach", "genieclust", "ggplot2", "glmet", "kernlab", 
"igraph", "intrinsicDimension", "loe", "missRanger", "mixKernel", "SNFtool", 
"parallel", "PerfMeas", "plyr", "randomForest", "ranger", "Rdimtools", "readr", 
"readxl", "ROCR", "R.utils", "rsvd", "StatMatch", "stringr", "tsne", "umap", 
"tableone", "writexl", "intRinsic", "cytominer"))

# Packages installed from Bioconductor:
BiocManager::install("impute")
BiocManager::install("limma")
BiocManager::install("mixOmics")
BiocManager::install("phyloseq")
BiocManager::install("MOFA2")
BiocManager::install("preprocessCore")
```

Additionally CRAN package [RMKL](https://cran.r-project.org/web/packages/RMKL/index.html) 
is currently archived. You can install it from source after downloading the 
file "RMKL_1.0.tar.gz":

```
install.packages("./RMKL_1.0.tar.gz", repos = NULL, type="source")
```

2. Download the datasets from this [link](https://drive.google.com/drive/folders/1zNxg-DBXKWsolag_4EAR26otTyUnPJTy?usp=sharing) and
move the files into the folder "code_DR+integration_pipelines/data".

3. In the file "main.R" you need to modify the settings to use required Python packages
through "reticulate" interface. Minimal instructions are available at the beginning 
of the file.

4. Load the code:

```
source('main.R')
```

5. You can now use the function main_best_models() to run the analyses. Intermediate 
results are automatically saved in the folder "code_DR+integration_pipelines/results".

Here there is also an examples on a toy matrix to run only 
the ID-estimation by block-analysis:
1. Load the code:

```
source("utilities.R")
source("compute_distance.R")
source("estimate_ID_twoNN.R")
source("blocking_ID.R")
```

2. Create a toy matrix (samples x features):

```
set.seed(123)

M <- matrix(rnorm(3000), nrow=30)
rownames(M) <- paste0("P", 1:nrow(M))
colnames(M) <- paste0("F", 1:ncol(M))

```

3. Compute the unbiased estimate of the intrinsic dimensionality:

```
ID_orig <- estimate_ID_twoNN(mat_data = M)$id
```

4. Run the blocking-id analysis:

```
factor = 2
args_ID = list(dist_fun_twoNN= 'canberra', perc_points = 0.9, maxit = 11, 
               ncores = min(8, detectCores()-1))

results <- blocking_ID(M, ID_orig = ID_orig, factor = factor, L0 = round(ID_orig*factor)) 
```

