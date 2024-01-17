#' Compute distance matrix
#' 
#' @description Compute distance matrix from a features matrix (samples x 
#' features).
#'
#' @param mat_data matrix. Matrix (samples x features).
#' @param method character. String with the name of the distance metric or 
#' correlation to use.
#' Possible choices are: 
#' - Distance metric: "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski"
#' - correlation: "pearson", "kendall", "spearman"
#' 
#' In case a correlation is used, the obtained correlation matrix M is converted
#' into a distance matrix by M = max(abs(M))-M.
#' 
#' @return Distance matrix.
#' @export
compute_distance <- function(mat_data, method = 'euclidean'){
    if (method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
        #cat('computing distance ', method, '\n')
        dm <- as.matrix(dist(mat_data,method = method))
    }else{
        #cat('using correlation based distance: ', method, '\n')
        if (method %in% c("pearson", "kendall", "spearman")){
            dm <- abs(cor(t(mat_data), method = method, use = "pairwise.complete.obs"))
            dm = max(dm)-dm
        }else{
            cat("distance function must be any metric distance [",
                c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski"), '] or correlation distance [',
                c("pearson", "kendall", "spearman"), '\n', 'I am using standard euclidean\n')
            dm <- as.matrix(dist(mat_data,method = "euclidean"))
        }
        
    }
    return(dm)
}
