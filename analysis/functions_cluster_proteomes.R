# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#             
#     clusterData
#
# Displays the hierarchically clustered data by "pheatmap", 
#  (default) the number of clusters along the markers (num_clusters_row) and/or samples (num_clusters_col) are estimated by (multivariate) metafeatures "cafr", or
#  (optional) 'num_clusters_row' and 'num_clusters_col' can be set by user, clusters are estimated by pair-wise analysis.
# Returns the hierarchical tree and the clusters of markers and samples.
#
# output: tree, cluster_IDs_of_rows, cluster_IDs_of_columns
#
# NOTE: This function has a dependency "cafr": https://github.com/weiyi-bitw/cafr, https://doi.org/10.1371/journal.pcbi.1002920.
# 
clusterData = function(data, annotation_col = NULL, annotation_row= NULL, annotation_colors = NULL, 
                       stringency_col = 6, stringency_row = 4, clustering_distance_cols = 'euclidean',
                       clustering_distance_rows = 'euclidean',
                       num_clusters_col = NULL, num_clusters_row = NULL, main = 'Hierarchical cluster',
                       cluster_cols = T, cluster_rows = T, annotate_new_clusters_col = F){
  
  if (!is.null(num_clusters_col) | !cluster_cols) {att_col = NULL}
  if (!is.null(num_clusters_row) | !cluster_rows) {att_row = NULL}
  
  if (is.null(num_clusters_col) & cluster_cols) {att_col = cafr::attractorScanning(as.matrix(t(data)), maxIter = 1E2, epsilon = 1E-14, a = stringency_col); num_clusters_col = dim(att_col)[1]}
  if (is.null(num_clusters_row) & cluster_rows) {att_row = cafr::attractorScanning(as.matrix(data), maxIter = 1E2, epsilon = 1E-14, a = stringency_row); num_clusters_row = dim(att_row)[1]}
  
  if (is.null(num_clusters_col)) {num_clusters_col = 1}
  if (is.null(num_clusters_row)) {num_clusters_row = 1}
  
  tree = pheatmap::pheatmap(data, annotation_col = annotation_col, annotation_colors = annotation_colors, 
                            clustering_distance_cols = clustering_distance_cols, 
                            clustering_distance_rows = clustering_distance_rows,
                            cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                            cutree_cols = num_clusters_col, show_colnames = T, 
                            cutree_rows = num_clusters_row, show_rownames = T, main = main)
  
  if (cluster_cols) {cluster_IDs_col = cutree(tree$tree_col, k = num_clusters_col)} else {cluster_IDs_col = NULL}
  if (cluster_rows) {cluster_IDs_row = cutree(tree$tree_row, k = num_clusters_row)} else {cluster_IDs_row = NULL}
  
  if (annotate_new_clusters_col) {
    annotation_col_new = cbind(annotation_col,data.frame(cluster_IDs_col)[rownames(annotation_col),])
    colnames(annotation_col_new)[length(annotation_col_new)] = 'Putative.subtype'
    tree = pheatmap::pheatmap(data, annotation_col = annotation_col_new, annotation_colors = annotation_colors, 
                              clustering_distance_cols = clustering_distance_cols, 
                              clustering_distance_rows = clustering_distance_rows,
                              cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                              cutree_cols = num_clusters_col, show_colnames = T, 
                              cutree_rows = num_clusters_row, show_rownames = T, main = main)
  }
  
  return(list(tree, cluster_IDs_row, cluster_IDs_col, att_row, att_col))
}

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#     
#    dropMarkers
#
# drops markers with >percent_NA % missing values
# drops markers with low variability expression rates (low-mean and std or low-variance)
#
# output: data_filtered, dropped_marker_names
#
dropMarkers = function(data, percent_NA = .2, low_mean_and_std = .5, q_low_var = .75, force_drop = NULL){
  low_var = quantile(apply(data, 1, 'var', na.rm = T), q_low_var)
  if (!is.null(force_drop)) {data = data[!(rownames(data) %in% force_drop),]}
  data = data[rowSums(is.na(data))/dim(data)[2] < percent_NA,]
  dropped_markers = ((apply(data, 1, 'mean', na.rm = T) < low_mean_and_std) & (apply(data, 1, 'sd', na.rm = T) < low_mean_and_std)) | (apply(data, 1, 'var', na.rm = T) < low_var)
  return(list(data[!dropped_markers,], rownames(data)[dropped_markers]))
}


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     qNorm 
#
# quantile normalization
#
qNorm = function(data, method = 'quantile') {
  data_norm = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(data)))
  rownames(data_norm) = rownames(data)
  colnames(data_norm) = colnames(data)
  return(data_norm)
}


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     listout
#
# returns multiple objects from a function 
# by Gabor Grothendieck, 2004
#
# examples
#   list[QR,,QRaux]  <- qr(c(1,1:3,3:1))
#   list[,Green,Blue]  <- col2rgb('aquamarine')
#
listOut <- structure(NA,class='result')
'[<-.result' <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# 
#     uqNorm
# 
# Upper quantile normalization
#
uqNorm = function(data, q = .75) {
  data_qexp <- apply(data, 2, function(x){quantile(x[x>0 & !is.na(x)], q)})
  return(t(t(data) / data_qexp))
}


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# 
#     diffExp
#
# Rank-orders differentially expressed markers for a given subtype
#
# NOTE: This function has a dependency "limma": 10.18129/B9.bioc.limma
#
diffExp = function(data, annotations, subtypes, sub = 1, ann_levels = NULL) {
  if (is.null(ann_levels)){
    if (is.factor(annotations)){ann_levels = levels(annotations[,subtypes])
    }else{ann_levels = levels(as.factor(annotations[,subtypes]))}
  }
  data = data[,colnames(data) %in% rownames(annotations)]
  model = as.factor((annotations[colnames(data), subtypes]==ann_levels[sub])*1)
  design = stats::model.matrix(~0+model)
  cont_mat = limma::makeContrasts(subtype1 = 'model1-model0', levels = design)
  res = limma::topTable(limma::eBayes(limma::contrasts.fit(limma::lmFit(data, design), cont_mat)), adjust="BH", sort.by = 'logFC', number = dim(data)[1])
  return(res)
} 


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     diffExpSubtypes  
# 
# Lists differentially expressed markers for all subtypes
#
diffExpSubtypes = function(data, annotations, subtypes, criteria = 'logFC') {
  if (is.factor(annotations)){ann_levels = levels(annotations[,subtypes])
  }else{ann_levels = levels(as.factor(annotations[,subtypes]))}
  output = data.frame(matrix(ncol = length(ann_levels), nrow = dim(data)[1])); 
  colnames(output) = ann_levels
  rownames(output) = rownames(data)
  for (i in 1:length(ann_levels)) {
    res = diffExp(data, annotations, subtypes, sub = i, ann_levels)
    output[rownames(res),i] = res[,criteria]
  }
  return(output)
}

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     compareFCs
#
# Compare the proteins' fold changes cross cancer subtypes 
# or, within cancer subtypes list all proteins with FC > thr_FC  
#
compareFCs = function(diff_exp_A, diff_exp_B = NULL, sim_FC = .5, thr_FC = .5) {
  if (is.null(diff_exp_B)){flag_within_cancer = T; sim_FC = 1E-6; diff_exp_B = diff_exp_A}
  else{flag_within_cancer = F}
  common_markers = intersect(rownames(diff_exp_A),rownames(diff_exp_B))
  diff_exp_A = diff_exp_A[common_markers,]
  diff_exp_B = diff_exp_B[common_markers,]
  for (i in 1:dim(diff_exp_A)[2]){
    for (j in 1:dim(diff_exp_B)[2]){
      dist_FC_markers = abs(diff_exp_A[,i] - diff_exp_B[,j])
      similarFC_markers = which((dist_FC_markers < sim_FC) & !is.na(dist_FC_markers))
      for (k in similarFC_markers){
        if ((diff_exp_A[k,i] > thr_FC & diff_exp_B[k,j] > thr_FC) | (diff_exp_A[k,i] < -thr_FC & diff_exp_B[k,j] < -thr_FC)){
          if (!flag_within_cancer){cat( rownames(diff_exp_A)[k], '\tsimilarly differentially-expressed for\t', colnames(diff_exp_A)[i], '\tand\t', colnames(diff_exp_B)[j], '\n')}
          else{cat( rownames(diff_exp_A)[k], '\tdifferentially-expressed for\t', colnames(diff_exp_A)[i], '\n')}
        }
      }
    }
  }
}


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     plotBICs
#
# Plot BICs for various GMMs
# NAs are allowed for imputaion
#
plotBICs = function(dat, impute = T){
  if (any(is.na(dat)) & impute) {dat = mice::complete(mice::mice(dat, m = 1))}
  bic1 = mclust::mclustBIC(dat)
  plot(bic1, xlab = 'Num Clusters in Markers') 
  bic2 = mclust::mclustBIC(t(dat))
  plot(bic2, xlab = 'Num Clusters in Samples')
  return(list(bic1, bic2))
} 


# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     writeFCTable
#
# Write tables of differential expression patterns, 
# appended by a column of log FC variation over the subtypes
#
writeFCTable = function(dat, main){
  out = cbind(dat, data.frame(var_logFC = apply(dat, 1, 'var', na.rm = T)))
  write.table(x = out, file = main, sep = '\t', quote = F, row.names = T, col.names = T)
  return(out)
}

