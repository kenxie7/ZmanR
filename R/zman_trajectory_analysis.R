# Smooth Zman AUC trajectory:
# Author: Ken Xie & Assaf Weiner

require(ComplexHeatmap)
require(grid)
require(circlize)
require(ggarchery)
require(ggnewscale)

new_scale <- function(new_aes) {
   structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

normalize_mm <- function (x, na.rm = TRUE) 
{
  return((x - min(x))/(max(x) - min(x)))
}
# Function for normalizing mc expression data
cl_compute_e_gc <- function (cl_mc, us, norm_by_mc_meansize = T) 
{
  f_g_cov = Matrix::rowSums(us) > 10
  e_gc = t(tgstat::tgs_matrix_tapply(us[f_g_cov, ], cl_mc, function(y) {
    exp(mean(log(1 + y))) - 1
  }))
  rownames(e_gc) = rownames(us)[f_g_cov]
  if (norm_by_mc_meansize) {
    mc_meansize = tapply(Matrix::colSums(us), cl_mc, mean)
    e_gc = t(t(e_gc)/as.vector(mc_meansize))
  }
  return(e_gc)
}

require(tgstat)
#tgs_matrix_tapply <- function (x, index, fun, ...) 
#{
#    if (missing(x) || missing(index) || missing(fun)) 
#        stop("Usage: tgs_matrix_tapply(x, index, fun)", call. = F)
#    args <- as.list(substitute(list(...)))[-1L]
#    if (!is.factor(index)) 
#        index <- factor(index)
#    .Call("tgs_matrix_tapply", x, index, fun, as.character(substitute(fun)), 
#        args, new.env(parent = parent.frame()))
#}

corr_dist <- function(x, y = NULL, method = NULL, use = "everything")
{
  x <- t(as.matrix(x))
  if(!is.null(y))
    y <- t(as.matrix(y))
  1 - (stats::cor(x, y) + 1) / 2
}

auc_order_path <- function(avg_points, nc){
  order_path <- order(avg_points$mean_pts[,3], decreasing = T)

  ordResult <- (match(c(1:nc), order_path))[avg_points$cl_res]
  tmp_dist <- sum(diff(avg_points$mean_pts[,3][order_path]))
  #print(tmp_dist)
  normResult <- normalize_mm(ordResult)
  return(list(raw_path=ordResult, norm_path = normResult, ord_path = order_path))
}

get_k_average_results <- function(cell_gene_exprs, df_mc_2d, k = 3, cl_res = NULL){
  if(is.null(cl_res)){
    #dist_mat <- redPATH:::corr_dist(cell_gene_exprs)
    dist_mat <- corr_dist(cell_gene_exprs)
    
    scaled_distance <- .Call(stats:::C_DoubleCentre, dist_mat^2)
    
    fit_hclust <- hclust(as.dist(as.matrix(dist(cell_gene_exprs))))
    cluster_result <- as.integer(cutree(fit_hclust, k))
    
  }else{
    cluster_result <- cl_res
  }
  #print(cluster_result)
  # Calculate average 2d for AUC:
  store_points <- c()
  tmp_df_mc_2d = df_mc_2d#[match(tmp_toshow, to_show), ]
  #cluster1 <- (as.numeric(factor(ab20_mf_annotations_v2$group[tmp_toshow], levels =c("Monocytes", "MonMAC1", "Cd72_Acp5_TAM", "MonMAC2", "Arg1_TAM", "Gpnmb_TAM"))))
  
  for(i in c(1:k)){
    avg_point <- apply(tmp_df_mc_2d[which(cluster_result==i),c(1:2,match("norm_auc", colnames(tmp_df_mc_2d)))], 2, mean)
    store_points <- rbind(store_points, avg_point)
    
  }
  
  # Calculate average expression:
  avg_exprs <- apply(cell_gene_exprs, 2, function(a) aggregate(a, list(cluster_result), FUN=mean)) 
  store_avg_exprs <- c()
  for(i in c(1:length(names(avg_exprs)))){
    store_avg_exprs <- cbind(store_avg_exprs, avg_exprs[[i]]$x)
  }
  #print(dim(store_avg_exprs))
  colnames(store_avg_exprs) <- names(avg_exprs)
  
  return(list(mean_exprs = store_avg_exprs, mean_pts = store_points, cl_res = as.numeric(as.factor(cluster_result))))
  
}

#' Computing the CDF / AUC value for each metacell
#'
#' This function calculates the AUC time value for each metacell using the time bin assignment information from FACS
#'
#' @param new_id id for loading the relevant metacell objects for scdb_mc and scdb_mat
#' @param time_per_cell This is the annotated dataframe for each single cell with the time group assignment 
#' @param select_mcs Take the relevant mc id for normalization. Here we take the T cells and NK cells together for normalization as an example
#' @param mc_annotations Dataframe with cell type annotations for each metacell id
#' @param time_points A vector containing the unique time points that you have in the experiment (ie c("12H", "24H", "36H"))
#' @param time_for_auc Corresponding numeric vector that provides the time values with an additional start of 0 for calculation of AUC (ie c(0,12,24,36))
#' @param min_cell_per_mc_mouse Threshold for minimum number of cells per mc of a mouse. Set to 10 for the NK example dataset.
#' @param Gate Use CD45 gating for normalization. Defaults to NULL if all cells have been gated as CD45. The column "Gating" is used for subseting
#'
#' @export
#' @examples
#' # Assuming `my_data` is your FACS data
#' new_id = "T_clean"
#' time_per_cell = well_fcs_GBM_T_time
#' mc_cdf <- compute_mc_cdf(new_id, time_per_cell, select_mcs = 1:37, mc_annotations = mc_annotations, time_points = c("12H","24H","36H"), time_for_auc = c(0,12,24,36), Gate="CD45_high")
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix
#' @import dplyr
#' @import metacell
compute_mc_cdf <- function(new_id, time_per_cell,
                           select_mcs = 1:37, mc_annotations=NULL, #to_show = c(13:37),
                           time_points = c("12H","24H","36H"), 
                           time_for_auc = c(0,12,24,36),min_cell_per_mc_mouse = 10,
                           Gate = NULL) {
  mat_all = scdb_mat(new_id)
  mc = scdb_mc(new_id)
  mc2d = scdb_mc2d(new_id)
  
  shared_names = intersect(names(mc@mc), rownames(time_per_cell))
  #print(length(shared_names))
  cell_metadata = mat_all@cell_metadata[shared_names,]
  cell_metadata[shared_names, "time_group"] = time_per_cell[shared_names, "groundtruth_group"]
  mouseID = paste0('TimeGroup-', cell_metadata$time_group, '-Treatment-', cell_metadata$Treatment, '-m', cell_metadata$Mouse)
  names(mouseID) <- rownames(cell_metadata)
  if(!is.null(Gate)){
    cell_metadata <- cell_metadata[which(cell_metadata$Gating == Gate),]
  }
  
  
  #print((mouseID))
  
  mouse_per_mc = table(mouseID[shared_names], mc@mc[shared_names])
  #print(mc@mc[shared_names])
  
  mouse_per_mc = mouse_per_mc[rowSums(mouse_per_mc) > min_cell_per_mc_mouse,]
  #print(mouse_per_mc)
  mouse_per_mc_norm = (0.01 + mouse_per_mc) / rowSums(0.01 + mouse_per_mc)
  mouse_per_mc_norm_smooth = mouse_per_mc_norm
  #print(head(mouse_per_mc_norm_smooth))
  mc2d_graph = mc2d@graph
  mc2d_df = data.frame(cbind(mc2d@mc_x, mc2d@mc_y))
  colnames(mc2d_df) = c("x", "y")
  for (j in c(1:3)) {
    mouse_per_mc_norm_pre = mouse_per_mc_norm_smooth
    for (i in 1:max(mc@mc)) {
      if (length(mc2d_graph$mc2[mc2d_graph$mc1 == i]) > 1) {
        neigh = mc2d_graph$mc2[mc2d_graph$mc1 == i]
        neigh2 = mc2d_graph$mc2[mc2d_graph$mc1 %in% neigh]
        mouse_per_mc_norm_smooth[, i] = 0.6 * mouse_per_mc_norm_pre[, i] + 0.3 * rowMeans(mouse_per_mc_norm_pre[, neigh]) + 0.1 * rowMeans(mouse_per_mc_norm_pre[, neigh2])
      }
    } 
  }
    
 
  #print(head(mouse_per_mc_norm_smooth))
  mouse_per_mc_norm_smooth_tumor = mouse_per_mc_norm_smooth[, select_mcs] / rowSums(mouse_per_mc_norm_smooth[, select_mcs])
  #print(mouse_per_mc_norm_smooth_tumor)
  df_mc = data.frame(mouse_per_mc_norm_smooth_tumor)
  #print(df_mc)
  
  ab_time = c()
  treatment = c()
  mouse = c()
  for (i in 1:dim(df_mc)[1]){
    ab_time[i] = paste(strsplit(as.character(df_mc$Var1[i]),'-')[[1]][2],collapse = '-')
    treatment[i] = paste(strsplit(as.character(df_mc$Var1[i]),'-')[[1]][4],collapse = '-')
    mouse[i] = paste(strsplit(as.character(df_mc$Var1[i]),'-')[[1]][5],collapse = '-')
  }
  
  df_mc$ab_time = ab_time
  df_mc$treatment = treatment
  df_mc$mouse = mouse
  
  time_cnt = table(cell_metadata$time_group)
  norm_constant <- 1/ (table(cell_metadata$time_group[match(names(mc@mc), rownames(cell_metadata))])[-length(time_cnt)] / sum(time_cnt[-length(time_cnt)]))
  norm_constant_mc <- rep(norm_constant[-length(time_cnt)], length(select_mcs))
  x1 = df_mc[df_mc$ab_time %in% time_points,]
  x1$Freq = x1$Freq * norm_constant_mc
  #print(x1)
  
  
  gd_nk <- x1 %>%
    select(Var2, Freq,ab_time) %>%
    group_by(Var2,ab_time) %>%
    summarise(Freq = mean(Freq)) %>%
    mutate(prop=Freq/sum(Freq))
  
  # Calc AUC
  x=time_for_auc
  mc_nk = unique(gd_nk$Var2)
  cm_nk_df = data.frame()
  for (i in mc_nk){
    c1 = cumsum=c(0,unlist(cumsum(gd_nk[gd_nk$Var2==i,'prop'])))
    tmp_df = data.frame(row.names=paste0(i,'_',x),time=x,mc=i,
                        cums=c1,auc = sum(c1*(36-x)))
    cm_nk_df= rbind(cm_nk_df,tmp_df)
  }
  
  
  df_mc_2d = data.frame(x=mc2d@mc_x,y=mc2d@mc_y,mc = names(mc2d@mc_y))
  cm_nk_df_0 = cm_nk_df[cm_nk_df$time==0,]
  
  df_mc_2d$auc = cm_nk_df_0$auc#[to_show]
  df_mc_2d$norm_auc = normalize_mm(df_mc_2d$auc) 
  if(!is.null(mc_annotations)){
    mc_annotations = mc_annotations[mc_annotations$mc_id %in% select_mcs, ]

    cm_nk_df$celltype <- mc_annotations$celltype[as.numeric(cm_nk_df$mc)]
    df_mc_2d$celltype <- mc_annotations$celltype[as.numeric(df_mc_2d$mc)]
    
  }
  #length(unique(cm_nk_df$mc))
  
  #norm_constant <- ... # (Rest of the code for computing norm_constant and other variables)
  
  #cm_nk_df = ... # (Rest of the code for computing cm_nk_df)
  
  return(list(mc_auc_time = cm_nk_df, mc2d_auc_time = df_mc_2d, 
              mc2d_graph = mc2d_graph, mc2d_df = mc2d_df))
}

#' Re-normalizing MC data for Zman trajectory
#'
#' This function prepares the data required for smoothing the Zman AUC derived trajectory
#'
#' @param new_id id for loading the relevant metacell objects for scdb_mc and scdb_mat
#' @param mc_annotations Dataframe with cell type annotations for each metacell id
#' @param mc2d_auc_time Dataframe containing the coordinates of the metacells and its corresponding AUC values and cell type information
#' @param select_celltypes Select cell types used for constructing the trajectory. i.e. we subset NKs or myeloid cells separately to look at their temporal profile.
#'
#' @export
#' @examples
#' # Normalize mc expression data for selected lineage of cell types
#' NK_mc_exprs = normalize_mc_exprs(new_id, mc_annotations, mc_cdf$mc2d_auc_time, select_celltypes = c("chemotactic", "cytotoxic", "intermediate", "dysfunctional"))
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix
#' @import dplyr
#' @import metacell
#' @import tgstat
normalize_mc_exprs <- function(new_id, mc_annotations, mc2d_auc_time, select_celltypes = c("chemotactic", "cytotoxic", "intermediate", "dysfunctional"), norm_by_mc_meansize = TRUE) {
  # Select markers and labels based on specified groups
  mc=scdb_mc(new_id)
  mat_all=scdb_mat(new_id)
  select_mcs <- mc_annotations$mc_id[which(mc_annotations$celltype %in% select_celltypes)]
  select_labels <- mc_annotations$celltype[select_mcs]
  
  # Create a subset DataFrame
  select_df <- mc2d_auc_time[match(select_mcs, mc2d_auc_time$mc),]
  
  # Adjust mc data
  tmp_mcr <- mc
  tmp_mcr@mc <- mc@mc[which(mc@mc %in% select_mcs)]
  tmp_mcr@cell_names <- mc@cell_names[match((names(mc@mc) %in% select_mcs), mc@cell_names)]
  
  # Create a subset of the matrix
  tmp_matr <- as.matrix(mat_all@mat)
  tmp_matr <- tmp_matr[,match(names(tmp_mcr@mc), colnames(tmp_matr))]
  
  # Compute expression values
  tmp_egc <- cl_compute_e_gc(tmp_mcr@mc, tmp_matr, norm_by_mc_meansize = norm_by_mc_meansize)
  select_exprs <- log2(1 + tmp_egc * 5000)
  
  # Return the relevant data
  return(list(select_df = select_df, select_exprs = select_exprs, select_labels = select_labels))
}

#' Smooth Zman Trajectory
#'
#' This function smoothes the AUC derived time of metacells to construct the Zman trajectory of the metacells
#'
#' @param disp_go_df Gene expression matrix after selection of GO genes and dispersion filter
#' @param normed_mc_res Result list from normalize_mc_exprs containing the raw AUC dataframe and corresponding metacell cell type labels
#' @param ref_k The number of k to build the reference path corresponding to the number of clusters, defaults to 4
#'
#' @export
#' @examples
#' # Smooth AUC trajectory:
#' # First we filter genes by dispersion to reduce the number of genes to more significant ones:
#' select_GO_mc <- get_GO_exp(select_exprs,gene_type = "gene_names", organism = "mouse", takelog=F)
#' # The threshold is set roughly at 90% quantile of the mc, since metacells have fewer cells and pooled expression, we expect in general lower dispersion values
#' filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threadnum = 1)
#' NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs, ref_k = 4)
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix
#' @import dplyr
smooth_zman_trajectory <- function(disp_go_df, normed_mc_res, ref_k = 4){
  select_df = normed_mc_res$select_df
  select_labels = normed_mc_res$select_labels
  #select_labels = as.numeric(as.factor(select_labels))
  avg_points <- get_k_average_results(disp_go_df, select_df, k = ref_k, cl_res = as.numeric(as.factor(select_labels)))
  order_path = auc_order_path(avg_points, ref_k)
  ref_path <- order_path$norm_path
  #print(ref_path)
  avg_points$ref_path = ref_path#order_path$raw_path
  avg_points$ord_path = order_path$ord_path
  avg_points$mean_pts = data.frame(avg_points$mean_pts)
  avg_points$mean_pts$ref_auc = as.vector(aggregate(ref_path, list(as.numeric(as.factor(select_labels))), FUN=mean)$x)
  avg_points$mean_pts$smoothed_auc = 1 - normalize_mm(avg_points$mean_pts$norm_auc)
  store_ai_paths <- c()
  #seq(11, 66, )
  for(i in c((ref_k+1):nrow(select_df))){
    tmp_avg <- get_k_average_results(disp_go_df, select_df, k = i, cl_res = NULL)
    #print(i)
    #tmp_path <- guided_paths_ai(tmp_avg,nc = i, base_path_reference = k_paths, alpha = 1, beta=2)
    tmp_path <- auc_order_path(tmp_avg, nc = i)
    store_ai_paths <- rbind(store_ai_paths, tmp_path$norm_path)
  }
  smoothed_trajectory <- colMeans(rbind(ref_path*1, 1*normalize_mm(colMeans(store_ai_paths[which(apply(store_ai_paths, 1, function(a) cor(ref_path, a, method = "spearman"))>0.75),]))))
  #smoothed_trajectory <- colMeans(rbind(ref_path*1, 1*normalize_mm(colMeans(store_ai_paths[which(apply(store_ai_paths, 1, function(a) cor(ref_path, a, method = "spearman"))>0.75),]))))
  select_df$smoothed_auc = smoothed_trajectory
  return(list(smoothed_auc = smoothed_trajectory, mc2d_auc_time = select_df, reference_path=avg_points))
}

#' Calculate correlated genes along time
#'
#' This function calculates the correlation between the gene expression and the time
#'
#' @param new_id ID from metacell object
#' @param smoothed_res Result object from smooth_zman_trajectory()
#' @param method Defaults to "spearman", but can be "pearson" or "kendall" as seen from cor.test function
#'
#' @export
#' @examples
#' # Smooth AUC trajectory:
#' # First we filter genes by dispersion to reduce the number of genes to more significant ones:
#' select_GO_mc <- get_GO_exp(select_exprs,gene_type = "gene_names", organism = "mouse", takelog=F)
#' # The threshold is set roughly at 90% quantile of the mc, since metacells have fewer cells and pooled expression, we expect in general lower dispersion values
#' filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threadnum = 1)
#' NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs, ref_k = 4)
#' NK_corr = calculate_corr_genes(new_id, NK_smoothed_res, "spearman") # list containing both correlation and pvalues vectors for each gene
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix

calculate_corr_genes <- function(new_id, smoothed_res, method = "spearman"){
    mc_fp = scdb_mc(new_id)@mc_fp
    mc_fp = mc_fp[, match(smoothed_res$mc2d_auc_time$mc, colnames(mc_fp))]
    
    store_cor_c1 <- c()
    store_cor_p1 <- c()
    defaultW <- getOption("warn")
    options(warn = -1)

    for(i in c(1:nrow(mc_fp))){
        ct = cor.test(x = sort(smoothed_res$mc2d_auc_time$smoothed_auc), y = mc_fp[i, order(smoothed_res$mc2d_auc_time$smoothed_auc)], method = "spearman")
        store_cor_p1 <- c(store_cor_p1, ct$p.value)
        store_cor_c1 <- c(store_cor_c1, ct$estimate)
    }
    options(warn = defaultW)
    names(store_cor_p1) <- rownames(mc_fp)
    names(store_cor_c1) <- rownames(mc_fp)
    return(list(correlation = store_cor_c1, pvalues = store_cor_p1))
}

#' Predict gene expression across time for insights and visualization
#'
#' This function utilizes the smoothed trajectory to generate predicted expression along time
#'
#' @param new_id id for loading the relevant metacell objects for scdb_mc and scdb_mat
#' @param smoothed_res Result object from smooth_zman_trajectory()
#' @param out_len The total of outcome observations along time, ie we usually set a number close to the total number of metacells to prevent over/under fitting
#' @param loess_degree Loess parameter to set the degree of fitting
#' @param loess_span Loess parameter to set spanning parameter
#' @param reverse_traj If the reverse coloring is preferred, set TRUE. Defaults to FALSE
#' @export
#' @examples
#' filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threadnum = 1)
#' NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs$select_df, ref_k = 4)
#' NK_predicted_res = predict_expression_along_time(new_id, NK_smoothed_res, out_len = 36, loess_degree=1, loess_span = .9)
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats                                                                                                   
#' @import Matrix
#' @import dplyr
predict_expression_along_time <- function(new_id, smoothed_res, out_len = 36, loess_degree=1, loess_span = .9, reverse_traj=FALSE) {
  tmp_df_mc_2d = smoothed_res$mc2d_auc_time
  mc_fp_nk = scdb_mc(new_id)@mc_fp
  tmp_fp_nk_only <- mc_fp_nk[,match(tmp_df_mc_2d$mc, colnames(mc_fp_nk))]
  tmp_df_mc_2d_treated <- tmp_df_mc_2d
  
  new_x <- seq(min(tmp_df_mc_2d_treated$smoothed_auc), max(tmp_df_mc_2d_treated$smoothed_auc), len = out_len)
  mc_predict_nk <- data.frame(new_x)
  if(reverse_traj){
    
    rev_auc = 1 - as.numeric(tmp_df_mc_2d_treated$smoothed_auc)
  }else{
    rev_auc = as.numeric(tmp_df_mc_2d_treated$smoothed_auc)
    
  }
  for (i in 1:(nrow(mc_fp_nk) - 1)) {
    tmp_dat <- data.frame(gene = tmp_fp_nk_only[i,],
                          auc = rev_auc)
    fit_i <- loess(gene ~ auc, data = tmp_dat, degree = loess_degree, span = loess_span)
    smooth_i <- predict(fit_i, newdata = data.frame(auc = new_x))
    mc_predict_nk <- cbind(mc_predict_nk, smooth_i)
    colnames(mc_predict_nk)[i + 1] <- rownames(mc_fp_nk)[i]
  }
  
  mc_predict_nk <- mc_predict_nk[, -1]
  
  normalized_predicted_exprs <- mc_predict_nk
  for (i in 1:ncol(normalized_predicted_exprs)) {
    normalized_predicted_exprs[, i] <- normalize_mm(normalized_predicted_exprs[, i])
  }
  
  return(normalized_predicted_exprs)
}


#' Plot gene expression across time for insights and visualization
#'
#' This function utilizes the smoothed trajectory to plot the gene expression trend across time using loess fitting
#'
#' @param predicted_exprs Result object from predict_expression_along_time()
#' @param up_regulated_genes Genes identified from correlation analysis that are up-regulated along time
#' @param down_regulated_genes Genes identified from correlation analysis that are up-regulated along time
#' @param k Number of gene groups for clustering separation, defaults to 2 for up- and down- regulated genes
#' @param row_size Size of text for genes (Defaults to 5.5)
#' @export
#' @examples
#' filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threadnum = 1)
#' NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs$select_df, ref_k = 4)
#' NK_predicted_res = predict_expression_along_time(new_id, NK_smoothed_res, out_len = 36, loess_degree=1, loess_span = .9)
#' plot_zman_genes_heatmap(NK_predicted_res, NK_smoothed_res, c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), c("Ccl3","Prf1","Gzma", "Gzmb","Nkg7"), k = 2)
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix
#' @import dplyr
#' @import circlize
#' @import ComplexHeatmap
#' @import grid
plot_zman_genes_heatmap <- function(predicted_exprs,smoothed_res, up_regulated_genes, down_regulated_genes, k=2, row_size=5.5){
  predicted_exprs = t(predicted_exprs)
  uniq_genes = unique(c(up_regulated_genes, down_regulated_genes))
  predicted_exprs = predicted_exprs[match(uniq_genes, rownames(predicted_exprs)), ]
  out_len = ncol(predicted_exprs)
  tmp_df_mc_2d = smoothed_res$mc2d_auc_time
  new_x <- seq(min(tmp_df_mc_2d$smoothed_auc),max(tmp_df_mc_2d$smoothed_auc), len=out_len)
  col_fun  <-  circlize::colorRamp2(seq(from = min(new_x), to = max(new_x), length.out = out_len), viridis::viridis(out_len, direction=-1))
  time_an  <-  ComplexHeatmap::HeatmapAnnotation(df = data.frame(time = new_x), border = FALSE, col = list(time = col_fun) )
  
  h1 <- ComplexHeatmap::Heatmap(predicted_exprs, row_split = factor(cutree(hclust(dist(predicted_exprs)), k)),
                                col = circlize::colorRamp2(seq(0,1, by=0.05), viridis::rocket(40, direction = -1)[1:21]),
                                heatmap_legend_param = list(title = "LFP"),#, at = seq(-3.5, 3.5, 1)),
                                bottom_annotation = time_an,
                                show_column_names=F,cluster_columns =F, #column_names_gp = gpar(fontsize = col_size),
                                row_names_gp = gpar(fontsize = row_size),#col = c("green", "orange"), fontsize = c(10, 14)),
  ) 
  return(h1)
}
#' Plot gene expression across time for insights and visualization
#'
#' This function utilizes the smoothed trajectory to plot the gene expression trend across time using loess fitting
#'
#' @param smoothed_res Result object from smooth_zman_trajectory()
#' @param mc_cdf Result object from compute_mc_cdf()
#' @param mc_only Defaults to FALSE, showing the trajectory arrows. If TRUE, it only shows the graph of metacells
#' @export
#' @examples
#' filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threadnum = 1)
#' NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs$select_df, ref_k = 4)
#' NK_predicted_res = predict_expression_along_time(new_id, NK_smoothed_res, out_len = 36, loess_degree=1, loess_span = .9)
#' plot_zman_genes_heatmap(NK_predicted_res, NK_smoothed_res, c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), c("Ccl3","Prf1","Gzma", "Gzmb","Nkg7"), k = 2)
#' plot_smoothed_trajectory(NK_smoothed_res, mc_cdf, mc_only=FALSE)
#' @section Trajectory Analysis
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import stats
#' @import Matrix
#' @import dplyr
#' @import circlize
#' @import ComplexHeatmap
#' @import grid
#' @import ggarchery
#' @import ggnewscale                                                                                        
plot_smoothed_trajectory <- function(smoothed_res, mc_cdf_res, mc_only = FALSE){
    
    ############ df for plotting 2d with connections:
    select_mcs = smoothed_res$mc2d_auc_time$mc
    full_mc2d = mc_cdf_res$mc2d_graph
    mc2d_df = mc_cdf_res$mc2d_df
    full_mc2d = full_mc2d[which( (full_mc2d[,1] %in% select_mcs) & (full_mc2d[,2] %in% select_mcs)),]
    full_mc2d = full_mc2d[-which(full_mc2d[,1]==full_mc2d[,2]),]
    full_sgtbl <- c()
    for(i in c(1:nrow(full_mc2d))){
        full_sgtbl <- rbind(full_sgtbl, cbind(mc2d_df[full_mc2d[i,1],1],mc2d_df[full_mc2d[i,1],2],
                                                 mc2d_df[full_mc2d[i,2],1],mc2d_df[full_mc2d[i,2],2]))
    }
    full_sgtbl <- data.frame(full_sgtbl)
    colnames(full_sgtbl) = c("x", "y", "xend", "yend")


    full_sgtbl <- tibble(xs=full_sgtbl$x,
                ys = full_sgtbl$y,
                x=full_sgtbl$x,
                y = full_sgtbl$y,
                xend = full_sgtbl$xend,
                yend = full_sgtbl$yend)
    full_sgtbl <- data.frame(full_sgtbl)
    ############ 
    sg_tbl <- c()
    store_points = smoothed_res$reference_path$mean_pts
    #print(head(store_points))
    order_path = smoothed_res$reference_path$ord_path
    
    for(i in c(1:(length(order_path) - 1))){
        sg_tbl <- rbind(sg_tbl, cbind(store_points[order_path[i], 1:2], store_points[order_path[i+1], 1:2]))
    }
    sg_tbl <- rbind(sg_tbl, cbind(store_points[order_path[i+1], 1:2], store_points[order_path[i+1], 1:2]))
    #print(head(sg_tbl))
    #print(sg_tbl)
    colnames(sg_tbl)[3:4] = c("xend", "yend")
    sg_tbl <- tibble(xs=sg_tbl$x,
                ys = sg_tbl$y,
                xend = sg_tbl$xend,
                yend = sg_tbl$yend)
    tmp_table <- data.frame(cbind(store_points, sg_tbl))

    #return(list(subset_tab = tmp_table, full_tab = full_sgtbl))
    if(!mc_only){
    traj_g01 <- 
        ggplot()+
        geom_segment(data = full_sgtbl,aes(x = xs, xend = xend, y = ys, yend = yend), 
                      #arrow_positions = c(0.5, 0.6), 
                    #arrows = arrow(length = unit(0.1, "inches")), 
                    size=.2, color='lightgrey', alpha = .7,
                    position = position_attractsegment(start_shave = 0.05, 
                                                       end_shave = 0.05, 
                                                       type_shave = "distance"))+
       geom_point(data=smoothed_res$mc2d_auc_time, aes(x = x, y = y, fill = smoothed_auc),size = 6.5, shape=21, color = 'black',alpha=0.5,)+# alpha=0.5,
       scale_fill_viridis_c(direction=-1, limits=c(0,1))+
        new_scale_fill()+
        new_scale_color()+

    geom_point(data=tmp_table, aes(x=x, y=y, fill=smoothed_auc, col=NULL), shape = 21, size = 10)+#,0.15,0.15
    labs(fill="smoothed_avg_auc")+
    scale_fill_viridis_c(direction=-1, limits=c(0,1))+
    new_scale_fill()+
    new_scale_color()+
    #[-c(2,3),]
    geom_arrowsegment(data = sg_tbl[1:(nrow(sg_tbl)-1),],aes(x = xs, xend = xend, y = ys, yend = yend), 
                      arrow_positions = c(.7), #color='black',
                    arrows = arrow(length = unit(0.2, "inches"), type = "closed"), 
                    size=.6,
                    position = position_attractsegment(start_shave = 0.05, 
                                                       end_shave = 0.05, 
                                                       type_shave = "distance"))+
     guides(alpha = "none")+
    theme(text = element_text(size = 20), axis.line=element_blank(), axis.ticks=element_blank(),
         axis.title = element_blank(), axis.text=element_blank())#ggplot(data=data.frame(store_points), aes(x=x, y=y, fill=auc))+geom_point(shape = 21)+scale_fill_viridis_c()
    }else{
    traj_g01 <- 
        ggplot()+
        geom_segment(data = full_sgtbl,aes(x = xs, xend = xend, y = ys, yend = yend), 
                      #arrow_positions = c(0.5, 0.6), 
                    #arrows = arrow(length = unit(0.1, "inches")), 
                    size=.2, color='lightgrey', alpha = .7,
                    position = position_attractsegment(start_shave = 0.05, 
                                                       end_shave = 0.05, 
                                                       type_shave = "distance"))+
       geom_point(data=smoothed_res$mc2d_auc_time, aes(x = x, y = y, fill = smoothed_auc),size = 6.5, shape=21, color = 'black',alpha=0.5,)+# alpha=0.5,
       scale_fill_viridis_c(direction=-1, limits=c(0,1))+
        guides(alpha = "none")+
    theme(text = element_text(size = 20), axis.line=element_blank(), axis.ticks=element_blank(),
         axis.title = element_blank(), axis.text=element_blank())#ggplot(data=data.frame(store_points), aes(x=x, y=y, fill=auc))+geom_point(shape = 21)+scale_fill_viridis_c()

    }
    return(traj_g01)
}