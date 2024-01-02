# Pre process functions for the gene expression matrix to obtain the input to redPATH

#load("preprocessed_GO/GO_differential_genes.RData")
#load("preprocessed_GO/GO_differential_genes_hs.RData")
#load("preprocessed_GO/processed_GO.RData")
require(doParallel)

#' Preprocess Gene Expression
#' 
#' This function reduces the gene expression matrix to selected Gene Ontology genes. 
#'@param gene_cell_matrix A genes (row) by cells (col) matrix.
#'@param gene_type The type of gene names. Default is "gene_names" - ie "HOPX", alternatively "entrez" is provided.
#'@param organism The type of species the dataset belongs to: defaults to "mouse". "human" option is available.
#'@param takelog Defaults to TRUE. Whether to perform log2(exprs + 1)
#'@keywords get_GO_exp
#'@section Preprocessing:
#'@export
#'@examples
#'subset_GO_data <- get_GO_exp(data, gene_type = "gene_names", organism = "human", takelog = F)
#'filtered_GO_data <- filter_exp(subset_GO_data, dispersion_threshold = 10, threads = 2)
#'
get_GO_exp <- function(gene_cell_matrix, gene_type = "gene_names", organism = "mouse", takelog = T){
  gene_names <- rownames(gene_cell_matrix)
  
  if(organism == "mouse"){
    if(gene_type == "gene_names"){
      matched_genes <- match(tolower(GO_unique_mgi), tolower(gene_names))
    }else if(gene_type == "entrez"){
      matched_genes <- match(tolower(GO_unique_entrez), tolower(gene_names))
    }
    m1_check <- which(is.na(matched_genes))
    if(sum(m1_check)==0){
      matched_genes = matched_genes
    }else{
      matched_genes <- matched_genes[-which(is.na(matched_genes))]   
    }
    
    diff_exp <- as.matrix(t(gene_cell_matrix[matched_genes, ]))
    
  }else if(organism == "human"){
    if(gene_type == "gene_names"){
      matched_genes <- match(tolower(GO_unique_hg), tolower(gene_names))
    }else if(gene_type == "entrez"){
      matched_genes <- match(tolower(GO_unique_entrez_hs), tolower(gene_names))
    }
    m1_check <- which(is.na(matched_genes))
    if(sum(m1_check)==0){
      matched_genes = matched_genes
    }else{
      matched_genes <- matched_genes[-which(is.na(matched_genes))]   
    }
 
    diff_exp <- as.matrix(t(gene_cell_matrix[matched_genes, ]))
    
  }

  if(takelog == T){
    diff_exp <- log2(diff_exp + 1)
  }
  return(diff_exp)
}


#' Filter genes
#' 
#' This function reduces the dimensionality by filtering out low dispersion ratio genes. 
#'@param cell_gene_matrix A cell (row) by gene (col) matrix. Typically the result returned from get_GO_exp.
#'@param dispersion_threshold Defaults to 10. Depending on the type of normalization / data, the ratio should be changed to retain a few hundred genes.
#'@param threads The number of threads to use. 
#'@keywords filter_exp
#'@section Preprocessing:
#'@export
#'@examples
#'subset_GO_data <- get_GO_exp(data, gene_type = "gene_names", organism = "human", takelog = F)
#'filtered_GO_data <- filter_exp(subset_GO_data, dispersion_threshold = 10, threads = 2)
#'
filter_exp <- function(cell_gene_matrix, dispersion_threshold = 10, threads = 2){
  require(doParallel)
  exp_dispersion <- dispersion_ratio(t(cell_gene_matrix), threads = threads)
  print("Dispersion Summary")
  print(summary(exp_dispersion))
  print("Dispersion Quantiles")
  print(quantile(exp_dispersion,0.95))
  print(quantile(exp_dispersion,0.9))
  print(quantile(exp_dispersion,0.85))

  retained_exp <- cell_gene_matrix[, which(exp_dispersion >= dispersion_threshold)]
  return(retained_exp)
  
}


# dispersion ratio
dispersion_ratio <- function(gene_cell_matrix, threads = 2){
  m_genes <- nrow(gene_cell_matrix)
  n_cells <- ncol(gene_cell_matrix)
  gene_cell_matrix <- gene_cell_matrix
  store_ratio <- c()
  require(doParallel)
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  resultLst <- foreach(i = 1: m_genes, .combine=c) %dopar%
  {
    #for(i in c(1:m_genes)){
    me1 <- sum(gene_cell_matrix[i, ]) / n_cells
    #print(me1)
    if(me1 > 0){
      ratio <- sum((gene_cell_matrix[i, ]-me1)^2) / (n_cells-1) * me1
      #print(ratio)
      #store_ratio <- c(store_ratio, ratio)
      return(ratio)
    }else{
      ratio <- 0
      #store_ratio <- c(store_ratio, ratio)
      return(ratio)
    }
    
    #ratio <- 
  }
  stopCluster(cl)
  #closeAllConnections()
  return(resultLst)
}

# Additional functions to convert gene names


#' Gene name conversion: MGI -> ENTREZ
#' 
#' Provides gene name conversion for ease of processing.
#'@param mgi_names The list of gene names in MGI. (Mouse)
#'@param gene_table mm_IDs (loaded from preprocessed table)
#'@section Preprocessing:
#'@export
#'@examples
#'mgi_to_entrez(mgi_list, gene_table = mm_IDs)

mgi_to_entrez <- function(mgi_names, gene_table){
  get_id <- match(mgi_names, gene_table$mgi_symbol)
  get_id <- get_id[-which(is.na(get_id))]
  
  entrez_id <- gene_table$entrez[get_id]
  return(entrez_id)
}

#' Gene name conversion: ENTREZ -> MGI
#' 
#' Provides gene name conversion for ease of processing.
#'@param entrez_names The list of gene names in ENTREZ. (Mouse)
#'@param gene_table mm_IDs (loaded from preprocessed table)
#'@section Preprocessing
#'@export
#'@examples
#'entrez_to_mgi(entrez_list, gene_table = mm_IDs)

entrez_to_mgi <- function(entrez_names, gene_table){
  get_id <- match(entrez_names, gene_table$entrez)
  if(is.na(sum(get_id))){
    get_id <- get_id[-which(is.na(get_id))]
  }
  
  
  mgi_id <- gene_table$mgi_symbol[get_id]
  return(mgi_id)
}

#' Gene name conversion: ENSEMBL -> MGI
#' 
#' Provides gene name conversion for ease of processing.
#'@param ensembl_names The list of gene names in ENSEMBL. (Mouse)
#'@param gene_table mm_IDs (loaded from preprocessed table)
#'@section Preprocessing:
#'@export
#'@examples
#'ensembl_to_mgi(ensembl_names, gene_table = mm_IDs)

ensembl_to_mgi <- function(ensembl_names, gene_table){
  get_id <- match(ensembl_names, gene_table$ensembl_gene_id)
  get_id <- get_id[-which(is.na(get_id))]
  
  entrez_id <- gene_table$mgi_symbol[get_id]
  return(entrez_id)
}
