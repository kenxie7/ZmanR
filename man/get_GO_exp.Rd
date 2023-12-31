% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_gene_expression.R
\name{get_GO_exp}
\alias{get_GO_exp}
\title{Preprocess Gene Expression}
\usage{
get_GO_exp(
  gene_cell_matrix,
  gene_type = "gene_names",
  organism = "mouse",
  takelog = T
)
}
\arguments{
\item{gene_cell_matrix}{A genes (row) by cells (col) matrix.}

\item{gene_type}{The type of gene names. Default is "gene_names" - ie "HOPX", alternatively "entrez" is provided.}

\item{organism}{The type of species the dataset belongs to: defaults to "mouse". "human" option is available.}

\item{takelog}{Defaults to TRUE. Whether to perform log2(exprs + 1)}
}
\description{
This function reduces the gene expression matrix to selected Gene Ontology genes.
}
\section{Preprocessing}{

}

\examples{
subset_GO_data <- get_GO_exp(data, gene_type = "gene_names", organism = "human", takelog = F)
filtered_GO_data <- filter_exp(subset_GO_data, dispersion_threshold = 10, threads = 2)

}
\keyword{get_GO_exp}
