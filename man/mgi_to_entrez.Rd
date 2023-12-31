% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_gene_expression.R
\name{mgi_to_entrez}
\alias{mgi_to_entrez}
\title{Gene name conversion: MGI -> ENTREZ}
\usage{
mgi_to_entrez(mgi_names, gene_table)
}
\arguments{
\item{mgi_names}{The list of gene names in MGI. (Mouse)}

\item{gene_table}{mm_IDs (loaded from preprocessed table)}
}
\description{
Provides gene name conversion for ease of processing.
}
\section{Preprocessing}{

}

\examples{
mgi_to_entrez(mgi_list, gene_table = mm_IDs)
}
