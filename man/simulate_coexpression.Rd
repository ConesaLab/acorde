% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coexpr_simulation.R
\name{simulate_coexpression}
\alias{simulate_coexpression}
\title{Simulate co-expression patterns on an already-simulated scRNA-seq matrix}
\usage{
simulate_coexpression(sim_data, feature_no, patterns, cluster_size)
}
\arguments{
\item{sim_data}{A \code{SingleCellExperiment} object containing an already-simulated
scRNA-seq matrix.}

\item{feature_no}{An integer indicating the number of high expression ("top")
and low expression ("bottom") features to be selected for co-expression
simulation. Note that \code{feature_no*2} features will be used in total.}

\item{patterns}{A \code{data.frame} or \code{tibble} containing cell types as
columns (ordered as in the \code{colData} slot in \code{sim_data}) and
co-expression clusters as rows. For each co-expression cluster,
a logical vector indicating the desired expression pattern must be
provided, in a row-wise manner. Insert \code{TRUE} if
high expression in that cell type is desired, \code{FALSE} if the opposite.}

\item{cluster_size}{An integer indicating the number of features to include
per cluster.}
}
\value{
A \code{list}, containing two objects:

\enumerate{

\item \code{sim_matrix}: a \code{tibble} containing the same number of cells
as in \code{sim_data} in the columns and \code{feature_no*2} in the rows.
Feature IDs are defined in the \code{feature} column

\item \code{sim_clusters}: a \code{list} with as many elements as simulated
clusters, where each element contains all feature IDs that were simulated to
follow the same co-expression pattern (that is, the clusters).

}
}
\description{
Using a existing simulated scRNA-seq matrix, this function
creates co-expression relationships between the features, following the
cell type-specific patterns of high/low expression supplied by the user.
In particular, the output of the \code{SymSim}
simulator is expected as input (see \code{SymSim} package documentation for details).
}
