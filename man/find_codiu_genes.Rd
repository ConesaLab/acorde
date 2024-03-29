% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coDIU_genes.R
\name{find_codiu_genes}
\alias{find_codiu_genes}
\title{Detect genes with co-Differential Isoform Usage}
\usage{
find_codiu_genes(cluster_list, gene_tr_table, parallel = TRUE, t = 4)
}
\arguments{
\item{cluster_list}{A list of character vectors containing isoform IDs.
Each element of the list represents a cluster of isoforms.}

\item{gene_tr_table}{A data.frame or tibble object containing two columns
named \code{transcript_id} and \code{gene_id}, indicating gene-isoform
correspondence.}

\item{parallel}{A logical. When \code{TRUE}, parallelization is enabled.
The \code{\link{future_map_lgl}} function in the \code{furrr} is used.}

\item{t}{An integer indicating the number of threads to be used for
parallelization. This will be passed to the \code{\link{plan}} function from
the \code{future} package via the \code{workers} argument.}
}
\value{
A matrix containing two rows and as many columns as potentially
co-DIU gene pairs detected, that is, genes co-expressing isoforms across
two or more clusters. Genes will be represented by the IDs provided in
\code{{gene_tr_table}}.
}
\description{
This function reports pairs of genes that present co-expressed
isoforms given a list of previously-detected isoform clusters.
The aim is to enable co-Differential Isoform Usage (co-DIU) analysis on the
returned gene pairs, which will be candidates for co-DIU.
}
\details{
We define coordinated splicing patterns as a situation where
post-transcriptional regulation, defined by isoform expression,
can be detected independently of transcriptional regulation,
i.e. gene-level expression. To detect splicing coordination, we defined
co-Differential Isoform Usage (co-DIU) as a pattern where a group of genes
shows co-expression of their isoforms, but no co-expression can be detected
when only gene expression is considered. In the context of our pipeline,
a set of potentially co-DIU genes will have at least two of their isoforms
assigned to the same clusters, therefore showing detectable isoform-level
co-expression, and suggesting coordinated splicing regulation in that group
of genes (see Arzalluz-Luque et al. 2021).
}
\references{
\insertRef{Arzalluz-Luque2021}{acorde}
}
