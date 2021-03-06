# acorde: isoform co-usage networks from single-cell RNA-seq data

### Introduction
The **acorde** R package contains an implementation of the pipeline showcased in
Arzalluz-Luque et al. 2021 [[1]](#1). *acorde* is an end-to-end pipeline designed 
for the study of isoform co-usage networks using single-cell RNA-seq data (scRNA-seq). 

[![DOI](https://zenodo.org/badge/365198830.svg)](https://zenodo.org/badge/latestdoi/365198830)

The pipeline includes three basic analysis blocks:

1. Single-cell isoform quantification and Differential Expression (DE) filtering. 

2. Detection of isoform co-expression relationships using percentile correlations
and semi-automated clustering.

3. Differential and co-Differential Isoform Usage analysis. To couple these 
analysis with a biologically interpretable readout, we incorporate functional 
annotations onto isoform models, and use 
[tappAS](https://github.com/ConesaLab/tappAS) for functional analysis.

Since both the long read-transcriptome definition/quantification procedure and the 
functional analyses in [[1]](#1) are based on external tools, the present R 
package does **not** include them. 

To reproduce our pipeline's **long read transcriptome building** strategy, 
please refer to our manuscript's dedicated Supplementary Note [[1]](#1).

Regarding **functional analyses**, we have included a section in the vignette 
providing instructions to perform them, including:

- How to obtain a functionally-annotated transcriptome using [isoAnnotLite](https://isoannot.tappas.org/isoannot-lite/).
- How to generate input files that are compatible with the [tappAS](https://github.com/ConesaLab/tappAS) application.

*acorde* contains the necessary functions and documentation to obtain 
a set of DIU and co-DIU genes using an single-cell, isoform-level expression 
matrix as input. In addition, we provide instructions to reproduce the figures 
and additional analyses included in Arzalluz-Luque et al. [[1]](#1).
The isoform expression matrix employed during the study is provided as internal 
data in the package.

![](images/acorde_pipeline-small.png)


### Installation
The *acorde* R package and all the necessary dependencies  can be installed 
from GitHub using `devtools`:

```
install.packages("devtools")
devtools::install_github("ConesaLab/acorde")
```

To access vignettes, you will need to force building with
`devtools::install_github(build_vignettes = TRUE)`. Please note that this will
also install all suggested packages required for vignette build and might 
increase install time. Alternatively, an HTML version of the vignette is
available under the [vignettes](https://github.com/ConesaLab/acorde/tree/master/vignettes)
folder.


### Getting started

In order to use acorde, you will need the following items:

- An isoform-level quantification matrix.
- Isoform-to-gene relationships and cell-level identity labels (i.e. cell type
or state).
- A functionally-annotated transcriptome (provided that you wish to perform
downstream functional analyses).

Please note that an example dataset (the one used in the acorde manuscript [[1]](#1))
including all required objects is provided along with the package. Details to
load and use them can be found in the package's vignette.

Importantly, acorde can be applied on both standard reference and long 
read-defined transcriptomes -the only limitation concerns the annotation of 
functional features for those isoforms.

For our study, we used IsoAnnotLite (v2.6) to transfer functional labels 
from tappAS??? pre-annotated mouse RefSeq 78 reference as well as from the mouse
neural transcriptome used in the tappAS publication, which was possible thanks 
to the fact that we used a compatible reference for long-read processing (i.e. 
an updated version of the mouse RefSeq transcriptome). Therefore, users will
need to take this into consideration and rely on pre-annotated references and 
isoAnnotLite until a *de novo* annotation tool is available (currently 
a [work in progress](https://github.com/ConesaLab/tappAS/issues/12)).

In the meantime, all available references can be viewed and downloaded [here](https://app.tappas.org/resources/downloads/gffs/). To use the mouse brain
transcriptome in our manuscript, download the 
`Mus_Musculus_GRCm38.p6_PacBioENCODE_RefSeq108.gff3` file.


### Contact
If you encounter a problem, please 
[open an issue](https://github.com/ConesaLab/acorde/issues) 
via GitHub or send an email to angeles.arzalluz [at] gmail.com.

  
### References
If you use *acorde* in your research, please cite the original publication:

<a id="1">[1]</a>
Arzalluz-Luque, A., Salguero, P., Tarazona, S. et al. acorde unravels functionally 
interpretable networks of isoform co-usage from single cell data. 
*Nat Commun* 13, 1828 (2022). [https://doi.org/10.1038/s41467-022-29497-w](https://doi.org/10.1038/s41467-022-29497-w)
