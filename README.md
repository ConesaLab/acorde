# acorde: isoform co-usage networks from single-cell RNA-seq data

### Introduction
The **acorde** R package contains an implementation of the pipeline showcased in
Arzalluz-Luque et al. 2021 [[1]](#1). *acorde* is an end-to-end pipeline designed 
for the study of isoform co-usage networks using single-cell RNA-seq data (scRNA-seq). 

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


### Contact
Note that *acorde* is currently under development. If you encounter a 
problem, please [open an issue](https://github.com/ConesaLab/acorde/issues) 
via GitHub or send an email to angeles.arzalluz [at] gmail.com.

  
### References
If you use *acorde* in your research, please cite:

<a id="1">[1]</a>
Angeles Arzalluz-Luque, Pedro Salguero, Sonia Tarazona, Ana Conesa:
*Acorde*: unraveling functionally-interpretable networks of isoform co-usage 
from single cell data. *bioRxiv* 2021.05.07.441841; 
doi: https://doi.org/10.1101/2021.05.07.441841
