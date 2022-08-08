# UXM fragment-level deconvolution algorithm
UXM is a computational fragment-level reference-based deconvolution algorithm for DNA methylation sequencing data.
It constructs a reference atlas where the percentages of unmethylated fragments is computed for every marker (row) in each cell type (column).
A non-negative least squares (NNLS) algorithm is then used to fit an input sample, and estimate its relative contributions.

<p align='center'>
    <img src="docs/img/Atlas.U25.l4.png" width="600" height="400" />
</p>
<p align='center'>
    <em>Visualization of the reference atlas published in [Loyfer et al.](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1)</em>
</p>

<!--![alt text](docs/img/Atlas.U25.l4.png "U25 atlas")-->
This project is developed by Netanel Loyfer in [Prof. Tommy Kaplan's lab](https://www.cs.huji.ac.il/~tommy/) at the Hebrew University, Jerusalem, Israel.<br>
If you are using the UXM deconvolution in a paper, please cite:

[Loyfer, N. *et al.* (2022)](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1) ‘A human DNA methylation atlas reveals principles of cell type-specific methylation and identifies thousands of cell type-specific regulatory elements’, *bioRxiv.* doi:10.1101/2022.01.24.477547.

