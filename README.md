# UXM fragment-level deconvolution algorithm
UXM is a computational fragment-level reference-based deconvolution algorithm for DNA methylation sequencing data.
It generates a reference atlas where the percentages of unmethylated fragments is computed for every marker (row) in each cell type (column).
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

[Loyfer, N. *et al.* (2022)] ‘A human DNA methylation atlas reveals principles of cell type-specific methylation and identifies thousands of cell type-specific regulatory elements’, *bioRxiv.* doi:10.1101/2022.01.24.477547.(https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1)



It converts data from standard formats (e.g., bam, bed) into tailored compact yet useful and intuitive formats ([pat](docs/pat_format.md), [beta](docs/beta_format.md)).
These can be visualized in terminal, or analyzed in different ways - subsample, merge, slice, mix, segment and more.
