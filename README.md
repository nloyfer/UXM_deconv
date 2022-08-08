# UXM fragment-level deconvolution algorithm
UXM is a computational fragment-level reference-based deconvolution algorithm for DNA methylation sequencing data.
It generates a reference atlas where the percentages of unmethylated fragments is computed for every marker (row) in each cell type (column).
A non-negative least squares (NNLS) algorithm is then used to fit an input sample, and estimate its relative contributions.


It converts data from standard formats (e.g., bam, bed) into tailored compact yet useful and intuitive formats ([pat](docs/pat_format.md), [beta](docs/beta_format.md)).
These can be visualized in terminal, or analyzed in different ways - subsample, merge, slice, mix, segment and more.
