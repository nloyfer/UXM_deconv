# UXM fragment-level deconvolution algorithm
UXM is a computational fragment-level reference-based deconvolution algorithm for DNA methylation sequencing data.
It constructs a reference atlas where the percentages of unmethylated fragments is computed for every marker (row) in each cell type (column).
A non-negative least squares (NNLS) algorithm is then used to fit an input sample, and estimate its relative contributions.

<p align='center'>
    <img src="docs/img/Atlas.U25.l4.png" width="600" height="400" />
</p>
<p align='center'>
    <em>Visualization of the reference atlas published in Loyfer et al.</em>
</p>


## Quick start
### Installation
- Make sure [`wgbstools`](https://github.com/nloyfer/wgbs_tools#installation) is installed and found in `$PATH`.<br>
- Clone UXM. 
- Optionally, add the "uxm" main scrip to your `$PATH`.
```bash
git clone https://github.com/nloyfer/UXM_deconv.git
cd UXM_deconv/
# Optionally, add to PATH (or link "./uxm" to some directory in $PATH)
export PATH=${PATH}:$PWD
```

#### Dependencies
- python 3+
- pandas, numpy, scipy
- wgbstools

### Usage examples
#### Running deconvolution
```bash
$ uxm deconv data/*pat.gz -o output.csv
dumped atlas to output.csv
```

#### Plotting deconvolution results
```bash
# Plot a stacked bar plot
$ uxm plot output.csv -o output.pdf
dumped figure to output.pdf
```

This project is developed by Netanel Loyfer in [Prof. Tommy Kaplan's lab](https://www.cs.huji.ac.il/~tommy/) at the Hebrew University, Jerusalem, Israel.<br>
If you are using the UXM deconvolution in a paper, please cite:

[Loyfer, N. *et al.* (2022)](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1) ‘A human DNA methylation atlas reveals principles of cell type-specific methylation and identifies thousands of cell type-specific regulatory elements’, *bioRxiv.* doi:10.1101/2022.01.24.477547.

### See the [tutorial](tutorial) for more information 
