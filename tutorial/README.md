# UXM deconvolution tutorial

## Installation and configuration
First make sure [`wgbstools`](https://github.com/nloyfer/wgbs_tools#installation) is installed and found in `$PATH`.<br>
Second, "install" (clone) UXM. Optionally, add the "uxm" main scrip to your `$PATH`.
```bash
git clone https://github.com/nloyfer/UXM_deconv.git
cd UXM_deconv/
# Optionally, add to PATH (or link "./uxm" to some directory in $PATH)
export PATH=${PATH}:$PWD
```

## Tutorial
### Data and regions
#### Input samples
For this short tutorial, we will use the following publicly available samples from the [Roadmap Epigenomic Project](https://www.nature.com/articles/nature14248). The fastq files were downloaded from [GEO](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE16256), and mapped to hg19 using [bwa-meth](https://github.com/brentp/bwa-meth). `pat` files were generated using `wgbstools bam2pat`, and finally, the large `pat` files were sliced (using `wgbstools view -L` and `wgbstools index`) into small pat files, containing only reads mapped to regions listed in the reference atlas. 


| SRX  | Tissue  |  Donor |
|---|---|---|
| [SRX388734](https://www.ncbi.nlm.nih.gov/sra?term=SRX388734) |  Lung cells          | STL001
| [SRX175350](https://www.ncbi.nlm.nih.gov/sra?term=SRX175350) |  Lung cells          | STL002
| [SRX388743](https://www.ncbi.nlm.nih.gov/sra?term=SRX388743) |  Pancreas cells      | STL002
| [SRX175354](https://www.ncbi.nlm.nih.gov/sra?term=SRX175354) |  Pancreas cells      | STL003
| [SRX175348](https://www.ncbi.nlm.nih.gov/sra?term=SRX175348) |  Sigmoid colon cells | STL001
| [SRX190161](https://www.ncbi.nlm.nih.gov/sra?term=SRX190161) |  Sigmoid colon cells | STL003
| [SRX213280](https://www.ncbi.nlm.nih.gov/sra?term=SRX213280) |  Liver cells         | STL011


#### Reference atlas
The reference atlas used for deconvolution is a table where the rows are regions (markers) and the columns are reference samples/cell types. We provide here the a reference atlas [supplemental/Atlas.U25.l4.hg19.tsv](../supplemental/Atlas.U25.l4.hg19.tsv), we published in [Loyfer *et al.* (2022)](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1). 
The markers (rows) are the top 25 specifically unmethylated blocks for each of ~40 healthy human cell types. The values are the "proportion of unmethylated reads", as computed from the pat files. <br>
The `uxm` tool allows you to build your own reference atlas, choose which cell types and samples to include, which markers and more. (see `uxm build --help`. Documentation for this feature will be added soon).


```bash
$ cd tutorial
# These are the samples we will be deconvolving
$ ls -1 data/*pat.gz
data/Liver_STL011.pat.gz
data/Lung_STL001.pat.gz
data/Lung_STL002.pat.gz
data/Pancreas_STL002.pat.gz
data/Pancreas_STL003.pat.gz
data/Sigmoid_Colon_STL001.pat.gz
data/Sigmoid_Colon_STL003.pat.gz
# This is the reference atlas, containing 900 markers
$ wc -l ../supplemental/Atlas.U25.l4.hg19.tsv 
901 ../supplemental/Atlas.U25.l4.hg19.tsv
```

#### Running deconvolution
```bash
$ uxm deconv data/*pat.gz -o output.csv
dumped atlas to output.csv
```

#### Plotting deconvolution results
```bash
# Plot stacked bar plots for the whole table
$ uxm plot output.csv -o output.pdf
dumped figure to output.pdf

# Plot bar plots for only one cell type (e.g., endothel)
$ uxm plot output.csv --stub Endothel
dumped figure to output.pdf
```

##### To be continued (last updated on October 10, 2022)
coming soon:  <br>
Multiple available atlases <br>
--ignore/include flags <br>
uxm build <br>

