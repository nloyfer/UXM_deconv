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
For this short tutorial, we will use the following publicly available samples from the [Roadmap Epigenomic Project](https://www.nature.com/articles/nature14248). The fastq files were downloaded from [GEO](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE16256), and mapped to hg19 using [bwa-meth](https://github.com/brentp/bwa-meth). `pat` files were generated using `wgbstools bam2pat`, and finally, the large `pat` files were sliced (using `wgbstools view -L` and `wgbstools index`) into small pat files, containing only reads mapped to regions listed in the bed file `Markers.U25.bed` (see below). 

#### Markers
We will be using a set of differentialy methylated regions / markers found and published in [Loyfer *et al.* (2022)](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1) (Supplemental table S4). They are the top (up to) 25 specifically unmethylated blocks for each of 39 healthy human cell types. They are listed here in file `Markers.U25.bed`.


| SRX  | Tissue  |  Donor |
|---|---|---|
| [SRX388734](https://www.ncbi.nlm.nih.gov/sra?term=SRX388734) |  Lung cells          | STL001
| [SRX175350](https://www.ncbi.nlm.nih.gov/sra?term=SRX175350) |  Lung cells          | STL002
| [SRX388743](https://www.ncbi.nlm.nih.gov/sra?term=SRX388743) |  Pancreas cells      | STL002
| [SRX175354](https://www.ncbi.nlm.nih.gov/sra?term=SRX175354) |  Pancreas cells      | STL003
| [SRX175348](https://www.ncbi.nlm.nih.gov/sra?term=SRX175348) |  Sigmoid colon cells | STL001
| [SRX190161](https://www.ncbi.nlm.nih.gov/sra?term=SRX190161) |  Sigmoid colon cells | STL003
| [SRX213280](https://www.ncbi.nlm.nih.gov/sra?term=SRX213280) |  Liver cells         | STL011

```bash
$ cd tutorial
$ ls -1 data/*pat.gz
data/Liver_STL011.pat.gz
data/Lung_STL001.pat.gz
data/Lung_STL002.pat.gz
data/Pancreas_STL002.pat.gz
data/Pancreas_STL003.pat.gz
data/Sigmoid_Colon_STL001.pat.gz
data/Sigmoid_Colon_STL003.pat.gz
```

##### To be continued (last updated on October 7, 2022)
