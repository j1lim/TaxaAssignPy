
# TaxaAssignPy

## Introduction

The TaxaAssignPy is a convenient and rapid pipeline for taxonomy assignments based on rRNA sequence. The TaxaAssignPy automatically defines Operational Taxonomic Units (OTUs) with raw sequence data and classify the OTUs based on the naive Bayesian classifier method. For both processes, TaxaAssignPy utilizes two other published tools, [VSEARCH](https://github.com/torognes/vsearch/) and [RDP Classifier](http://rdp.cme.msu.edu/classifier/classifier.jsp).

Another goal of this project is to create a pipeline for the person who is not familiar with computation by using cloud services. For that reason, two cloud services of Google, [Colaboratory](https://colab.research.google.com/notebooks/intro.ipynb) and [Drive](https://www.google.com/drive/), were recommended. Google Colaboratory service was chosen to run the pipeline and Google drive was chosen to import input data and store output results.

## Install

For the users who have own Linux system, they just need to install two other preprocess tools, VSEARCH and RDP Classifier. The TaxaAssignPy is a python library so pipeline itself does not need extra installation steps.

However, for the users who are going to utilize cloud service, they can run the pipeline just with run the 'ipynb' on the Google Colaboratory.

VSEARCH can be installed easily by using 'apt',
```
apt install vsearch -y
```
or 'anaconda',
```
conda install -c bioconda vsearch -y
```
or install manually with [VSEARCH github](https://github.com/torognes/vsearch/).


RDP Classifier can be installed by using 'anaconda',
```
conda install -c bioconda rdp_classifier -y
```
or install manually with [RDP Classifier github](https://github.com/rdpstaff/classifier).

## Usage
In the TaxaAssignPy pipeline, there are seven different modules in detail. If you want to see the results of each step, you can use an individual module as below.

```
import TaxaAssign

path = "directory of the folder which contains sequence data"
TaxaAssign.taxaassign( path )
TaxaAssign.taxaassign.fastq_merge()
TaxaAssign.taxaassign.fastq_filtering()
TaxaAssign.taxaassign.OTU_defining()
TaxaAssign.taxaassign.taxonomy_assigning()
TaxaAssign.taxaassign.analysis_statistic()
TaxaAssign.taxaassign.result_summary()
```

or if you want to run all the processes at once, the module as below can be used.

```
impot TaxaAssign

path = "directory of the folder which contains sequence data"
TaxaAssign.autorun_all( path )
```

## References

* Blaxter, M.; Mann, J.; Chapman, T.; Thomas, F.; Whitton, C.; Floyd, R.; Abebe, E. (2005)
**Defining operational taxonomic units using DNA barcode data** *Philos Trans R Soc Lond B Biol Sci*. 360 (1462): 1935–43. DOI: [10.1098/rstb.2005.1725](https://royalsocietypublishing.org/doi/10.1098/rstb.2005.1725)

* Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016)
**VSEARCH: a versatile open source tool for metagenomics.** *PeerJ* 4:e2584. DOI: [10.7717/peerj.2584](https://peerj.com/articles/2584/)

* Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. (2007)
**Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy.** *Appl Environ Microbiol.* 73(16):5261-7. DOI: [10.1128/AEM.00062-07](https://aem.asm.org/content/73/16/5261.short)

