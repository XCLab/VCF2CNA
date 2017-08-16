# Proto-Oncogene / Tumor Supressor Gene Finder

This repository is comprises a suite of tools to aid with copy number analysis and tumor purity estimation. These tools are designed for use with data in the VCF format, Mutation Annotation format (MAF), or Bambino format (High20).  The types of data accepted are 1) Whole Exome Sequncing Data, 2) Paired Tumor/Normal Whole Genome Sequencing Data, 3) Germline (unpaired analysis) 4) RNA-Seq Analysis - to determine allelic expression.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Linux Operating System
Perl version 5.10.1 or greater
R    version 3.0.1 or greater
R    Packages: Tree, gplots, RColorBrewer, ggplot2, reshape2, grid
```

### Installing

Download the repository and place it in desired directory
The repository will be in a folder called VCF2CNA

From the VCF2CNA Folder run the extraction command

```
bash extract.sh
```

### Running the application

To run the application use the run.sh bash script

```
bash run.sh [FILENAME][FILE_DIR][SEQ_TYPE (EXOME or WHOLEGENOME)
bash run.sh TEST.high_20.out test EXOME
```

### Additional Information

## Authors

* **Daniel K. Putnam, Xiaotu Ma, Stephen V. Rice, Yu Liu, Jinghui Zhang, Xiang Chen** 
