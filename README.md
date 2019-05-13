# VCF2CNA

This repository is comprises a suite of tools to aid with copy number analysis and tumor purity estimation. These tools are designed for use with data in the VCF format, Mutation Annotation format (MAF), or Bambino format (High20).  The types of data accepted are 1) Whole Exome Sequncing Data, 2) Paired Tumor/Normal Whole Genome Sequencing Data,  3) RNA-Seq Analysis - to determine allelic expression.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Linux Operating System
Perl version 5.10.1 or greater
R    version 3.4.0 
R    Packages: Tree, gplots, RColorBrewer, ggplot2, reshape2, grid
Python version 2.7.12
```

### Installing

Download the repository and place it in desired directory
The repository will be in a folder called VCF2CNA

From the VCF2CNA Folder run the extraction command

```
bash extract.sh
```

### Running the application

To run the application use the execute.py python script

## Mandatory Arguments
|Argument       | Description               |
|---------------|---------------------------|
|FileName       | File to Analyze           |
|FilePath       | Path to FileName          |
|OutputDirectory| Path to store results     |


```
python execute.py [FILENAME][FILEPATH][OUTPUT_DIRECTORY]
```

### Additional Parameters

|Argument       | Description                      | Default|
|---------------|----------------------------------|--------|
|VcfOrder       | Tumor/Normal TN or NT            | TN     |
|Diploid Chr    | Chromosome used for normalization| -1     |      
|StartChr       | Start location on Diploid Chr    | -1     |
|EndChr         | End location on Diploid Chr      | -1     |
|Median         | Median Normal Coverage           | -1     |
|MinSf          | Minimum Scale Factor             | 0.50   | 
|MaxSf          | Maximum Scale Factor             | 1.50   | 
|XminSf         | Minimum X Scale Factor           | 0.25   |
|XmaxSf         | Maximum X Scale Factor           | 1.50   |
|R_Loc          | Path to R Program                | ""     |
|R_ScriptLoc    | Path to RScript Program          | ""     |

## Authors

* **Daniel K. Putnam, Xiaotu Ma, Stephen V. Rice, Yu Liu, Scott Newman, Jinghui Zhang, Xiang Chen** 
