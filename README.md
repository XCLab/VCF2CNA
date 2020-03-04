# VCF2CNA

This repository is comprises a suite of tools to aid with copy number analysis and tumor purity estimation. These tools are designed for use with data in the VCF format, Mutation Annotation format (MAF), or Bambino format (High20).  The types of data accepted are 1) Whole Exome Sequncing Data, 2) Paired Tumor/Normal Whole Genome Sequencing Data,  3) RNA-Seq Analysis - to determine allelic expression.

https://github.com/XCLab/VCF2CNA/blob/master/Figure1_flowchart.pdf

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

VCF2CNA relies on germline and somatic SNVs.  Do not filter germline SNVs from the input data.
VCF2CNA expects paired tumor/normal read count data for heterozygous sites from all chromosomes.

To run the application use the execute.py python script

## Mandatory Arguments
|Argument       | Description               |
|---------------|---------------------------|
|FileName       | File to Analyze           |
|FilePath       | Path to FileName          |
|OutputDirectory| Path to store results     |

### Simple Commandline

```
python execute.py [FILENAME][FILEPATH][OUTPUT_DIRECTORY]
```

### Run Example Data

```
Navigate to test data directory

Extract individual xz files and combine them into one file
sh extract_example.sh

Run example data with default parameters
python execute.py test.high_20.out ~/test_data ~/example_output

Run example data with user specified Diploid Chromosome
python execute.py test.high_20.out ~/test_data ~/example_output -DiploidChr=13

Run example data with user specified Region inside Chromosome
python execute.py test.high_20.out ~/test_data ~/example_output -DiploidChr=13 -StartChr=19020701 -EndChr=26516137 
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

### Output Files

|Filename                                 | Description                      |
|-----------------------------------------|----------------------------------|
|Samplename_CONSERTING_Mapability_100.txt | Copy Number Segments             | 
|Samplename_LOH_RegTree.txt               | Allelic Imbalance by segment     |
|Samplename.jpg                           | Visualization of Copy Number     |
|purity_Samplename.txt                    | Overall Tumor Purity Estimation  |
|purity_Samplename.csv                    | Tumor Purity by Segment and Type | 
|purity_Samplename.png                    | Visualization of Purity Segments | 

### Samplename_CONSERTING_Mapability_100.txt

|Column      | Description                                                                            | 
|------------|----------------------------------------------------------------------------------------|
|chrom       | Chromosome (X=23 and Y=24)                                                             | 
|loc.start   | Starting coordinate of the segments                                                    |
|loc.end     | Ending coordinates of the segments                                                     |
|num.mark    | Number of windows in the segment (Sequencing gaps and low mapability scores excluded)  |
|seg.mean    | The estimated GC corrected difference signal (2 copy gain will have a seg.mean of 1)   | 
|GMEAN       | The mean coverage in the germline sample (A value of 1 represents diploid)             | 
|DMEAN       | The mean coverage in the tumor sample                                                  | 
|LogRatio    | The log2ratio of the coverage signal between tumor and germline samples                | 



## Authors

* **Daniel K. Putnam, Xiaotu Ma, Stephen V. Rice, Yu Liu, Scott Newman, Jinghui Zhang, Xiang Chen** 
