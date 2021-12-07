# Single Sample Somatic Variant Calling Pipeline using Snakemake
## Overview
This Snakemake pipeline provides an ensemble method for calling somatic SNVs in a single sample (i.e., w/out a normal match) 
consisting of either RNA-seq or DNA-seq reads. 

*The pipeline is currently only configured for use with the Slurm scheduler (i.e., usage on Canopus only).*

## Dependencies
* Python 3
* Anaconda (or Miniconda)
* Snakemake
* Pandas

## Getting Started
### Installation
The pipeline currently uses conda environments to automatically deploy required software for each process (i.e., mapping, variant calling, etc.). 
Ensure `conda` is available for use on your system before proceeding. To create a new environment with the above dependencies, run the following:
```
# specify any environment name after -n flag
$ conda create -n snakemake python=3 pandas snakemake
```

Clone the pipeline in the current working directory using git:
```
$ git clone https://github.com/bwinnacott/SNP_pipeline.git
```

**Note**: At the moment, Canopus is not configured for use with git. SSH onto Rustang and clone the repository in the desired working directory. 
The file transfer system in place will permit access through Canopus.

To verify proper configuration, run the following:

```
# output should provide a laundry list of Snakemake parameters
$ snakemake -h
```

### Sample Data Preparation
The directory structure of the pipeline is set up to read input fastq files from a specified directory in the [Data Folder](data/) folder

### Configuration
Before running the pipeline,