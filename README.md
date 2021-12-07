# Single Sample Somatic Variant Calling Pipeline using Snakemake
## Overview
This Snakemake pipeline provides an ensemble method for calling somatic SNVs in a single sample (i.e., w/out a normal match) 
consisting of either RNA-seq or DNA-seq reads. The core workflow is shown in the example below.



## Dependencies
* Python 3
* Anaconda (or Miniconda)
* Snakemake
* Pandas

## Installation
The pipeline currently uses conda environments to automatically deploy required software for each process (i.e., mapping, variant calling, etc.). 
Ensure `conda` is available for use on your system before proceeding. To create a new environment with the above dependencies, run the following:
```
# specify any environment name after -n flag
$ conda create -n snakemake python=3 pandas snakemake
# activate the previously created environment
$ conda activate snakemake
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

## Getting Started
### Sample Data Preparation
The pipeline is set up to read input fastq files from a user specified directory in [data](data/). The idea is to create a new directory 
for each project and add all associated sample fastq files to it. Once fastq files are added to the created directory, the specific sample 
input details for the project need to be added to [this table](config/samples.tsv). More info [here](config/README.md).

Additionally, add the reference assembly file and, if running in "RNA" mode, a gene annotation file to a user specified directory 
in the [resources](resources/) folder. Again, create a directory that represents the project in question (i.e., if working with mouse samples, 
a simple directory named "mouse_GRCm39" would hold the file(s)). Once the reference indexes for a given organism are initially created, they 
can be reused with future executions of the pipeline, without the need to recreate them (Snakemake will auto detect). 

### General Pipeline Configuration
Before running the pipeline, configure general settings in [this config file](config/config.yaml). Details for parameters are provided in 
the file. 

### Cluster Configuration
Currently, the pipeline is built for submission to a cluster using the slurm workload manager (i.e., Canopus). To configure default resources 
used by all rules such as memory, nodes, and time, modify entries [here](config/slurm/config.yaml). For advanced configuration, specifying resources 
at rule definition will override the default parameters in the cluster configuration file.

### Example Pipeline Run
Given all input files are organized according to the following structure (defined above):

```
├── README.md
├── workflow
├── config
├── resources
│   ├── organism1
│   │   └── reference.fa
│   └── organism2
├── data
│   ├── project1
│   │   ├── sample1_R1.fq
│   │   ├── sample1_R2.fq
│   │   ├── sample2_R1.fq
│   │   └── sample2_R2.fq
│   └── project2
└── results
```

...and [samples.tsv](config/samples.tsv) is populated in the following manner:

| Sample_name | R1 | R2 |
| :---: | :---: | :---: |
| Sample1 | sample1_R1.fq | sample1_R1.fq |
| Sample2 | sample2_R1.fq | sample2_R1.fq |

...and the workflow is properly configured for the current run, one can submit a pipeline instance to Canopus (slurm) with the following commands:

```
# move to directory where Snakefile is present
cd workflow/
# submit jobs to slurm scheduler
snakemake --profile ../config/slurm --config reference=organism1 data=project1 mode=DNA
```

Snakemake automatically generates, submits, and monitors job submissions to the slurm scheduler. Log files for all rules will 
be output in a newly created "logs/" folder within the [workflow](workflow/) directory, named according to the sample name and 
slurm job ID. All pipeline output for each sample will go to the [results](results/) folder under a newly created directory 
named after the sample. Final variant callsets will be located in the *final_calls* directory.

## Current Limitations
