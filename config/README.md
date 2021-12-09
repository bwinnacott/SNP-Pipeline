### Cluster Configuration
See the config file in *slurm* for details on specifying cluster submission parameters, as well as 
some general Snakemake parameters. For now, the pipeline is only configured to run on clusters using 
the slurm scheduling manager.

### Pipeline Configuration
See the config file for details on specifying parameters specific to individual tools used throughout 
the pipeline, as well as some general pipeline parameters. Default parameters should work for most cases; 
however, some need to be changed more frequently such as ploidy, number of vcf files to intersect, etc.

### Sample Configuration
Update the samples.tsv file before running an instance of the pipeline. Make sure the entries are tab separated 
after modifying. Columns are defined as:

* **Sample_name**: Name of a single sample; will be used as the directory name for storing pipeline output in results folder
* **R1**: File name of the forward fastq file; needs to be present in the data directory specified at the command line
* **R2**: File name of the reverse fastq file; needs to be present in the data directory specified at the command line

Single end data is accepted and the file name needs to be specified under R1 column. In this case, the R2 column may be left 
blank (or input with NA). 