# Performance Benchmarking Progress

## Introduction

Below is a summary of the methods and initial results for the evaluation of DNA- and RNA-seq pipeline performance on 
a mix of gold standard and simulated SNVs. The gold standard set of variants used in this analysis were derived from the 
Genome In A Bottle (GIAB) consortium<sup>1</sup>, which were then used to benchmark calls generated from the ensemble pipeline. 
Benchmarking in this scenario was performed using the ```vcfeval``` function provided by RTG Tools<sup>2</sup>, which enables 
more sophisticated variant comparisons between VCF files. For variant simulation studies, the tool Somatosim<sup>3</sup> was used 
to induce SNVs in experimental BAMs of Scen DOE152z at various VAF's across a range of coverages. Precision, recall, 
and accuracy (measure as F1) were used to quantify overall performance.

## Datasets

### Gold Standard -- DNA and RNA

GIAB's NA12878 (HG001) human gold standard reference dataset was selected for the evaluation of pipeline performance. One 
of the benefits of using this sample for performance evaluation is the availability of both DNA- and RNA-seq data.
The reference Fasta was obtained from the following location:

```https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz```

DNA Fastq files were downloaded from the following location:

```https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/```

Fastq files from two of the project folders in the above directory were combined into two final Fastq files (R1 and R2), 
for a total coverage of ~50x. High confidence regions used for benchmarking purposes were pulled from the following:

```https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed```

The reference callset:

```https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz```

RNA Fastq files:

```https://www.ncbi.nlm.nih.gov/sra/?term=SRR5665260```

For RNA variant calling evaluation, the same confident regions can be subset to only include the coding regions. Genome 
stratification files are used to accomplish this:

```https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz```

Lastly, the gene annotation file:

```https://ftp-trace.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz```

### Scen DOE152Z -- DNA and RNA

The Scen DOE152Z Fastq files were used as experimental datasets to validate the ensemble SNV pipeline, through SNV simulation. 
As the GIAB reference consists of human variants that are considered low hanging fruit, it is necessary to evaluate performance 
on samples more representative of those the team consistently works with. The following files are located on Canopus.

Reference Fasta file:

```/panfs/biopan04/scratch-krusec/Scen152Z_RNAseq/Hap1Genome/scenedesmus_obliquus_doe0152z_complete_genome.fasta```

DNA Fastq files were pulled from:

```/panfs/biopan03/home/hovdebt/Scenedesmus/Scenedsums_raw_data_2017/A1-ScenedesmusDOE0152-Z_S5_L002_R*_001.fastq.gz```

RNA Fastq files were pulled from the following location and combined into two files (R1 and R2):

```/panfs/biopan04/scratch-krusec/Scen152Z_RNAseq/2152_*_R*_001.fastq.gz```

Gene annotation file:

```/panfs/biopan04/scratch-krusec/Scen152Z_RNAseq/Final152ZGaplessAssembly/FINAL_scenedesmus_Complete.gff3```

## Analysis

The following diagram details the methods applied to benchmarking the ensemble SNV pipeline against the GIAB gold 
standard truth set (both DNA- and RNA-seq input):

![GIAB methods](../../media/GIAB_methods.png?raw=true)

The RTG ```vcfeval``` tool was used to perform the comparison of the experimental callset to the baseline callset:

```rtg vcfeval -b <baseline callset VCF> -c <experimental callset VCF> -e <confident regions BED file> -t <SDF of reference> -o <output directory>```

This command generates an output table with the following metrics:

| Threshold | True-pos-baseline | True-pos-call | False-pos | False-neg | Precision | Sensitivity | F-measure |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 6.000 | 3192886 | 3192899 | 9026 | 62491 | 0.9972 | 0.9808 | 0.9889 |
| None | 3217514 | 3217527 | 12925 | 37863 | 0.9960 | 0.9884 | 0.9922 |

**True-pos-baseline:** Number of variants called in baseline VCF matching those found in experimental VCF
**True-pos-call:** Number of variants called in experimental VCF matching those found in baseline VCF
**False-pos:** Number of variants called in experimental VCF that do not match the baseline VCF calls (w/in confident regions only)
**False-neg:** Number of variants missed in experimental VCF that are found in the baseline VCF callset (w/in confident regions only)
**Precision:** TP<sub>call</sub>/(TP<sub>call</sub> + FP); where TP == true positives and FP == false positives
**Sensitivity:** TP<sub>baseline</sub>/(TP<sub>baseline</sub> + FN); where TP == true positives and FN == false negatives
**F-measure:** Harmonic mean of precision and sensitivity metrics; accuracy metric accounting for class imbalance











