#!/bin/bash

### NOTES ###
### This script uses the bedops tool (v2.4.39)
### This script needs to run from the directory in which the original bed file is contained! ###

# first argument to specify at command line represents input bed file
bed=$1
# iterate over each chr/contig and get a total count of genomic regions listed in bed file;
# add counts to new file
for chr in `bedextract --list-chr $bed`; do
    bedextract $chr $bed > ${bed/%.bed/.$chr.bed}
    count=`wc -l ${bed/%.bed/.$chr.bed} | cut -d ' ' -f 1`
    echo -e "$chr\t$count"
done > count.txt
# second argument to specify at command line represents total number of genomic regions to output
sample_size=$2
# get total number of regions in bed file
sum=`cat count.txt | awk '{sum+=$2} END {print sum}'`
# use previous variables to get proportion of regions for each chr/contig in original bed file; 
# given sample size, output number of regions for each chr/contig to select for output bed file
awk -v sum=$sum -v sample_size=$sample_size '{print $0"\t" ($2/sum)"\t"int(sample_size*$2/sum)}' count.txt > proportions.txt

# remove file if already present; the following appends so will add to files already existing for the same bed file
if test -f ${bed/%.bed/.subsampled.bed}; then
    rm ${bed/%.bed/.subsampled.bed}
fi
# iterate over each separate chr/contig file; shuffle each and extract the number provided in 'proportions.txt'
while read chr regions proportion count; do
    shuf ${bed/%.bed/.$chr.bed} | head -n $count >> ${bed/%.bed/.subsampled.bed}
done < proportions.txt

# clean up intermediate files; can easily restore by rerunning script
rm proportions.txt
rm count.txt
for chr in `bedextract --list-chr $bed`; do
    rm ${bed/%.bed/.$chr.bed}
done