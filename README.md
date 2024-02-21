# pRMNE
The probability of Not Excluding a Random Man based on alleles

A novel bioinformatics tool for calculating the probability of Random Man Not Excluded (pRMNE) based on multiple loci alleles and popualtion frequency.
The tool is programmed in Python3.

## Installation

The tool is prograqmmed in Pyhton3.8. No installation is needed. It can be run after downloading from Github using the following command in a terminal.

` git clone https://github.com/XuewenWangUGA/pRMNE`


or click the download button in the Github page to download the .zip file , then uncompressed the file.

## How to use

Usage: 

`python AllelExclusionProb.py [options]`

e.g., 

`python3 AllelExclusionProb.py -d databaseAllele.txt -e evidenceAllele.txt -s suspectAllele.txt -t 2 > outResults.txt`

Options:
    
        -h, --help: show this help message and exit
        -d, --databaseAlle: string, required input file of dataBase with Allele count data, tab delimited
        -e, --evidenceAlle: string, required input file of observed evidence Alleles, tab delimited
        -s, --suspectAlle: string, required input file of suspect Alleles, tab delimited
        -l, --log10: flag only, no value,output the probability after log10 conversion if -l is given.
        -t, --threads: int, the number of parallelized computing threads, default 2
         the result will be output to standard output/screen. use > to redirect to a file
    
## Support
Version: 1.1.0, Dec,5th,2023
Support: xwang.kib@gmail.com
