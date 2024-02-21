# pRMNE
The probability of Not Excluding a Random Man  (pMNE)

A novel bioinformatics and statiatical tool for calculating the probability of Random Man Not Excluded (pRMNE) based on multiple loci alleles and popualtion frequency.
The tool is programmed in Python3.

## Installation

The tool is prograqmmed in Pyhton v3.9. No installation is needed. It can be run after downloading from Github using the following command in a terminal.

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
    
## Output

The output is a tab seperatedly text file. The data is the exclusion probability of a random man in each super population at each STR locus site and all combined sites for  26 populations (FIN,CDX,IBS,MXL,CHB,CHS,CEU,JPT,ESN,KHV,TSI,CLM,YRI,GBR,PEL,STU,BEB,GIH,PJL,MSL,ITU,GWD,LWK,ASW,PUR,ACB), 5 supperpopulations (EAS, EUR, AFR,SAS, AMR).

e.g. for AFR (Africa popupation): 20 CODIS loci and  All loci combined at a p value 8.23071958151542135332E-36

    Population	Locus	Probability_of_Exclusion(PE)	Probability_of_Inclusion(PI)	Adjusted_PE²	Adjusted_PI²	Adjusted_PE³	Adjusted_PI³
    AFR	CSF1PO	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93217121942747780534E-1	6.78287805725221946563E-3	9.93217121942747780534E-1	6.78287805725221946563E-3
    AFR	D10S1248	9.99777233236801069280E-1	2.22766763198930719537E-4	9.93316271208893620440E-1	6.68372879110637956039E-3	9.93316271208893620440E-1	6.68372879110637956039E-3
    AFR	D12S391	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93266787826630159409E-1	6.73321217336984059149E-3	9.93266787826630159409E-1	6.73321217336984059149E-3
    AFR	D13S317	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93217121942747780534E-1	6.78287805725221946563E-3	9.93217121942747780534E-1	6.78287805725221946563E-3
    AFR	D16S539	9.99777233236801069280E-1	2.22766763198930719537E-4	9.84961610220010645989E-1	1.50383897799893540109E-2	9.84961610220010645989E-1	1.50383897799893540109E-2
    AFR	D18S51	9.99943467691785855616E-1	5.65323082141443835152E-5	9.84738524371182506202E-1	1.52614756288174937977E-2	9.84738524371182506202E-1	1.52614756288174937977E-2
    AFR	D19S433	9.99777233236801069280E-1	2.22766763198930719537E-4	9.73265084835574481758E-1	2.67349151644255182415E-2	9.73265084835574481758E-1	2.67349151644255182415E-2
    AFR	D1S1656	9.99943467691785855616E-1	5.65323082141443835152E-5	9.98304280485686945134E-1	1.69571951431305486641E-3	9.98304280485686945134E-1	1.69571951431305486641E-3
    AFR	D21S11	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93217121942747780534E-1	6.78287805725221946563E-3	9.93217121942747780534E-1	6.78287805725221946563E-3
    AFR	D22S1045	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93217121942747780534E-1	6.78287805725221946563E-3	9.93217121942747780534E-1	6.78287805725221946563E-3
    AFR	D2S1338	9.99777233236801069280E-1	2.22766763198930719537E-4	9.73265084835574481758E-1	2.67349151644255182415E-2	9.73265084835574481758E-1	2.67349151644255182415E-2
    AFR	D2S441	9.99777233236801069280E-1	2.22766763198930719537E-4	9.73265084835574481758E-1	2.67349151644255182415E-2	9.73265084835574481758E-1	2.67349151644255182415E-2
    AFR	D3S1358	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93266787826630159409E-1	6.73321217336984059149E-3	9.93266787826630159409E-1	6.73321217336984059149E-3
    AFR	D5S818	9.99943467691785855616E-1	5.65323082141443835152E-5	9.93217121942747780534E-1	6.78287805725221946563E-3	9.93217121942747780534E-1	6.78287805725221946563E-3
    AFR	D7S820	9.99943467691785855616E-1	5.65323082141443835152E-5	9.84738524371182506202E-1	1.52614756288174937977E-2	9.84738524371182506202E-1	1.52614756288174937977E-2
    AFR	D8S1179	9.99943467691785855616E-1	5.65323082141443835152E-5	9.84738524371182506202E-1	1.52614756288174937977E-2	9.84738524371182506202E-1	1.52614756288174937977E-2
    AFR	FGA	9.99943467691785855616E-1	5.65323082141443835152E-5	9.98304280485686945134E-1	1.69571951431305486641E-3	9.98304280485686945134E-1	1.69571951431305486641E-3
    AFR	TH01	9.99777233236801069280E-1	2.22766763198930719537E-4	9.84961610220010645989E-1	1.50383897799893540109E-2	9.84961610220010645989E-1	1.50383897799893540109E-2
    AFR	TPOX	9.99777233236801069280E-1	2.22766763198930719537E-4	9.84961610220010645989E-1	1.50383897799893540109E-2	9.84961610220010645989E-1	1.50383897799893540109E-2
    AFR	vWA	9.99943467691785855616E-1	5.65323082141443835152E-5	9.84813068822410803974E-1	1.51869311775891960257E-2	9.84813068822410803974E-1	1.51869311775891960257E-2
    AFR	All_loci	1.00000000000000000000E+0	8.23071958151542135332E-36	1.00000000000000000000E+0	4.38875716538970662812E-41	1.00000000000000000000E+0	4.38875716538970662812E-41



## Support
Version: 1.1.0, Dec,5th,2023
Support: xwang.kib@gmail.com
