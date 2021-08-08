# IPAFinder



## Description

IPAFinder performs *de novo* identification and quantification of dynamic IpA events using standard RNA-seq, regardless of any prior poly(A) site annotation. Assuming there is an intronic poly(A) site used in a given intron, IPAFinder models the normalized single-nucleotide resolution RNA-seq read coverage profiles and identifies profound drop in coverage to infer the used poly(A) site. To detect skipped IpA (or splicing-coupled IpA), IPAFinder recognized cryptic 3′ splice site by junction-spanning reads and concatenated preceding exon to potential terminal exon. IPAFinder also has the ability to exclude alternative splicing events such as alternative 5′ splice site and cryptic exon activation by recognizing junction-spanning reads.



##  Diagram depicts the IPAFinder algorithm. 

<img src="https://github.com/ZhaozzReal/IPAFinder/blob/master/Diagram.jpg" width="600" height="400"/>

## Installation

IPAFinder consists of both Python (3.5+) and R scripts:

1. Install the following software pre-requisites:

   i. python (required packages HTSeq, pysam, itertools, numpy, collections, multiprocessing, scipy, argparse, copy, and subprocess)

   ii. R (required packages optparse, dplyr, stringr, and DEXSeq)

2. Clone the lastest development version of IPAFinder and change directory:

 ```
  git clone https://github.com/ZhaozzReal/IPAFinder.git
  cd IPAFinder
 ```



## Usage 

IPAFinder has three sub-commands:

1.```IPAFinder_GetAnno.py```: Generate annotation file containing intron and exon information

2.```IPAFinder_DetectIPA.py```: Detect and quantify IPA sites, and calculate read counts of all exons

3.```Infer_DUIPA.R```: Infer differential usage of IPA sites



## Generate specialized annotation file

RefSeq GTF file could be downloaded from the UCSC website: https://genome.ucsc.edu/.

The UCSC tool ```gtfToGenePred``` is required here.

**Command**

```
python IPAFinder_GetAnno.py -gtf hg38refGene.gtf -output IPAFinder_anno_hg38.txt
```

We have generated annotation file for hg19, hg38 and mm10, and ones could download and use it directly.



## Detect and quantify IPA sites, and calculate read counts of all exons

**Command** 

```
python IPAFinder_DetectIPA.py -b allbamfiles.txt -anno IPAFinder_anno_hg38.txt -p 5 -o IPAFinder_IPUI.txt
```

allbamfiles.txt contains all filename of bamfile between two conditions, as shown below:

```
condition1=ctrl1.bam,ctrl2.bam 
condition2=case1.bam,case2.bam
```

Following counting reads mapped to all exons, IPAFinder expects the results to be located inside its own sub-directory. For example, new generated results may appear with the following directory structure:

```
project/
  |-- ctrl1_exoncount.txt
  |-- ctrl2_exoncount.txt
  |-- case1_exoncount.txt
  |-- case2_exoncount.txt
```



## Infer differential usage of IPA sites

DEXSeq, which is widely used for differential exon usage analysis on RNA-seq data, was applied to detect differential usage of IPA sites. This statistical framework could account for biological variability between replicates and is robust to changes in isoform abundance between conditions.

**Command**

```
Rscript Infer_DUIPA.R -b allbamfiles.txt -I IPAFinder_IPUI.txt -d project -o IPAFinder_DUIPA.txt
```

Final results will be saved in the file ```IPAFinder_DUIPA.txt```.

The final output format is as follows:

| Column        | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| SYMBOL        | gene symbol                                                  |
| Intron_rank   | rank number of intron contains IpA event                     |
| Terminal_exon | genomic location of corresponding terminal exon of IpA isoform |
| IPAtype       | type of terminal exon (Skipped or Composite)                 |
| *ctrl1*       | IPUI estimate for sample *ctrl1*                             |
| *ctrl2*       | IPUI estimate for sample *ctrl2*                             |
| *case1*       | IPUI estimate for sample *case1*                             |
| *case2*       | IPUI estimate for sample *case2*                             |
| IPUI_diff     | difference of mean IPUI between conditions                   |
| pvalue        | *P* value for testing differential usage of terminal exon    |
| padj          | adjusted *P* value for testing differential usage of terminal exon |
| change        | Direction of changed IpA event (UP, DOWN, or NOT)            |





### IPAFinder analysis on paired samples without replicates 

**Option 1: Infer differential used IPA sites using Fisher's exact test**

 ```
 python IPAFinder_PS_FET.py -b1 ctrl.bam -b2 case.bam -anno IPAFinder_anno_hg38.txt -p 5 -o IPAFinder_output.txt
 ```



**Option 2: Infer differential used IPA sites using bootstrapping-based method**

```
python IPAFinder_PS_FDR.py -b1 ctrl.bam -b2 case.bam -anno IPAFinder_anno_hg38.txt -p 5 -o IPAFinder_output.txt
```

