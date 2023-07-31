# IPAFinder



## Description

IPAFinder performs *de novo* identification and quantification of IPA events using standard RNA-seq data based on "change point" model, which has been widely used for tandem 3' UTR APA analysis. IPAFinder could exclude the interference of alternative splicing events such as alternative 5' splice site and cryptic exon activation by recognizing junction-spanning reads.



##  Diagram illuminates the IPAFinder algorithm. 

Panel A indicates the schematic diagram of IPAFinder in detecting composite terminal exon IPA event.

Panel B indicates the schematic diagram of IPAFinder in quantifying the usage of composite IPA site.

![image](https://github.com/ZhaozzReal/IPAFinder/blob/master/Diagram.jpg)

## Installation

IPAFinder consists of both Python (3.5+) and R scripts:

1. Install the following software pre-requisites:

   i. python (required packages HTSeq, itertools, numpy, collections, multiprocessing, scipy, argparse, os, warnings, and subprocess)

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
python IPAFinder_GetAnno.py -gtf /path/to/hg38refGene.gtf -output /path/to/IPAFinder_anno_hg38.txt
```

**We have generated annotation file for hg19, hg38 and mm10, and we suggest users utilize it directly.**



## Detect IPA sites and quantify their usage

**Command** 

```
python IPAFinder_DetectIPA.py -b /path/to/allbamfiles.txt -anno /path/to/IPAFinder_anno_hg38.txt -p 10 -o /path/to/IPAFinder_IPUI.txt
```

allbamfiles.txt contains all filename of bamfile between two conditions, as shown below:

```
condition1=/path/to/ctrl1.bam,/path/to/ctrl2.bam 
condition2=/path/to/case1.bam,/path/to/case2.bam
```

Following counting reads mapped to all exons, IPAFinder expects the results to be located inside its own sub-directory. For example, new generated results may appear with the following directory structure:

```
project/
  |-- ctrl1_exoncount.txt
  |-- ctrl2_exoncount.txt
  |-- case1_exoncount.txt
  |-- case2_exoncount.txt
```



## Infer statistically differential usage of IPA sites

DEXSeq, the model for differential exon usage analysis based on standard RNA-seq data, was applied to detect differential usage of IPA terminal exon. This statistical framework could account for biological variability between replicates and is robust to changes in isoform abundance between conditions.

**Command**

```
Rscript Infer_DUIPA.R -b /path/to/allbamfiles.txt -I /path/to/IPAFinder_IPUI.txt -d /path/to/project -o /path/to/IPAFinder_DUIPA.txt
```

Final results will be saved in the file ```IPAFinder_DUIPA.txt```.

The final output format is as follows:

| Column        | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| SYMBOL        | gene symbol                                                  |
| Intron_rank   | rank number of intron contains IPA event                     |
| Terminal_exon | genomic location of corresponding terminal exon of IPA isoform |
| IPAtype       | type of terminal exon (Skipped or Composite)                 |
| *ctrl1*       | IPUI estimate for sample *ctrl1*                             |
| *ctrl2*       | IPUI estimate for sample *ctrl2*                             |
| *case1*       | IPUI estimate for sample *case1*                             |
| *case2*       | IPUI estimate for sample *case2*                             |
| IPUI_diff     | difference of mean IPUI between conditions                   |
| pvalue        | *P* value for testing differential usage of terminal exon    |
| padj          | adjusted *P* value for testing differential usage of terminal exon |
| change        | Direction of changed IPA event (UP, DOWN, or NOT)            |



## IPAFinder analysis on paired samples without replicates 



**Option 1: Infer differentially used IPA sites using Fisher's exact test-based method**

 ```
 python IPAFinder_PS_FET.py -b1 /path/to/ctrl.bam -b2 /path/to/case.bam -anno /path/to/IPAFinder_anno_hg38.txt -p 10 -o /path/to/IPAFinder_DUIPA.txt
 ```



**Option 2: Infer differentially used IPA sites using bootstrapping-based method**

```
python IPAFinder_PS_FDR.py -b1 /path/to/ctrl.bam -b2 /path/to/case.bam -anno /path/to/IPAFinder_anno_hg38.txt -p 10 -o /path/to/IPAFinder_DUIPA.txt
```

## Population-level IPA detection and quantification
Given that a large amount of publicly available standard RNA-seq data has been accumulated, we update the pipeline for population-level IPA analysis, which could be used on datasets from TCGA or GTEx project.

```
python IPAFinder_Population.py -b /path/to/allbamfiles.txt -anno /path/to/IPAFinder_anno_hg38.txt -p 10 -o /path/to/Allsamples_IPUI.txt
```
allbamfiles.txt contains all filename of bamfiles, as shown below:

```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

**We also provide substeps for population-level IPA analysis**

***Step1: Detect IPA events from single RNA-seq sample***
```
python IPAFinder_DetectIPA.py -b /path/to/sample1.bam -anno /path/to/IPAFinder_anno_hg38.txt -p 10 -o /path/to/sample1_IPA.txt
```

***Step2: Merge and ontain all non-redundant IPA events from all samples***
```
python IPAFinder_MergeIPA.py -path /path of sampleN_IPA.txt/ -o /path/to/IPA_merge.txt
```

***Step3: Quantify the usage of IPA events across all samples***
```
python IPAFinder_QuanIPA.py -b /path/to/allbamfiles.txt -IPA /path/to/IPA_merge.txt -p 10 -o /path/to/IPA_merge_IPUI.txt
```



## Citation

*Please cite the following articles if you use IPAFinder in your research:*

* Zhao Z, Xu Q, Wei R, Wang W, Ding D, Yang Y, Yao J, Zhang L, Hu YQ, Wei G, Ni T. Cancer-associated dynamics and potential regulators of intronic polyadenylation revealed by IPAFinder using standard RNA-seq data. Genome Res. 2021 Sep 2. doi: 10.1101/gr.271627.120. PMID: 34475268.

* Zhao Z, Xu Q, Wei R, Huang L, Wang W, Wei G, Ni T. Comprehensive characterization of somatic variants associated with intronic polyadenylation in human cancers. Nucleic Acids Res. 2021 Oct 11. doi: 10.1093/nar/gkab772. PMID: 34508351.



## Contact

If you have any comments, suggestions, questions, bug reports, etc, feel free to contact Zhaozhao Zhao (zz_zhaobioinfo@126.com). And PLEASE attach your command line and log messages if possible.

