# IPAFinder

## Description
IPAFinder performs *de novo* identification and quantification of dynamic IpA events using standard RNA-seq, regardless of any prior poly(A) site annotation. Assuming there is an intronic poly(A) site used in a given intron, IPAFinder models the normalized single-nucleotide resolution RNA-seq read coverage profiles and identifies profound drop in coverage, which can be interpreted as evidence of polyadenylation processing, to infer the used poly(A) site, as reflected by the lowest ratio of the sum of mean squared error (MSE) in the upstream and downstream segments splited by break point and the MSE computed for the entire intron region (Ratiomse). To detect skipped IpA (or splicing-coupled IpA), IPAFinder recognized cryptic 3′ splice site by junction-spanning reads and concatenated preceding exon to potential terminal exon. IPAFinder also has the ability to exclude alternative splicing events such as alternative 5′ splice site and cryptic exon activation by recognizing junction-spanning reads.

##  Diagram depicts the IPAFinder algorithm. 
![Sketch](https://github.com/ZhaozzReal/IPAFinder/blob/master/IPAFinder_diagram.jpg)

## Manual

### Step 1: 
### Get specialized annotation file contain intron and exon information from RefSeq GTF file(could be downloaded from the UCSC website: https://genome.ucsc.edu/).


**Command Example:**

```python IPAFinder_GetAnno.py -gtf hg38refGene.gtf -output IPAFinder_anno_hg38.txt```

### Step 2:
### Detect dynamic IpA events from samples with replicates between two conditions using standard RNA-seq

**Command Example:**

```python IPAFinder_multisamples.py -cfg cfg.txt -anno_txt IPAFinder_anno_hg38.txt -processors 5 -output IPAFinder_output.txt```

#### example of cfg.txt:

```
condition1=ctrl1.bam,ctrl2.bam 
condition2=case1.bam,case2.bam
```


### Step 2*:
### Detect dynamic IpA events from samples without replicates between two conditions using standard RNA-seq
###  input files
 1、 Mapped reads .bam file, Should be sorted and index using samtools
 2、 Specialized annotation file generated from Step 1
 
 **Command Example:**
 
 ```python IPAFinder_Pairsample.py -N normal.bam -T tumor.bam -anno_txt IPAFinder_anno_hg38.txt -processors 5 -output IPAFinder_output.txt```

## The following python packages are necessary
HTSeq、pysam、itertools、numpy、collections、multiprocessing、scipy、argparse、copy、subprocess
