# IPAFinder

## Description
IPAFinder performs de novo identification and quantification of dynamic IpA events using standard RNA-seq, regardless of any prior polyA (pA) site annotation. Assuming there is an intronic pA site used in a given intron, IPAFinder models the normalized single-nucleotide resolution RNA-seq read coverage profiles and identifies profound drop in coverage, which can be interpreted as evidence of PAS processing, to infer the used pA site, as reflected by the lowest ratio of the sum of mean squared error (MSE) in the upstream and downstream segments splited by break point and the MSE computed for the entire intron region (Ratiomse).

##  Diagram depicts the IPAFinder algorithm. 
![Sketch](https://github.com/ZhaozzReal/IPAFinder/blob/master/GitHub_manual.jpg)
