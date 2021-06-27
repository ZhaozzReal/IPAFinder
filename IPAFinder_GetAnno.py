import HTSeq
import collections
import copy
import subprocess
from argparse import ArgumentParser,ArgumentTypeError
parser = ArgumentParser(description="De novo analysis of Intronic Polyadenylation from standard RNA-seq (convert gtf file to annotation file recognized by IPAFinder)")
parser.add_argument("-gtf",dest='RefSeq',action="store",type=str,help="Input RefSeq gtf file")
parser.add_argument("-output",dest="Anno",action="store",type=str,help="Output specialized annotation files contain intron and exon information")
args = parser.parse_args()
GTF2GPD = 'gtfToGenePred -genePredExt -ignoreGroupsWithoutExons -geneNameAsName2 %s %s'
gpdfile = (args.RefSeq).split(".")[0]+".gpd"
cmd = GTF2GPD % (args.RefSeq,gpdfile)
status = subprocess.call(cmd,shell=1)

def getOverlapFeatures(key, val):
    assert isinstance(val, list)
    ls_sorted = sorted(val, key=lambda x:(x[3], x[4]))
    n = 1
    ls_sorted[0].append(key + "/" + str(n))
    for i in range(1, len(ls_sorted)):
        if ls_sorted[i-1][4] > ls_sorted[i][3]:
            ls_sorted[i].append(key + "/" + str(n))
        else:
            n += 1
            ls_sorted[i].append(key + "/" + str(n))
    return ls_sorted


def check_overlap(iv, gas):
    assert isinstance(iv, HTSeq.GenomicInterval)
    assert isinstance(gas, HTSeq.GenomicArrayOfSets)
    features = set()
    for iv2, val in gas[iv].steps():
        features |= val
    if len(features) == 1:
        return 'Clean'
    elif len(features) > 1:
        return 'Overlap'
    else:
        return 'NA'


gene_label_dict = collections.OrderedDict()
for line in open(gpdfile,"r"):
    transcript_id,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,a,genename,b,c,d = line.strip().split()
    label = chrom+":"+genename+"|"+strand
    fields = [transcript_id, chrom, strand, int(txStart), int(txEnd), cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, a, genename, b,c,d]
    if "_" in chrom:
        pass
    else:
        if label not in gene_label_dict:
            gene_label_dict[label] = [fields]
        else:
            gene_label_dict[label].append(fields)


converted_list = []
for label in gene_label_dict:
    gene_label_dict[label] = getOverlapFeatures(label, gene_label_dict[label])
    for ls in gene_label_dict[label]:
        converted_list.append(tuple(ls))


gene_dict = collections.OrderedDict()
for line in converted_list:
    label = line[-1]
    gene_dict.setdefault(label, []).append(line)



shared_feature_dict = collections.OrderedDict()
gas_global = HTSeq.GenomicArrayOfSets("auto", stranded=False)#if two gene from different strand has overlap,discard them 
for label in gene_dict:
    gas_local = HTSeq.GenomicArrayOfSets("auto",stranded=True)
    last_exon_dict = {}
    for line in gene_dict[label]:
        transcript_id,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,a,genename,b,c,d,label = line
        exonStartsList = list(map(int, exonStarts.strip().rstrip(',').split(',')))
        exonEndsList = list(map(int, exonEnds.strip().rstrip(',').split(',')))
        intronStartList = copy.deepcopy(exonEndsList)[0:len(exonEndsList)-1]
        intronEndList = copy.deepcopy(exonStartsList)[1:len(exonStartsList)]
        if strand == "+":
            exonStartsList = exonStartsList[:-1]
            exonEndsList = exonEndsList[:-1]
        else:
            exonStartsList = exonStartsList[1:]
            exonEndsList = exonEndsList[1:]
        for s,e in zip(exonStartsList,exonEndsList):
            iv = HTSeq.GenomicInterval(chrom,s,e,strand)
            gas_local[iv] += (label,"exon")
            gas_global[iv] += (label,"exon")
        for s,e in zip(intronStartList,intronEndList):
            if e-s < 1:
                continue
            iv = HTSeq.GenomicInterval(chrom,s,e,strand)
            gas_local[iv] += (label,"intron")
            gas_global[iv] += (label,"intron")
    boundary_left,boundary_right = min([int(x[3]) for x in gene_dict[label]]),max([int(x[4]) for x in gene_dict[label]])
    for iv,val in gas_local[HTSeq.GenomicInterval(chrom,boundary_left,boundary_right,strand)].steps():
        if len(val) == 1:
            shared_feature_dict.setdefault(tuple(val)[0][0],[]).append((iv,tuple(val)[0][1]))


file_annot = open(args.Anno, "w")
for label in shared_feature_dict:
    SYMBOL = label.split(":")[1].split("|")[0]
    features_filtered = [item for item in shared_feature_dict[label] if check_overlap(item[0],gas_global) == 'Clean']
    for iv,feature in features_filtered:
        rank = features_filtered.index((iv,feature)) + 1
        pos = ('%s:%d-%d:%s') %(iv.chrom,iv.start,iv.end,iv.strand)
        exon_rank_left = -1
        exon_rank_right = -1
        intron_rank_left = -1
        intron_rank_right = -1
        if feature == 'intron':
            rank_exon = [features_filtered.index(x)+1 for x in features_filtered if x[1]=="exon"]
            try:
                exon_rank_left = max([x for x in rank_exon if x < rank])
            except:
                exon_rank_left = -1
            try:
                exon_rank_right = min([x for x in rank_exon if x > rank])
            except:
                exon_rank_right = -1
            file_annot.write('\t'.join(map(str,[label,feature,rank,pos,iv.end-iv.start,exon_rank_left,exon_rank_right]))+'\n')
        if feature == 'exon':
            rank_intron = [features_filtered.index(x)+1 for x in features_filtered if x[1]=="intron"]
            try:
                intron_rank_left = max([x for x in rank_intron if x < rank])
            except:
                intron_rank_left = -1
            try:
                intron_rank_right = min([x for x in rank_intron if x > rank])
            except:
                intron_rank_right = -1
            file_annot.write('\t'.join(map(str,[label,feature,rank,pos,iv.end-iv.start,intron_rank_left,intron_rank_right]))+'\n')

file_annot.close()