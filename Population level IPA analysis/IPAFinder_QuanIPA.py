import HTSeq,itertools, warnings, collections,sys
import numpy as np
from multiprocessing import Pool
from argparse import ArgumentParser,ArgumentTypeError
parser = ArgumentParser(description = "IPA detection from single RNA-seq sample")
parser.add_argument("-b",dest = 'bamfiles',action = "store",type = str,help = "Input file contains all RNA-seq bamfiles")
parser.add_argument('-IPA',dest = 'IPAevents',action = "store",type = str,help = "Input file contains all IPA events")
parser.add_argument("-p",dest = "processors",action = "store",default = 10,type = int,help = "<INT> Number of processors used [default: 10]")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all IPA usage")
args = parser.parse_args()



def Get_region_cvg_list(region,bam_reader):
    warnings.simplefilter("ignore")
    chrom,start_end = region.split(":")
    start,end = start_end.split("-")
    ga = HTSeq.GenomicArray("auto", stranded = False, typecode = "i")
    read_seq = bam_reader.fetch(region = region)
    for a in read_seq:
        iv_seq = (cigop.ref_iv for cigop in a.cigar if cigop.type == "M" and cigop.size >0)
        for iv in iv_seq:
            ga[iv] += 1
    region_iv = HTSeq.GenomicInterval(chrom,int(start),int(end),"-")
    cvg_list = list(ga[region_iv])
    return cvg_list


def Get_Skipend_num(chrom,Skipped_pos,bam_reader,strand):
    region_fetch = chrom + ":" + str(Skipped_pos - 10) + "-" + str(Skipped_pos + 10)
    read_seq = bam_reader.fetch(region = region_fetch)
    skip_list = []
    for a in read_seq:
        if strand == "+":
            skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
    skip_dict = dict(collections.Counter(skip_list))
    num_list = [value for key,value in skip_dict.items() if Skipped_pos -2 < key < Skipped_pos + 2]
    if num_list == []:
        num = 0.1
    else:
        num = max(num_list)
    return num



def Cal_IPUI(input):
    bamfile,infor = input
    SYMBOL,intron,strand,IPAtype,IPA_region,exon_region = infor.split("\t")
    chrom = IPA_region.split(":")[0]
    bam_reader = HTSeq.BAM_Reader(bamfile)
    exon_cvg = Get_region_cvg_list(exon_region,bam_reader)
    exonabundance = round(np.mean(sorted(exon_cvg,reverse = True)[:30]),3)
    if "Skipped" in IPAtype:
        Skipped_pos = int(IPAtype.split("_")[1])
        IPAabundance = Get_Skipend_num(chrom,Skipped_pos,bam_reader,strand)
    else:
        IPA_cvg = Get_region_cvg_list(IPA_region,bam_reader)
        IPAabundance = np.mean(IPA_cvg)
        ratio = len(list(filter(lambda x:x>IPAabundance/3,IPA_cvg)))/len(IPA_cvg)
        if ratio < 0.8:
            IPAabundance = 0
    if exonabundance > 25 and IPAabundance != 0:
        IPUI = round(IPAabundance/exonabundance + 0.001,3)
        if IPUI > 1:
            IPUI = None
    else:
        IPUI = None
    return IPUI



bamfiles = open(args.bamfiles).readline().strip("\n").split(",")
out = open(args.outfile,"w")
out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Genename","Intron","Strand","IPAtype","IPA_region","exon_region","\t".join(bamfiles)))


for line in open(args.IPAevents,"r"):
    if "Genename" in line:
        continue
    input_tuple = list(zip(bamfiles,[line.strip()]*len(bamfiles)))
    pool = Pool(args.processors)
    result_list = pool.map(Cal_IPUI,input_tuple)
    out.write("{}\t{}\n".format(line.strip(),"\t".join(map(str,result_list))))


out.close()



