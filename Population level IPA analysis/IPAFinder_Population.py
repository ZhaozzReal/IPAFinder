import HTSeq,itertools, warnings, collections
import numpy as np
from multiprocessing import Pool
from argparse import ArgumentParser,ArgumentTypeError
parser = ArgumentParser(description = "IPA detection and quantification from population-level RNA-seq samples")
parser.add_argument("-b",dest = 'bamfiles',action = "store",type = str,help = "Input file contains all RNA-seq bamfiles")
parser.add_argument('-anno',dest = 'anno_txt',action = "store",type = str,help = "Input annotation file contains intron and flanking exons information")
parser.add_argument("-p",dest = "processors",action = "store",default = 10,type = int,help = "<INT> Number of processors used [default: 10]")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all inferred intronic poly(A) sites and their IPUI values")
args = parser.parse_args()
annot = collections.OrderedDict()
for line in open(args.anno_txt,"r"):
    gene_label,feature,rank,position,length,exon_rank_left,exon_rank_right = line.strip().split('\t')
    chrom,iv_str,strand = position.strip().split(':')
    start,end = map(int,iv_str.strip().split('-'))
    SYMBOL = gene_label.split(":")[1].split("|")[0]
    annot.setdefault(SYMBOL,[]).append((feature,int(rank),chrom,start,end,strand,int(length),exon_rank_left,exon_rank_right))


def Get_region_cvg_list(region,bam_reader):
    warnings.simplefilter("ignore")
    chrom,start_end = region.split(":")
    start,end = start_end.split("-")
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    read_seq = bam_reader.fetch(region = region)
    for a in read_seq:
        iv_seq = (cigop.ref_iv for cigop in a.cigar if cigop.type == "M" and cigop.size >0)
        for iv in iv_seq:
            ga[iv] += 1
    region_iv = HTSeq.GenomicInterval(chrom,int(start),int(end),"-")
    cvg_list = list(ga[region_iv])
    return cvg_list


def Get_region_cvg(region,bam_reader):
    chrom,start_end = region.split(":")
    start,end = start_end.split("-")
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    read_seq = bam_reader.fetch(region = region)
    for a in read_seq:
        iv_seq = (cigop.ref_iv for cigop in a.cigar if cigop.type == "M" and cigop.size >0)
        for iv in iv_seq:
            ga[iv] += 1
    region_iv = HTSeq.GenomicInterval(chrom,int(start),int(end),"-")
    return round(np.mean(sorted(list(ga[region_iv]),reverse = True)[:30]),3)


def Get_Skipend_dict(region_fetch,bamfile,strand):
    warnings.simplefilter("ignore")
    bam_reader = HTSeq.BAM_Reader(bamfile)
    read_seq = bam_reader.fetch(region = region_fetch)
    skip_list = []
    for a in read_seq:
        if strand == "+":
            skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
    skip_dict = dict(collections.Counter(skip_list))
    return skip_dict


def Get_Skipstart_dict(region_fetch,bamfile,strand):
    warnings.simplefilter("ignore")
    skip_list = []
    bam_reader = HTSeq.BAM_Reader(bamfile)
    read_seq = bam_reader.fetch(region=region_fetch)
    for a in read_seq:
        if strand == "+":
            skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
    skip_dict = dict(collections.Counter(skip_list))
    return skip_dict


def Get_SS_skip(skipstart_dict,intron_start,intron_end):
    SS_list = []
    for key,value in skipstart_dict.items():
        if int(intron_start) -10 < int(key) < int(intron_start) + 10 or int(intron_end) -10 < int(key) < int(intron_end) + 10:
            SS_list.append(value)
    SS_value = 100
    if SS_list != []:
        SS_value = max(SS_list)
    return SS_value


def Get_max_skip(skipend_dict,intron_start,intron_end):
    Intron_list = []
    for key,value in skipend_dict.items():
        if int(intron_start) + 100 < int(key) < int(intron_end) - 100:
            Intron_list.append(value)
    if Intron_list != []:
        num = max(Intron_list)
        pos = [key for key,value in skipend_dict.items() if value == num][0]
    else:
        num = 0
        pos = 0
    return pos,num


def Estimation_abundance(Region_Coverage,break_point):
    downstream_cov_mean = np.mean(Region_Coverage[break_point:])
    upstream_cov_mean = np.mean(Region_Coverage[:break_point])
    Coverage_diff = Region_Coverage[:break_point]-upstream_cov_mean
    Coverage_diff = np.append(Coverage_diff,Region_Coverage[break_point:]-downstream_cov_mean)
    Mean_Squared_error = np.mean(Coverage_diff**2)
    return Mean_Squared_error


def Get_min_mseratio(cvg_region):
    search_result = [[],[]]
    global_mse = np.mean((cvg_region - np.mean(cvg_region))**2)
    for curr_search_point in range(50,len(cvg_region),int(len(cvg_region)/50)):
        Mean_Squared_error = Estimation_abundance(cvg_region, curr_search_point)
        mse_ratio = Mean_Squared_error/global_mse
        search_result[0].append(mse_ratio)
        search_result[1].append(curr_search_point)
    min_mse = min(search_result[0])
    min_mse_index = search_result[0].index(min_mse)
    min_mse_point = search_result[1][min_mse_index]
    return min_mse,min_mse_point


def Get_IPAevent(input_tuple):
    warnings.simplefilter("ignore")
    label,bamfile = input_tuple
    IPA_result = []
    intron_list = [region for region in annot[label] if region[0] == 'intron' and int(region[6]) > 250]
    if intron_list != []:
        for intron in intron_list:
            feature,rank,chrom,intron_start,intron_end,strand,length,exon_rank_left,exon_rank_right = intron
            intron_region = chrom + ":" + str(intron_start) + "-" + str(intron_end)
            bam_reader = HTSeq.BAM_Reader(bamfile)
            cvg = Get_region_cvg_list(intron_region,bam_reader)
            if strand == "-":
                cvg = cvg[::-1]
            IPAtype = "Composite"
            coverage_threshold = 20
            if max(cvg) > coverage_threshold:
                skipend_dict = Get_Skipend_dict(intron_region,bamfile,strand)
                skipstart_dict = Get_Skipstart_dict(intron_region,bamfile,strand)
                max_skip_pos,max_skip_num = Get_max_skip(skipend_dict,intron_start,intron_end)
                SS_skip_num = Get_SS_skip(skipstart_dict,intron_start,intron_end)
                if  max_skip_num > SS_skip_num/5 or max_skip_num > 10:
                    if strand == "+":
                        skip_position = int(max_skip_pos) - int(intron_start)
                        intron_start = int(max_skip_pos)
                    else:
                        skip_position = int(intron_end) - int(max_skip_pos)
                        intron_end = int(max_skip_pos)
                    cvg = cvg[skip_position:]
                    IPAtype = "Skipped" + "_" + str(max_skip_pos)
                if len(cvg) > 200:
                    min_mseratio,min_mse_point = Get_min_mseratio(cvg)
                    if min_mseratio < 0.5:
                        IPA_point = int(min_mse_point)
                        IPAabundance = np.mean(cvg[:IPA_point])
                        upstream_cov = len(list(filter(lambda x:x>IPAabundance/2,cvg[:IPA_point])))/IPA_point
                        downstream_cov = len(list(filter(lambda x:x>IPAabundance/2,cvg[IPA_point:])))/(len(cvg)-IPA_point)
                        if upstream_cov > 0.5 and downstream_cov < 0.2:
                            if strand == "+":
                                IPA_location = int(intron_start) + IPA_point
                                IPA_region = chrom + ":" + str(intron_start) + "-" + str(IPA_location)
                                exon_rank = int(exon_rank_left)
                            else:
                                IPA_location = int(intron_end) - IPA_point
                                IPA_region = chrom + ":" + str(IPA_location) + "-" + str(intron_end)
                                exon_rank = int(exon_rank_right)
                            for key,value in skipstart_dict.items():
                                if IPA_location - 50 < int(key) < IPA_location + 50  and int(intron_start) + 50 < int(key) < int(intron_end) - 50 and int(value) >  IPAabundance/2:
                                    break
                            else:
                                if len([feature for feature in annot[label] if feature[1] == exon_rank]) == 1:
                                    exon_inf = [feature for feature in annot[label] if feature[1] == exon_rank][0]
                                    exon_region = exon_inf[2] + ':' + str(exon_inf[3]) + "-" + str(exon_inf[4])
                                    exon_abundance = Get_region_cvg(exon_region,bam_reader)
                                    IPUI = round(IPAabundance/exon_abundance,3)
                                    if IPUI < 1:
                                        IPA_inf = label + ";" + feature + "_" + str(rank) + ";" + strand + ";" +  IPAtype + ";" + IPA_region + ";" + exon_region
                                        IPA_result.append(IPA_inf)
    return IPA_result




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
    SYMBOL,intron,strand,IPAtype,IPA_region,exon_region = infor.split(";")
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


bamfile_list = open(args.bamfiles).readline().strip("\n").split(",")
out = open(args.outfile,"w")
out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Genename","Intron","Strand","IPAtype","IPA_region","exon_region","\t".join(bamfile_list)))


mergeIPA_list = []
IPAevents_list = []
for bamfile in bamfile_list:
    input_tuple = list(zip(annot.keys(),[bamfile]*len(annot)))
    pool = Pool(args.processors)
    result_list = pool.map(Get_IPAevent,input_tuple)
    result_list = [ i for i in result_list if i != []]
    for gene_list in result_list:
        for IPA in gene_list:
            Genename,Intron,Strand,IPAtype,IPA_terminal_region,Upstream_exon_region = IPA.split(";")
            IPAevents = Genename + ":" + Intron + ":" + IPAtype
            if IPAevents not in IPAevents_list:
                mergeIPA_list.append(IPA)
                IPAevents_list.append(IPAevents)


for IPA in mergeIPA_list:
    input_tuple = list(zip(bamfile_list,[IPA]*len(bamfile_list)))
    pool = Pool(args.processors)
    result_list = pool.map(Cal_IPUI,input_tuple)
    out.write("{}\t{}\n".format("\t".join(IPA.split(";")),"\t".join(map(str,result_list))))


out.close()