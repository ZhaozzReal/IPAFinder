import os
from argparse import ArgumentParser,ArgumentTypeError
parser = ArgumentParser(description = "Merge and obtain all non-redundant IPA events")
parser.add_argument("-path",dest = 'path',action = "store",type = str,help = "Path containing files of IPA events detected from single sample")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all non-redundant IPA events")
args = parser.parse_args()
IPA_files = os.listdir(args.path)

outfile = open(args.outfile,"w")
IPAevents_list = []
for filename in IPA_files:
    for line in open(args.path + filename,"r"):
        Genename,Intron,strand,IPAtype,IPA_terminal_region,Upstream_exon_region,IPUI,Samplename = line.strip().split("\t")
        IPAevents = Genename + ":" + Intron + ":" + IPAtype
        if IPAevents not in IPAevents_list:
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(Genename,Intron,strand,IPAtype,IPA_terminal_region,Upstream_exon_region))
            IPAevents_list.append(IPAevents)

outfile.close()
