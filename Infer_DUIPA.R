suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
warning("off")
options(warn = -1)
args = commandArgs(TRUE)
option.list = list(
    make_option(c("-b","--bam"),action = "store_true",default = FALSE,help = "input file contains all filenames of bamfile [%default]"),
    make_option(c("-I","--IPUI"),action = "store_true",default = FALSE,help = "input file generated in the previous step which contains all IPUI value [%default]"),
    make_option(c("-d","--dir"),action = "store_true",default = FALSE,help = "input directory name contains all exoncount files [%default]"),
    make_option(c("-o","--output"),action = "store_true",default = FALSE,help = "output final results[%default]")
    )
desc = "Infer differentially used IPA sites using DEXSeq model"
parser = OptionParser(option_list = option.list, description = desc)
opt = parse_args(parser, args = args, positional_arguments = TRUE)
file_name = opt$args[1]
cfg_file = read.table(opt$args[1])
filenames = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist()) %>% as.character()
filename_prefix = c()
for (filename in filenames){
  file_prefix = strsplit(filename,".",fixed=T) %>% unlist() %>% .[1]
  filename_prefix = c(filename_prefix,file_prefix)
}
condition_length = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist() %>% length())
sampleTable = data.frame(row.names = filename_prefix,condition = c(rep("condition1",condition_length[1]),rep("condition2",condition_length[2])))
countFiles = list.files(opt$args[3], pattern="exoncount.txt", full.names=TRUE)
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData=sampleTable,design= ~ sample + exon + condition:exon,flattenedfile=NULL)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dexresult = DEXSeqResults( dxd )
dexresult = dexresult %>% as.data.frame() %>% dplyr::filter(grepl("intron",featureID)) %>% dplyr::select("groupID","featureID","pvalue","padj") %>% 
  dplyr::mutate(featureID=stringr::str_sub(featureID,2))
colnames(dexresult) = c("SYMBOL","intron_rank","pvalue","padj")
IPUI_table = read.table(opt$args[2],header = T)
dt = inner_join(IPUI_table,dexresult)
final = dt %>% dplyr::mutate(change=ifelse(padj<0.05&abs(IPUI_diff)>0.1,ifelse(IPUI_diff>0,"UP","DOWN"),"NOT"))
write.table(final,opt$args[4],quote = FALSE,sep = "\t",row.names = FALSE)
