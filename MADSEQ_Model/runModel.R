library("MADSEQ")
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
options(stringsAsFactors = F)


######################################################
####### 1. Run Model for Each Chromosome #############
######################################################
hetero = paste0("heterozygous/",ID,".raw.variants.vcf.gz_filtered_heterozygous.txt")  ## specify the location you put processed heterozygous data
cov = paste0("normalized_cov/",ID,"_normed_depth.txt") ## specify the location you put normalized coverage data

print(ID)
pdf(file=paste0("plot/",ID,"_chr",chr,"_AAF.pdf"),width=6,height=4)  ## if you want to produce a AAF figure for each chromosome
res = runMadSeq(hetero = hetero, ## heterozygous file
                coverage = cov,  ## coverage file
                target_chr = 1,  ## chromosome number you want to test, please note you should use "chr1" if your assembly contig ID is "chr1"
                adapt=5000,  ## adapt steps
                burnin = 5000)  ## burn in steps
dev.off()
saveRDS(res,file=paste0("result/",ID,"_chr1.RDS")) ## save result



#################################################
####### 2. Examine Result using BIC #############
#################################################
## 1). Load result for each chromosome, and plot model selection figures using BIC values
source("plotBIC.R")  ## functions to plot BIC 

CHR = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
        "chr18","chr19","chr20","chr21","chr22","chrX")

BIC = as.data.frame(t(sapply(CHR,function(x)readBIC(paste0("result/",ID,"_",x,".RDS"))))) ## read BIC for chr1 to chrX
saveRDS(BIC,file=paste0("BIC/",indv,"_BIC.RDS")) ## save BIC data

## transform BIC value to make it plot friendly
BIC_net = BIC
BIC_net[,2:5] = BIC_net[,2:5] + 10  ## this step is optional, if you only want to keep high confidence detection, add 10 (threshold) to other models but not null model
BIC_net = BIC_net - apply(BIC_net,1,min)  ## calculate delta BIC
BIC_net_log = -log(BIC_net)  ## log transformation
BIC_log = apply(BIC_net_log,2,function(x)ifelse(x==Inf,1,x))  ## convert infinite values to 1 to make it plottable
rownames(BIC_log) = CHR
BIC_log = data.frame(BIC_log,chr=CHR)
BIC_melt = melt(BIC_log,id.vars = "chr")
colnames(BIC_melt) = c("chr","chromosome_status","deltaBIC")
BIC_melt$chromosome_status = gsub("BIC_","",BIC_melt$chromosome_status)
## plot model selection figure
plotBIC(BIC_melt,ID)

## from the BIC plot you can know which chromosome is not normal


#################################################
####### 3. Plot Whole Genome Data   #############
#################################################
source("plotGenome.R")  ## functions to plot genome
## plot alternative allele frequency (AAF) for whole genome
g_aaf = plotAAF(hetero)  ## input is the position of processed heterozygous sites file
g_aaf = g_aaf + labs(title=ID)
## plot coverage for whole genome
g_cov = plotCov(cov) + labs(title=ID) ## input is the position of normalized coverage file

grid.arrange(g_aaf,g_cov,nrow=2)