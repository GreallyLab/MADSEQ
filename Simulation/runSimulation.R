library("GenomicRanges")
library("MADSEQ")
library("BSgenome")
options(stringsAsFactors = F)

######################################
####### 1. Simulate Data #############
######################################
## For simulation, we run each condition for 500 times
set.seed(100) ## set seed to create reproducible simulation
source("simulateData.R")  ## load functions to generate simulated data
data = generateData(nSNP=1000,  ## number of SNPs
                    nCov=1000,  ## number of coverage points
                    coverage=100,  ## sequencing depth
                    vartom=20,  ## variation level of coverage
                    m=0.47,  ## midpoint
                    fraction=0.3,  ## fraction of mosaic aneuploidy (0,1)
                    type="LOH",  ## type of aneuploidy  ("normal","trisomy","meiosis","monosomy","LOH")
                    noise=0.1) ## noise level

################################################################
####### 2. Convert Simulated Data to MADSEQ Format #############
################################################################
# 1. generate AAF table
AAF = data.frame(seqnames="chr1",
                 start = seq(1,data$nSNP),
                 end = seq(1,data$nSNP),
                 width = 1,
                 strand = "*",
                 REF = "C",
                 ALT = "A",
                 QUAL = 1,
                 FILTER= ".",
                 GT = "0/1",
                 DP = data$N,
                 Ref_D = data$N-data$z,
                 Alt_D = data$z)
write.table(AAF,file=paste0("simulation/LOH_F0.3_DP100_seed100_AAF.txt"),quote = F,row.names = F,col.names = T,sep="\t")

# 2. generate coverage table
cov = data.frame(seqnames = "chr1",
                 start = seq(1,data$nSites*100000,100000),
                 end = seq(100000,data$nSites*100000,100000),
                 width = 100000,
                 strand = "*",
                 depth = NA,
                 quantiled_depth = NA,
                 GC = runif(data$nSites,0,1),
                 normed_depth = data$N_cov,
                 ref_depth = data$m0)
write.table(cov,file=paste0("simulation/LOH_F0.3_DP100_seed100_Cov.txt"),quote = F,row.names = F,col.names = T,sep="\t")

###########################################################
####### 3. Run MADSEQ Model on Simulated Data #############
###########################################################
AAF = "simulation/LOH_F0.3_DP100_seed100_AAF.txt"
cov = "simulation/LOH_F0.3_DP100_seed100_Cov.txt"
res = runMadSeq(hetero = AAF,
                coverage = cov,
                target_chr = "chr1",
                adapt = 5000,
                burnin = 5000)
saveRDS(res,file=paste0("simulation/result/LOH_F0.3_DP100_seed100.RDS"))
