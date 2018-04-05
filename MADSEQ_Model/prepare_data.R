library("GenomicRanges")
library("MADSEQ")
library("BSgenome")
library("BSgenome.Hsapiens.1000genomes.hs37d5")
options(stringsAsFactors = F)
#####################################


######################################################
####### 1. Prepare Coverage from Bam Files ###########
######################################################
## ID/id here is the samples you have

## specify the position of your bam file and the bed files containing captured regions
bam.file = paste0("bam/",ID,".merged.processed.bam")
target.file = "MadSeq_capture_targets.bed"

## prepare coverage 
cov = prepareCoverageGC(target_bed = target.file,bam=bam.file)
saveRDS(cov,file=paste0("coverage/",ID,"_cov.RDS"))


#####################################################################
####### 2. Normalize Coverage from Cohorts of Individuals ###########
#####################################################################
## load the coverage prepared before
for (id in all_id){
    print(id)
    tmp = readRDS(paste0("coverage/",id,"_cov.RDS"))
    assign(id,tmp)
}

## normalize female and male separately, and plot the normalization figures
## 1. female

pdf(file="plot/normalize_coverage_female.pdf",width=8,height=12)
normalizeCoverage(id1, ## here id1/id2/id3 is the name you assigned for your coverage file, each row is one sample 
                  id2,
                  id3,
                  writeToFile = T,
                  destination = "normalized_cov/")  ## where to put normalized coverage data
dev.off()

## 2. male
pdf(file="plot/normalize_coverage_male.pdf",width=8,height=12)
normalizeCoverage(id4,
                  id5,
                  id6,
                  writeToFile = T,
                  destination = "normalized_cov/")
dev.off()

###########################################################################
####### 3. Prepare Heterozygous Sites from genotyping VCF files ###########
###########################################################################
## Note: the vcf file here is one vcf for each individual, the vcf file has to contain GT and DP field
vcf.file = paste0("vcf/",ID,".raw.variants.vcf.gz")  ## location of your vcf file
prepareHetero(vcffile = vcf.file,
              target_bed = target.file,  ## location of your bed file, the same as step 1
              genome = "BSgenome.Hsapiens.1000genomes.hs37d5",  ## which genome assembly your bam/vcf file used
              writeToFile = T,
              destination = "heterozygous/",  ## where to put processed heterozygous data
              plot = T)  ## if you want to generate a plot of alternative allele frequency before and after processing, set it to TRUE


## Now you have all the data ready, please go to runModel.R to find the script to run the model.