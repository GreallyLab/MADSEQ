library("ggplot2")
library("gridExtra")
library("reshape2")
library("grid")
library("dplyr")
options(stringsAsFactors = F)

##############################################################################
#########1.Prepare Coverage and Number of Heterozygous Sites Info#############
##############################################################################
id = read.table("1000G.id.txt",sep="\t",header=F)[,1]  ## read all the IDs of 1000G (or any cohort of your samples)
info = read.table("1000G.sample.info.txt",sep="\t",header=T) ## read population and sex information of 1000G samples (or any cohort of your samples)
chr = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
## 1. prepare number of heterozygous for each chromosome
hetero = mclapply(id,function(x){
    res = read.table(paste0("heterozygous/",x,".raw.variants.vcf.gz_filtered_heterozygous.txt"),sep="\t",header=T)  ## read processed heterozygous files
    tmp_chr=table(res$seqnames)[chr]  ## count number of heterozygous sites for each chromosome
    names(tmp_chr) = chr
    tmp_chr},
    mc.cores = 10)

## heterozygous matrix: each row is one individual, each column is one chromosome
hetero_mat = do.call(rbind,hetero) ## convert list to matrix
rownames(hetero_mat) = id

## 2. prepare average sequencing depth for each chromosome
cov = mclapply(id,function(x){
    tmp_cov = read.table(paste0("normalized_cov/",x,"_normed_depth.txt"),sep="\t",header=T) ## read normalized coverage files
    tmp_cov_chr = tmp_cov %>% group_by(seqnames) %>% summarise(depth_mean=mean(depth))  ## calculate average depth for each chromosomes
    tmp_cov_chr = as.data.frame(tmp_cov_chr)
    rownames(tmp_cov_chr) = tmp_cov_chr$seqnames
    tmp_cov_chr = tmp_cov_chr[chr,"depth_mean"]
    names(tmp_cov_chr) = chr
    tmp_cov_chr},
    mc.cores = 10)
## coverage matrix: each row is one individual, each column is one chromosome
cov_mat = do.call(rbind,cov)
rownames(cov_mat) = id

## save matrix
saveRDS(cov_mat,file="coverage_chr.RDS")
saveRDS(hetero_mat,file="hetero_count_chr.RDS")

################################################
#########2.Prepare Matrix for Plot #############
################################################
# heterozygous data
hetero_mat = ifelse(is.na(hetero_mat),0,hetero_mat)
colnames(hetero_mat) = seq(1,24)
hetero_df = melt(hetero_mat)
colnames(hetero_df) = c("ID","Chr","Count")
hetero_df$Sex = info[match(hetero_df$ID,info$ID),"sex"]  ## get sex info
hetero_df$Pop = info[match(hetero_df$ID,info$ID),"population"]  ## get population info

pdat = hetero_df%>%group_by(Sex,Chr)%>%do(data.frame(loc=density(.$Count)$x,dens=density(.$Count)$y))
pdat_female = pdat[pdat$Sex=="female",]
pdat_male = pdat[pdat$Sex=="male",]
female_max = sapply(split(pdat_female,f=pdat_female$Chr),function(x)max(x$dens))
male_max = sapply(split(pdat_male,f=pdat_male$Chr),function(x)max(x$dens))
female_scale = 0.7/female_max
male_scale = 0.7/male_max
for (i in 1:24){
    pdat$dens= ifelse(pdat$Sex=="female"&pdat$Chr==i,pdat$dens*female_scale[i],pdat$dens)
    pdat$dens= ifelse(pdat$Sex=="male"&pdat$Chr==i,pdat$dens*male_scale[i],pdat$dens)
}
pdat$dens = ifelse(pdat$Sex=="female",pdat$dens*-1,pdat$dens)
pdat$dens = pdat$dens+2*pdat$Chr
pdat_df = as.data.frame(pdat)

# coverage data
cov_mat = ifelse(is.na(cov_mat),0,cov_mat)
colnames(cov_mat) = seq(1,24)
cov_df = melt(cov_mat)
colnames(cov_df) = c("ID","Chr","Depth")
cov_df$Sex = info[match(cov_df$ID,info$ID),"sex"]  ## get sex info
cov_df$Pop = info[match(cov_df$ID,info$ID),"population"]  ## get population info

pdat_cov = cov_df%>%group_by(Sex,Chr)%>%do(data.frame(loc=density(.$Depth)$x,dens=density(.$Depth)$y))
pdat_cov_female = pdat_cov[pdat_cov$Sex=="female",]
pdat_cov_male = pdat_cov[pdat_cov$Sex=="male",]
female_max = sapply(split(pdat_cov_female,f=pdat_cov_female$Chr),function(x)max(x$dens))
male_max = sapply(split(pdat_cov_male,f=pdat_cov_male$Chr),function(x)max(x$dens))
female_scale = 0.7/female_max
male_scale = 0.7/male_max
for (i in 1:24){
    pdat_cov$dens= ifelse(pdat_cov$Sex=="female"&pdat_cov$Chr==i,pdat_cov$dens*female_scale[i],pdat_cov$dens)
    pdat_cov$dens= ifelse(pdat_cov$Sex=="male"&pdat_cov$Chr==i,pdat_cov$dens*male_scale[i],pdat_cov$dens)
}
pdat_cov$dens = ifelse(pdat_cov$Sex=="female",pdat_cov$dens*-1,pdat_cov$dens)
pdat_cov$dens = pdat_cov$dens+2*pdat_cov$Chr
pdat_cov_df = as.data.frame(pdat_cov)

#######################################################
#########3.Plot Split Violin Plot (Fig 5) #############
#######################################################
## 1. Number of Hetero
split_g = ggplot(pdat_df, aes(dens, loc,fill=Sex,group=interaction(Sex,Chr)) )+geom_polygon(alpha=0.35,colour="gray40") 
split_g = split_g + geom_point(data=hetero_df,aes(x=Chr*2,y=Count,colour=Pop),size=0.1,alpha=0.5,position=position_jitter(width = 0.1))
split_g = split_g + scale_colour_manual(values = c(EUR="blue",EAS="green4",AMR="red",AFR="orange",SAS="magenta","NA"="gray26"))
split_g = split_g + theme_classic()
split_g = split_g + guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))
split_g = split_g + labs(y="number of heterozygous sites",title="1000G WXS Data",x="Chromosome")
split_g = split_g + scale_y_continuous(breaks = c(0,2000,4000,6000),labels = c(0,2000,4000,6000),limits=c(0,10000))
split_g = split_g + scale_x_continuous(breaks = seq(2,48,2),labels = CHR,limits=c(1,49))
split_g = split_g + theme(legend.position = "bottom",
                          axis.title = element_text(size=16,colour="black"),
                          axis.text = element_text(size=14,colour="black"),
                          title = element_text(size=16))
split_g = split_g + geom_hline(yintercept = 1000,colour="red",linetype="dashed")
#split_g

## Ave. Depth
split_g2 = ggplot(pdat_cov, aes(dens,loc,fill=Sex,group=interaction(Sex,Chr)))+geom_polygon(alpha=0.35,colour="gray40")
split_g2 = split_g2 + geom_point(data=cov_df,aes(x=2*Chr,y=Depth,colour=Affected),size=0.1,alpha=0.5,position=position_jitter(width = 0.1))
split_g2 = split_g2 + theme_bw()
split_g2 = split_g2 + guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))
split_g2 = split_g2 + scale_colour_manual(values = c(EUR="blue",EAS="green4",AMR="red",AFR="orange",SAS="magenta","NA"="gray26"))
split_g2 = split_g2 + scale_x_continuous(breaks = seq(2,48,2),labels = CHR,limits=c(1,49))
split_g2 = split_g2 + scale_y_continuous(breaks = c(50,100),labels = c(50,100),limits=c(-400,200))
split_g2 = split_g2 + theme(legend.position = "bottom",
                            axis.title = element_text(size=16,colour="black"),
                            axis.text = element_text(size=14,colour="black"),
                            #axis.ticks.y = element_blank(),
                            title = element_text(size=16))
split_g2 = split_g2 + theme(panel.grid=element_blank()) + theme(panel.background = element_rect(fill = NA))
#split_g2
source("Double_Axis_Graph.R")
plot(double_axis_graph(split_g,split_g2))