options(stringsAsFactors = F)

################################################################################
#########1.Prepare Number of Heterozygous Sites for each Population#############
################################################################################



CHR=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
      "chr18","chr19","chr20","chr21","chr22","chrX","chrY")
#-----------DATA-------------------
info = read.table("1000G.sample.info.txt",sep="\t",header=T) ## read population and sex information of 1000G samples (or any cohort of your samples)
hetero_mat = readRDS(file="hetero_count_chr.RDS")  ## read the heterozygous count matrix we generated before (Violin(Fig5).R)

#---------------ALL CHROM----------------
ALL = NULL
for (pop in POP){
    print(pop)
    ## subset IDs in this population
    tmp_pop_mat = hetero_mat[info[info$population==pop,"ID"],]
    all = as.vector(tmp_pop_mat)
    assign(paste0(pop,"_all"),all)
    ALL = c(ALL,all)
}
#----------------PLOT-----------------
plot(1:26,y=rep(50000,26),type="n",xaxt="n",ylim=c(0,3500),xlab="",ylab="number of heterozygous sites/chromosome",yaxt="n",bty="l",main="Common SNPs Design")
axis(2,at=seq(0,3500,500),cex.axis=0.8,hadj=0.5,tck=0.02)

points(seq(0.7,1.3,length.out = length(FIN_all)),sort(FIN_all),pch=3,cex=0.5,col="cyan4")
points(1,mean(FIN_all),pch=20,col="black",lwd=2,cex=0.7)
text(1,max(FIN_all)+50,"FIN",cex=0.8)

points(seq(1.7,2.3,length.out = length(GBR_all)),sort(GBR_all),pch=3,cex=0.5,col="cyan2")
points(1.95,mean(GBR_all),pch=20,col="black",lwd=2,cex=0.7)
text(2,max(GBR_all)+50,"GBR",cex=0.8)

points(seq(2.7,3.3,length.out = length(CEU_all)),sort(CEU_all),pch=3,cex=0.5,col="dodgerblue3")
points(3,mean(CEU_all),pch=20,col="black",lwd=2,cex=0.7)
text(3,max(CEU_all)+50,"CEU",cex=0.8)

points(seq(3.7,4.3,length.out = length(IBS_all)),sort(IBS_all),pch=3,cex=0.5,col="cornflowerblue")
points(3.98,mean(IBS_all),pch=20,col="black",lwd=2,cex=0.7)
text(4,max(IBS_all)+50,"IBS",cex=0.8)

points(seq(4.7,5.3,length.out = length(TSI_all)),sort(TSI_all),pch=3,cex=0.5,col="blue4")
points(4.95,mean(TSI_all),pch=20,col="black",lwd=2,cex=0.7)
text(5,max(TSI_all)+50,"TSI",cex=0.8)

points(seq(5.7,6.3,length.out = length(CHS_all)),sort(CHS_all),pch=3,cex=0.5,col="limegreen")
points(6,mean(CHS_all),pch=20,col="black",lwd=2,cex=0.7)
text(6,max(CHS_all)+50,"CHS",cex=0.8)

points(seq(6.7,7.3,length.out = length(CDX_all)),sort(CDX_all),pch=3,cex=0.5,col="green4")
points(7,mean(CDX_all),pch=20,col="black",lwd=2,cex=0.7)
text(7,max(CDX_all)+50,"CDX",cex=0.8)

points(seq(7.7,8.3,length.out = length(CHB_all)),sort(CHB_all),pch=3,cex=0.5,col="olivedrab3")
points(8,mean(CHB_all),pch=20,col="black",lwd=2,cex=0.7)
text(8,max(CHB_all)+50,"CHB",cex=0.8)

points(seq(8.7,9.3,length.out = length(JPT_all)),sort(JPT_all),pch=3,cex=0.5,col="seagreen4")
points(9,mean(JPT_all),pch=20,col="black",lwd=2,cex=0.7)
text(9,max(JPT_all)+50,"JPT",cex=0.8)

points(seq(9.7,10.3,length.out = length(KHV_all)),sort(KHV_all),pch=3,cex=0.5,col="springgreen2")
points(10,mean(KHV_all),pch=20,col="black",lwd=2,cex=0.7)
text(10,max(KHV_all)+50,"KHV",cex=0.8)

points(seq(10.7,11.3,length.out = length(GIH_all)),sort(GIH_all),pch=3,cex=0.5,col="purple4")
points(11,mean(GIH_all),pch=20,col="black",lwd=2,cex=0.7)
text(11,max(GIH_all)+50,"GIH",cex=0.8)

points(seq(11.7,12.3,length.out = length(STU_all)),sort(STU_all),pch=3,cex=0.5,col="orchid4")
points(12,mean(STU_all),pch=20,col="black",lwd=2,cex=0.7)
text(12,max(STU_all)+50,"STU",cex=0.8)

points(seq(12.7,13.3,length.out = length(PJL_all)),sort(PJL_all),pch=3,cex=0.5,col="violetred")
points(13,mean(PJL_all),pch=20,col="black",lwd=2,cex=0.7)
text(13,max(PJL_all)+50,"PJL",cex=0.8)

points(seq(13.7,14.3,length.out = length(ITU_all)),sort(ITU_all),pch=3,cex=0.5,col="maroon4")
points(14,mean(ITU_all),pch=20,col="black",lwd=2,cex=0.7)
text(14,max(ITU_all)+50,"ITU",cex=0.8)

points(seq(14.7,15.3,length.out = length(BEB_all)),sort(BEB_all),pch=3,cex=0.5,col="darkmagenta")
points(15,mean(BEB_all),pch=20,col="black",lwd=2,cex=0.7)
text(15,max(BEB_all)+50,"BEB",cex=0.8)

points(seq(15.7,16.3,length.out = length(PEL_all)),sort(PEL_all),pch=3,cex=0.5,col="firebrick1")
points(16,mean(PEL_all),pch=20,col="black",lwd=2,cex=0.7)
text(16,max(PEL_all)+50,"PEL",cex=0.8)

points(seq(16.7,17.3,length.out = length(MXL_all)),sort(MXL_all),pch=3,cex=0.5,col="red")
points(17,mean(MXL_all),pch=20,col="black",lwd=2,cex=0.7)
text(17,max(MXL_all)+50,"MXL",cex=0.8)

points(seq(17.7,18.3,length.out = length(CLM_all)),sort(CLM_all),pch=3,cex=0.5,col="firebrick")
points(18,mean(CLM_all),pch=20,col="black",lwd=2,cex=0.7)
text(18,max(CLM_all)+50,"CLM",cex=0.8)

points(seq(18.7,19.3,length.out = length(PUR_all)),sort(PUR_all),pch=3,cex=0.5,col="orangered3")
points(19,mean(PUR_all),pch=20,col="black",lwd=2,cex=0.7)
text(19,max(PUR_all)+50,"PUR",cex=0.8)

points(seq(19.7,20.3,length.out = length(ASW_all)),sort(ASW_all),pch=3,cex=0.5,col="orangered")
points(20,mean(ASW_all),pch=20,col="black",lwd=2,cex=0.7)
text(20,max(ASW_all)+50,"ASW",cex=0.8)

points(seq(20.7,21.3,length.out = length(ACB_all)),sort(ACB_all),pch=3,cex=0.5,col="orange")
points(21,mean(ACB_all),pch=20,col="black",lwd=2,cex=0.7)
text(21,max(ACB_all)+50,"ACB",cex=0.8)

points(seq(21.7,22.3,length.out = length(GWD_all)),sort(GWD_all),pch=3,cex=0.5,col="darkorange")
points(22,mean(GWD_all),pch=20,col="black",lwd=2,cex=0.7)
text(22,max(GWD_all)+50,"GWD",cex=0.8)

points(seq(22.7,23.3,length.out = length(YRI_all)),sort(YRI_all),pch=3,cex=0.5,col="orange2")
points(23,mean(YRI_all),pch=20,col="black",lwd=2,cex=0.7)
text(23,max(YRI_all)+50,"YRI",cex=0.8)

points(seq(23.7,24.3,length.out = length(LWK_all)),sort(LWK_all),pch=3,cex=0.5,col="yellow4")
points(24,mean(LWK_all),pch=20,col="black",lwd=2,cex=0.7)
text(24,max(LWK_all)+50,"LWK",cex=0.8)

points(seq(24.7,25.3,length.out = length(ESN_all)),sort(ESN_all),pch=3,cex=0.5,col="gold")
points(25,mean(ESN_all),pch=20,col="black",lwd=2,cex=0.7)
text(25,max(ESN_all)+50,"ESN",cex=0.8)

points(seq(25.7,26.3,length.out = length(MSL_all)),sort(MSL_all),pch=3,cex=0.5,col="gold3")
points(26,mean(MSL_all),pch=20,col="black",lwd=2,cex=0.7)
text(26,max(MSL_all)+50,"MSL",cex=0.8)

grid(nx=NA,ny=NULL)
abline(h=1000,col="red",lty=3,lwd=2)
mtext(paste0("total chromosomes: ",length(ALL)), side = 1)
mtext(paste0("chromosomes with heterozygous sites <1000: ",sum(ALL<1000), " (",sum(ALL<1000)/length(ALL),"%)"),side=1,line=1)