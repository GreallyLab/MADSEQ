even = c("2","4","6","8","10","12","14","16","18","20","22")
CHR = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

plotAAF = function(hetero){
    tmp_aaf = read.table(hetero,sep="\t",header=T)
    tmp_add = data.frame(seqnames=1,start=-30000000,end=1,width=1,strand="*", paramRangeID=NA ,REF="A", ALT="G",
                         QUAL="190", FILTER=".",GT="0/1", DP=10, Ref_D=5 ,Alt_D=5,AAF=0.5, binned_AAF1=NA, binned_AAF2=NA)
    tmp_aaf = rbind(tmp_add,tmp_aaf)
    tmp_aaf$seqnames = gsub("chr","",tmp_aaf$seqnames)
    tmp_aaf$group = ifelse(tmp_aaf$seqnames%in%even,2,1)
    tmp_aaf$group = factor(tmp_aaf$group)
    tmp_aaf$seqnames = factor(tmp_aaf$seqnames,levels=CHR)
    
    g = ggplot(tmp_aaf,aes(x=start,y=AAF,colour=group)) + geom_point(size=0.5) + ylim(0,1)
    g = g + facet_wrap(~seqnames,switch = "x",scales="free_x",nrow=1) + theme_classic()
    g = g + labs(x="Chromosome",y="Alternative Allele Frequency")
    g = g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
    g = g + theme(strip.text = element_text(size=14),
                  axis.title.y = element_text(size=16),
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size=14),
                  plot.title = element_text(size=16,face="bold"),
                  legend.position = "none",
                  strip.background  = element_blank(),
                  strip.placement = "outside")
    g
}


plotCov = function(cov,bin_size=20){
    tmp_cov = read.table(cov,sep="\t",header=T)
    tmp_add_cov = data.frame(seqnames=1,start=-30000000,end=1,width=1,strand="*", depth=10,quantiled_depth=10,GC=0.5, normed_depth=tmp_cov$ref_depth[1], ref_depth=51)
    tmp_cov = rbind(tmp_add_cov,tmp_cov)
    tmp_cov$seqnames = gsub("chr","",tmp_cov$seqnames)
    ## bin tmp_cov
    data_bin = NULL
    for (j in 1:24){
        chr = j
        if(j==23) chr="X"
        if(j==24) chr="Y"
        #print(chr)
        tmp = tmp_cov[tmp_cov$seqnames==chr,]
        if(nrow(tmp)==0)next
        for (i in seq(1,nrow(tmp),bin_size)){
            if((i+bin_size-1)<nrow(tmp)){
                start = tmp[i,"start"]
                end = tmp[i+bin_size-1,"end"]
                bin_normed_depth = mean(tmp[i:(i+bin_size-1),"normed_depth"])
            }
            else{
                start = tmp[i,"start"]
                end = tmp[nrow(tmp),"end"]
                bin_normed_depth = mean(tmp[i:nrow(tmp),"normed_depth"])
            }
            data_bin = rbind(data_bin,data.frame(chr,start,end,bin_normed_depth))
        }
    }
    names(data_bin) = c("seqnames","start","end","normed_depth")
    data_bin = data_bin[!is.na(data_bin$normed_depth),]
    data_bin$group = ifelse(data_bin$seqnames%in%even,2,1)
    data_bin$group = factor(data_bin$group)
    data_bin$seqnames = factor(data_bin$seqnames,levels=CHR)
    data_mean = data_bin%>%group_by(seqnames)%>%summarise(mean=mean(normed_depth))
    data_mean = as.data.frame(data_mean)
    
    g_cov = ggplot(data_bin,aes(x=start,y=normed_depth,colour=group)) + geom_point(size=0.5) + ylim(-50,200)
    g_cov = g_cov + facet_wrap(~seqnames,switch = "x",scales="free_x",nrow=1) + theme_classic()
    g_cov = g_cov + geom_hline(data=data_mean,aes(yintercept = mean),colour="red",size=1)
    g_cov = g_cov + labs(x="Chromosome",y="Normalized Depth",title="")
    g_cov = g_cov + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
    g_cov = g_cov + theme(strip.text = element_text(size=14),
                          axis.title = element_text(size=16),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_blank(),
                          plot.title=element_blank(),
                          axis.text.y = element_text(size=14),
                          legend.position = "none",
                          strip.background  = element_blank(),
                          strip.placement = "outside")
    g_cov
}

