readBIC = function(file){
    if(file.exists(file)){
        tmp_chr = readRDS(file)
        if(is.null(tmp_chr)){
            tmp_BIC = c("BIC_normal"=0,"BIC_mitotic_trisomy"=1,"BIC_monosomy"=1,"BIC_LOH"=1,"BIC_meiotic_trisomy"=1)
        } 
        else{
            tmp_BIC = deltaBIC(tmp_chr)
            tmp_BIC = tmp_BIC[c("BIC_normal","BIC_mitotic_trisomy","BIC_monosomy","BIC_LOH","BIC_meiotic_trisomy")]
        }
    }
    else{
        tmp_BIC = c("BIC_normal"=0,"BIC_mitotic_trisomy"=1,"BIC_monosomy"=1,"BIC_LOH"=1,"BIC_meiotic_trisomy"=1)
    }
    tmp_BIC
}


plotBIC = function(data,ID){
    g = ggplot(data,aes(x=chr,y=deltaBIC,colour=chromosome_status,group=chromosome_status)) + geom_point(alpha=0.8) + geom_line()
    g = g + theme_bw() + scale_x_discrete(limits=CHR)
    g = g + scale_color_manual(values = c(normal="black",meiotic_trisomy="orange",monosomy="green4",LOH="magenta",mitotic_trisomy="blue"))
    g = g + labs(x="Chromosome",y="log(delta BIC)",title=ID)
    g = g + theme(legend.position="bottom",
                  axis.title = element_text(size=14,colour="black"),
                  axis.text.x = element_text(size=10,colour="black",angle=45,hjust = 1,vjust = 1),
                  axis.text.y = element_text(size=10,colour="black"),
                  title = element_text(size=14),
                  legend.text = element_text(size=12),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.major.x = element_line(size=0.3))
    g = g + geom_hline(yintercept = -log(10), linetype=2, colour="black",alpha=0.5)
    g
}