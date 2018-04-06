options(stringsAsFactors = F)
library(ggplot2)
library(grid)
library(gridExtra)
source("../Simulation/simulateData.R")
levels = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99) 
even_levels = c(0,0.2,0.4,0.6,0.8,0.99)## level of mosaicism we will simulate

################################################################
#########1.Plot Monosomy at Different Mosaic Levels#############
################################################################
## 1. simulate data
monosomy_AAF_df = NULL
monosomy_cov_df = NULL
for (tmp_level in levels){
    tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                            vartom=5,m=0.5,fraction=tmp_level,
                            type="monosomy",noise = 0.001)  ## simulate data at 1000 sites and 100X coverage for different mosaic levels
    tmp_AAF_df = data.frame(start=seq(1:tmp_data$nSNP),level=tmp_level,AAF=tmp_data$z/tmp_data$N) ## generate df for AAF
    monosomy_AAF_df = rbind(monosomy_AAF_df,tmp_AAF_df)
    
    tmp_cov_df = data.frame(start=seq(1:tmp_data$nSites),level=tmp_level,cov=mean(tmp_data$N_cov/tmp_data$m0)*(tmp_data$N_cov))
    monosomy_cov_df = rbind(monosomy_cov_df,tmp_cov_df)
}

## 2. convert the data to df to plot
## AAF
tmp_add = data.frame(start=-200,level=0,AAF=0.5)
monosomy_AAF_df = rbind(tmp_add,monosomy_AAF_df)
monosomy_AAF_df$group = ifelse(monosomy_AAF_df$level%in%even_levels,1,2)
monosomy_AAF_df$group = factor(monosomy_AAF_df$group)
monosomy_AAF_df$level = factor(monosomy_AAF_df$level,levels=levels)
levels(monosomy_AAF_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## coverage
tmp_add_cov = data.frame(start=-200,level=0,cov=100)
monosomy_cov_df = rbind(tmp_add_cov,monosomy_cov_df)
monosomy_cov_df$group = ifelse(monosomy_cov_df$level%in%even_levels,1,2)
monosomy_cov_df$group = factor(monosomy_cov_df$group)
monosomy_cov_df$level = factor(monosomy_cov_df$level,levels=levels)
levels(monosomy_cov_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## 3. plot AAF distribution
monosomy_AAF_g = ggplot(monosomy_AAF_df,aes(x=start,y=AAF,colour=group)) + geom_point(size=0.5) + ylim(0,1)
monosomy_AAF_g = monosomy_AAF_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
monosomy_AAF_g = monosomy_AAF_g + labs(x="",y="Alternative Allele Frequency")
monosomy_AAF_g = monosomy_AAF_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
monosomy_AAF_g = monosomy_AAF_g + theme(strip.text = element_text(size=10),
                                        axis.title.y = element_text(size=12),
                                        axis.title.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.text.y = element_text(size=10),
                                        plot.title = element_text(size=12,face="bold"),
                                        legend.position = "none",
                                        strip.background  = element_blank(),
                                        strip.placement = "outside")
## 4. plot Coverage distribution
monosomy_cov_g = ggplot(monosomy_cov_df,aes(x=start,y=cov,colour=group)) + geom_point(size=0.5) + ylim(-200,400)
monosomy_cov_g = monosomy_cov_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
monosomy_cov_g = monosomy_cov_g + geom_hline(yintercept = 100,colour="red",linetype=2)
monosomy_cov_g = monosomy_cov_g + labs(x="",y="Relative Sequencing Depth")
monosomy_cov_g = monosomy_cov_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
monosomy_cov_g = monosomy_cov_g + scale_y_continuous(limits=c(-200,400),breaks=c(50,100,150),labels = c(50,100,150))
monosomy_cov_g = monosomy_cov_g + theme(strip.text = element_text(size=10),
                      axis.title = element_text(size=12),
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      plot.title=element_blank(),
                      axis.text.y = element_text(size=10),
                      legend.position = "none",
                      strip.background  = element_blank(),
                      strip.placement = "outside")


#grid.arrange(monosomy_AAF_g,monosomy_cov_g,nrow=1)


#######################################################################
#########2.Plot Mitotic Trisomy at Different Mosaic Levels#############
#######################################################################
## 1. simulate data
trisomy_AAF_df = NULL
trisomy_cov_df = NULL
for (tmp_level in levels){
    tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                            vartom=5,m=0.5,fraction=tmp_level,
                            type="trisomy",noise = 0.001)  ## simulate data at 1000 sites and 100X coverage for different mosaic levels
    tmp_AAF_df = data.frame(start=seq(1:tmp_data$nSNP),level=tmp_level,AAF=tmp_data$z/tmp_data$N) ## generate df for AAF
    trisomy_AAF_df = rbind(trisomy_AAF_df,tmp_AAF_df)
    
    tmp_cov_df = data.frame(start=seq(1:tmp_data$nSites),level=tmp_level,cov=mean(tmp_data$N_cov/tmp_data$m0)*(tmp_data$N_cov))
    trisomy_cov_df = rbind(trisomy_cov_df,tmp_cov_df)
}

## 2. convert the data to df to plot
## AAF
tmp_add = data.frame(start=-200,level=0,AAF=0.5)
trisomy_AAF_df = rbind(tmp_add,trisomy_AAF_df)
trisomy_AAF_df$group = ifelse(trisomy_AAF_df$level%in%even_levels,1,2)
trisomy_AAF_df$group = factor(trisomy_AAF_df$group)
trisomy_AAF_df$level = factor(trisomy_AAF_df$level,levels=levels)
levels(trisomy_AAF_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## coverage
tmp_add_cov = data.frame(start=-200,level=0,cov=100)
trisomy_cov_df = rbind(tmp_add_cov,trisomy_cov_df)
trisomy_cov_df$group = ifelse(trisomy_cov_df$level%in%even_levels,1,2)
trisomy_cov_df$group = factor(trisomy_cov_df$group)
trisomy_cov_df$level = factor(trisomy_cov_df$level,levels=levels)
levels(trisomy_cov_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## 3. plot AAF distribution
trisomy_AAF_g = ggplot(trisomy_AAF_df,aes(x=start,y=AAF,colour=group)) + geom_point(size=0.5) + ylim(0,1)
trisomy_AAF_g = trisomy_AAF_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
trisomy_AAF_g = trisomy_AAF_g + labs(x="",y="Alternative Allele Frequency")
trisomy_AAF_g = trisomy_AAF_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
trisomy_AAF_g = trisomy_AAF_g + theme(strip.text = element_text(size=10),
                                        axis.title.y = element_text(size=12),
                                        axis.title.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.text.y = element_text(size=10),
                                        plot.title = element_text(size=12,face="bold"),
                                        legend.position = "none",
                                        strip.background  = element_blank(),
                                        strip.placement = "outside")
## 4. plot Coverage distribution
trisomy_cov_g = ggplot(trisomy_cov_df,aes(x=start,y=cov,colour=group)) + geom_point(size=0.5) 
trisomy_cov_g = trisomy_cov_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
trisomy_cov_g = trisomy_cov_g + geom_hline(yintercept = 100,colour="red",linetype=2)
trisomy_cov_g = trisomy_cov_g + labs(x="",y="Relative Sequencing Depth")
trisomy_cov_g = trisomy_cov_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
trisomy_cov_g = trisomy_cov_g + scale_y_continuous(limits=c(-200,400),breaks=c(50,100,150),labels = c(50,100,150))
trisomy_cov_g = trisomy_cov_g + theme(strip.text = element_text(size=10),
                                        axis.title = element_text(size=12),
                                        axis.ticks.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        plot.title=element_blank(),
                                        axis.text.y = element_text(size=10),
                                        legend.position = "none",
                                        strip.background  = element_blank(),
                                        strip.placement = "outside")

#grid.arrange(trisomy_AAF_g,trisomy_cov_g,nrow=1)

#######################################################################
#########3.Plot Meiotic Trisomy at Different Mosaic Levels#############
#######################################################################
## 1. simulate data
meiosis_AAF_df = NULL
meiosis_cov_df = NULL
for (tmp_level in levels){
    if(tmp_level==0){
        tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                                vartom=5,m=0.5,fraction=tmp_level,
                                type="normal",noise = 0.001)
    }
    else{
        tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                                vartom=5,m=0.5,fraction=tmp_level,
                                type="meiosis",noise = 0.001)  ## simulate data at 1000 sites and 100X coverage for different mosaic levels
    }
    tmp_AAF_df = data.frame(start=seq(1:tmp_data$nSNP),level=tmp_level,AAF=tmp_data$z/tmp_data$N) ## generate df for AAF
    meiosis_AAF_df = rbind(meiosis_AAF_df,tmp_AAF_df)
    
    tmp_cov_df = data.frame(start=seq(1:tmp_data$nSites),level=tmp_level,cov=mean(tmp_data$N_cov/tmp_data$m0)*(tmp_data$N_cov))
    meiosis_cov_df = rbind(meiosis_cov_df,tmp_cov_df)
}

## 2. convert the data to df to plot
## AAF
tmp_add = data.frame(start=-200,level=0,AAF=0.5)
meiosis_AAF_df = rbind(tmp_add,meiosis_AAF_df)
meiosis_AAF_df$group = ifelse(meiosis_AAF_df$level%in%even_levels,1,2)
meiosis_AAF_df$group = factor(meiosis_AAF_df$group)
meiosis_AAF_df$level = factor(meiosis_AAF_df$level,levels=levels)
levels(meiosis_AAF_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## coverage
tmp_add_cov = data.frame(start=-200,level=0,cov=100)
meiosis_cov_df = rbind(tmp_add_cov,meiosis_cov_df)
meiosis_cov_df$group = ifelse(meiosis_cov_df$level%in%even_levels,1,2)
meiosis_cov_df$group = factor(meiosis_cov_df$group)
meiosis_cov_df$level = factor(meiosis_cov_df$level,levels=levels)
levels(meiosis_cov_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## 3. plot AAF distribution
meiosis_AAF_g = ggplot(meiosis_AAF_df,aes(x=start,y=AAF,colour=group)) + geom_point(size=0.5) + ylim(0,1)
meiosis_AAF_g = meiosis_AAF_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
meiosis_AAF_g = meiosis_AAF_g + labs(x="",y="Alternative Allele Frequency")
meiosis_AAF_g = meiosis_AAF_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
meiosis_AAF_g = meiosis_AAF_g + theme(strip.text = element_text(size=10),
                                      axis.title.y = element_text(size=12),
                                      axis.title.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_text(size=10),
                                      plot.title = element_text(size=12,face="bold"),
                                      legend.position = "none",
                                      strip.background  = element_blank(),
                                      strip.placement = "outside")
## 4. plot Coverage distribution
meiosis_cov_g = ggplot(meiosis_cov_df,aes(x=start,y=cov,colour=group)) + geom_point(size=0.5) 
meiosis_cov_g = meiosis_cov_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
meiosis_cov_g = meiosis_cov_g + geom_hline(yintercept = 100,colour="red",linetype=2)
meiosis_cov_g = meiosis_cov_g + labs(x="",y="Relative Sequencing Depth")
meiosis_cov_g = meiosis_cov_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
meiosis_cov_g = meiosis_cov_g + scale_y_continuous(limits=c(-200,400),breaks=c(50,100,150),labels = c(50,100,150))
meiosis_cov_g = meiosis_cov_g + theme(strip.text = element_text(size=10),
                                      axis.title = element_text(size=12),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      plot.title=element_blank(),
                                      axis.text.y = element_text(size=10),
                                      legend.position = "none",
                                      strip.background  = element_blank(),
                                      strip.placement = "outside")

#grid.arrange(meiosis_AAF_g,meiosis_cov_g,nrow=1)


###########################################################
#########4.Plot LOH at Different Mosaic Levels#############
###########################################################
## 1. simulate data
LOH_AAF_df = NULL
LOH_cov_df = NULL
for (tmp_level in levels){
    if(tmp_level==0){
        tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                                vartom=5,m=0.5,fraction=tmp_level,
                                type="normal",noise = 0.001)
    }
    else{
        tmp_data = generateData(nSNP=1000,nCov=1000,coverage=100,
                                vartom=5,m=0.5,fraction=tmp_level,
                                type="LOH",noise = 0.001)  ## simulate data at 1000 sites and 100X coverage for different mosaic levels
    }
    tmp_AAF_df = data.frame(start=seq(1:tmp_data$nSNP),level=tmp_level,AAF=tmp_data$z/tmp_data$N) ## generate df for AAF
    LOH_AAF_df = rbind(LOH_AAF_df,tmp_AAF_df)
    
    tmp_cov_df = data.frame(start=seq(1:tmp_data$nSites),level=tmp_level,cov=mean(tmp_data$N_cov/tmp_data$m0)*(tmp_data$N_cov))
    LOH_cov_df = rbind(LOH_cov_df,tmp_cov_df)
}

## 2. convert the data to df to plot
## AAF
tmp_add = data.frame(start=-200,level=0,AAF=0.5)
LOH_AAF_df = rbind(tmp_add,LOH_AAF_df)
LOH_AAF_df$group = ifelse(LOH_AAF_df$level%in%even_levels,1,2)
LOH_AAF_df$group = factor(LOH_AAF_df$group)
LOH_AAF_df$level = factor(LOH_AAF_df$level,levels=levels)
levels(LOH_AAF_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## coverage
tmp_add_cov = data.frame(start=-200,level=0,cov=100)
LOH_cov_df = rbind(tmp_add_cov,LOH_cov_df)
LOH_cov_df$group = ifelse(LOH_cov_df$level%in%even_levels,1,2)
LOH_cov_df$group = factor(LOH_cov_df$group)
LOH_cov_df$level = factor(LOH_cov_df$level,levels=levels)
levels(LOH_cov_df$level) = c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")

## 3. plot AAF distribution
LOH_AAF_g = ggplot(LOH_AAF_df,aes(x=start,y=AAF,colour=group)) + geom_point(size=0.5) + ylim(0,1)
LOH_AAF_g = LOH_AAF_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
LOH_AAF_g = LOH_AAF_g + labs(x="",y="Alternative Allele Frequency")
LOH_AAF_g = LOH_AAF_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
LOH_AAF_g = LOH_AAF_g + theme(strip.text = element_text(size=10),
                                      axis.title.y = element_text(size=12),
                                      axis.title.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_text(size=10),
                                      plot.title = element_text(size=12,face="bold"),
                                      legend.position = "none",
                                      strip.background  = element_blank(),
                                      strip.placement = "outside")
## 4. plot Coverage distribution
LOH_cov_g = ggplot(LOH_cov_df,aes(x=start,y=cov,colour=group)) + geom_point(size=0.5) 
LOH_cov_g = LOH_cov_g + facet_wrap(~level,switch = "x",scales="free_x",nrow=1) + theme_classic()
LOH_cov_g = LOH_cov_g + geom_hline(yintercept = 100,colour="red",linetype=2)
LOH_cov_g = LOH_cov_g + labs(x="",y="Relative Sequencing Depth")
LOH_cov_g = LOH_cov_g + scale_colour_manual(values=c("1"="blue","2"="cornflowerblue"))
LOH_cov_g = LOH_cov_g + scale_y_continuous(limits=c(-200,400),breaks=c(50,100,150),labels = c(50,100,150))
LOH_cov_g = LOH_cov_g + theme(strip.text = element_text(size=10),
                                      axis.title = element_text(size=12),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      plot.title=element_blank(),
                                      axis.text.y = element_text(size=10),
                                      legend.position = "none",
                                      strip.background  = element_blank(),
                                      strip.placement = "outside")
quartz(w=14,h=10)
grid.arrange(meiosis_AAF_g,meiosis_cov_g,
             trisomy_AAF_g,trisomy_cov_g,
             monosomy_AAF_g,monosomy_cov_g,
             LOH_AAF_g,LOH_cov_g,nrow=4,ncol=2)
