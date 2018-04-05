options(stringsAsFactors = F)

########################## function1: simulate data ##########################
generateData = function(nSNP,nCov,coverage,vartom,m,fraction,type="trisomy",noise=0.01){
    ## first case: mitotic trisomy
    if (type == "trisomy"){
        
        # generate the coverage for trisomy chromosome, the mean coverage should be (1+f)*coverage
        # and the coverage will be described by a negative binomial distribution
        cov_mean = coverage
        cov_var = cov_mean*vartom
        # use the mean and variance of coverage to calculate r and p of negative binomial model
        p = cov_mean/cov_var
        r = cov_mean^2/(cov_var-cov_mean)
        
        # the prob (mean BAF) for the two mixtures: theta1 and theta2
        # which could be calculated by the m and fraction of aneuploidy cells
        theta1 = m*(1+fraction)/(1+m*fraction)
        theta2 = m/(1+fraction-m*fraction)
        # generate the depth for B allele according to the coverage and BAF
        alt_reads = NULL
        depth = NULL
        for (i in 1:((1-noise)*nSNP)){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            # make the weight for the two mixture equals to 0.5
            random = runif(1)
            theta = ifelse(random>=0.5,theta1,theta2)
            # generate the depth for alternative allele by binomial distribution
            sub_alt = rbinom(1,sub_depth,theta)
            # add this while loop to make sure all the data points are heterozygous site
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        # 1% of the points are errors
        for (i in 1:(noise*nSNP)){
            sub_depth = rnbinom(1,size=r,prob=p)
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            while(sub_alt==sub_depth){
                sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        # randomly generate coverage information besides the heterozygous sites
        N_cov = rnbinom((1-noise)*nCov,size=r,prob=p)
        # generate noise% of the coverage as outliers
        N_cov_out = round(runif(noise*nCov,min=1,max=max(N_cov)))
        N_cov = c(N_cov[N_cov!=0],N_cov_out)
        # randomly generate m0 from a poisson distribution to account for the randomness of average coverage for chrom, however the variation is small
        coverage_mean = mean(N_cov)
        print(coverage_mean)
        m0 = round(rnorm(1,2*coverage_mean/(fraction+2),0.05*coverage),2)
        while(m0>coverage_mean|(2*coverage_mean/m0-2)>1){
            m0 = round(rnorm(1,2*coverage_mean/(fraction+2),0.05*coverage),2)
        }
        if (fraction==0) m0 = rpois(1,coverage_mean)
        print(m0)
        dataList = cbind(alt_reads,depth)
    }
    
    ## second case: meiosis aneuploidy
    if (type == "meiosis"){
        # generate the coverage for trisomy chromosome, the mean coverage should be (1+f)*coverage
        # and the coverage will be described by a negative binomial distribution
        cov_mean = coverage
        cov_var = cov_mean*vartom
        # use the mean and variance of coverage to calculate r and p of negative binomial model
        p = cov_mean/cov_var
        r = cov_mean^2/(cov_var-cov_mean)
        
        # meiosis aneuploidy will have 4 mixtures, plus one outlier component
        theta1 = m*(1+fraction)/(1+m*fraction)
        theta2 = m/(1+fraction-m*fraction)
        theta3 = m*fraction/(2-2*m+m*fraction)
        theta4 = 2*m/(2*m+fraction-m*fraction)
        alt_reads = NULL
        depth = NULL
        
        for (i in 1:((1-noise)*nSNP)){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            
            # set the weight for each mixture as cut1 and cut2
            random = runif(1)
            if (random < 0.165) theta = theta3
            else if (random<0.5 & random>=0.165) theta = theta2
            else if (random<0.835 & random >= 0.5) theta = theta1
            else theta = theta4
            sub_alt = rbinom(1,sub_depth,theta)
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        for (i in 1:(noise*nSNP)){
            sub_depth = rnbinom(1,size=r,prob=p)
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            while(sub_alt==sub_depth){
                sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        
        # randomly generate coverage information besides the heterozygous sites
        N_cov = rnbinom((1-noise)*nCov,size=r,prob=p)
        # generate noise% of the coverage as outliers
        N_cov_out = round(runif(noise*nCov,min=1,max=max(N_cov)))
        N_cov = c(N_cov[N_cov!=0],N_cov_out)
        coverage_mean = mean(N_cov)
        print(coverage_mean)
        m0 = round(rnorm(1,2*coverage_mean/(fraction+2),0.05*coverage),2)
        while(m0>coverage_mean|(2*coverage_mean/m0-2)>1){
            m0 = round(rnorm(1,2*coverage_mean/(fraction+2),0.05*coverage),2)
        }
        if (fraction==0) m0 = rpois(1,coverage_mean)
        print(m0)
        dataList = cbind(alt_reads,depth)
    }
    
    ## third case: mosaic monosomy
    if (type == "monosomy"){
        cov_mean = coverage
        cov_var = cov_mean*vartom
        # use the mean and variance of coverage to calculate r and p of negative binomial model
        p = cov_mean/cov_var
        r = cov_mean^2/(cov_var-cov_mean)
        
        # there are two mixtures in mitotic monosomy, 
        # and one more mixture is the outlier component which account for 0.1% of all sites
        
        # the prob (mean BAF) for the two mixtures: theta1 and theta2
        # which could be calculated by the m and fraction of aneuploidy cells
        theta1 = m/(1-fraction+m*fraction)
        theta2 = m*(1-fraction)/(1-m*fraction)
        # generate the depth for B allele according to the coverage and BAF
        alt_reads = NULL
        depth = NULL
        for (i in 1:((1-noise)*nSNP)){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            # also depth have to greater than 1 to make sure site is heterozygous
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            
            # make the weight for the two mixture equals to 0.5
            random = runif(1)
            theta = ifelse(random>=0.5,theta1,theta2)
            # generate the depth for alternative allele by binomial distribution
            sub_alt = rbinom(1,sub_depth,theta)
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        for (i in 1:(noise*nSNP)){
            sub_depth = rnbinom(1,size=r,prob=p)
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            while(sub_alt==sub_depth){
                sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        
        # randomly generate coverage information besides the heterozygous sites
        N_cov = rnbinom((1-noise)*nCov,size=r,prob=p)
        # generate noise% of the coverage as outliers
        N_cov_out = round(runif(noise*nCov,min=1,max=max(N_cov)))
        N_cov = c(N_cov[N_cov!=0],N_cov_out)
        coverage_mean = mean(N_cov)
        print(coverage_mean)
        m0 = round(rnorm(1,2*coverage_mean/(2-fraction),0.05*coverage),2)
        while(m0<coverage_mean|(2-2*coverage_mean/m0)<0){
            m0 = round(rnorm(1,2*coverage_mean/(2-fraction),0.05*coverage),2)
        }
        if (fraction==0) m0 = rpois(1,coverage_mean)
        print(m0)
        dataList = cbind(alt_reads,depth)
    }
    
    #4. LOH
    if (type=="LOH"){
        cov_mean = coverage
        cov_var = cov_mean*vartom
        # use the mean and variance of coverage to calculate r and p of negative binomial model
        p = cov_mean/cov_var
        r = cov_mean^2/(cov_var-cov_mean)
        
        ######## AAF ##########
        # the prob (mean BAF) for the two mixtures: theta1 and theta2
        # which could be calculated by the m and fraction of aneuploidy cells
        theta1 = m*(1+fraction)/(1-fraction+2*m*fraction)
        theta2 = m*(1-fraction)/(fraction-2*m*fraction+1)
        
        ### only part of the reads are UPD, randomly generate the percentage of region being UPD
        region = runif(1,min=0.1,max=1)
        
        # 1. generate B allele for the sites that located in UPD region
        ab_alt_reads = NULL
        ab_depth = NULL
        alt_reads = NULL
        depth = NULL
        for (i in 1:(round(region*nSNP))){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            ab_depth = c(ab_depth,sub_depth)
            # make the weight for the two mixture equals to 0.5
            random = runif(1)
            theta = ifelse(random>=0.5,theta1,theta2)
            # generate the depth for alternative allele by binomial distribution
            sub_alt = rbinom(1,sub_depth,theta)
            # add this while loop to make sure all the data points are heterozygous site
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            ab_alt_reads= c(ab_alt_reads,sub_alt)
        }
        
        # 2. generate B allele for the rest of the regions, except % of noise
        for (i in 1:(round((1-region-noise)*nSNP))){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            # for normal region, theta equals to m
            theta = m
            # generate the depth for alternative allele by binomial distribution
            sub_alt = rbinom(1,sub_depth,theta)
            # add this while loop to make sure all the data points are heterozygous site
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        
        # 3. genarate noise
        for (i in 1:(round(noise*nSNP))){
            sub_depth = rnbinom(1,size=r,prob=p)
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            while(sub_alt==sub_depth){
                sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        ######### COVERAGE ###############
        # randomly generate coverage information besides the heterozygous sites
        N_cov = rnbinom((1-noise)*nCov,size=r,prob=p)
        # generate noise% of the coverage as outliers
        N_cov_out = round(runif(noise*nCov,min=1,max=max(N_cov)))
        N_cov = c(N_cov[N_cov!=0],N_cov_out)
        # randomly generate m0 from a poisson distribution to account for the randomness of average coverage for chrom, however the variation is small
        m0 = round(rnorm(1,mean(N_cov),0.05*coverage),2)
        dataList = cbind(alt_reads,depth)
        ab_dataList = cbind(alt_reads=ab_alt_reads,depth=ab_depth)
        ########## Process the dataList ##############
        data = as.data.frame(dataList)
        data = data[data$alt_reads>=2&(data$depth-data$alt_reads)>=2,]
        ab_data = as.data.frame(ab_dataList)
        ab_data = ab_data[ab_data$alt_reads>=2&(ab_data$depth-ab_data$alt_reads)>=2,]
        print(dim(data))
        print(dim(ab_data))
        # random arrange normal data
        dat = data[sample(dim(data)[1]),]
        
        # randomly insert abnormal region into normal data
        pos = sample(seq(1:nrow(dat)),1)
        final_data = rbind(dat[1:pos,],ab_data,dat[(pos+1):nrow(dat),])
        final_data = final_data[!is.na(final_data[,1])&!is.na(final_data[,2]),]
        
        
        z = final_data$alt_reads
        N = final_data$depth
        BAF = z/N
        m = mean(BAF)
        nSNP = length(z)
        datalist = list(nSNP=nSNP,
                        nSites = length(N_cov),
                        N_cov = N_cov,
                        m = mean(BAF),
                        m0 = m0,
                        z = z,
                        N = N)
        return(datalist)
    }
    
    #5. normal case
    if (type == "normal"){
        cov_mean = coverage
        cov_var = cov_mean*vartom
        # use the mean and variance of coverage to calculate r and p of negative binomial model
        p = cov_mean/cov_var
        r = cov_mean^2/(cov_var-cov_mean)
        
        # random generate the noise level between 0-10%
        noise = runif(1,min=0,max=0.1)
        print(noise)
        # generate the depth for B allele according to the coverage and BAF
        alt_reads = NULL
        depth = NULL
        for (i in 1:((1-noise)*nSNP)){
            # generate depth by a negative binomial distribution
            sub_depth = rnbinom(1,size=r,prob=p)
            # add this while loop to make sure that the coverage for all the sites are greater than zero
            # also depth have to greater than 1 to make sure site is heterozygous
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            
            # at normal case, theta equals m
            theta = m
            # generate the depth for alternative allele by binomial distribution
            sub_alt = rbinom(1,sub_depth,theta)
            while(sub_alt==sub_depth|sub_alt/sub_depth<0.02|sub_alt/sub_depth>0.98){
                sub_alt = rbinom(1,sub_depth,theta)
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        for (i in 1:(noise*nSNP)){
            sub_depth = rnbinom(1,size=r,prob=p)
            while(sub_depth<=1){
                sub_depth = rnbinom(1,size=r,prob=p)
            }
            depth = c(depth,sub_depth)
            sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            while(sub_alt==sub_depth){
                sub_alt = rbinom(1,sub_depth,rbeta(1,1,1))
            }
            alt_reads= c(alt_reads,sub_alt)
        }
        
        
        N_cov = rnbinom((1-noise)*nCov,size=r,prob=p)
        # generate noise% of the coverage as outliers
        N_cov_out = round(runif(noise*nCov,min=1,max=max(N_cov)))
        N_cov = c(N_cov[N_cov!=0],N_cov_out)
        m0 = round(rnorm(1,mean(N_cov),0.05*coverage),2)
        dataList = cbind(alt_reads,depth)
    }
    
    data = as.data.frame(dataList)
    data = data[data$alt_reads>=2&(data$depth-data$alt_reads)>=2,]
    dat = data[sample(dim(data)[1]),]
    nSNP = dim(dat)[1]
    z = dat$alt_reads
    N = dat$depth
    BAF = z/N
    m = mean(BAF)
    datalist = list(nSNP=nSNP,
                    nSites = length(N_cov),
                    N_cov = N_cov,
                    m = mean(BAF),
                    m0 = m0,
                    z = z,
                    N = N)
    print(mean(N_cov))
    print(m0)
    return(datalist)
}
