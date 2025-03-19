library(loo)
library(rjags)
library(mvtnorm)

options(scipen=999)
rm(list=ls())

Niter = 60000
burnIn = 30000

pp <- 3

tau<-0.5
beta <- c(6,2,1,0,1)


#-------Load Data-----------------------------------------------------------------
for(N in c(100)){
  name <- paste('y.','pp',pp,'.N',N,'.txt',sep='') 
  #name1	<-	paste('MCAR1.','pp',pp,'.N',N,'.txt',sep='')
  #name2	<-	paste('MCAR2.','pp',pp,'.N',N,'.txt',sep='')
  #name3	<-	paste('MCAR3.','pp',pp,'.N',N,'.txt',sep='')
  name4	<-	paste('MAR1.','pp',pp,'.N',N,'.txt',sep='')
  name5	<-	paste('MAR2.','pp',pp,'.N',N,'.txt',sep='')
  name6	<-	paste('MAR3.','pp',pp,'.N',N,'.txt',sep='')
  name7	<-	paste('MNAR1.','pp',pp,'.N',N,'.txt',sep='')
  name8	<-	paste('MNAR2.','pp',pp,'.N',N,'.txt',sep='')
  name9	<-	paste('MNAR3.','pp',pp,'.N',N,'.txt',sep='')
  for(dataname in c(name,name4,name5,name6,name7,name8,name9)){
    # name1,name2,name3,name4,name5,name6
    data <- read.table(dataname)
    colnames(data) <- c('id','V1','V2','V3','V4')
    
    #data <- read.table('y.pp1.N1002.txt')
    #colnames(data) <- c('id','V1','V2','V3','V4','V5')
    
    #load(file = paste0('data/avg_tri_1819_n1000.Rdata'))
    
    sum1 <- 0
    sum2 <- 0
    ESE.sum1 <- 0   
    num1 <- 0
    
    for(rep2 in 1:500){
      y <- data[which(data$id == rep2),-1] 
      N <- nrow(y)
      Time <- ncol(y)
      Lambda <- matrix(c(rep(1,Time), (c(1:Time)-1)), nrow = Time)
      
      
      #------------------------------------------
      # GCM
      dat <- list("N" = N, "y" = y, "tau" = tau, "Time" = ncol(y))
      
      # Set initial values
      initial=list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = 20231210)
      
      jags.gmm.md1 <- jags.model( file = "gcm_median_MNAR.txt", 
                                  data=dat, inits= initial, n.chains=1, n.adapt=1000 )
      params <- c("par")
      samps.gmm.md1 <- coda.samples(jags.gmm.md1, params, n.iter = Niter)
      
      smp.md1 <- window(samps.gmm.md1[[1]], start = burnIn)
      geweke.md1 <- apply(smp.md1, 2, function(x) geweke.diag(x)$z)
      
      if(all(geweke.md1[1:5] > -2 & geweke.md1[1:5] < 2)){
        res <- summary(smp.md1)$quantiles[1:5,3]
        sum1 <- sum1 + res
        ESE.sum1 <- ESE.sum1 + (res- beta)^2
        num1 <- num1  + 1
        
        res1 <- summary(smp.md1)$quantiles[6:9,3]
        sum2 <- sum2 + res1
      } #end if
    } # end for
    
    # est, average se
    est.mean1 <- sum1/num1
    est.mean2 <- sum2/num1
    # est bias
    bias1 <- est.mean1-beta
    # empiricl standard error (ESE)
    ESE1 <-  sqrt(ESE.sum1/(num1-1))
    # MSE
    MSE1 <- bias1^2 + ESE1^2
    # relative bias
    relative.bias1 <- 100*bias1/beta
    #model coverage rate
    model.coverage <- 100*num1/500
    #est, bias, ESE, MSE, relative bias, model coverage rate
    est.GCM <- 	cbind(dataname, est.mean1, bias1, ESE1*100, MSE1*1000, relative.bias1,model.coverage)
    name.GCM <- paste('GCMmedianMNAR100.','pp',pp,'.txt', sep='')
    write.table(est.GCM,file= name.GCM,append=TRUE,row.names=F,col.names=F)
    #
   # est.GCM2 <- 	cbind(dataname, est.mean2)
   # name.GCM2 <- paste('GCMmedianPar6_9.','pp',pp,'.txt', sep='')
   # write.table(est.GCM2,file= name.GCM2,append=TRUE,row.names=F,col.names=F)
    
  }
}
