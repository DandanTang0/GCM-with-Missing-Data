
pp <- 3


library(MASS)
library(rsem)
library(missForest)
library(mice)
library(lattice)
library(rpart)
library(VIM)#KNN
library(magrittr) # %>%
library(mvtnorm)
library(dplyr) # if_else
library(lavaan)

gcmodel<-'i =~ 1*V1 + 1*V2 + 1*V3  + 1*V4
	s =~ 0*V1 + 1*V2 + 2*V3 + 3*V4'


beta <- c(1,1,0,6,2)

num <- 500

# two stage robust method
for(N in c(100,200,500)){
  #we don't need recompute the data without missing data
  name <- paste('y.','pp',pp,'.N',N,'.txt',sep='') 
  name1	<-	paste('MAR1.','pp',pp,'.N',N,'.txt',sep='')
  name2	<-	paste('MAR2.','pp',pp,'.N',N,'.txt',sep='')
  name3	<-	paste('MAR3.','pp',pp,'.N',N,'.txt',sep='')
  name4	<-	paste('MNAR1.','pp',pp,'.N',N,'.txt',sep='')
  name5	<-	paste('MNAR2.','pp',pp,'.N',N,'.txt',sep='')
  name6	<-	paste('MNAR3.','pp',pp,'.N',N,'.txt',sep='')
  for(dataname in c(name,name1,name2,name3,name4,name5,name6)){
    
    data <- read.table(dataname)
    colnames(data) <- c('id','V1','V2','V3','V4')
    # two-robust stage 
    sum2 <- 0
    ESE.sum2 <- 0
    coverage.beta2.1 <- 0
    #coverage.beta2.2 <- 0
    res22 <- list()
    robust.se2 <- list()
    for(rep2 in 1:num){
      data.rep2 <- data[which(data$id == rep2),] 
      #two-stage robust method
      pat <- rsem.pattern(data.rep2[,-1])
      
      musig <- rsem.emmusig(pat, max.it = 100000)
      res2 <- growth(gcmodel, sample.cov=musig$sigma, sample.mean=musig$mu, sample.nobs=N,mimic='EQS')
      
      ascov <- rsem.Ascov(pat, musig)
      robust.se <- rsem.se(res2, ascov$Gamma)
      res.est <- cbind(coef(res2), robust.se$se[[1]])[5:9,]
      sum2 <- sum2 + res.est
      ESE.sum2 <- ESE.sum2 + (res.est[,1] - beta)^2
      
      
      # coverage
      coverage.beta2.1 <- coverage.beta2.1 + if_else(( (res.est[,1] - 1.96*res.est[,2]) < beta & (res.est[,1] + 1.96*res.est[,2]) > beta), 1, 0)
      # coverage.beta2.2 <- coverage.beta2.2 + if_else(( (res.est[2,1] - 1.96*res.est[2,2]) < 2 & (res.est[2,1] + 1.96*res.est[2,2]) >2), 1, 0)
    }
    ##two-stage robust 
    # est, average se
    est.mean2 <- sum2/num 
    # est bias
    bias2 <- est.mean2[,1]-beta
    # empiricl standard error (ESE)
    ESE2 <-  sqrt(ESE.sum2/(num-1))
    # MSE
    MSE2 <- bias2^2 + ESE2^2
    # relative bias
    relative.bias2 <- bias2/beta
    
    # coverage rate
    coverage2 <- as.data.frame(100*coverage.beta2.1/num)
    # time
    time.average2 <- time.sum2/num
    #est, bias, ASE, ESE, MSE, relative bias, confidence interval,coverage rate,time 
    est.robust <- cbind(dataname, est.mean2[,1], bias2, est.mean2[,2]*100,ESE2*100,MSE2*1000,
                        relative.bias2*100,coverage2)
    name.Robust <- paste('Robust.','pp',pp,'.txt',sep='')
    
    write.table(est.robust,file= name.Robust,append=TRUE,row.names=F,col.names=F)

  }
}
