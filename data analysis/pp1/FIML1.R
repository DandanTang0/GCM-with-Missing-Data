

pp <- 1

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
set.seed(20231210)
gcmodel<-'i =~ 1*V1 + 1*V2 + 1*V3  + 1*V4
	s =~ 0*V1 + 1*V2 + 2*V3 + 3*V4'


beta <- c(1,1,0,6,2)


#FIML
num <- 500
for(N in c(100,200,500)){
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
    #FIML
    sum1 <- 0
    ESE.sum1 <- 0
    time.sum1 <- 0
    coverage.beta1.1 <- 0
    
    # coverage.beta1.2 <- 0
    for(rep2 in 1:num){
      data.rep2 <- data[which(data$id == rep2),] 
      # FIML
      
      res1 <- growth(gcmodel, data=data.rep2, missing = "fiml")
      
      est1 <- parameterEstimates(res1)
      sum1 <- sum1 + est1[c(13:15,20:21),4:9]
      ESE.sum1 <- ESE.sum1 + (est1[c(13:15,20:21),4] - beta)^2
      coverage.beta1.1 <- coverage.beta1.1 + if_else((est1[c(13:15,20:21),8] < beta & est1[c(13:15,20:21),9] > beta), 1, 0)
      
      # coverage.beta1.2 <- coverage.beta1.2 + if_else((est1[c(13:15,20:21),8] < 2 & est1[21,9] >2), 1, 0)
      # coverage.beta1.1 <- coverage.beta1.1 + if_else((est1[13,8] < 1 & est1[13,9] >1), 1, 0)
      # coverage.beta1.2 <- coverage.beta1.2 + if_else((est1[14,8] < 2 & est1[21,9] >2), 1, 0)
    }
    ## FIML
    # est, average se
    est.mean1 <- sum1/num
    # est bias
    bias1 <- est.mean1[,1]-beta
    # empiricl standard error (ESE)
    ESE1 <-  sqrt(ESE.sum1/(num-1))
    # MSE
    MSE1 <- bias1^2 + ESE1^2
    # relative bias
    relative.bias1 <- bias1/beta
    # coverage rate
    coverage1 <- 100*coverage.beta1.1/num
    
    #est, bias, ASE, ESE, MSE, relative bias, confidence interval,coverage rate,time 
    est.FIML <- 	cbind(dataname,est.mean1[,1],bias1,est.mean1[,2]*100, ESE1*100,MSE1*1000,
                       relative.bias1*100,est.mean1[,5:6],coverage1)
    name.FIML <- paste('FIML.','pp',pp,'.txt', sep='')
    
    write.table(est.FIML,file= name.FIML,append=TRUE,row.names=F,col.names=F)
  
    #name.res1 <- paste(dataname,'FIML','.rda', sep='')
    
    #save(res11,file =  name.res1)
  }
}



