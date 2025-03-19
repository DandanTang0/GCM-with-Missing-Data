# lognormal
pp <- 4


library(magrittr) 
library(mvtnorm)
library(dplyr) 

set.seed(20231210)


T	<-	4

miss.rate0 <- 0
miss.rate1 <- 0.05
miss.rate2 <- 0.15
miss.rate3	<-	0.3


beta		<-	c(6,2)
lambda	<-	matrix(c(1,1,1,1, 0,1,2,3),T,2)
mu		<-	lambda%*%beta

slognorm <- function(n,nsd){
  u <- rnorm(n)
  z <- exp(u)
  sz <- (z-exp(0.5))/sqrt(exp(1)*(exp(1)-1))*nsd		##"standardization"
  sz}


for(N in c(100,200,500)){
  NT <- N*T
  diaE <- 1
  varE <-	diag(diaE,T)
  covD <- 0
  D		<-	array(c(1,covD,covD,1), dim=c(2,2))
  
  for(rep1 in 1:500){
    id <- rep(rep1,N)
    e	<- matrix(slognorm(NT,sqrt(diaE)),N,T) 
    u 	<-	rmvnorm(N,c(0,0),D)
    y	<-	matrix(rep(mu,N),N,T,byrow=TRUE)+t(lambda%*%t(u)) + e
    yout <- cbind(id,y)
    name <- paste('y.','pp',pp,'.N',N,'.txt',sep='') 
    write.table(yout,file=name,append=TRUE,row.names=F,col.names=F)
    
    #MAR
    yMAR1	<-	y
    yMAR2	<-	y
    yMAR3	<-	y
    miss1 <- round(2*N*miss.rate1/(T-1))
    miss2 <- round(2*N*miss.rate2/(T-1))
    miss3 <- round(2*N*miss.rate3/(T-1))
    for (t in 1:(T-1)){
      yMAR1	<-	yMAR1[order(yMAR1[,t], decreasing=T),]
      yMAR1[(N-t*miss1+1):N,(t+1)]	<-	-99
      yMAR2	<-	yMAR2[order(yMAR2[,t], decreasing=T),]
      yMAR2[(N-t*miss2+1):N,(t+1)]	<-	-99
      yMAR3	<-	yMAR3[order(yMAR3[,t], decreasing=T),]
      yMAR3[(N-t*miss3+1):N,(t+1)]	<-	-99
    }
    for(i in 1:N){
      for(j in 2:T){
        if(yMAR1[i,j]==-99){yMAR1[i,j]<-NA}
        if(yMAR2[i,j]==-99){yMAR2[i,j]<-NA}
        if(yMAR3[i,j]==-99){yMAR3[i,j]<-NA}
      }
    }
    yMAR1out <- cbind(id,yMAR1)
    yMAR2out <- cbind(id,yMAR2)
    yMAR3out <- cbind(id,yMAR3)
    
    name4	<-	paste('MAR1.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR1out,file=name4,row.names=F,col.names=F,append=TRUE)
    name5	<-	paste('MAR2.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR2out,file=name5,row.names=F,col.names=F,append=TRUE)
    name6	<-	paste('MAR3.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMAR3out,file=name6,row.names=F,col.names=F,append=TRUE)
    
    
    #MNAR
    b1	<-	u[,2]+beta[2]
    corAb <- 0.8
    a <- corAb/sqrt(1-corAb^2)
    Aux <- a*b1+rnorm(N,0,1)
    
    yMNAR1	<-	y
    yMNAR2	<-	y
    yMNAR3	<-	y
    miss1 <- 2*miss.rate1/(T-1)
    miss2 <- 2*miss.rate2/(T-1)
    miss3 <- 2*miss.rate3/(T-1)
    for(j in 2:T){
      crit1 <- qnorm((1-(j-1)*miss1),a*beta[2],sqrt(a^2+1))
      yMNAR1[which(Aux>crit1),j] <- NA
      crit2 <- qnorm((1-(j-1)*miss2),a*beta[2],sqrt(a^2+1))
      yMNAR2[which(Aux>crit2),j] <- NA
      crit3 <- qnorm((1-(j-1)*miss3),a*beta[2],sqrt(a^2+1))
      yMNAR3[which(Aux>crit3),j] <- NA
    }
    yMNAR1out <- cbind(id,yMNAR1)
    yMNAR2out <- cbind(id,yMNAR2)
    yMNAR3out <- cbind(id,yMNAR3)
    
    name7	<-	paste('MNAR1.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR1out,file=name7,row.names=F,col.names=F,append=TRUE)
    name8	<-	paste('MNAR2.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR2out,file=name8,row.names=F,col.names=F,append=TRUE)
    name9	<-	paste('MNAR3.','pp',pp,'.N',N,'.txt',sep='')
    write.table(yMNAR3out,file=name9,row.names=F,col.names=F,append=TRUE)
  }   
  
}

