##################################################################
##################################################################
#Apply multilevel functional principal component analysis to the combined dataset
#between MENU and RFH studies
##################################################################
##################################################################

#prep work: for both Menu and RFH datasets, find days with >720min consecutive wear time and discard the rest of the records. Keep only the first 720min of data points from each day. This results in a total of 3413 daily records, each having 720 minute-level data points. 
#combine menu and rfh data
load("/Users/Selene/Desktop/MENU_Acc_DataTransfer_Selene_20160922/complete consecutive profiles.rdata")
gooddf1=gooddf

load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles.rdata")
gooddf2=gooddf

gooddf=rbind(gooddf2,gooddf1)
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata')

#visualize mean
n=nrow(gooddf)/720
x=gooddf$activity
y=matrix(x,nrow=n,byrow=TRUE)
t=1:720
mu=apply(y, 2,mean)
plot(mu,type='l')

#plot mu +/- pointwise se
sd=apply(y,2,sd)
sd=sd/sqrt(n)
mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=50)
sd.sm=ksmooth(1:720,sd,kernel='normal',bandwidth=50)
plus=mu.sm$y+sd.sm$y
minus=mu.sm$y-sd.sm$y
plot(mu.sm,type='l',ylim=c(-500,1200),lwd=2,ylab='activity counts',xlab='time')
lines(plus,type='l',lty=2,col='orange')
lines(minus, type='l',lty=2,col='orange')

#perform multilevel functional PCA (code reference: Di's paper)
#1. find total covariance matrix, between subject covariance matrix, and 
#within subject covariance matrix
resd=matrix(0, nrow=n, ncol=720) 
resd=t( t(y) - mu ) 

index=function(J){
  col.1=gl(J,J-1)
  v=1:J
  col.2=rep(0,J*(J-1))
  for (i in 1:J){
    col.2[((i-1)*(J-1)+1):((i-1)*(J-1)+J-1)]=v[-i]
  }
  cbind(col.1,col.2)
}

gooddf$identifier=as.character(gooddf$identifier)

M=length(unique(gooddf$identifier))
J=rep(0,M)
for (i in 1:M){
  J[i]=nrow(gooddf[gooddf$identifier==unique(gooddf$identifier)[i],])/720
}
Jsum=J*(J-1)
SUM=sum(J*(J-1))

mat1 <- matrix(0, SUM, ncol=720)
mat2 <- matrix(0, SUM, ncol=720)
for (m in 1:M){
  if (J[m] > 1){
  for (k in 1:Jsum[m]){
    mat1[ ifelse(m==1,0,sum(Jsum[1:(m-1)]))+k, ]=resd[ifelse(m==1,0,sum(J[1:(m-1)])) + index(J[m])[k,1], ]
     mat2[ ifelse(m==1,0,sum(Jsum[1:(m-1)]))+k, ]=resd[ifelse(m==1,0,sum(J[1:(m-1)])) + index(J[m])[k,2], ]
  }
  }
}

N=720
G <- matrix(0, N, N)
Gb <- matrix(0, N, N)
Gw <- matrix(0, N, N)
for(i in 1:N){ 
  for(j in i:N) {
    G[i,j] <- cov(resd[,i],resd[,j])
    G[j,i] <- G[i,j]
    Gb[i,j] <- cov(mat1[,i],mat2[,j]) 
    Gb[j,i] <- Gb[i,j]
  }
}
Gw = G - Gb   

setwd("/Users/Selene/Desktop/menu+rfh")
save(Gb,file='Gb.rdata')
save(G,file='G.rdata')
save(Gw,file='Gw.rdata')

#2. eigen-decomposition of the covariance matrices
e1 <- eigen(Gb)
e2 <- eigen(Gw)
fpca1.value <- e1$values 
fpca2.value <- e2$values 
fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
percent1 <- (fpca1.value)/sum(fpca1.value)
percent2 <- (fpca2.value)/sum(fpca2.value)
K1 <- max( which(cumsum(percent1) < 0.9 | percent1 > 1/N ) + 1)
K2 <- max( which(cumsum(percent2) < 0.9 | percent2 > 1/N ) + 1)
rho=sum(fpca1.value)/(sum(fpca1.value)+sum(fpca2.value))
#K1=87 K2=313 rho=0.1125758

fpca1.vectors <- e1$vectors[, 1:K1]
fpca2.vectors <- e2$vectors[, 1:K2]

for(i in 1:K1) {
  v2 <- fpca1.vectors[,i]
  tempsign <- sum(v2)
  fpca1.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
for(i in 1:K2) {
  v2 <- fpca2.vectors[,i]
  tempsign <- sum(v2)
  fpca2.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
}
save(fpca1.vectors,file='fpca1.vectors.rdata')
save(fpca2.vectors,file='fpca2.vectors.rdata')

#3. calculate principal component scores using the projection method
#introduced in Di's paper
cross.integral <- matrix(0, K1, K2)
for(i in 1:K1)
  for(j in 1:K2) 
    cross.integral[i,j] <- sum(fpca1.vectors[,i]* fpca2.vectors[,j]) 

n=nrow(gooddf)/720
int1 <- matrix(0, n, K1)
int2 <- matrix(0, n, K2)
for(i in 1:n)   {
  for(j in 1:K1)  {
    int1[ i ,j] <- sum( resd[i,] * fpca1.vectors[,j] ) 
  }
  for(j in 1:K2) {
    int2[ i ,j] <- sum( resd[i,] * fpca2.vectors[,j] )    
  }
}

s1 <- matrix(0, n, K1)
s2 <- matrix(0, n, K2)
library(MASS)
design.xi <- ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
for(m in 1:M) {
  resid <- rep(0, K1)
  for(j in 1:J[m]) {
    index <-  ifelse(m==1,0,sum(J[1:(m-1)])) + j
    resid <- resid + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/J[m]
  }
  index.m <- ( ifelse(m==1,0,sum(J[1:(m-1)])) + 1 ) : (sum(J[1:m]))
  xi.temp <- design.xi %*% resid
  s1[index.m,] <- matrix(rep(xi.temp, each=J[m]), nrow=J[m])
  for(j in 1:J[m]) {
    index <- ifelse(m==1,0,sum(J[1:(m-1)])) + j
    s2[index,] <- int2[index,] - drop( t(cross.integral) %*% xi.temp )
  }
}

coef.1=s1
save(coef.1,file='coef.1.rdata')
coef.2=s2
save(coef.2,file='coef.2.rdata')




