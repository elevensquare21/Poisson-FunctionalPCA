##################################################################
##################################################################
#Account for over-dispersion via two models (detailed in thesis paper)
#1. compound Poisson model
#2. quasi Poisson model 
#The following code computes the parameters for each method
##################################################################
##################################################################
library('pracma')
#library('elliptic')
load("/Users/Selene/Desktop/menu+rfh/coef.1.rdata")
load("/Users/Selene/Desktop/menu+rfh/fpca1.vectors.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")

#get intensity curve for each person and each day
inten.d=coef.1 %*% t(fpca1.vectors)
M=length(unique(gooddf$identifier))
inten.p=matrix(nrow=M,ncol=720)
for(i in 1:M){
	row=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	inten.p[i,]=inten.d[row,]
}
inten.p=t(inten.p)
inten.d=t(inten.d)

n=nrow(gooddf)/720
x=gooddf$activity
y=matrix(x,nrow=n,byrow=TRUE)
t=1:720
mu=apply(y, 2,mean)
intensity.p=inten.p+mu
intensity=inten.d+mu
for(i in 1:ncol(intensity)){
	intensity[,i]=ifelse(intensity[,i]>0,intensity[,i],0)
}
for(i in 1:ncol(intensity.p)){
	intensity.p[,i]=ifelse(intensity.p[,i]>0,intensity.p[,i],0)
}
save(intensity, file='/Users/Selene/Desktop/menu+rfh/intensity.rdata')
save(intensity.p, file='/Users/Selene/Desktop/menu+rfh/intensity.p.rdata')

#find cumulative integral matrix based on intensity.p
intemat=matrix(0,720,M)
for (i in 1:ncol(intemat)){
	p.1=intensity.p[,i]
	p.1.t=c(p.1[1],p.1)
	inte=rep(0,720)
	for(j in 1:720){
		a=1:(j+1)
		b=p.1.t[1:(j+1)]
		inte[j]=trapz(a,b)
	}
	intemat[,i]=inte
}
save(intemat, file='/Users/Selene/Desktop/menu+rfh/intemat.rdata')

#group data into panels/bins (refer to paper for details)
#rough criteria for how big bins are: on average, want 100 bins per day
panel=max(c(mean(intemat[720,])/100,max(intemat[1,]))) #panel = 2129.545

#record where to set bins
ind.list=list()
for (i in 1:length(unique(gooddf$identifier))){
	inte=intemat[,i]
	n.pan=floor(inte[720]/panel)
	ind=rep(0,n.pan)
	for (j in 1:n.pan){
		p=j*panel
		ind[j]=min(which(inte>p))-1
	}	
	ind.list[[i]]=ind
}
save(ind.list, file='/Users/Selene/Desktop/menu+rfh/ind.list.rdata')

#get sum of accelerometer counts in each bin
Sp=list()
for(i in 1:length(unique(gooddf$identifier))){
	id=unique(gooddf$identifier)[i]
	ma=gooddf[gooddf$identifier==id,]
	n.day=nrow(ma)/720
	ind=ind.list[[i]]
	n.p=length(ind)
	S.mat=matrix(0,n.day,n.p)
	for(j in 1:n.day){
		act=ma[(j-1)*720+1:720,3]
		vec=rep(0,n.p)
		for(k in 1:n.p){			
			vec[k]=sum(act[(ifelse(k==1,0,ind[k-1])+1): ind[k]])
		}
		S.mat[j,]=vec
	}
	Sp[[i]]=S.mat
}
save(Sp,file='/Users/Selene/Desktop/menu+rfh/Sp.rdata')
v.Sp=lapply(Sp,as.vector)
Sn=unlist(v.Sp)
save(Sn,file='/Users/Selene/Desktop/menu+rfh/Sn.rdata')


##################################################################
#Model 1: compound Poisson 
##################################################################

#model mark distribution using negative binomial/ quasi Poisson distribution
#estimate same mark distribution for all subjects using Sn
rqpois = function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}
e.Sn=mean(Sn)
v.Sn=var(Sn)
e.mark=e.Sn/panel
v.mark=v.Sn/panel-e.mark^2
theta=v.mark/e.mark

#estimate same mark distribution for all subjects using X(t)
e.X=mean(gooddf$activity)
v.X=var(gooddf$activity)
e.mark=e.X/mean(intensity)
v.mark=v.X/mean(intensity)-e.mark^2
theta=v.mark/e.mark

#estimate different mark distribution for each person using Sn 
e.mark = rep(NA, length(Sp))
v.mark = rep(NA, length(Sp))
theta = rep(NA, length(Sp))
for (i in 1: length(Sp)){
	e.Sn=mean(as.vector(Sp[[i]]))
	v.Sn=var(as.vector(Sp[[i]]))
	e.mark[i]=e.Sn/panel
	v.mark[i]=v.Sn/panel-e.mark[i]^2
	theta[i]=v.mark[i]/e.mark[i]

}
markdist = data.frame(emark = rep(e.mark,J), vmark = rep(v.mark, J), thetamark = rep(theta, J))
save(markdist,file='/Users/Selene/Desktop/menu+rfh/markdist.rdata')

#estimate different mark distribution for each person using X(t)
e.mark = rep(NA, length(J))
v.mark = rep(NA, length(J))
theta = rep(NA, length(J))
for (i in 1: length(J)){
	start = (ifelse(i==1,0,sum(J[1:(i-1)])))*720+1
	end = sum(J[1:i]) * 720
	col = (ifelse(i==1,0,sum(J[1:(i-1)])))+1
	tempdf = gooddf[start:end,]
	e.X=mean(tempdf$activity)
	v.X=var(tempdf$activity)
	e.mark[i]=e.X/mean(intensity[,col])
	v.mark[i]=v.X/mean(intensity[,col])-e.mark[i]^2
	theta[i]=v.mark[i]/e.mark[i]

}
markdist2 = data.frame(emark = rep(e.mark,J), vmark = rep(v.mark, J), thetamark = rep(theta, J))
save(markdist2,file='/Users/Selene/Desktop/menu+rfh/markdist2.rdata')

#estimate different mark distribution for each person using Sn with mixed model
#still use the estimate for e.mark from before, but do mixed model for v.mark
e.mark = rep(NA, length(Sp))
v.mark = rep(NA, length(Sp))
theta = rep(NA, length(Sp))
for (i in 1: length(Sp)){
	e.Sn=mean(as.vector(Sp[[i]]))
	v.Sn=var(as.vector(Sp[[i]]))
	e.mark[i]=e.Sn/panel
	v.mark[i]=v.Sn/panel-e.mark[i]^2
	theta[i]=v.mark[i]/e.mark[i]

}
vlist = list()
for (i in 1:length(Sp)){
	a = e.mark[i]
	Sn.v = as.vector(Sp[[i]])
	response = ((Sn.v-panel)^2)/panel - (a^2)
	vlist[[i]]=data.frame(Y = response, id = i)
	
}
library(plyr)
mixeddf = ldply(vlist, data.frame)
library(nlme)
mod = lme(Y~1, random = ~1|id, data= mixeddf)
fix = mod$coef$fixed
v.mark = fix + mod$coef$random[[1]]
theta = v.mark/e.mark
markdist3 = data.frame(emark = rep(e.mark,J), vmark = rep(v.mark, J), thetamark = rep(theta, J))
save(markdist3,file='/Users/Selene/Desktop/menu+rfh/markdist3.rdata')

#estimate different mark distribution for each person using X(t) with mixed model
#still use the estimate for e.mark from before, but do mixed model for v.mark
e.mark = rep(NA, length(J))
v.mark = rep(NA, length(J))
theta = rep(NA, length(J))
for (i in 1: length(J)){
	start = (ifelse(i==1,0,sum(J[1:(i-1)])))*720+1
	end = sum(J[1:i]) * 720
	col = (ifelse(i==1,0,sum(J[1:(i-1)])))+1
	tempdf = gooddf[start:end,]
	e.X=mean(tempdf$activity)
	v.X=var(tempdf$activity)
	e.mark[i]=e.X/mean(intensity[,col])
	v.mark[i]=v.X/mean(intensity[,col])-e.mark[i]^2
	theta[i]=v.mark[i]/e.mark[i]

}
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
intensity = intensity +1
gooddf$intensity = as.vector(intensity)
gooddf$emark = rep(e.mark, J*720)
gooddf$Y = ((gooddf$activity - gooddf$intensity)^2)/gooddf$intensity - (gooddf$emark ^2)
mod = lme(Y~1, random = ~1|identifier, data= gooddf)
fix = mod$coef$fixed
v.mark = fix + mod$coef$random[[1]]
theta = v.mark/e.mark
markdist4 = data.frame(emark = rep(e.mark,J), vmark = rep(v.mark, J), thetamark = rep(theta, J))
save(markdist4,file='/Users/Selene/Desktop/menu+rfh/markdist4.rdata')


##################################################################
#Model 2: quasi Poisson 
##################################################################

#use an over-dispersed poisson/negative binomial distribution instead of regular poisson for intensity curve
#estimate same theta for all subjects using Sn
theta.poi = var(Sn)/panel #4118.763

#estimate same theta for all subjects using X(t)
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
intensity = intensity +1
gooddf$intensity = as.vector(intensity)
gooddf$Y = ((gooddf$activity - gooddf$intensity)^2)/gooddf$intensity 
theta.poi = mean(gooddf$Y) #4755.991

#estimate different theta distribution for each person using Sn 
theta = rep(NA, length(Sp))
for (i in 1: length(Sp)){
	v.Sn=var(as.vector(Sp[[i]]))
	theta[i]=v.Sn/panel

}
thetadist = data.frame(theta = rep(theta, J))
save(thetadist,file='/Users/Selene/Desktop/menu+rfh/thetadist.rdata')

#estimate different theta distribution for each person using X(t)
theta=aggregate(gooddf$Y, list(gooddf$identifier),mean)[[2]]
thetadist2 = data.frame(theta = rep(theta, J))
save(thetadist2,file='/Users/Selene/Desktop/menu+rfh/thetadist2.rdata')

#estimate different theta distribution for each person using Sn with mixed model
vlist = list()
for (i in 1:length(Sp)){
	Sn.v = as.vector(Sp[[i]])
	response = ((Sn.v-panel)^2)/panel
	vlist[[i]]=data.frame(Y = response, id = i)
	
}
mixeddf = ldply(vlist, data.frame)
library(nlme)
mod = lme(Y~1, random = ~1|id, data= mixeddf)
fix = mod$coef$fixed
theta = fix + mod$coef$random[[1]]
thetadist3 = data.frame(theta = rep(theta, J))
save(thetadist3,file='/Users/Selene/Desktop/menu+rfh/thetadist3.rdata')


#estimate different theta distribution for each person using X(t) with mixed model
mod = lme(Y~1, random = ~1|identifier, data= gooddf)
fix = mod$coef$fixed
theta = fix + mod$coef$random[[1]]
thetadist4 = data.frame(theta = rep(theta, J))
save(thetadist4,file='/Users/Selene/Desktop/menu+rfh/thetadist4.rdata')














