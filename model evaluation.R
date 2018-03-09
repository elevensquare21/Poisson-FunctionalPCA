##################################################################
##################################################################
#Evaluate quasi Poisson model under two metrics:
#1. correlation between simulated data and original data
#2. L2 distance in frequencies (within preset intervals) between simulated data and original data
##################################################################
##################################################################

##################################################################
#use correlation to validate results 
##################################################################
library('Hmisc')

#same theta for all subjects using Sn
theta.poi = var(Sn)/panel #4118.763
no=nrow(gooddf)/720
NN=no
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long1=long

#same theta for all subjects using X(t)
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
intensity = intensity +1
gooddf$intensity = as.vector(intensity)
gooddf$Y = ((gooddf$activity - gooddf$intensity)^2)/gooddf$intensity 
theta.poi = mean(gooddf$Y) #4755.991
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long4=long 

#different theta distribution for each person using Sn
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	theta.poi=thetadist[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long2 = long

#different theta distribution for each person using X(t)
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	theta.poi=thetadist2[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long5=long

#different theta distribution for each person using Sn with mixed model
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	theta.poi=thetadist3[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long3 = long

#different theta distribution for each person using X(t) with mixed model
long=rep(0,NN)
for(i in 1:NN){
	day.sim=intensity[,i]
	theta.poi=thetadist4[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	day.true=y[i,]
	long[i]=rcorr(sim,day.true,type='pearson')$P[1,2]
	
}
long6 = long

#regular Poisson model 
q=rep(0,no)
for(i in 1:no){
	q[i]=rcorr(intensity[,i],y[i,],type='pearson')$P[1,2]
}

#boxplot results
v1=c(long1,long2, long3,long4,long5,long6, q)
v2.1=rep('method 1',NN)
v2.2=rep('method 2',NN)
v2.3=rep('method 3',NN)
v2.4=rep('intensity curves',no)
v2=c(v2.1,v2.2, v2.3,v2.1,v2.2, v2.3, v2.4)
v3.1=rep('Nm',3*NN)
v3.2=rep('Xij(t)',3*NN)
v3.3=rep('n/a',NN)
v3=c(v3.1,v3.2,v3.3)
df=as.data.frame(cbind(v1,v2,v3))
df$v1=as.numeric(levels(df$v1))[df$v1]
names(df)=c('correlation','method','tool')
save(df,file='/Users/Selene/Desktop/menu+rfh/metric1result.rdata')
#boxplot(correlation~tool*method,data=df,main='Poisson Intensity Curves VS Quasi-Poisson Intensity Curves',xlab='Model Types',ylab='correlation')
library(ggplot2)
ggplot(df, aes(method, correlation, group=interaction(method, tool), colour=factor(tool)))+geom_boxplot(position=position_dodge(width=0.4), width=0.4)+labs(colour = "Quantities used:")


##################################################################
#compare frequencites in preset intervals 
##################################################################

#use increments of 50 from 0 to 2000 as intervals 
breaks=seq(from=0,to=2000,by=50)
nn=length(breaks)-1

#true data frequencies
true.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	act=gooddf[((i-1)*720+1):(i*720),3]
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(act>=breaks[j] & act<breaks[j+1])
	}
	true.data[i,]=v
}

#regular Poisson model 
pca.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	act=intensity[,i]
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(act>=breaks[j] & act<breaks[j+1])
	}
	pca.data[i,]=v

}
pca.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=pca.data[i,]
	pca.eval[i]=dist(rbind(a, b))
}
q=pca.eval

#same theta for all subjects using Sn
theta.poi = var(Sn)/panel
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v

}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long1=sim.eval

#same theta for all subjects using X(t)
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
intensity = intensity +1
gooddf$intensity = as.vector(intensity)
gooddf$Y = ((gooddf$activity - gooddf$intensity)^2)/gooddf$intensity 
theta.poi = mean(gooddf$Y) #4755.991
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v

}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long4=sim.eval

#different theta distribution for each person using Sn
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long2=sim.eval

#different theta distribution for each person using X(t)
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist2[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long5=sim.eval

#different theta distribution for each person using Sn with mixed model
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist3[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long3=sim.eval

#different theta distribution for each person using X(t) with mixed model
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist4[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long6=sim.eval

#boxplot results
v1=c(long1,long2, long3,long4,long5,long6, q)
v2.1=rep('method 1',NN)
v2.2=rep('method 2',NN)
v2.3=rep('method 3',NN)
v2.4=rep('intensity curves',no)
v2=c(v2.1,v2.2, v2.3,v2.1,v2.2, v2.3, v2.4)
v3.1=rep('Nm',3*NN)
v3.2=rep('Xij(t)',3*NN)
v3.3=rep('n/a',NN)
v3=c(v3.1,v3.2,v3.3)
df=as.data.frame(cbind(v1,v2,v3))
df$v1=as.numeric(levels(df$v1))[df$v1]
names(df)=c('FrenquencyDistance','method','tool')
save(df,file='/Users/Selene/Desktop/menu+rfh/metric2result.rdata')
#boxplot(correlation~tool*model,data=df,main='Poisson Intensity Curves VS Quasi-Poisson Intensity Curves',xlab='Model Types',ylab='correlation')
library(ggplot2)
ggplot(df, aes(method, FrenquencyDistance, group=interaction(method, tool), colour=factor(tool)))+geom_boxplot(position=position_dodge(width=0.4), width=0.4)+labs(colour = "Quantities used:")


#use cutoff points for sedentary activity and MVPA to set intervals
breaks=c(0,100,1952,max(gooddf$activity))
nn=length(breaks)-1

#true data frequencies
true.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	act=gooddf[((i-1)*720+1):(i*720),3]
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(act>=breaks[j] & act<breaks[j+1])
	}
	true.data[i,]=v
}

#regular Poisson model
pca.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	act=intensity[,i]
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(act>=breaks[j] & act<breaks[j+1])
	}
	pca.data[i,]=v

}

pca.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=pca.data[i,]
	pca.eval[i]=dist(rbind(a, b))
}
q=pca.eval

#same theta for all subjects using Sn
theta.poi = var(Sn)/panel
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v

}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long1=sim.eval

#same theta for all subjects using X(t)
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
intensity = intensity +1
gooddf$intensity = as.vector(intensity)
gooddf$Y = ((gooddf$activity - gooddf$intensity)^2)/gooddf$intensity 
theta.poi = mean(gooddf$Y) #4755.991
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v

}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long4=sim.eval

#different theta distribution for each person using Sn
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long2=sim.eval

#different theta distribution for each person using X(t)
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist2[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long5=sim.eval

#different theta distribution for each person using Sn with mixed model
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist3[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long3=sim.eval

#different theta distribution for each person using X(t) with mixed model
sim.data=matrix(nrow=n,ncol=nn)
for(i in 1:n){
	day.sim=intensity[,i]
	theta.poi=thetadist4[i,1]
	sim=rep(0,720)
	for(j in 1:720){
		lam=day.sim[j]
		sim[j]=rqpois(1,lam, theta.poi)-1
		sim[j]=ifelse(sim[j]<0,0,sim[j])
	}
	v=rep(NA,nn)
	for(j in 1:nn){
		v[j]=sum(sim>=breaks[j] & sim<breaks[j+1])
	}
	sim.data[i,]=v
}
sim.eval=rep(NA,n)
for(i in 1:n){
	a=true.data[i,]
	b=sim.data[i,]
	sim.eval[i]=dist(rbind(a, b))
}
long6=sim.eval

#boxplot results
v1=c(long1,long2, long3,long4,long5,long6, q)
v2.1=rep('method 1',NN)
v2.2=rep('method 2',NN)
v2.3=rep('method 3',NN)
v2.4=rep('intensity curves',no)
v2=c(v2.1,v2.2, v2.3,v2.1,v2.2, v2.3, v2.4)
v3.1=rep('Nm',3*NN)
v3.2=rep('Xij(t)',3*NN)
v3.3=rep('n/a',NN)
v3=c(v3.1,v3.2,v3.3)
df=as.data.frame(cbind(v1,v2,v3))
df$v1=as.numeric(levels(df$v1))[df$v1]
names(df)=c('FrenquencyDistance','method','tool')
save(df,file='/Users/Selene/Desktop/menu+rfh/metric2result2.rdata')
#boxplot(correlation~tool*model,data=df,main='Poisson Intensity Curves VS Quasi-Poisson Intensity Curves',xlab='Model Types',ylab='correlation')
library(ggplot2)
ggplot(df, aes(method, FrenquencyDistance, group=interaction(method, tool), colour=factor(tool)))+geom_boxplot(position=position_dodge(width=0.4), width=0.4)+labs(colour = "Quantities used:")












