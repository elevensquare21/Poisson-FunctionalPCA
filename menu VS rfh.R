##################################################################
##################################################################
#Compare activity patterns between the Menu study(no cancer) and the RFH study(cancer)
##################################################################
##################################################################

#mean comparison plot
nc.x=gooddf1$activity
nc.y=matrix(nc.x,nrow=nrow(gooddf1)/720,byrow=TRUE)
nc.mu=apply(nc.y, 2,mean)
nc.sm=ksmooth(1:720,nc.mu,kernel='normal',bandwidth=30)
c.x=gooddf2$activity
c.y=matrix(c.x,nrow=nrow(gooddf2)/720,byrow=TRUE)
c.mu=apply(c.y, 2,mean)
c.sm=ksmooth(1:720,c.mu,kernel='normal',bandwidth=30)
plot(c.sm,type='l',xlab='time',ylab='activity counts',col='red',lty=2,lwd=3)
lines(nc.sm,lwd=3)
legend("topright",col=c('black','red'),legend=c('non-cancer','cancer'),text.col=c('black','red'),lty=c(1,2))

#t test
n=nrow(gooddf2)/720
c.m=coef.1[1:n,]
nc.m=coef.1[(n+1):nrow(coef.1),]
c.m=unique(c.m)
nc.m=unique(nc.m)
t.test(c.m[,1],nc.m[,1]) #p=0.5169
t.test(c.m[,2],nc.m[,2]) #p=0.04673
t.test(c.m[,3],nc.m[,3]) #p=0.8382
t.test(c.m[,4],nc.m[,4]) #p=0.03685

c.l2=rep(0,nrow(c.m))
nc.l2=rep(0,nrow(nc.m))
for(i in 1:nrow(c.m)){
	v=c.m[i,]
	c.l2[i]=sqrt(sum(v^2))
}
for(i in 1:nrow(nc.m)){
	v=nc.m[i,]
	nc.l2[i]=sqrt(sum(v^2))
}
t.test(c.l2,nc.l2) #p=0.2279

#test for mvpa, sed and PA
#first get daily mvpa, sedentary time for the two studies
load("/Users/Selene/Desktop/MENU_Acc_DataTransfer_Selene_20160922/extended data with activity and date.rdata")
edata=data.e
agg=aggregate(edata$activity,list(edata$dt,edata$identifier),mean,na.rm=TRUE)
names(agg)=c('dt','identifier','activity')
mvpa.0=list()
sed.0=list()
miss=list()
for(i in 1:length(unique(agg$identifier))){
	x=edata[edata$identifier==unique(agg$identifier)[i],]
	s=rep(0,length(unique(x$dt)))
	t=rep(0,length(unique(x$dt)))
	m=rep(0,length(unique(x$dt)))
	for(j in 1:length(unique(x$dt))){
		s[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<100,3])
		t[j]=length(x[x$dt==unique(x$dt)[j] & x$activity>1951,3])
		m[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<0,3])

	}
	sed.0[[i]]=s
	mvpa.0[[i]]=t
	miss[[i]]=m
}
mvpa.0=unlist(mvpa.0)
sed.0=unlist(sed.0)
miss=unlist(miss)
agg$miss=miss
agg$mvpa.0=mvpa.0
agg$sed.0=sed.0
agg$mvpa=agg$mvpa.0-agg$miss
agg$sed=agg$sed.0-agg$miss
agg$wt=1440-agg$miss
agg$mvpa.adj=agg$mvpa/agg$wt*720
agg$sed.adj=agg$sed/agg$wt*720
ag1=na.omit(agg)

load("/Users/Selene/Desktop/multi level pca/extended data with activity and date.rdata")
edata=data.e
agg=aggregate(edata$activity,list(edata$dt,edata$identifier),mean,na.rm=TRUE)
names(agg)=c('dt','identifier','activity')
mvpa.0=list()
sed.0=list()
miss=list()
for(i in 1:length(unique(agg$identifier))){
	x=edata[edata$identifier==unique(agg$identifier)[i],]
	s=rep(0,length(unique(x$dt)))
	t=rep(0,length(unique(x$dt)))
	m=rep(0,length(unique(x$dt)))
	for(j in 1:length(unique(x$dt))){
		s[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<100,3])
		t[j]=length(x[x$dt==unique(x$dt)[j] & x$activity>1951,3])
		m[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<0,3])

	}
	sed.0[[i]]=s
	mvpa.0[[i]]=t
	miss[[i]]=m
}
mvpa.0=unlist(mvpa.0)
sed.0=unlist(sed.0)
miss=unlist(miss)
agg$miss=miss
agg$mvpa.0=mvpa.0
agg$sed.0=sed.0
agg$mvpa=agg$mvpa.0-agg$miss
agg$sed=agg$sed.0-agg$miss
agg$wt=1440-agg$miss
agg$mvpa.adj=agg$mvpa/agg$wt*720
agg$sed.adj=agg$sed/agg$wt*720
ag2=na.omit(agg)

t.test(ag1$activity, ag2$activity) #0.003355
t.test(ag1$mvpa.adj, ag2$mvpa.adj) #0.4173
t.test(ag1$sed.adj, ag2$sed.adj) #0.0004068













