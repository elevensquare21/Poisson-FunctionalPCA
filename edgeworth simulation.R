##################################################################
##################################################################
#Perform simulations to show edgeworth expansion out-performs CLT
##################################################################
##################################################################
qnm=qnorm(0.975)
n.rep=1000
size=10
no=10000
lam=2000
q=rep(NA, n.rep)
p=rep(NA,n.rep)
for (j in 1:n.rep){
	data=rep(NA,no)
	for(i in 1:no){
		random=rpois(size,lam)
		std=((random-lam)^2-lam)/sqrt(lam*(1+2*lam))
		data[i]=sum(std)/sqrt(size)
	}
	q[j]=quantile(data,0.975)
	p[j]=sum(data<qnm)/no
}
sim.10=cbind(p,q)
save(sim.10,file="/Users/Selene/Desktop/sim.10.rdata")

sim.10=as.data.frame(sim.10)
qnm=qnorm(0.975) #1.959964
k3=(lam+22*lam^2+8*lam^3)/sqrt(size)/(lam+2*lam^2)^1.5
qcorr=qnm+(qnm^2-1)/6*k3
pcorr=pnorm(qnm)-dnorm(qnm)*k3*(qnm*2-1)/6
sim.10$p.diff=sim.10$p-pnorm(qnm) 
sim.10$pcorr.diff=sim.10$p-pcorr
sim.10$q.diff=sim.10$q-qnm 
sim.10$qcorr.diff=sim.10$q-qcorr
boxplot(sim.10$p.diff,sim.10$pcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated p-value and theoretical p-value')
abline(a=0,b=0)
boxplot(sim.10$q.diff,sim.10$qcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated quantile value and theoretical quantile value')
abline(a=0,b=0)

qnm=qnorm(0.975)
n.rep=1000
size=100
no=10000
lam=2000
q=rep(NA, n.rep)
p=rep(NA,n.rep)
for (j in 1:n.rep){
	data=rep(NA,no)
	for(i in 1:no){
		random=rpois(size,lam)
		std=((random-lam)^2-lam)/sqrt(lam*(1+2*lam))
		data[i]=sum(std)/sqrt(size)
	}
	q[j]=quantile(data,0.975)
	p[j]=sum(data<qnm)/no
}
sim.100=cbind(p,q)
save(sim.100,file="/Users/Selene/Desktop/sim.100.rdata")

sim.100=as.data.frame(sim.100)
qnm=qnorm(0.975) #1.959964
k3=(lam+22*lam^2+8*lam^3)/sqrt(size)/(lam+2*lam^2)^1.5
qcorr=qnm+(qnm^2-1)/6*k3
pcorr=pnorm(qnm)-dnorm(qnm)*k3*(qnm*2-1)/6
sim.100$p.diff=sim.100$p-pnorm(qnm) 
sim.100$pcorr.diff=sim.100$p-pcorr
sim.100$q.diff=sim.100$q-qnm 
sim.100$qcorr.diff=sim.100$q-qcorr
boxplot(sim.100$p.diff,sim.100$pcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated p-value and theoretical p-value')
abline(a=0,b=0)
boxplot(sim.100$q.diff,sim.100$qcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated quantile value and theoretical quantile value')
abline(a=0,b=0)

qnm=qnorm(0.975)
n.rep=1000
size=1000
no=10000
lam=2000
q=rep(NA, n.rep)
p=rep(NA,n.rep)
for (j in 1:n.rep){
	data=rep(NA,no)
	for(i in 1:no){
		random=rpois(size,lam)
		std=((random-lam)^2-lam)/sqrt(lam*(1+2*lam))
		data[i]=sum(std)/sqrt(size)
	}
	q[j]=quantile(data,0.975)
	p[j]=sum(data<qnm)/no
}
sim.1000=cbind(p,q)
save(sim.1000,file="/Users/Selene/Desktop/sim.1000.rdata")

sim.1000=as.data.frame(sim.1000)
qnm=qnorm(0.975) #1.959964
k3=(lam+22*lam^2+8*lam^3)/sqrt(size)/(lam+2*lam^2)^1.5
qcorr=qnm+(qnm^2-1)/6*k3
pcorr=pnorm(qnm)-dnorm(qnm)*k3*(qnm*2-1)/6
sim.1000$p.diff=sim.1000$p-pnorm(qnm) 
sim.1000$pcorr.diff=sim.1000$p-pcorr
sim.1000$q.diff=sim.1000$q-qnm 
sim.1000$qcorr.diff=sim.1000$q-qcorr
boxplot(sim.1000$p.diff,sim.1000$pcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated p-value and theoretical p-value')
abline(a=0,b=0)
boxplot(sim.1000$q.diff,sim.1000$qcorr.diff,names=c('CLT','Edgeworth Correction'),main='Difference between simulated quantile value and theoretical quantile value')
abline(a=0,b=0)



























