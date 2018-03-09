##################################################################
##################################################################
#Compare true data against simulated data under different models for a random day
##################################################################
##################################################################
load("/Users/Selene/Desktop/menu+rfh/coef.1.rdata")
load("/Users/Selene/Desktop/menu+rfh/fpca1.vectors.rdata")
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/intensity.rdata")
load("/Users/Selene/Desktop/menu+rfh/thetadist3.rdata")
nn=nrow(gooddf)/720

#get true data for a random day
set.seed(3)
sam=sample(1:nn,1) #sam=1798
df=gooddf[(720*(sam-1)+1):(720*sam),]
y=df[[3]]
ysmooth=ksmooth(1:720,y,kernel='normal',bandwidth=20)

#simulate data using regular Poisson model 
z=intensity[,sam]
zsim=rep(0,720)
for(j in 1:720){
	lam=day.sim[j]
	zsim[j]=rpois(1,lam)
}
zsm=ksmooth(1:720,z,kernel='normal',bandwidth=20)

#simulate data using quasi-Poisson model
rqpois = function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}
day.sim=z+1
theta.poi=thetadist3[sam,1]
sim=rep(0,720)
for(j in 1:720){
	lam=day.sim[j]
	sim[j]=rqpois(1,lam, theta.poi)-1
	sim[j]=ifelse(sim[j]<0,0,sim[j])
}
predsmooth=ksmooth(1:720,sim,kernel='normal',bandwidth=20)

#plot
plot(ysmooth$y,type='l',ylab='activity counts',col='grey', xlab='time',lwd=2)
lines(zsm$y,col='orange',lwd=2)
lines(predsmooth$y,col=' green',lwd=2)
legend("topright",col=c('grey','orange','green'),legend=c('True data','Simulated data with regular Poisson intensity', 'Simulated data with quasi-Poisson intensity'),text.col=c('grey','orange','green'), lty=1,lwd=2)













