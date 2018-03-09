##################################################################
##################################################################
#Show what information can be captured by the L2 norm of level 1 PC scores and 
#the over-dispersion parameter theta through plots
##################################################################
##################################################################
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")
load("/Users/Selene/Desktop/menu+rfh/coef.1.rdata")
load("/Users/Selene/Desktop/menu+rfh/df.rdata")
M=length(unique(gooddf$identifier))
com=matrix(nrow=M,ncol=ncol(coef.1))
for(i in 1:M){
	row=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	com[i,]=coef.1[row,]
}
com.l2=rep(0,M)
for(i in 1:M){
	v=com[i,]
	com.l2[i]=sqrt(sum(v^2))
}
df$l2=com.l2/1000

#find person with largest L2 norm VS smallest L2 norm
large= which(df$l2==max(df$l2))
small = which(df$l2==min(df$l2))
order(df$l2)
small = 372 
large = 520

begin = (sum(J[1:(large-1)]))*720+1
end = sum(J[1:large]) * 720
sm.large = ksmooth(1:720,gooddf[(begin+720*6):(begin+720*7-1),3],kernel='normal',bandwidth=20)
#plot(sm.large)
begin = (sum(J[1:(small-1)]))*720+1
end = sum(J[1:small]) * 720
sm.small = ksmooth(1:720,gooddf[(begin+720*0):(begin+720*1-1),3],kernel='normal',bandwidth=20)
#plot(sm.small)
plot(sm.small,type='l',xlab='time',ylab='activity counts',col='red',lwd=3,ylim=c(0,max(unlist(sm.large))))
lines(sm.large,lwd=3)
legend("topright",col=c('black','red'),legend=c('person with the largest L2 norm','person with the smallest L2 norm'),text.col=c('black','red'),lty=c(1,2))


#find person with largest theta VS smallest theta
load("/Users/Selene/Desktop/menu+rfh/df_theta.rdata")
large= which(df$theta3==max(df$theta3))
small = which(df$theta3==min(df$theta3))
order(df$theta3)
small = 396

begin = (sum(J[1:(large-1)]))*720+1
end = sum(J[1:large]) * 720
sm.large = ksmooth(1:720,gooddf[(begin+720*2):(begin+720*3-1),3],kernel='normal',bandwidth=20)
#plot(sm.large)
begin = (sum(J[1:(small-1)]))*720+1
end = sum(J[1:small]) * 720
sm.small = ksmooth(1:720,gooddf[(begin+720*3):(begin+720*4-1),3],kernel='normal',bandwidth=20)
#plot(sm.small)
plot(sm.small,type='l',xlab='time',ylab='activity counts',col='red',lwd=3,ylim=c(0,max(unlist(sm.large))))
lines(sm.large,lwd=3)
legend("topright",col=c('black','red'),legend=c('person with the largest theta','person with the smallest theta'),text.col=c('black','red'),lty=c(1,2))






















