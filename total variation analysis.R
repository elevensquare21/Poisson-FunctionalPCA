##################################################################
##################################################################
#Regression analysis on total variation of the derivative of smoothed intensity curve
##################################################################
##################################################################

load("/Users/Selene/Desktop/menu+rfh/fpca1.vectors.rdata")
load("/Users/Selene/Desktop/menu+rfh/coef.1.rdata")
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")

#get intensity curve
pred=fpca1.vectors %*% t(coef.1)
M=length(unique(gooddf$identifier))
com=matrix(nrow=720,ncol=M)
for(i in 1:M){
	col=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	com[,i]=pred[,col]
}

#get total variation of the derivative of the smoothed intensity curve
totv=rep(NA,M)
for (i in 1: M){
	sm=ksmooth(1:720,com[,i],kernel='normal',bandwidth=100)
	v1=sm$y[1:719]
	v2=sm$y[2:720]
	totv[i]=sum(abs(v2-v1))
}

#can also get total variantion on raw data Xij(t)
totalv = rep(NA, sum(J))
for (i in 1:sum(J)){
	start = ifelse(i==1,0,i-1) * 720+1
	end = i*720
	sm=ksmooth(1:720, gooddf[start:end,3],kernel='normal',bandwidth=100)
	v1=sm$y[1:719]
	v2=sm$y[2:720]
	totalv[i]=sum(abs(v2-v1))
#	id[i] = gooddf[start, 1]
}
id = rep(1:M, J)
tempdf = data.frame(id = id, totalv = totalv)
redf = aggregate(tempdf$totalv, by= list(tempdf$id), mean)
#test its correlation with the L2 norm of level 1 PC scores
cor(redf$x, df$l2)
cor.test(redf$x, df$l2)
cor(redf$x, com.l2)
cor.test(redf$x, com.l2)

#perform regression
load("/Users/Selene/Desktop/menu+rfh/df.rdata")
df$totv=totv
mod5.2=lm(insulin~age+education+bmi+smoke+cancer+totv, data=df)
mod6.2=lm(CRP~age+education+bmi+smoke+cancer+totv, data=df)
mod7.2=lm(QOLp~age+education+bmi+smoke+cancer+totv, data=df)
mod8.2=lm(QOLm~age+education+bmi+smoke+cancer+totv, data=df)


