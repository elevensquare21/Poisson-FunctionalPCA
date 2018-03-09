##################################################################
##################################################################
#Plot eigen functions, PC scores, covariance matrix, and more
##################################################################
################################################################

#plot first 4 level 1 pca (smoothed)
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca1.vectors[,1]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plus1=mu.sm$y+pc1.sm$y
minus1=mu.sm$y-pc1.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (30.5%)',ylim=c(100,500))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca1.vectors[,2]*1000
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plus2=mu.sm$y+pc2.sm$y
minus2=mu.sm$y-pc2.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (10.2%)',ylim=c(100,500))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca1.vectors[,3]*1000
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plus3=mu.sm$y+pc3.sm$y
minus3=mu.sm$y-pc3.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (7.6%)',ylim=c(100,500))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca1.vectors[,4]*1000
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plus4=mu.sm$y+pc4.sm$y
minus4=mu.sm$y-pc4.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (4.7%)',ylim=c(100,500))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')

#plot level 1 PC by itself (without mu)
par(mfrow=c(2,2))

pc1=fpca1.vectors[,1]
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plot(pc1.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (30.5%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc2=fpca1.vectors[,2]
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plot(pc2.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (10.2%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc3=fpca1.vectors[,3]
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plot(pc3.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (7.6%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc4=fpca1.vectors[,4]
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plot(pc4.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (4.7%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)


#plot first 4 level 2 pca (smoothed)
par(mfrow=c(2,2))

mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=100)
pc1=fpca2.vectors[,1]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plus1=mu.sm$y+pc1.sm$y
minus1=mu.sm$y-pc1.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (6.0%)',ylim=c(200,450))
lines(plus1,lty=2,col='red')
lines(minus1,lty=2,col='blue')

pc2=fpca2.vectors[,2]*1000
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plus2=mu.sm$y+pc2.sm$y
minus2=mu.sm$y-pc2.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (4.9%)',ylim=c(200,450))
lines(plus2,lty=2,col='red')
lines(minus2,lty=2,col='blue')

pc3=fpca2.vectors[,3]*1000
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plus3=mu.sm$y+pc3.sm$y
minus3=mu.sm$y-pc3.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (3.6%)',ylim=c(200,450))
lines(plus3,lty=2,col='red')
lines(minus3,lty=2,col='blue')

pc4=fpca2.vectors[,4]*1000
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plus4=mu.sm$y+pc4.sm$y
minus4=mu.sm$y-pc4.sm$y
plot(mu.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (3.1%)',ylim=c(200,450))
lines(plus4,lty=2,col='red')
lines(minus4,lty=2,col='blue')

#plot level 2 PC by itself (without mu)
par(mfrow=c(2,2))

pc1=fpca2.vectors[,1]
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=100)
plot(pc1.sm,type='l',xlab='time',ylab='intensity function',main='PC 1 (6.0%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc2=fpca2.vectors[,2]
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=100)
plot(pc2.sm,type='l',xlab='time',ylab='intensity function',main='PC 2 (4.9%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc3=fpca2.vectors[,3]
pc3.sm=ksmooth(1:720,pc3,kernel='normal',bandwidth=100)
plot(pc3.sm,type='l',xlab='time',ylab='intensity function',main='PC 3 (3.6%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)

pc4=fpca2.vectors[,4]
pc4.sm=ksmooth(1:720,pc4,kernel='normal',bandwidth=100)
plot(pc4.sm,type='l',xlab='time',ylab='intensity function',main='PC 4 (3.1%)',ylim=c(-0.1,0.1))
abline(a=0,b=0)


#simulate comparison between two profiles with the same avg PA but different L2 norm
mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=50)
pc1=fpca1.vectors[,1]*1000
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=50)
plus1=mu.sm$y+pc1.sm$y

pc2=fpca1.vectors[,1]*1000+fpca1.vectors[,2]*500+fpca1.vectors[,3]*500-fpca1.vectors[,4]*3300
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=50)
plus2=mu.sm$y+pc2.sm$y

plot(plus1,type='l',xlab='time',ylab='intensity function',ylim=c(50,600), col='blue', lwd=3)
lines(plus2,type='l',xlab='time',ylab='intensity function',col ='green', lwd=3)
legend("topright",col=c('green','blue'),legend=c('person with larger L2 norm','person with smaller L2 norm'),text.col=c('green', 'blue'), lty=c(1,1))


#simulate comparison between two profiles with contrast between pc1 & pc2
par(mfrow=c(1,2))
ex=coef.1[7,]
mu1=fpca1.vectors %*% ex
mu.sm=ksmooth(1:720,mu,kernel='normal',bandwidth=50)

pc1=fpca1.vectors[,1]*1000+fpca1.vectors[,2]*100
pc1.sm=ksmooth(1:720,pc1,kernel='normal',bandwidth=50)
plus1=mu.sm$y+pc1.sm$y
plot(plus1,type='l',xlab='time',ylab='intensity function',main='Daily Profile with Large PC1 and Small PC2',ylim=c(50,520))
lines(mu.sm,lty=2)

pc2=fpca1.vectors[,2]*1000+fpca1.vectors[,1]*100
pc2.sm=ksmooth(1:720,pc2,kernel='normal',bandwidth=50)
plus2=mu.sm$y+pc2.sm$y
plot(plus2,type='l',xlab='time',ylab='intensity function',main='Daily Profile with Small PC1 and Large PC2',ylim=c(50,520))
lines(mu.sm,lty=2)


#3D plot of total covariance matrix
load("/Users/Selene/Desktop/menu+rfh/G.rdata")
library('plotly')

m=cov2cor(G)
p <- plot_ly(z=m, type="surface",showscale=TRUE)
p

#add some smoothness
df <- data.frame(x = rep(seq_len(ncol(m)), each = nrow(m)),
                 y = rep(seq_len(nrow(m)), times = ncol(m)),
                 z = c(m))
require("mgcv")
mod <- gam(z ~ te(x, y), data = df)
m2 <- matrix(fitted(mod), ncol = 720)
require("lattice")
#wireframe(m2)
correlation_matrix=m2
p <- plot_ly(z=correlation_matrix, type="surface",showscale=TRUE)
p

#try heat map with specified colors
p <- plot_ly(z=m2, type="heatmap",showscale=TRUE,colorscale = list(c(0, "rgb(0, 0, 0)"), list(1, "rgb(0, 255, 0)")))
library(RColorBrewer)
palette <- colorRampPalette(c("blue","green","yellow", "red", "darkred"))
p <- plot_ly(z=m2, type="heatmap",showscale=TRUE,colors = palette(50))
p

#try a different package
library(fields)
library(gplots)
m=cov2cor(G)
heatmap.2(m2, col=tim.colors(64), scale="none", Rowv=NULL, Colv=NULL,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=1.2, dendrogram ='none')


#summary and boxplot of PC1-4 (and their interaction with cancer status)
com=matrix(nrow=M,ncol=ncol(coef.1))
for(i in 1:M){
	row=ifelse(i==1,0,sum(J[1:(i-1)]))+1
	com[i,]=coef.1[row,]
}
com[,1]=com[,1]/1000
com[,2]=com[,2]/1000
com[,3]=com[,3]/1000
com[,4]=com[,4]/1000
dat=data.frame(PC1=com[,1],PC2=com[,2],PC3=com[,3],PC4=com[,4])

n=length(unique(gooddf2$identifier))
dat$cancer=c(rep('cancer',n),rep('non-cancer',nrow(dat)-n))

par(mfrow=c(2,2))
boxplot(dat$PC1, main='PC1',ylim=c(-13,15))
boxplot(dat$PC2, main='PC2',ylim=c(-13,15))
boxplot(dat$PC3, main='PC3',ylim=c(-13,15))
boxplot(dat$PC4, main='PC4',ylim=c(-13,15))

par(mfrow=c(2,2))
boxplot(PC1~cancer,data=dat, main='PC1')
boxplot(PC2~cancer,data=dat, main='PC2')
boxplot(PC3~cancer,data=dat, main='PC3')
boxplot(PC4~cancer,data=dat, main='PC4')






