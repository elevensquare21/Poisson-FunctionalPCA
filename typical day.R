##################################################################
##################################################################
#Plot a typical day of accelerometer data
##################################################################
##################################################################
load("/Users/Selene/Dropbox/FPCA/multi level pca/extended data with activity and date.rdata")
y=data.e[1:1440,3]
sm=ksmooth(1:863,y[!is.na(y)],kernel='normal',bandwidth=15)
sm=c(rep(NA,501),sm$y,rep(NA,76))
plot(time, y,type='l',ylab='activity counts',col='grey')
abline(h=100,col='red',lty=2)
abline(h=1951,col='green',lty=2)
points(time,  y, pch=16, cex=0.5, col=ifelse(y<100,'red',ifelse(y>1951,'green','grey')))
legend("topright",col=c('green','red'),legend=c('MVPA','sedentary'),lty=c(2,2))
lines(time,sm,lwd=1.5)
legend("topleft",col='black',legend='smoothed activity curve',lwd=1.5,cex=0.7)

sum(y[!is.na(y)]<100) #591
sum(y[!is.na(y)]>1951) #7
sum(is.na(y)) #577 
