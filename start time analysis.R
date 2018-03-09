##################################################################
##################################################################
#Regression analysis on start time
##################################################################
##################################################################
library(nlme)
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")
load("/Users/Selene/Desktop/menu+rfh/df.rdata")

n=nrow(gooddf)/720
start=rep(NA,n)
for (i in 1:n){
	start[i]=gooddf[(i-1)*720+1,'time']
}
id=rep(unique(gooddf$identifier),J)
QOLp=rep(df$QOLp,J)
QOLm=rep(df$QOLm,J)
insulin=rep(df$insulin,J)
CRP=rep(df$CRP,J)
age=rep(df$age,J)
education=rep(df$education,J)
bmi=rep(df$bmi,J)
smoke=rep(df$smoke,J)
data=data.frame(id=id,QOLp=QOLp,QOLm=QOLm,insulin=insulin,CRP=CRP,age=age,education=education,bmi=bmi,smoke=smoke,start=start)
data=na.omit(data)

mod1=lme(start~age+education+bmi+smoke+QOLp,random=~1|id,data,method='REML')
#significant 
mod2=lme(start~age+education+bmi+smoke+QOLm,random=~1|id,data,method='REML')
mod3=lme(start~age+education+bmi+smoke+insulin,random=~1|id,data,method='REML')
#significant
mod4=lme(start~age+education+bmi+smoke+CRP,random=~1|id,data,method='REML')

agg=aggregate(start,list(id),mean)
meanstart=agg[[2]]
df$meanstart=meanstart
mod1=lm(QOLp~age+education+bmi+smoke+meanstart,df)
#significant
mod2=lm(QOLm~age+education+bmi+smoke+meanstart,df)
mod3=lm(insulin~age+education+bmi+smoke+meanstart,df)
#significant
mod4=lm(CRP~age+education+bmi+smoke+meanstart,df)

#test whether start time is different between menu and rfh 
load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles.rdata")
gooddf2=gooddf
n=nrow(gooddf2)/720

c.start=start[1:n]
nc.start=start[(n+1):length(start)]
mean(c.start)
sd(c.start)
mean(nc.start)
sd(nc.start)
t.test(c.start,nc.start)







