##################################################################
##################################################################
#Perform regression analysis on level 2 elements
##################################################################
##################################################################

#get first level 2 score and the L2 norm of L2 scores for each day
load("/Users/Selene/Desktop/menu+rfh/coef.2.rdata")
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/menu+rfh/J.rdata")
l2=rep(NA,nrow(coef.2))
for (i in 1:nrow(coef.2)){
	v=coef.2[i,]/1000
	l2[i]=sqrt(sum(v^2))
}

#aggregate into subject level summary statistics
ind=rep(1:length(J),J)
dt=data.frame(ind=ind,l2=l2,pc1=coef.2[,1])
dt.m=aggregate(dt$l2,list(dt$ind),mean)
l2.p=dt.m[[2]]
dt.r=aggregate(dt$pc1,list(dt$ind),function(x) max(x)-min(x))
PC1.p=dt.r[[2]]/1000

#perform regression
load("/Users/Selene/Desktop/menu+rfh/df.rdata")
df$l2=l2.p
mod1=lm(insulin~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+l2, data=df)
mod2=lm(CRP~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+l2, data=df)
mod3=lm(QOLp~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+l2, data=df)
mod4=lm(QOLm~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+l2, data=df)

df$lev2.PCrange=PC1.p
df=df[df$lev2.PCrange!=0,]
mod1.2=lm(insulin~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+lev2.PCrange, data=df)
mod2.2=lm(CRP~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+lev2.PCrange, data=df)
mod3.2=lm(QOLp~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+lev2.PCrange, data=df)
mod4.2=lm(QOLm~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+lev2.PCrange, data=df)

#set level 2 pc scores to be the dependent variable, lme on start time, weekday, etc
agg=aggregate(gooddf$activity, list(gooddf$dt, gooddf$identifier),mean)
names(agg)=c('date','id','activity')
agg$date=as.Date(agg$date)
agg$dow=weekdays(agg$date)
agg$ind=ifelse(agg$dow=='Saturday'|agg$dow=='Sunday',1,0)
agg$ind=as.factor(agg$ind)
n=nrow(gooddf)/720
start=rep(NA,n)
for (i in 1:n){
	start[i]=gooddf[(i-1)*720+1,'time']
}
agg$start=as.numeric(start)
agg$QOLp=rep(df$QOLp,J)
agg$QOLm=rep(df$QOLm,J)
agg$insulin=rep(df$insulin,J)
agg$CRP=rep(df$CRP,J)
agg$age=rep(df$age,J)
agg$education=rep(df$education,J)
agg$bmi=rep(df$bmi,J)
agg$smoke=rep(df$smoke,J)
agg$PC1=coef.2[,1]
agg$PC2=coef.2[,2]
agg$PC3=coef.2[,3]
agg$PC4=coef.2[,4]
agg$l2=l2
agg=na.omit(agg)

mod1=lme(PC1~age+education+bmi+smoke+QOLp+insulin+CRP+start+ind,random=~1|id,agg,method='REML')
mod2=lme(PC2~age+education+bmi+smoke+QOLp+insulin+CRP+start+ind,random=~1|id,agg,method='REML')
mod3=lme(PC3~age+education+bmi+smoke+QOLp+insulin+CRP+start+ind,random=~1|id,agg,method='REML')
mod4=lme(PC4~age+education+bmi+smoke+QOLp+insulin+CRP+start+ind,random=~1|id,agg,method='REML')
mod5=lme(l2~age+education+bmi+smoke+QOLp+insulin+CRP+start+ind,random=~1|id,agg,method='REML')













