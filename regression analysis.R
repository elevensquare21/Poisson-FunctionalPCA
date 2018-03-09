##################################################################
##################################################################
#Perform regression analysis between health variables (biomarkers and QoL) and 
#physical activity variables (PC scores and over-dispersion parameters)
##################################################################
##################################################################
setwd('/Users/Selene/Desktop/menu+rfh')

#load data
load("/Users/Selene/Desktop/menu+rfh/complete consecutive profiles.rdata")
M=length(unique(gooddf$identifier))
J=rep(0,M)
for (i in 1:M){
  J[i]=nrow(gooddf[gooddf$identifier==unique(gooddf$identifier)[i],])/720
}
save(J,file='J.rdata')

#load level 1 PC scores and compute L2 norm 
load("/Users/Selene/Desktop/menu+rfh/coef.1.rdata")
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

#scale the first 4 PC scores
com[,1]=com[,1]/1000
com[,2]=com[,2]/1000
com[,3]=com[,3]/1000
com[,4]=com[,4]/1000

#load health variables (df1 and df2 were prepared under RFH and MENU separately)
df=rbind(df1,df2)
df$PC1=com[,1]
df$PC2=com[,2]
df$PC3=com[,3]
df$PC4=com[,4]
save(df,file='df.rdata')
df$l2=com.l2/1000

#add theta info
load("/Users/Selene/Desktop/menu+rfh/thetadist.rdata")
load("/Users/Selene/Desktop/menu+rfh/thetadist2.rdata")
load("/Users/Selene/Desktop/menu+rfh/thetadist3.rdata")
load("/Users/Selene/Desktop/menu+rfh/thetadist4.rdata")
for (i in 1:length(J)){
	j=(ifelse(i==1,0,sum(J[1:(i-1)])))+1
	df[i,'theta1']=thetadist[j,1]
	df[i,'theta2']=thetadist2[j,1]
	df[i,'theta3']=thetadist3[j,1]
	df[i,'theta4']=thetadist4[j,1]
	
}
save(df,file='df_theta.rdata')

#add in mark info
load("/Users/Selene/Desktop/menu+rfh/markdist.rdata")
load("/Users/Selene/Desktop/menu+rfh/markdist2.rdata")
load("/Users/Selene/Desktop/menu+rfh/markdist3.rdata")
load("/Users/Selene/Desktop/menu+rfh/markdist4.rdata")
for (i in 1:length(J)){
	j=(ifelse(i==1,0,sum(J[1:(i-1)])))+1
	df[i,'emark1']=markdist[j,1]
	df[i,'vmark1']=markdist[j,2]
	df[i,'emark2']=markdist2[j,1]
	df[i,'vmark2']=markdist2[j,2]
	df[i,'emark3']=markdist3[j,1]
	df[i,'vmark3']=markdist3[j,2]
	df[i,'emark4']=markdist4[j,1]
	df[i,'vmark4']=markdist4[j,2]
	
}
save(df,file='df_mark.rdata')


#perform regression
mod1=lm(insulin~age+education+bmi+smoke+PC1*cancer+PC2*cancer+PC3*cancer+PC4*cancer, data=df)
mod1.2=lm(insulin~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df)
mod2=lm(CRP~age+education+bmi+smoke+PC1*cancer+PC2*cancer+PC3*cancer+PC4*cancer, data=df)
mod2.2=lm(CRP~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df)
mod3=lm(QOLp~age+education+bmi+smoke+PC1*cancer+PC2*cancer+PC3*cancer+PC4*cancer, data=df)
mod3.2=lm(QOLp~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df)
mod4=lm(QOLm~age+education+bmi+smoke+PC1*cancer+PC2*cancer+PC3*cancer+PC4*cancer, data=df)
mod4.2=lm(QOLm~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df)
mod5=lm(insulin~age+education+bmi+smoke+l2*cancer, data=df)
mod5.2=lm(insulin~age+education+bmi+smoke+cancer+l2, data=df)
mod6=lm(CRP~age+education+bmi+smoke+l2*cancer, data=df)
mod6.2=lm(CRP~age+education+bmi+smoke+cancer+l2, data=df)
mod7=lm(QOLp~age+education+bmi+smoke+l2*cancer, data=df)
mod7.2=lm(QOLp~age+education+bmi+smoke+cancer+l2, data=df)
mod8=lm(QOLm~age+education+bmi+smoke+l2*cancer, data=df)
mod8.2=lm(QOLm~age+education+bmi+smoke+cancer+l2, data=df)

library(lmtest)
lrtest(mod1.2, mod1)
lrtest(mod2.2, mod2)
lrtest(mod5.2, mod5)
lrtest(mod6.2, mod6)
lrtest(mod3.2, mod3)
lrtest(mod7.2, mod7)
lrtest(mod4.2, mod4)
lrtest(mod8.2, mod8)

#regression of HOMA 
df$insulinstatus=ifelse(df$HOMA>3,1,0)
df$insulinstatus=as.factor(df$insulinstatus)
df$cancer=as.factor(df$cancer)

mod1.HOMA=lm(HOMA~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df)
mod2.HOMA=lm(HOMA~age+education+bmi+smoke+cancer+l2, data=df)
mod3.HOMA=lm(HOMA~scale(age)+education+scale(bmi)+smoke+cancer+scale(PC1)+scale(PC2)+scale(PC3)+scale(PC4)+scale(theta3), data=df)
mod4.HOMA=lm(HOMA~scale(age)+education+scale(bmi)+smoke+cancer+scale(l2)+scale(theta3), data=df)
mod5.HOMA=glm(insulinstatus~scale(age)+education+scale(bmi)+smoke+cancer+scale(PC1)+scale(PC2)+scale(PC3)+scale(PC4)+scale(theta3), data=df, family = "binomial")
mod6.HOMA=glm(insulinstatus~scale(age)+education+scale(bmi)+smoke+cancer+scale(l2)+scale(theta3), data=df, family = "binomial")
mod7.HOMA=lm(HOMA~scale(age)+education+scale(bmi)+smoke+cancer+scale(PC1)+scale(PC2)+scale(PC3)+scale(PC4)+scale(emark4)+scale(vmark4), data=df)
mod8.HOMA=lm(scale(HOMA)~scale(age)+education+scale(bmi)+smoke+cancer+scale(l2)+scale(emark4)+scale(vmark4), data=df)
mod9.HOMA=glm(insulinstatus~scale(age)+education+scale(bmi)+smoke+cancer+scale(PC1)+scale(PC2)+scale(PC3)+scale(PC4)+scale(emark4)+scale(vmark4), data=df, family = "binomial")
mod10.HOMA=glm(insulinstatus~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4, data=df, family = "binomial")
mod11.HOMA=glm(insulinstatus~age+education+bmi+smoke+cancer+l2,  data=df, family = "binomial")

#regress PC on cancer and age, (opt: add bmi)
mod1=lm(PC1~cancer+age,data=df)
mod1.2=lm(PC1~cancer+age+bmi, data=df)
mod2=lm(PC2~cancer+age,data=df)
mod2.2=lm(PC2~cancer+age+bmi, data=df)
mod3=lm(PC3~cancer+age,data=df)
mod3.2=lm(PC3~cancer+age+bmi, data=df)
mod4=lm(PC4~cancer+age,data=df)
mod4.2=lm(PC4~cancer+age+bmi, data=df)

#regression analysis with mvpa and sed
gooddf$mvpa=ifelse(gooddf$activity>1951,1,0)
gooddf$sed=ifelse(gooddf$activity<100,1,0)
ag.mvpa=aggregate(gooddf$mvpa,list(gooddf$identifier),sum)
ag.sed=aggregate(gooddf$sed,list(gooddf$identifier),sum)
ag.pa=aggregate(gooddf$activity,list(gooddf$identifier),mean)
df$mvpa=ag.mvpa[[2]]
df$sed=ag.sed[[2]]
df$pa=ag.pa[[2]]

mod1.2=lm(insulin~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+mvpa, data=df)
mod2.2=lm(CRP~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+mvpa, data=df)
mod3.2=lm(QOLp~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+mvpa, data=df)
mod4.2=lm(QOLm~age+education+bmi+smoke+cancer+PC1+PC2+PC3+PC4+mvpa, data=df)
mod5.2=lm(insulin~age+education+bmi+smoke+cancer+l2+mvpa, data=df)
mod6.2=lm(CRP~age+education+bmi+smoke+cancer+l2+mvpa, data=df)
mod7.2=lm(QOLp~age+education+bmi+smoke+cancer+l2+mvpa, data=df)
mod8.2=lm(QOLm~age+education+bmi+smoke+cancer+l2+mvpa, data=df)


#correlation analysis between mvpa, sed, pa and pc scores
cor(df$PC1,df$mvpa)
cor.test(df$PC1,df$mvpa) #p-value < 2.2e-16, r=0.6793644

mat=data.frame(PC1=df$PC1,PC2=df$PC2,PC3=df$PC3,PC4=df$PC4,L2=df$l2,pa=df$pa,mvpa=df$mvpa,sed=df$sed)
cor(mat)[1:5,6:8]

pmat=matrix(NA,nrow=8,ncol=8)
for(i in 1:8){
	for(j in i:8){
		pmat[i,j]=cor.test(mat[[i]],mat[[j]])$p.value
		pmat[j,i]=pmat[i,j]
	}
}
pmat[1:5,6:8]








