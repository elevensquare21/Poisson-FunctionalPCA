##################################################################
##################################################################
#Count how many minutes were lost due to truncating daily record to 720 min
##################################################################
##################################################################

load("/Users/Selene/Desktop/MENU_Acc_DataTransfer_Selene_20160922/complete consecutive profiles.rdata")
load("/Users/Selene/Desktop/MENU_Acc_DataTransfer_Selene_20160922/aggregate on mean.rdata")
sub=ag[ag$prop<0.5,]
sub$min=1440-sub$prop*1440
agg=aggregate(gooddf$activity,list(gooddf$dt,gooddf$identifier),mean,na.rm=TRUE)
names(agg)=c('dt','identifier','activity')
for(i in 1:nrow(agg)){
	agg$min[i]=sub[sub$identifier==agg$identifier[i] & sub$dt==agg$dt[i],5]
}

load("/Users/Selene/Desktop/multi level pca/aggregate on mean.rdata")
load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles.rdata")
sub=ag[ag$prop<0.5,]
sub$min=1440-sub$prop*1440
agg=aggregate(gooddf$activity,list(gooddf$dt,gooddf$identifier),mean,na.rm=TRUE)
names(agg)=c('dt','identifier','activity')
for(i in 1:nrow(agg)){
	agg$min[i]=sub[sub$identifier==agg$identifier[i] & sub$dt==agg$dt[i],5]
}
agg2=agg

min=c(agg$min, agg2$min)
min=min-720
sum(min)/length(min) #137
median(min) #129





