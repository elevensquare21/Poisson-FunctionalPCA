##################################################################
##################################################################
#Improve curve registration by not only aligning start time but also aligning end time
##################################################################
##################################################################

#from the "extended data with activity and date.rdata" for each study, do following for 
#each daily record, say v
v=v[!is.na(v)]
l=length(v)
interval=(l-1)/719
newv=rep(NA,720)
newv[1]=v[1]
newv[720]=v[l]
for(i in 1:718){
	f=floor(i*interval)
	c=ceiling(i*interval)
	w2=i*interval-f
	w1=c-i*interval
	ind1=f+1
	ind2=c+1
	newv[i+1]=v[ind1]*w1+v[ind2]*w2
}
#results are saved under "complete consecutive profiles2.rdata" under the directory of each study


load("/Users/Selene/Desktop/MENU_Acc_DataTransfer_Selene_20160922/complete consecutive profiles2.rdata")
gooddf1=gooddf

load("/Users/Selene/Desktop/multi level pca/complete consecutive profiles2.rdata")
gooddf2=gooddf

gooddf=rbind(gooddf2,gooddf1)
save(gooddf, file='/Users/Selene/Desktop/menu+rfh/complete consecutive profiles2.rdata')

