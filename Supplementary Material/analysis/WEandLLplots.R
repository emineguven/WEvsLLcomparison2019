load('myWEandLLEnvironment.RData')
#sd_fit_Gomp<-unlist(sd_fit_Gomp)
#mean_fit_Gomp<-unlist(mean_fit_Gomp)

sd_fit_Weib<-unlist(sd_fit_Weib)
mean_fit_Weib<-unlist(mean_fit_Weib)

sd_fit_LLogis<-unlist(sd_fit_LLogis)
mean_fit_LLogis<-unlist(mean_fit_LLogis)

sample.lifespan.BY4742<-c(lifespansTemp_list[[10]])
sample.lifespan.BY4741<-c(lifespansTemp_list[[44]])
#i=10 is BY4742
#i =44 is BY4741
pdf(paste("plots/", "singleLifespanBY4742.pdf", sep=''))
#hist(unlist(R_squared_Gomp),breaks=length(R_squared_Gomp)/5,xlim=c(0.75,1),ylim=c(0,150))
#R.Llogis<-unlist(R_squared_LLogis)
hist(sample.lifespan.BY4742,breaks=15,xlim=c(0,50),ylim=c(0,15),col="blue",main="RLS distribution of a WT BY4742",
         xlab="Lifespan,Minutes",probability = F)
#lines(density(x,bw=1), col='red', lwd=3)
box(lty = 'solid', col = 'black')
dev.off()


pdf(paste("plots/", "singleLifespanBY4741.pdf", sep=''))
#hist(unlist(R_squared_Gomp),breaks=length(R_squared_Gomp)/5,xlim=c(0.75,1),ylim=c(0,150))
#R.Llogis<-unlist(R_squared_LLogis)
hist(sample.lifespan.BY4741,breaks=17,xlim=c(0,60),ylim=c(0,15),col="red",main="RLS distribution of a WT BY4741",
     xlab="Lifespan,Minutes",probability = F)
box(lty = 'solid', col = 'black')
dev.off()
#pdf(paste("plots/", "Histogram_fitted_data_summarize.pdf", sep=''))
#par(mfrow=c(1,2)) 
#hist(sd_fit_Gomp,breaks=length(sd_fit_Gomp)/5)
#hist(mean_fit_Gomp,breaks=length(mean_fit_Gomp)/5)

#hist(sd_fit_Weib,breaks=length(sd_fit_Weib)/5,xlim=c(0,0.16))
#hist(mean_fit_Weib,breaks=length(mean_fit_Weib)/5,)

#hist(sd_fit_LLogis,breaks=length(sd_fit_LLogis)/5)
#hist(mean_fit_LLogis,breaks=length(mean_fit_LLogis)/5)


#calculate CV's
#CV_Gomp=sd_fit_Gomp/mean_fit_Gomp
CV_Weib=sd_fit_Weib/mean_fit_Weib
CV_LLogis=sd_fit_LLogis/mean_fit_LLogis

#pdf(paste("plots/", "Histogram_empirical_vs_fitted_data_CV.pdf", sep=''))

#hist(CV_Gomp,main="CV of fitted data",breaks=length(CV_Gomp)/5,xlim=c(0,1),ylim=c(0,120),col="green")
#hist(new_fit$CV_exp,main="CV of fitted data",breaks=length(CV_exp_fit)/5,xlim=c(0.2,0.6),ylim=c(0,120),col="blue")

hist(CV_Weib,main="CV of fitted data",breaks=length(CV_Weib)/5,xlim=c(0,1),ylim=c(0,120),col="blue")

hist(CV_LLogis,main="CV of fitted data",breaks=length(CV_LLogis)/5,xlim=c(0,1),ylim=c(0,120),col="blue")

#pdf(paste("plots/", "Histogram_empirical_data_model_params.pdf", sep=''))

#par(mfrow=c(3,1)) 
hist(new_fit$AvgLS,xlab="Mean Lifespan",main="Empirical Data",
     breaks=length(new_fit$AvgLS)/10,col="gray",xlim=c(10,40),ylim=c(0,70))
mtext('Total experiments=5303',line=-1.5,cex=0.5)
box(lty='1373',col='blue')

#tor1<- new_fit[new_fit$genotype=="tor1",]
#sir2<- new_fit[new_fit$genotype=="sir2",]

#BY4742 plots
WT.BY4742<- new_fit[new_fit$genotype=="BY4742",]
#WT.BY4742.temp<-WT.BY4742[WT.BY4742$temp==30,]
WT.BY4742.media<-WT.BY4742[WT.BY4742$media=="YPD",]
WT.BY4742.media<-WT.BY4742.media[WT.BY4742.media$mat=="MATalpha",]
WT.BY4742.media= WT.BY4742.media[!is.na(WT.BY4742.media[,1]), ]

#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4742_model_params_AvglS.pdf", sep=''))
par(mfrow=c(3,1)) 

hist(as.numeric(WT.BY4742.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4742",
     breaks=length(as.numeric(WT.BY4742.media$AvgLS))/5,col="gray",xlim=c(10,40),ylim=c(0,25))
mtext('Total experiments=2108',line=-3,cex=2)
#dev.off()
#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4742_model_params_G.pdf", sep=''))
hist(as.numeric(WT.BY4742.media$weibScale),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$weibScale))/5,col="blue")

#dev.off()

#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4742_model_params_R.pdf", sep='')
hist(as.numeric(WT.BY4742.media$weibShape),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$weibShape))/10,col="green")

#dev.off()

#pdf(paste("plots/", "Histogram_Delta_LL_empirical_data_WT_BY4742.pdf", sep=''))

fifty.percent2<-quantile(WT.BY4742.media$Delta_LL, c(.50)) 


hist(WT.BY4742.media$Delta_LL,breaks=length(WT.BY4742.media$Delta_LL)/5,col="blue",xlim=c(-50,120),ylim=c(0,600),xlab=expression(paste( ~delta[LL])),main="Histogram of Delta_LL for WT BY4742")


mtext('Total experiments=2108',line=-3,at=70,side=3,cex=0.8)
mtext('median=4.28',line=-4,at=70,side=3,cex=0.8)
mtext('mean=9.98',line=-5,at=70,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent2, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-20.5,at=500,cex=0.8)

#dev.off()

#BY4741 plots
WT.BY4741<- new_fit[new_fit$genotype=="BY4741",]



#WT.BY4742.temp<-WT.BY4742[WT.BY4742$temp==30,]
WT.BY4741.media<-WT.BY4741[WT.BY4741$media=="YPD",]
WT.BY4741.media<-WT.BY4741.media[WT.BY4741.media$mat=="MATa",]
WT.BY4741.media= WT.BY4741.media[!is.na(WT.BY4741.media[,1]), ]


hist(as.numeric(WT.BY4742.media$AvgLS))
hist(as.numeric(WT.BY4742.media$weibScale))
hist(as.numeric(WT.BY4742.media$weibShape))
hist(as.numeric(WT.BY4742.media$LlogisScale))
#Error in hist.default(as.numeric(WT.BY4742.media$LlogisScale)) : 
#  invalid number of 'breaks'
 hist(as.numeric(WT.BY4742.media$LLogisScale))
hist(as.numeric(WT.BY4742.media$LLogisShape))
summary(WT.BY4742.media$weibShape)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7815  2.5210  2.8410  2.8140  3.1290  8.7220 
 summary(WT.BY4742.media$weibScale)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.30   26.85   29.64   29.20   31.98   44.88 
summary(WT.BY4742.media$LLogisScale)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.247  22.570  25.260  24.600  27.220  37.060 
 summary(WT.BY4742.media$LLogisShape)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.029   3.631   4.130   4.092   4.571  12.640 