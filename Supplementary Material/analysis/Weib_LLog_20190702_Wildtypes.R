

library('varhandle')
library('flexsurv')
library('stringr')
library('plot3D')

##load original data for analyzing
new_rls<-read.csv("rls.csv")
my.data=new_rls

# uniquely determine the ref_lifespans of each experiments individually
big_data<-my.data[!duplicated(my.data$ref_lifespans),] 
big_data<-big_data[-c(114),]
#row 94 ref_lifespan is blank
ref_lifespan_mean<-list()

for (k in 1:(length(sort(unique(big_data$ref_lifespans))))){
#k=94
f1<-big_data$ref_lifespans
f<-unfactor(f1[k])
ref_lifespan_single<-as.numeric(unlist(str_split(f, ",")))
ref_lifespan_single_mean<-mean(ref_lifespan_single)
ref_lifespan_mean[[length(ref_lifespan_mean)+1]]<-ref_lifespan_single_mean
}


big_data$ref_lifespan_mean<-unlist(ref_lifespan_mean)

fit_names = c( 'popsize','genotype','media','mat','sd.ls','medianLS','AvgLS','weibAIC','weibLogLik','weibScale','weibShape', 'LLogisAIC','LLogisLogLik','LLogisScale','LLogisShape')
fit = data.frame(matrix(, nrow=length(fit_names), ncol=length(fit_names))) #set up a skeleton table
names(fit) = c( "popsize","genotype","media","mat","sd.ls","medianLS","AvgLS","weibAIC","weibLogLik","weibScale","weibShape", "LLogisAIC","LLogisLogLik","LLogisScale","LLogisShape")


lifespansChar_set<-big_data$set_lifespans
lifespansChar_ref<-big_data$ref_lifespans




#initilizations for CV calculations
new_fit<-read.csv("Rls_fitting_emine_13022017.csv")
#new_fit$CV_exp<-new_fit$sd.ls/new_fit$AvgLS



#R_squared_Gomp<-list()
R_squared_Weib<-list()
R_squared_LLogis<-list()


lifespansTemp_list<-list()

for (genotype_in in c(1:length(big_data$ref_name))){
  
  #genotype=117
  
  
  lifespansChar_set1<-unfactor(lifespansChar_set[genotype_in])
  #if lifespansChar_set1=="NA"
  lifespansTemp_set =as.numeric(unlist(str_split(lifespansChar_set1, ",")))

 
  lifespansChar_ref1<-unfactor(lifespansChar_ref[genotype_in])
  lifespansTemp_ref =as.numeric(unlist(str_split(lifespansChar_ref1, ",")))
  
  #lifespansTemp<-c(lifespansTemp_set,lifespansTemp_ref)  
  lifespansTemp<-c(lifespansTemp_ref,lifespansTemp_set)
  
                       
  
  lifespansTemp[lifespansTemp < 0] <- 0
  lifespansTemp <-floor(lifespansTemp+0.5)
  lifespansTemp <- lifespansTemp[ lifespansTemp != 0 ]
  
  lifespansTemp_list[[length(lifespansTemp_list)+1]]<-lifespansTemp
  
  pop_size<-length(lifespansTemp)

  
  ##lifespanGomp = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
  lifespanWeib = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'weibull')  
  lifespanLlog = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'llogis')  
    
   
  
## Fill in added columns in the controlConditions table with data from gompertz and weibull ##fitting. Columns are named according to respective variables
  media<-big_data$ref_media
  media=unfactor(media[genotype_in])
  mat<-big_data$ref_mating_type
  mat=unfactor(mat[genotype_in])
  avgLS = mean(lifespansTemp)
  StddevLS = sd(lifespansTemp)
  medianLS = median(lifespansTemp)
  #gompShape = lifespanGomp$res[1,1]
  #gompRate = lifespanGomp$res[2,1]
  #gompLogLik = lifespanGomp$loglik
  #gompAIC = lifespanGomp$AIC
  
  weibShape = lifespanWeib$res[1,1]
  weibScale = lifespanWeib$res[2,1]
  weibLogLik = lifespanWeib$loglik
  weibAIC = lifespanWeib$AIC  
  
  LLogisShape = lifespanLlog$res[1,1]
  LLogisScale = lifespanLlog$res[2,1]
  LLogisLogLik = lifespanLlog$loglik
  LLogisAIC = lifespanLlog$AIC   
  
  
  
  genotype<-big_data$ref_name
  
  genotype=unfactor(genotype[genotype_in])
  
  

  fit = rbind(fit,c(pop_size,genotype,media,mat,StddevLS,medianLS,avgLS,weibAIC,weibLogLik,weibScale,weibShape,LLogisAIC,LLogisLogLik,LLogisScale,LLogisShape))
 
  if (genotype=="BY4742"){
    pdf(paste("plots/BY4742/",genotype_in, "_BY4742_surv.pdf", sep='')) 
    plot(lifespanWeib, ylim=c(0,60),ci=FALSE, conf.int=FALSE, ylab="Survival Fraction", xlab="Lifespan,Minutes")
    lines(lifespanLlog, col="blue", ci=FALSE,ylim=c(0,60))    
    legend("topright", lty=c(1,1), lwd=c(2,2),
           col=c( "blue","red"),
           c("Weibull","Log-logistic"))
    dev.off()
    
    
  }
  if (genotype=="BY4741"){
    pdf(paste("plots/BY4741/",genotype_in, "_BY4741_surv.pdf", sep=''))  
    #plot(lifespanWeib,mark.time=TRUE,xlab="lifespan,minutes",ylab="survival fraction",main="Survival Plot of Weibull Model fit")
    plot(lifespanWeib, ylim=c(0,60),ci=FALSE, conf.int=FALSE, ylab="Survival Fraction", xlab="Lifespan,Minutes")
    lines(lifespanLlog, col="blue", ci=FALSE,ylim=c(0,60))    
    legend("topright", lty=c(1,1), lwd=c(2,2),
           col=c( "blue","red"),
           c("Weibull","Log-logistic"))
    dev.off()
    
  }
  
  
}




new_fit<-fit[(length(fit_names)+1):length(fit[,1]),]

write.csv(new_fit,"analysis_output_07022019.csv")


###read the analyzed data

new_fit<-read.csv("analysis_output_07022019.csv")

new_fit<-new_fit[-c(1),-c(21)]




#calculate R^2 manually

r.square = function( lifespan.observed, lifespan.fitted){

  # compute coefficient of determination
  if (class(lifespan.fitted) == "numeric") {
    return(cor(lifespan.observed, lifespan.fitted)^2)
  } else {
    R2 = double(ncol(lifespan.observed))
    for (ic in 1:ncol(lifespan.observed))
      R2[ic] = cor(lifespan.observed[,ic], lifespan.fitted[,ic])^2
    return(R2)
  }
}




R_squared<-list()
R_squared_Weib<-list()
R_squared_LLogis<-list()

sd_fit_Weib<-list()
mean_fit_Weib<-list()

sd_fit_LLogis<-list()
mean_fit_LLogis<-list()


#lifespansTemp_list


for (i in 1:length(new_fit[,1])){
  
  #lifespan_fit_Gomp<- new_fit$gompRate[i]* exp(new_fit$gompShape[i] *lifespansTemp_list[[i]])
  lifespan_fit_Weib<- new_fit$weibShape[i]*(lifespansTemp_list[[i]]^new_fit$weibShape[i])
  ####TODO next 04262018
  lifespan_fit_LLogis<-(new_fit$LLogisShape[i]/new_fit$LLogisScale[i])*(lifespansTemp_list[[i]]/new_fit$LLogisScale[i])^(new_fit$LLogisShape[i]-1)*1/(1+(lifespansTemp_list[[i]]/new_fit$LLogisScale[i])^(new_fit$LLogisShape[i]))
  lifespan.observed<-lifespansTemp_list[[i]]
  
  

  
  #R_squ_Gomp<- r.square(lifespan.observed,lifespan_fit_Gomp)
  R_squ_Weib<- r.square(lifespan.observed,lifespan_fit_Weib)
  R_squ_LLogis<- r.square(lifespan.observed,lifespan_fit_LLogis)
  
  #R_squared_Gomp[[length(R_squared_Gomp)+1]]<-R_squ_Gomp
  R_squared_Weib[[length(R_squared_Weib)+1]]<-R_squ_Weib
  R_squared_LLogis[[length(R_squared_LLogis)+1]]<-R_squ_LLogis
  
  
  
  #sd_fit_Gomp[[length(sd_fit_Gomp)+1]]<-sd(lifespan_fit_Gomp)
  #mean_fit_Gomp[[length(mean_fit_Gomp)+1]]<-mean(lifespan_fit_Gomp)
  
  sd_fit_Weib[[length(sd_fit_Weib)+1]]<-sd(lifespan_fit_Weib)
  mean_fit_Weib[[length(mean_fit_Weib)+1]]<-mean(lifespan_fit_Weib)
  
  sd_fit_LLogis[[length(sd_fit_LLogis)+1]]<-sd(lifespan_fit_LLogis)
  mean_fit_LLogis[[length(mean_fit_LLogis)+1]]<-mean(lifespan_fit_LLogis)
  
  
  
}


#R_squared_Gomp<-unlist(R_squared_Gomp)
#R_squared_Gomp<-R_squared_Gomp[!is.na(R_squared_Gomp)]
#summary(R_squared_Gomp)

R_squared_Weib<-unlist(R_squared_Weib)
R_squared_Weib<-R_squared_Weib[!is.na(R_squared_Weib)]
summary(R_squared_Weib)

R_squared_LLogis<-unlist(R_squared_LLogis)
R_squared_LLogis<-R_squared_LLogis[!is.na(R_squared_LLogis)]
summary(R_squared_LLogis)



pdf(paste("plots/", "Weib_Histogram_empirical_vs_fitted_R_squared.pdf", sep=''))
#hist(unlist(R_squared_Gomp),breaks=length(R_squared_Gomp)/5,xlim=c(0.75,1),ylim=c(0,150))
h<-hist(unlist(R_squared_Weib),breaks=length(R_squared_Weib)/7,xlim=c(0.7,1),ylim=c(0,60),col="cornflowerblue",main=TeX('R^2 relation of Weibull Model Fit and Observed Data'),
        xlab=TeX('R^2'),ylab="Lifespan")
xfit <- seq(min(unlist(R_squared_Weib)), max(unlist(R_squared_Weib)), length =length(unlist(R_squared_Weib)) ) 
yfit <- dnorm(xfit, mean = mean(unlist(R_squared_Weib)), sd = sd(unlist(R_squared_Weib))) 
yfit <- yfit * diff(h$mids[1:2]) * length(unlist(R_squared_Weib)) 

lines(xfit, yfit, col = "red", lwd = 3)
legend("topright", lty=c(1,1), lwd=c(2,2),
       col=c( "cornflowerblue","red"),
       c("Weibull","Gaussian"))
dev.off()




pdf(paste("plots/", "Llog_Histogram_empirical_vs_fitted_R_squared.pdf", sep=''))
#hist(unlist(R_squared_Gomp),breaks=length(R_squared_Gomp)/5,xlim=c(0.75,1),ylim=c(0,150))
R.Llogis<-unlist(R_squared_LLogis)
h2<-hist(unlist(R_squared_LLogis),breaks=length(R_squared_LLogis)/10,xlim=c(0.2,1),ylim=c(0,50),col="blue",main=TeX('R^2 relation of Log-logis Model Fit and Observed Data'),
        xlab=TeX('R^2'),ylab="Lifespan",probability = F)

xfit <- seq( min(unlist(R_squared_LLogis)), max(unlist(R_squared_LLogis)), length = length(unlist(R_squared_LLogis))) 
yfit <- dnorm(xfit, mean = mean(unlist(R_squared_LLogis)), sd = sd(unlist(R_squared_LLogis))) 
yfit <- yfit * diff(h2$mids[1:2]) * length(unlist(R_squared_LLogis)) 

lines(xfit, yfit, col = "red", lwd = 3)
legend("topright", lty=c(1,1), lwd=c(2,2),
       col=c( "blue","red"),
       c("Log-logis","Gaussian"))
dev.off()


load('myWEandLLEnvironment.RData')
#sd_fit_Gomp<-unlist(sd_fit_Gomp)
#mean_fit_Gomp<-unlist(mean_fit_Gomp)

sd_fit_Weib<-unlist(sd_fit_Weib)
mean_fit_Weib<-unlist(mean_fit_Weib)

sd_fit_LLogis<-unlist(sd_fit_LLogis)
mean_fit_LLogis<-unlist(mean_fit_LLogis)

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



hist(CV_Weib,main="CV of fitted data",breaks=length(CV_Weib)/5,xlim=c(0,1),ylim=c(0,120),col="green")

hist(CV_LLogis,main="CV of fitted data",breaks=length(CV_LLogis)/5,xlim=c(0,1),ylim=c(0,120),col="green")

#pdf(paste("plots/", "Histogram_empirical_data_model_params.pdf", sep=''))

#par(mfrow=c(3,1)) 
hist(new_fit$AvgLS,xlab="Mean Lifespan",main="Empirical Data",
     breaks=length(new_fit$AvgLS)/10,col="gray",xlim=c(10,40),ylim=c(0,70))
mtext('Total experiments=5303',line=-2,cex=0.5)
box(lty='1373',col='blue')

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
hist(as.numeric(WT.BY4742.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$gompShape))/5,col="blue",xlim=c(0,0.2),ylim=c(0,65))

#dev.off()

#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4742_model_params_R.pdf", sep='')
hist(log10(as.numeric(WT.BY4742.media$gompRate)),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$gompRate))/10,col="green",xlim=c(-3,-0.8),ylim=c(0,130))

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


#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4741_model_params_AvglS.pdf", sep=''))
par(mfrow=c(3,1)) 
hist(as.numeric(WT.BY4741.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4741",
     breaks=length(as.numeric(WT.BY4741.media$AvgLS))/5,col="gray",xlim=c(10,40),ylim=c(0,35))

mtext('Total experiments=381',line=-2,cex=2)
#dev.off()
#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4741_model_params_G.pdf", sep=''))
hist(as.numeric(WT.BY4741.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4741.media$gompShape))/5,col="blue",xlim=c(0.02,0.16),ylim=c(0,30))
#dev.off()

#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4741_model_params_R.pdf", sep=''))
hist(log10(as.numeric(WT.BY4741.media$gompRate)),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(WT.BY4741.media$gompRate))/5,col="green",xlim=c(-3,-1),ylim=c(0,25))
#dev.off()


#pdf(paste("plots/", "Histogram_Delta_LL_empirical_data_WT_BY4741.pdf", sep=''))
fifty.percent1<-quantile(WT.BY4741.media$Delta_LL, c(.50)) 
hist(WT.BY4741.media$Delta_LL,breaks=length(WT.BY4741.media$Delta_LL)/5,col="blue",xlim=c(-10,60),ylim=c(0,80),xlab=expression(paste( ~delta[LL])),main="Histogram of Delta_LL for WT BY4741 ")


mtext('Total experiments=381',line=-3,at=30,side=3,cex=0.8)
mtext('median=3.6',line=-4,at=30,side=3,cex=0.8)
mtext('mean=6.05',line=-5,at=30,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent1, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-24,at=60,cex=0.8)


#BY4743 plots
WT.BY4743<- new_fit[new_fit$genotype=="BY4743",]
WT.BY4743.media<-WT.BY4743[WT.BY4743$media=="YPD",]
WT.BY4743.media<-WT.BY4743.media[WT.BY4743.media$mat=="diploid",]
WT.BY4743.media= WT.BY4743.media[!is.na(WT.BY4743.media[,1]), ]

#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4743_model_params.pdf", sep=''))
par(mfrow=c(3,1)) 
hist(as.numeric(WT.BY4743.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4743",
     breaks=length(as.numeric(WT.BY4743.media$AvgLS)),col="gray",xlim=c(25,45),ylim=c(0,30))
hist(as.numeric(WT.BY4743.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4743.media$gompShape)),col="blue",xlim=c(0.05,0.25),ylim=c(0,30))
hist(as.numeric(WT.BY4743.media$gompRate),xlab="R rate parameters",main="",
     breaks=length(as.numeric(WT.BY4743.media$gompRate)),col="green",xlim=c(0,0.01),ylim=c(0,40))
#dev.off()



#pdf(paste("plots/", "WT-By4741-3D-plots-Delta_LL-G-R.pdf", sep=''))
WT.BY4741.media.neg<-WT.BY4741.media[WT.BY4741.media$Delta_LL<0,]

scatter3D(WT.BY4741.media.neg$gompShape,WT.BY4741.media.neg$gompRate,WT.BY4741.media.neg$Delta_LL,xlab="G",ylab="R",zlab="delta_LL",phi = 0, bty = "g",
          pch = 20, cex = 1,xlim=c(0.04,0.25),zlim=c(-10,60), ticktype = "detailed",main="WT BY4741",col= adjustcolor( "red", alpha.f = 0.2),colkey=F)

#par(new=TRUE)

WT.BY4741.media.pos<-WT.BY4741.media[WT.BY4741.media$Delta_LL>0,]

scatter3D(WT.BY4741.media.pos$gompShape,WT.BY4741.media.pos$gompRate,WT.BY4741.media.pos$Delta_LL,xlab="G",ylab="R",zlab="delta_LL",phi = 0, bty = "g",
          pch = 20, cex = 1, ticktype = "detailed",main="WT BY4741",col= adjustcolor( "blue", alpha.f = 0.2),add=TRUE,colkey=T)


#dev.off()

#pdf(paste("plots/", "WT-By4742-3D-plots-Delta_LL-G-R.pdf", sep=''))

WT.BY4742.media.neg<-WT.BY4742.media[WT.BY4742.media$Delta_LL<0,]

scatter3D(WT.BY4742.media.neg$gompShape,WT.BY4742.media.neg$gompRate,WT.BY4742.media.neg$Delta_LL,xlab="G",ylab="R",zlab="delta_LL",phi = 0, bty = "g",
           pch=20,cex = 1,xlim=c(0.04,0.35),zlim=c(-30,805), ticktype = "detailed",main="WT BY4742",col= adjustcolor( "red", alpha.f = 0.2),colkey=F)

#par(new=TRUE)

WT.BY4742.media.pos<-WT.BY4742.media[WT.BY4742.media$Delta_LL>0,]

scatter3D(WT.BY4742.media.pos$gompShape,WT.BY4742.media.pos$gompRate,WT.BY4742.media.pos$Delta_LL,xlab="G",ylab="R",zlab="delta_LL",phi = 0, bty = "g",
           pch=20,cex = 1, ticktype = "detailed",main="WT BY4742",col= adjustcolor( "blue", alpha.f = 0.2),add=TRUE,colkey=T)




#dev.off()







# We expect that larger CV (more noisy) 
#data will bring down LLH of both Gomeprtz and Weibull models. 
#However, Weibull’s decreasing LLH maybe slower than Gompertz. 










