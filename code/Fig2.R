######## 
# Name:			Fig2.R, a script to plot Figure 2 in the paper
# Description: 	A numerical simulation of phytoplankton growth and iron dynamics in EDTA buffered media
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18

########

require(deSolve)
require(sfsmisc)
source("ModelFunc.R")

# Numerically estimate phytoplankton growth, run Monte Carlo simulations and plot the results.


## estimate growth at high resolution
params <-c(	R		=	5.6E-6, 	# m
			Qmax	=	1.9E-5,
			Krho	=	5.1E-10, 	# mol * L^-1
			Vmax	=	1276, 		#nmol m^-2 d^-1 from sunda and huntsman 1997
			kfFeY	=	17,			# mol^-1 * sec^-1
			kdFeY	=	1.72E-6,		# sec^-1
			kfAFO1	=	4.1E+7,			# mol^-1 * sec^-1
			kdAFO1	=	2.8E-2,		# sec^-1
			kfAFO2	=	2E-6,
			kdAFO2	=	2.12E-6,
			koxFeII	=	6.2E-3,
			khvFeY 	=	4.3E-6,
			iue		=	0.98)#0.5712963 		#iron use efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 uE, 14:10 	)

# note starting values are from a kinetic model of media allowed to equilibrate for 24 hours. That model had the starting values: state <-c(Fe=0,FeII=0,FeY=3E-9,Y=1E-5,cells=1E3,Q=0,AFO1=0,AFO2=0)
# state parameters of T weiss 1E-5, Y 3E-9 Fe
state<-c(Fe=1.017e-10,FeII=2.067e-12,FeY=2.866e-09,Y=1e-5, cells=1E3, Q=1.245e-16, AFO1=1.778e-11, AFO2=1.447e-11)
times <-seq(0,2000000,1000)
out1<- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
names(out1)[10:15]<-c("rho","mu","FeYflux", "FeIIflux", "AFO1flux", "AFO2flux")

## Run uncertainty analysis
for (i in 1:1000){
sc.R<- rnorm(1,mean=5.6E-6,sd=5.6E-7)
while(sc.R<=0){sc.R<- rnorm(1,mean=5.6E-6,sd=5.6E-7)}

sc.Qmax<- rnorm(1,mean=1.9E-5,sd=1.9E-6)
while(sc.Qmax<=0){sc.R<- rnorm(1,mean=1.9E-5,sd=1.9E-6)}

sc.Krho<-rnorm(1,mean=5.1E-10, sd=1.275E-10) #25%cv guess from data in from Sunda and Huntsman 1997
while(sc.Krho<=0){sc.Krho<- rnorm(1,mean=5.1E-10, sd=1.275E-10)}

sc.Vmax<-rnorm(1,mean=1276,sd=128)
while(sc.Vmax<=0){sc.Vmax<-rnorm(1,mean=1276,sd=128)}

sc.kfFeY<-rnorm(1,mean=17,sd=1.7) # from hudson 1992
while(sc.kfFeY<=0){sc.kfFeY<-rnorm(1,mean=17,sd=1.7)}

sc.kdFeY<-rnorm(1,mean=1.72E-6,sd=1.72E-7)# sunda 2003
while(sc.kdFeY<=0){sc.kdFeY<-rnorm(1,mean=1.72E-6,sd=1.72E-7)}#from sunda 2003

sc.kfAFO1<- rnorm(1,mean=4.1E+7,sd=1.1E+7)#cv=measured from rose and waite 2003
while(sc.kfAFO1 <= 0) {sc.kfAFO1<- rnorm(1,mean=4.1E+7,sd=1.1E+7)}

sc.kdAFO1<- rnorm(1,mean=2.8E-2, sd= 1.4E-2)#cv=50%
while(sc.kdAFO1 <= 0){sc.kdAFO1<- rnorm(1,mean=2.8E-2, sd= 1.4E-2)}

sc.kfAFO2<-rnorm(1,mean=2E-6,sd=2E-7)
while(sc.kfAFO2<=0){sc.kfAFO2<-rnorm(1,mean=2E-6,sd=2E-7)}

sc.kdAFO2<-rnorm(1,mean=2.12E-6,sd=2.12E-7) 
while(sc.kdAFO2<=0){sc.kdAFO2<-rnorm(1,mean=2.12E-6,sd=2.12E-7)} #cv10%

sc.koxFeII<-rnorm(1,mean=6.2E-3, sd=3.1E-3)# from sunda 2003 ph 8.15 cv= 50%
while(sc.koxFeII<=0){sc.koxFeII<-rnorm(1,mean=6.2E-3, sd=3.1E-3)}

sc.khvFeY<-rnorm(1,mean=4.3e-6,sd=4.3e-7)#from sunda 2003 ph 8.15 cv=10%
while(sc.khvFeY<=0){sc.khvFeY<-rnorm(1,mean=4.3e-6,sd=4.3e-7)}

sc.iue<-rnorm(1,mean=0.98,sd=0.098)#from sunda 1997  cv=10%
while(sc.iue<=0){sc.iue<-rnorm(1,mean=0.98,sd=0.098)}

	
params <-c(	R		=	sc.R,		#5.6E-6, 	# m
			Qmax	=	sc.Qmax,
			Krho	=	sc.Krho,	#0.51E-9, 	# mol * L^-1
			Vmax	=	sc.Vmax,	#1276, 		#nmol m^-2 d^-1 from sunda and huntsman 1997
			kfFeY	=	sc.kfFeY,	#17,			# mol^-1 * sec^-1
			kdFeY	=	sc.kdFeY,	#1.76E-6,		# sec^-1
			kfAFO1	=	sc.kfAFO1,	#4.1E+7,			# mol^-1 * sec^-1
			kdAFO1	=	sc.kdAFO1,	#2.8E-2,		# sec^-1
			kfAFO2	=	sc.kfAFO2,	#2E-6,
			kdAFO2	=	sc.kdAFO2,	#2.12E-6,
			koxFeII	=	sc.koxFeII,	#6.2E-3,
			khvFeY 	=	sc.khvFeY,	#4.3E-6		#s^-1
			iue	 	=	sc.iue		#0.5712963 		#iron use efficiency mol C mol Fe-1 s-1,  corresponds to 500 uE, 14:10 		
					)

times <-seq(0,2000000,10000)
times <-times[-1]
times<- c(0,1000,5000,times)

if(i==1){data<- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))} else {temp <- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
	data <- rbind(data,temp)}
}
# process the model output
names(data)[10:15]<-c("rho","mu","FeYflux", "FeIIflux", "AFO1flux", "AFO2flux")
cn<- colnames(data)
times<-levels(factor(data$time))
summ<-as.data.frame(times)
names(summ)[1]<-"time"
for (i in 2:length(cn)){
col<-paste("data$",cn[i],sep="")
col25name<-paste(cn[i],"25",sep="")
col75name<-paste(cn[i],"75",sep="")
col25<- as.data.frame(tapply(data[,i],data$time,quantile,prob=0.25))
names(col25)[1]<-col25name
summ<- cbind(summ,col25)
col75<- as.data.frame(tapply(data[,i],data$time,quantile,prob=0.75))
names(col75)[1]<-col75name
summ<-cbind(summ,col75)
}
summ$time <- as.numeric(as.vector(summ$time))



##plotting the results

#data cleanup
all<- out1
summ<-summ[summ$time>0,]

#calculating buffer failure

all$blown <- .90 >= all$Fe / max(all$Fe)
allblown<-all[-1:-100,]
allblown <-allblown[allblown$blown==T,]
blowntime <- min(allblown$time)/86400

#cut data with a cell density of more than 10,000 cells /mL
all<-all[all$cells < 2e7, ]
rono <-dim(all)[1]
summ<-summ[summ$time < all$time[rono],]


pdf(file="fig2.pdf",width=4, height=9, family="Helvetica")
par(mfrow=c(3,1))
par(mar=c(2.1, 7.1, 4.1, 2.1))
##CELLS
plot(c(all$time/86400, summ$time/86400,summ$time/86400),c(all$cells/1000,summ$cells25/1000,summ$cells75/1000),log="y", ylab="", xlab="", type="n", bty="l",xaxt="n",ylim=c(1,1E5), axes=FALSE,frame=TRUE)
 aY <- axTicks(2)
 axis(2, at=aY, label= axTexpr(1, aY),las=2)
 #axis(1)
 mtext(expression(Cells~mL^-1), side=2, line=5)
x<- c(summ$time/86400,rev(summ$time/86400))
y<- c(summ$cells25/1000,rev(summ$cells75/1000))
polygon(x,y,col="grey80",border=NA)
abline(v=blowntime,lty="dotted")
text(0,5E4,label="A.", pos=4,cex=1.5)
lines(all$time/86400,all$cells/1000)



 ##Fe Concentration plots
 par(mar=c(2.1, 7.1, 0, 2.1))
plot(c(all$time/86400,summ$time/86400,summ$time/86400),c(all$Fe*1e12,summ$Fe25*1e12,summ$Fe75*1e12),log ="y",ylab="",xlab="", type="n",bty="l",xaxt="n",ylim=c(2,2E2),axes=FALSE,frame=TRUE)
x<- c(summ$time/86400,rev(summ$time/86400))
y<- c(summ$Fe25*1e12,rev(summ$Fe75*1E12))
polygon(x,y,col="grey80",border=NA)
abline(v=blowntime,lty="dotted")
abline(h=7e-10, lty="dashed") # milero
lines(all$time/86400,all$Fe*1E12)
#aY <- axTicks(2)
axis(2,las=1)#, at=aY, label= axTexpr(1, aY),las=2)
#axis(1)  
 #abline(h=, lty="dashed") # calculated fron sunda 2003
#text(7,1E-9, label="Measured solubility limit")
text(0,140,label="B.", pos=4,cex=1.5)
mtext(expression("Fe',"~pmol~L^-1), side=2, line=5,cex=0.8)

###FLUX PLOTS
 par(mar=c(4.1, 7.1, 0, 2.1))
 
all$tfluxall	<- all$FeYflux + all$FeIIflux #+ all$AFO2flux +all$AFO1flux
all$afoflux 	<- all$AFO2flux +all$AFO1flux
summ$tflux25 	<- summ$FeYflux25 + summ$FeIIflux25 #+ summ$AFO2flux25 +summ$AFO1flux25
summ$tflux75	<- summ$FeYflux75 + summ$FeIIflux75 #+ summ$AFO2flux75 +summ$AFO1flux75
 

 plot(c(all$time/86400,all$time/86400,all$time/86400, all$time/86400,summ$time/86400),
 c(all$tfluxall*1E+12,all$rho*all$cells*1E+12,all$FeYflux*1E+12,all$FeIIflux*1E+12,summ$tflux75*1E+12),
  ylab="",
  xlab="Time, d",
  type="n",
  bty="l",
  log="y",
  axes=FALSE,
  frame=TRUE)
x<- c(summ$time/86400,rev(summ$time/86400))
y<- c(summ$rho25*summ$cells25*1E+12,rev(summ$rho75*summ$cells75*1E+12))
polygon(x,y,col="grey80",border=NA)
 x<- c(summ$time/86400,rev(summ$time/86400))
 y<- c(summ$tflux25*1E+12,rev(summ$tflux75*1E+12))
 polygon(x,y,col="grey50",border=NA)
abline(v=blowntime,lty="dotted")
 lines(all$time/86400,all$rho*all$cells*1E+12, col="black")
 lines(all$time/86400,all$tfluxall*1E+12, col="black", lty=2)
 lines(all$time/86400,all$FeYflux*1E+12, col="black",lty=3)
 lines(all$time/86400,all$FeIIflux*1E+12 ,col="black",lty=4)
aY <- axTicks(2)
axis(2, at=aY, label= axTexpr(1, aY),las=2)
axis(1) 

mtext(expression("Flux, pmol"~s^-1),side=2,line=5,cex=0.8)
text(0,8E-3,label="C.", pos=4,cex=1.5)
legend("bottomright",legend=c("Cellular uptake", "Total Fe flux","EDTA dissociation", "Fe(II) oxidation"),lty=c(1,2,3,4),lwd=1,bty="n", cex=0.8)

dev.off()

