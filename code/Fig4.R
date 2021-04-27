######## 
# Name:			Fig4.R, a script to plot Figure 4 in the paper
# Description: 	A comparison of the growth rate predictions of the numerical model to growth rate data from Sunda and Huntsman (1995).
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18
########

require(deSolve)
source("ModelFunc.R")
sunda<-read.csv("sundahuntsman1995data.csv")
sunda$radius<-(0.75*sunda$cellvol.um3/pi)^(1/3)

weiss<- sunda[sunda$species=="T. weissflogii",]
pseudo<-sunda[sunda$species=="T. pseudonana",]

times <-seq(0,2000000,45000)

weiss.df<-NULL	
for( i in 1:200){
radius<- 1e-6*mean(weiss$radius,na.rm=T)
Qmax<-1.9E-5 # T. weissflogii
iron<-i*2e-09
params <-c(	R		=	radius, 	# m
			D		=	9E-10, 		# m^2*s^-1
			Qmax	=	Qmax,
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
			iue		=	.57)  		#iron use efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 µE m^-2 s^-1, 14h:10h light:dark cycle 	)
state<-c(Fe=0,FeII=0,FeY=iron,Y=1e-4, cells=1E3, Q=0, AFO1=0, AFO2=0)

out1 <- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
names(out1)[10:11]<-c("rho","mu")
out1<-out1[-1:-2,]
weiss.df<-rbind(weiss.df,data.frame(max(out1$Fe),max(out1$mu)))
}


pseudo.df<-NULL
for( i in 1:200){
radius<- 1e-6*mean(pseudo$radius,na.rm=T)
Qmax<-3.5E-5 #T. pseudanana

iron<-i*2e-09
params <-c(	R		=	radius, 	# m
			D		=	9E-10, 		# m^2*s^-1
			Qmax	=	Qmax,
			Krho	=	5.1E-10, 	# mol * L^-1
			Vmax	=	1276, 		#nmol m^-2 d^-1 from Sunda and Huntsman (1997)
			kfFeY	=	17,			# mol^-1 * sec^-1
			kdFeY	=	1.72E-6,	# sec^-1
			kfAFO1	=	4.1E+7,		# mol^-1 * sec^-1
			kdAFO1	=	2.8E-2,		# sec^-1
			kfAFO2	=	2E-6,
			kdAFO2	=	2.12E-6,
			koxFeII	=	6.2E-3,
			khvFeY 	=	4.3E-6,
			iue		=	.57)#0.5712963 		#iron use Efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 uE, 14:10 	)
state<-c(Fe=0,FeII=0,FeY=iron,Y=1e-4, cells=1E3, Q=0, AFO1=0, AFO2=0)
out1 <- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
names(out1)[10:11]<-c("rho","mu")
out1<-out1[-1:-2,]
pseudo.df<-rbind(pseudo.df,data.frame(max(out1$Fe),max(out1$mu)))
}

pdf("fig4.pdf", width=5,height=4.5,family="Helvetica")
par(ps=10)
plot(pseudo.df$max.out1.Fe.*1e12,pseudo.df$max.out1.mu.*86400,log="x", bty="l",las=1,type="n",ylim=c(0,2),xlab=expression("Fe',"~"pmol"~L^-1),ylab= expression("Growth rate,"~d^-1))
lines(weiss.df$max.out1.Fe.*1e12,weiss.df$max.out1.mu.*86400,lty=1,lwd=2)
lines(pseudo.df$max.out1.Fe.*1e12,pseudo.df$max.out1.mu.*86400,lty=1,lwd=2,col="grey70")
points(weiss$freeFe.pm,weiss$mu.perD, pch=1)
points(pseudo$freeFe.pm,pseudo$mu.perD, pch=2, col="grey70")
legend("topleft",legend=c(expression(paste(italic(T.~weissflogii)," dynamic model")),expression(paste(italic(T.~pseudonana)," dynamic model")),expression(paste(italic(T.~weissflogii)," data")),expression(paste(italic(T.~pseudonana)," data"))) ,lty=c(1,1,0,0),col=c("black",col="grey70","black",col="grey70"),pch=c(NA,NA,1,2),bty="n")
dev.off()
