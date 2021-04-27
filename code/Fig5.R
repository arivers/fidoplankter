######## 
# Name:			Fig5.R, a script to plot Figure 5 in the paper
# Description: 	A comparison of the buffer failure predictions from the numerical and analytical models.
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18
########

source("ModelFunc.R")
require(sfsmisc)
require(deSolve)

###Analysis

##1. Thalassiosira weissflogii:

#1.1 dynamic model
params <-c(	R		=	5.6E-6, 	# m
			D		=	9E-10, 		# m^2*s^-1
			Qmax	=	3.5E-5,
			Krho	=	5.1E-10, 	# mol * L^-1
			Vmax	=	1276, 		#nmol m^-2 d^-1 from sunda and huntsman 1997
			kfFeY	=	17,			# mol^-1 * sec^-1
			kdFeY	=	1.72E-6,	# sec^-1
			kfAFO1	=	4.1E+7,		# mol^-1 * sec^-1
			kdAFO1	=	2.8E-2,		# sec^-1
			kfAFO2	=	2E-6,
			kdAFO2	=	2.12E-6,
			koxFeII	=	6.2E-3,
			khvFeY 	=	4.3E-6,
			iue		=	0.98)		#iron use efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 uE, constant light	)
#state paramters, chemical concentrations are for chemically equilibrated media
state<-c(Fe=1.014213e-10,FeII=1.987063e-12,FeY=2.865067e-09,Y=1.000013e-05, cells=1E3, Q=7.5e-18, AFO1=1.768756e-11, AFO2=1.383751e-11)
times <-seq(0,2000000,1000)
outweiss<- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
names(outweiss)[10:15]<-c("rho","mu","FeYflux", "FeIIflux", "AFO1flux", "AFO2flux")

#1.2 uptake equilibrium model
dweiss			<- 11.2E-6						# m					Cell diameter; T. weissflogii = 2 x 5.6E-6, from Popp et al. (1998)
edta			<- 1e-5 						# mol L^-1			EDTA not chelated to Fe ~ all EDTA
feedta			<- 3e-9							# mol L^-1			Fe chelated to EDTA	~ all Fe
maxbiomass		<- 0.0002 						# mol L^-1			The maximum biomass to be plotted in
light			<- 500							#µE m^-2 s^-1

weiss<-staticbb(dweiss,edta=edta,Fetot=feedta,maxbiomass=maxbiomass,light=light)



##2. Thalassiosira pseudonana

#2.1 Dynamic model
params <-c(	R		=	2E-6, 		# m
			D		=	9E-10, 		# m^2*s^-1
			Qmax	=	1.9E-5,
			Krho	=	5.1E-10, 	# mol * L^-1
			Vmax	=	1276, 		#nmol m^-2 d^-1 from Sunda and Huntsman 1997
			kfFeY	=	17,			# mol^-1 * sec^-1
			kdFeY	=	1.72E-6,	# sec^-1
			kfAFO1	=	4.1E+7,		# mol^-1 * sec^-1
			kdAFO1	=	2.8E-2,		# sec^-1
			kfAFO2	=	2E-6,
			kdAFO2	=	2.12E-6,
			koxFeII	=	6.2E-3,
			khvFeY 	=	4.3E-6,
			iue	=	0.98)	 		#iron use efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 uE, constant light	)
#state paramters, chemical concentrations are for chemically equilibrated media
state<-c(Fe=1.014213e-10,FeII=1.987063e-12,FeY=2.865067e-09,Y=1.000013e-05, cells=1E3, Q=7.5e-18, AFO1=1.768756e-11, AFO2=1.383751e-11)
times <-seq(0,2000000,1000)
outpseudo<- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19))
names(outpseudo)[10:15]<-c("rho","mu","FeYflux", "FeIIflux", "AFO1flux", "AFO2flux")

#2.2 Uptake equilibrium model
dpseudo				<- 4E-6						# m					Cell diameter; T. pseudonana 

pseudo<-staticbb(dpseudo,edta=edta,Fetot=feedta,maxbiomass=maxbiomass,light=light)


#3 Calculate buffer failure points
pseudofail<-staticfailpoint(d=dpseudo,edta=edta,Fetot=feedta,maxbiomass=maxbiomass,light=light,failpct=0.1)
weissfail<-staticfailpoint(dweiss,edta=edta,Fetot=feedta,maxbiomass=maxbiomass,light=light,failpct=0.1)
wfd<-dynamicfailpoint(outweiss,0.1)
pfd<-dynamicfailpoint(outpseudo,0.1)



###PLOT

#plot Cells/mL against picomoles of Fe'
pdf("fig5.pdf", width=5,height=4.5,family="Helvetica")
par(ps=10)
plot(pseudo$Nlog/1000,pseudo$Feprimelog*1e12,
 log="x",
 las=1, 
 type="l",
 ylim=c(0,110),
 xlim=c(1,1e7),
ylab=expression("Fe',"~"pmol"~L^-1),
xlab=expression("Cells"~mL^-1),
col="grey70",
xaxt="n",
asp=0.75,
bty="l",
lwd=2)
aX <- axTicks(1); axis(1, at=aX, label= axTexpr(1, aX))
lines(outpseudo$cells/1000,outpseudo$Fe*1e12,lty=2,col="grey70",lwd=2)
lines(weiss$Nlog/1000,weiss$Feprimelog*1e12,lwd=2)
lines(outweiss$cells/1000,outweiss$Fe*1e12,lty=2,lwd=2)
abline(v=pseudofail[1], lty=1, col="grey70")
abline(v=pfd[2], lty=2, col="grey70")
abline(v=weissfail[1], lty=1)
abline(v=wfd[2],lty=2)
legend("bottomleft",legend=c(expression(paste(italic(T.~pseudonana)," equilibrium")),expression(paste(italic(T.~pseudonana)," dynamic")),expression(paste(italic(T.~weissflogii)," equilibrium")), expression(paste(italic(T.~weissflogii)," dynamic"))), lwd=2, lty=c(1,2,1,2), col=c("grey70","grey70","black","black"),bty="n",cex=.7)
dev.off()
