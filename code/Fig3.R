######## 
# Name:			Fig3.R, a script to plot Figure 3 in the paper
# Description: 	A numerical simulation of phytoplankton growth and iron dynamics in EDTA buffered media
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18
########
require(deSolve)
require(sfsmisc)
source("ModelFunc.R")


## Estimate the relationship between diameter and buffer failure
fail<-NULL # create variable to store cell denity at buffer failure point
failtime<-NULL # create variable to store time at buffer failure point
Rvec<-seq(0.2E-6, 10E-6,length.out=40) #create a vector of the radii to test
state <-c(Fe=0,FeII=0,FeY=3E-8,Y=1E-4,cells=1E3,Q=0,AFO1=0,AFO2=0) # initital state variables

for (i in 1:40){ # loop through all the Radii in Rvec
#set parameters for the model runs
params <-c(	R		=	Rvec[i], 	# m
			Qmax	=	1.9E-5,
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
			iue		=	0.98)		#0.5712963 		#iron use efficiency mol C (mol Fe)^-1 s^-1,  corresponds to 500 uE, 14:10 	)

times <-seq(0,2000000,1000) # times to record data
times <-times[-1]	# remove first timepoint
times<- c(0,100,500,times) # add three closly spaced timepints to the beginning of the simulation
out1<- as.data.frame(ode(state,times,dynamicbb,params,atol=1e-19)) # solve the ode's
names(out1)[10:15]<-c("rho","mu","FeYflux", "FeIIflux", "AFO1flux", "AFO2flux") # name extra columns in the ODE solver output
fail[i]<-dynamicfailpoint(out1,0.1)[2] # record the cell density at the failure point
failtime[i]<-dynamicfailpoint(out1,0.1)[3] # record the time of failure
}

#function to calculate biomass when given vectors of cell radius and cell number
biomassvect<-function(Rvec,cells){
	C.mol.per.vol	<- 15							# mol L of cell vol An estimate of moles Carbon per cel volume values from Sunda 2005
	vol.cell		<- (4/3) * pi *Rvec^3 * 1000		# L cell^-1			Volume of an individual cell (1000 converts m^3 to L)
	biomass.mol.est	<- vol.cell * C.mol.per.vol 	# mol cell^-1 		Moles of carbon per individual cell
	biomass.mol.est *cells*1000 # biomass in mol C per L of culture (1000 converts cells per mL into cells per L)
}

#calculate biomass at buffer failure
biomass<-biomassvect(Rvec,fail)

## Plot the results
pdf(file="fig3.pdf",width=6, height=4.5, family="Helvetica") #open pdf writer
par(mar=c(6.1, 6.1, 2.1, 6.1)) # expand the margins

plot(Rvec*2e6,fail,type="l", lty=1,lwd=2, log="y", axes= FALSE, frame=TRUE,ylab="",xlab="",ylim=c(2e3,2e7)) # create the initial plot of Fe' vs Diameter
aY1 <- axTicks(2,axp=c(3,7,1))  # define the first y axis
axis(2, at=aY1, label= axTexpr(2, aY1),las=2) # plot the first y axis
axis(1) # plot the first x axis
par(new=TRUE) # continue adding objects to the same plot

plot(Rvec*2e6,biomass,type="l",lty=2,,lwd=2,xaxt="n",yaxt="n",xlab="",ylab="",log="y",ylim=c(1e-6,1e-2)) # add a second plot of Biomass vs Diamter
aY2<-axTicks(4) # define the second y axis
axis(4, at=aY2, label= axTexpr(4, aY2),las=2) # plot the second y axis

mtext(expression("Biomass, mol carbon"~L^-1),side=4,line=4) # label the second y axis
mtext(expression("Cell density, cells"~mL^-1),side=2,line=4) # label the first y axis
mtext(expression(paste("Cell diameter, ",mu,"m")),side=1,line=3) # label the x axis
legend("topright",lty=c(1,2),lwd=2,legend=c("Cell density","Biomass"), bty="n") # print a legend
dev.off() # close pdf writer
