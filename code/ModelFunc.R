######## 
# Name:			ModelFunc.R, Blown buffer calculator functions
# Description: 	This is the master file for both the static and dynamic buffer failure calculation functions
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18
########


######## Dynamic model solved numerically with an ODE solver
#usage: dynamicbb(t,state,parameters)
#   	   t:  A vector of timepoints (in seconds) at which to evaluate the system
#      state:  Initial values for the state variables in form c(Fe=#,FeII=#,FeY=#,Y=#,cells=#,Q=#,AFO1=#,AFO2=#)
# parameters:  A numeric vector in the form: 
	#     c(R		=	,# m,
	#		Qmax	=	,# mol cell^-1,
	#		Krho	=	,# mol * L^-1,
	#		Vmax	=	,# nmol m^-2 d^-1,
	#		kfFeY	=	,# L * mol^-1 * sec^-1
	#		kdFeY	=	,# sec^-1
	#		kfAFO1	=	,# L * mol^-1 * sec^-1
	#		kdAFO1	=	,# sec^-1
	#		kfAFO2	=	,# sec^-1
	#		kdAFO2	=	,# sec^-1
	#		koxFeII	=	,# sec^-1
	#		khvFeY 	=	,# sec^-1
	#		iue		=	)# mol C (mol Fe)^-1 s^-1
########	
require(deSolve)
dynamicbb <-function(t,state,parameters)
{ with(as.list(c(state,parameters)),{
		
		##variables
		cellarea 	<- 	4*pi*R^2
		rhomax  	<-  cellarea*Vmax/(86400*1E+9) #vmax calculation L* s^-1 (6.33e-11)
		rho			<-  rhomax * Fe/(Fe + Krho)
		Ccell   	<-  15 * (4/3) * pi* R^3 * 1000 # estimated from Sunda et al. (2005), 15 molC / L cell culture * vol	# mol cell-1
		femu		<- 	max(iue*((Q/Ccell) - 3E-6),0)
		mu			<-  min(femu, iue*Qmax)
		Fei			<-  AFO1 + Fe
		FeYflux		<-	kdFeY*FeY
		FeIIflux 	<-	koxFeII*FeII
		AFO1flux	<-	kdAFO1*AFO1		
		AFO2flux	<-	kdAFO2*AFO2
		
		##State values	
		dFe 	<-	(koxFeII*FeII
					+kdFeY*FeY 
					+kdAFO1*AFO1
					+kdAFO2*AFO2
					-kfFeY*Fe*Y
					-kfAFO1*Fe*Fei
					-cells*rho)
		
		dFeII	<-	(khvFeY*FeY
					-koxFeII*FeII)
							
		dFeY 	<- 	(-kdFeY*FeY
					-khvFeY*FeY 
					+kfFeY*Fe*Y)
					
		dY		<- 	(kdFeY*FeY
					+khvFeY*FeY 
					-kfFeY*Fe*Y)
		
		dcells 	<-	cells*mu
		
		dQ		<-	rho-mu*Q
					
		dAFO1	<-	(kfAFO1*Fe*Fei
					-kdAFO1*AFO1
					-kfAFO2*AFO1)
					
		dAFO2	<-	(kfAFO2*AFO1
					-kdAFO2*AFO2)
	 	
	 		 		
		list(c(dFe,dFeII,dFeY,dY,dcells,dQ,dAFO1,dAFO2),c(rho,mu,FeYflux, FeIIflux, AFO1flux, AFO2flux))  
		})
}

## This function estimates the point of buffer failue in Cells per mL for the dynamic model
dynamicfailpoint<-function(out,failpct){
out$blown <- (1-failpct) >= out$Fe / max(out$Fe)
out<-out[-1:-100,]
out<-out[out$blown==T,]
c(out[1,2],out[1,6]/1000,out[1,1])
}


#usage individual values for: the cell diameter (m),EDTA concentration, FeEDTA concetration (M), the maximum biomass to be modeled (µmol/L), light intensity (µE * m^-2 * s^-1)

#####
#static solver
#usage staticbb(
	#	d			= 			,# Cell diameter (m)
	#	edta		=			,# EDTA concentration
	#	Fetot		=			,# Fe concentration (M)
	#	light		= 500  		,# light intensity (µE * m^-2 * s^-1)
	#	Krho		= 5.1E-10	,# mol * L^-1,
	#	Vmax		= 1276		,# nmol m^-2 d^-1,
	#	kdFeY		= 17		,# sec^-1
	#	kfFeY		= 1.72e-6	,# L * mol^-1 * sec^-1
	#	khvFeY 		= 4.3e-6  	,# sec^-1
	#	failpct		= 0.1	,# change from original that is considered buffer failure, e.g. 0.1 is 10%
	#	maxbiomass	= 0.0002	)# The maximum biomass to be modeled (mol/L)
# 	returns a data frame, the first three columns are associated, as are the last two
# 	N,	biomass,	Feprime,	Nlog,	Feprimelog

staticbb<-function(d, edta, Fetot, light=500, Krho= 5.1E-10,Vmax=1276,kfFeY=17, kdFeY=1.72e-6,khvFeY=4.3e-6, failpct=0.1, maxbiomass=0.0002)
{
##variables
		kd			<-	kdFeY # alias for dissociation rate constant of FeY complex, sec^-1 
		kf			<-  kfFeY # allias for association rate constant of Fe  and Y L * mol^-1 * sec^-1
		khv			<-	khvFeY # alias for photoreduction rate constant of FeY, sec^-1
		kdprime		<- kdFeY + (khv * light / 500)		# s^-1				The total FeEDTA dissociation rate
		r			<-	d/2
		cellarea 	<- 	4*pi*r^2
		rhomax  	<-  cellarea*Vmax/(86400*1E+9) #vmax calculation L* s^-1 (6.33e-11)
		#rho			<-  rhomax * Fe/(Fe + Krho)
		feedta		<-  Fetot # we make the assumption that FeEDTA ~ total Fe
		Ccell   	<-  15 * (4/3) * pi* r^3 * 1000 # estimated from sunda 2005, 15 molC / L cell culture * vol	# mol cell-1
# Estimating the blown buffer point	
		abiotic <- (kdprime*feedta)/(kf*edta) #
		fefail	<- ((1-failpct)*kdprime*feedta)/(kf*edta) # 90% of the abiotic Fe concentration
		failN   <-  -(edta*fefail^2*kf - fefail*feedta*kdprime + (edta*fefail*kf -feedta*kdprime)*Krho)/(fefail*rhomax)#
		failNmL	<- failN/1000 #
		failC	<- failN*Ccell*1e6;
		
#Generate a regularly spaced vector of values from 1 cell/L to the cell density equal to the maximum desired biomass value
N				<- seq(from=1, to=(maxbiomass/Ccell), length.out=200) 	# cells L^-1
biomass			<- Ccell * N
Feprime			<- -1/2*(Krho*edta*kf - feedta*kdprime + N*rhomax -sqrt(Krho^2*edta^2*kf^2 + 2*Krho*edta*feedta*kdprime*kf +feedta^2*kdprime^2 + N^2*rhomax^2 + 2*(Krho*edta*kf - feedta*kdprime)*N*rhomax))/(edta*kf)

#generate a log spaced vector of values from 1 cell/L to the cell density equal to the maximum desired biomass value
Nlog			<- 10^(seq(from=log10(1), to=log10(maxbiomass/Ccell),length.out=200)) # cells L^-1
#Calculate the Fe' for each cell concentration
Feprimelog<- -1/2*(Krho*edta*kf - feedta*kdprime + Nlog*rhomax -sqrt(Krho^2*edta^2*kf^2 + 2*Krho*edta*feedta*kdprime*kf +feedta^2*kdprime^2 + Nlog^2*rhomax^2 + 2*(Krho*edta*kf - feedta*kdprime)*Nlog*rhomax))/(edta*kf)
		
data.frame(N,biomass,Feprime,Nlog,Feprimelog)

}

#####
#static solver
#usage staticfailpoint(
	#	d			= 			,# Cell diameter (m)
	#	edta		=			,# EDTA concentration
	#	Fetot		=			,# Fe concentration (M)
	#	light		= 500  		,# light intensity (µE * m^-2 * s^-1)
	#	Krho		= 5.1E-10	,# mol * L^-1,
	#	Vmax		= 1276		,# nmol m^-2 d^-1,
	#	kdFeY		= 17		,# sec^-1
	#	kfFeY		= 1.72e-6	,# L * mol^-1 * sec^-1
	#	khvFeY 		= 4.3e-6  	,# sec^-1
	#	failpct		= 0.1	,# change from original that is considered buffer failure, e.g. 0.1 is 10%
	#	maxbiomass	= 0.0002	)# The maximum biomass to be modeled (mol/L)
# 	returns a vector with the cells per mL and biomass that buffer failure occurs at,  c(failNmL,failC) 
staticfailpoint<-function(d, edta, Fetot, light=500, Krho= 5.1E-10,Vmax=1276, kdFeY=1.72e-6,kfFeY=17,khvFeY=4.3e-6, failpct=0.1, maxbiomass=0.0002)
{
##variables
		kd			<-	kdFeY # alias for dissociation rate constant of FeY complex, sec^-1 
		kf			<-  kfFeY # allias for association rate constant of Fe  and Y L * mol^-1 * sec^-1
		khv			<-	khvFeY # alias for photoreduction rate constant of FeY, sec^-1
		kdprime		<- kdFeY + (khv * light / 500)		# s^-1				The total FeEDTA dissociation rate
		r			<-	d/2
		cellarea 	<- 	4*pi*r^2
		rhomax  	<-  cellarea*Vmax/(86400*1E+9) #vmax calculation L* s^-1 (6.33e-11)
		#rho			<-  rhomax * Fe/(Fe + Krho)
		feedta		<-  Fetot # we make the assumption that FeEDTA ~ total Fe
		Ccell   	<-  15 * (4/3) * pi* r^3 * 1000 # estimated from sunda 2005, 15 molC / L cell culture * vol	# mol cell-1
# Estimating the blown buffer point	
		abiotic <- (kdprime*feedta)/(kf*edta) #
		fefail	<- ((1-failpct)*kdprime*feedta)/(kf*edta) # 90% of the abiotic Fe concentration
		failN   <-  -(edta*fefail^2*kf - fefail*feedta*kdprime + (edta*fefail*kf -feedta*kdprime)*Krho)/(fefail*rhomax)#
		failNmL	<- failN/1000 #
		failC	<- failN*Ccell*1e6;
		
#Generate a regularly spaced vector of values from 1 cell/L to the cell density equal to the maximum desired biomass value
N				<- seq(from=1, to=(maxbiomass/Ccell), length.out=200) 	# cells L^-1
biomass			<- Ccell * N
Feprime			<- -1/2*(Krho*edta*kf - feedta*kdprime + N*rhomax -sqrt(Krho^2*edta^2*kf^2 + 2*Krho*edta*feedta*kdprime*kf +feedta^2*kdprime^2 + N^2*rhomax^2 + 2*(Krho*edta*kf - feedta*kdprime)*N*rhomax))/(edta*kf)

#Generate a log spaced vector of values from 1 cell/L to the cell density equal to the maximum desired biomass value
Nlog			<- 10^(seq(from=log10(1), to=log10(maxbiomass/Ccell),length.out=200)) # cells L^-1
#Calculate the Fe' for each cell concentration
Feprimelog<- -1/2*(Krho*edta*kf - feedta*kdprime + Nlog*rhomax -sqrt(Krho^2*edta^2*kf^2 + 2*Krho*edta*feedta*kdprime*kf +feedta^2*kdprime^2 + Nlog^2*rhomax^2 + 2*(Krho*edta*kf - feedta*kdprime)*Nlog*rhomax))/(edta*kf)
		
c(failNmL,failC)

}


#This function estimates cell biomass in Mol C when given the radius of a cell in meters.
biomassmol<-function(r){
	vol.cell		<- (4/3) * pi *r^3 * 1000		# L cell^-1			Volume of an individual cell
	C.mol.per.vol	<- 15						# mol L of cell vol An estimate of moles Carbon per cel volume values from Sunda et al. (2005)
	vol.cell * C.mol.per.vol }
	
	
