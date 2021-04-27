######## 
# Name:			SuppFig1.R, a script to plot Supp. Figure 1 in the paper.
# Description: 	Three estimates of Fe area-dependent iron update compared to data from Sunda and Huntsman (1995)
# Citation:		Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
# URL:			http://fidoplankter.uga.edu
# License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
# Modified:		2014-08-18
########

require(lattice)
require(sfsmisc)
data.o<-read.csv("sundahuntsman1995data.csv") # read data from Sunda and Huntsman (1995)
#data<-data.o[data.o$experiment %in% c("126","127a","127b","98","123"),] #select only a subset of experiments
data<-data.o 	# select all experiments
data$radius <- ((3*data$cellvol.um3)/(4*pi))^(1/3) # calculate cell radius from the cell volume in the table
data$up.percell<-data$uptakerate*(data$cellvol.um3/1e15) 	#convert their uptake rate to µmol Fe per cell per day
data$cellsurf<- 4*pi*(data$radius/1e6)^2  	#calculate the cell surface area (m^2 per cell)
data$uptakesurf <- data$up.percell/data$cellsurf 	#uptake rate in µmol m^-2 d^-1

#The estimation of cellular uptake based on Husdon and Morel (1993) that rate=kf * Mt * [Fe']
At <- 1.66E-17 			#area of an iron transporter from as assumed in Hudson and Morel (1990), m^2 
#fecA<- pi*20E-10^2 	#area of the FecA transporter from: Wyatt W. Yue, Sylvestre Grizot, Susan K. Buchanan. 2003. Journal of Molecular Biology. 332(2) 353-368. m^2 
kf<-2E6 				#the Fe-Mt formation rate constant, L mol^-1 s^-1
kuprime<-(.15*kf)/(At*6.022e23)*0.001#uptake rate constant, m  s^-1
conv<- 86.4 			#conversion factor for pM* m s^-1 -> µmol m^-2 d^-1
kuprime.mud<-kuprime*conv # kuprime times a scaling constant for picomoles Fe versus µmoles m-2 d-1
s05ku<-0.29*conv 		#m s^-1 from Sunda et al. (2005)

#Sunda and Huntsman (1997) uptake parameters
Vmax<-1276 # nmol m-^2 d^-1
Ks<- 0.51 # nM


#Plot data and add regression line and line from theoretical Fe uptake rates
pdf("suppFig1.pdf", family="Helvetica")
par(mar=c(5, 5, 4, 2) + 0.1)
xlim=c(1,1000)
ylim=c(0.01,10)
plot(uptakesurf ~ freeFe.pm,log="xy",data=data, bty="l",type="n", xlim=xlim, ylim=ylim, las=1, xlab=expression(paste(Fe*minute,",",~~pmol~L^-1)),ylab=expression(paste(  "Iron uptake per cell surface, ", mu * mol~m^-2~d^-1)))
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="106" & data$method=="total",],type="b", pch=18, col="grey70")#E. huxleyi
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="126" & data$method=="intracellular",],type="b", pch=6, col="grey70")#E. huxleyi
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="127a" & data$method=="total",],type="b", pch=19, col="grey70")#E. huxleyi
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="127b" & data$method=="total",],type="b", pch=20, col="grey70")#E. huxleyi
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="98" & data$method=="intracellular",],type="b", pch=1, col="grey70")#T. oceanica
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="123" & data$method=="intracellular",],type="b", pch=0, col="black")#T. pseudonana
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="113" & data$method=="intracellular",],type="b", pch=2, col="black")#T. weissflogii
lines(uptakesurf ~ freeFe.pm,data=data[data$experiment=="137" & data$method=="intracellular",],type="b", pch=5, col="grey70")#P. minimum
abline(a=0,b=kuprime.mud, col="blue",untf=T) # Sunda and Huntsman (1995)
abline(a=0,b=(s05ku*0.0864),col="green",lty=2,untf=T) # Sunda et al. (2005) 
curve((Vmax*x)/(1E3*(x+Ks*1E3)), from=1, to=1000, col="darkred", add=T) # Sunda and Huntsman (1997)
legend("topright", cex=0.8,
legend=c(expression(paste("Sunda et al. (2005), ", k[u]*minute == 2.9 %*% 10^-1 ~ m ~s^-1)),
		expression(paste("Sunda and Huntsman (1995), ", k[u]*minute ==3.0 %*% 10^-5~m~s^-1)),
		"Sunda and Huntsman (1997), Michaelis-Menten type",
		expression(paste(italic(E.~huxleyi), " exp. 106")),
		expression(paste(italic(E.~huxleyi), " exp. 126")),
		expression(paste(italic(E.~huxleyi), " exp. 127a")),
		expression(paste(italic(E.~huxleyi), " exp. 127b")),
		expression(paste(italic(T.~oceanica), " exp. 98")),
		expression(paste(italic(T.~pseudonana), " exp. 123")),
		expression(paste(italic(T.~weissflogii), " exp. 113")),
		expression(paste(italic(P.~minimum), " exp. 137"))
		),
col=c("green","blue","darkred","grey70","grey70","grey70","grey70","grey70","black","black","grey70"),
lty=c(2,1,1,1,1,1,1,1,1,1,1),
pch=c(NA,NA,NA,18,6,19,20,1,0,2,5),
bty="n")
dev.off()
