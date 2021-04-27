// Name:		calculator.js
// Description: A calculator to estimate the iron buffering capacity of a FeEDTA medium for growing phytoplankton
// Citation:	Rivers, Adam R., Rose, Andrew L., Webb, Eric A.(2013) An online calculator for marine phytoplankton iron culturing experiments. Journal of Phycology. 49(5) 1017-1021.
// URL:			http://www.marsci.uga.edu/fidoplankter
// License:		©2012 University of Georgia Research Foundation, Inc. Non-commercial use only.
// Modified:	2013-10-09
// Dependency:	requires jquery, flot.js (http://www.flotcharts.org/)

     var feplotting =   function() {
		//form input variables
		var d = ($('#celld').val())/1e6;
		var r = d/2;		// m			phytoplankton diameter in micrometers
		var edta = $('#EDTA').val();		// mol L^-1		EDTA not chelated to Fe ~ all EDTA
		var feedta = $('#Fe').val();		// mol L^-1		Fe chelated to EDTA	~ all Fe
		var light = $('#light').val();		// uE m^-2 s^-1 Light intensity	
		var hlight = $('#hlight').val();
		//user input variables
		var maxbiomass= 0.0002; 			// mol C L^-1	The maximum biomass to be plotted in moles C per liter of culture
		var resolution = 50;				// #			The number of points to generate on the plot
					
		//numerical constants
		var kd = 1.72e-6;						// s^-1 				FeEDTA dark dissociation rate constant at 20C (Sunda2003)
		var khv = 4.3e-6;						// s^-1 				FeEDTA photooxydation rate constant at 20C, 500 uE m^-2 s^-1 (Sunda2003)
		var kdprime= kd + (khv * (light / 500) * (hlight / 24));  // s^-1					The FeEDTA dissociation rate
		var kf=17;  							// L mol^-1 s^-1		The FeEDTA formation rate
		var Cmolpervol=15;						// mol L of cell vol 	An estimate of moles Carbon per cel volume values from Sunda 2005
		var Vmax=1276; 							// nmol m^-2 d^-1 		from sunda and huntsman 1997
		var Krho=5.1E-10;						// mol * L^-1
					
		//calculated model values
		var volcell = (4/3) * Math.PI *Math.pow(r,3) * 1000;	// L cell^-1		volume of an individual cell
		var biomassmolest =  volcell * Cmolpervol; 				// mol cell^-1 		Moles of carbon per individual cell
		var rhomax =  	4 * Math.PI *Math.pow(r,2) * Vmax/(86400*1E+9);		// vmax calculation L* s^-1 
			
		//Calculated values for plot d1 and d2
		var maxvect = maxbiomass/biomassmolest;
		var stepsize = maxvect/resolution;
		var lnmaxvect = Math.log(maxvect);
		var lnstepsize = lnmaxvect/resolution;

		//calculate buffer failure
		var failpct =  0.1; //change from original that is considered buffer failure, e.g. 0.1 is 10%
		var abiotic = (kdprime*feedta)/(kf*edta);
		var fefail	= ((1-failpct)*kdprime*feedta)/(kf*edta); // 90% of the abiotic Fe concentration
		var failN   =  -(edta*Math.pow(fefail,2)*kf - fefail*feedta*kdprime + (edta*fefail*kf -feedta*kdprime)*Krho)/(fefail*rhomax);
		var failNmL	= failN/1000;
		var failC	= failN*biomassmolest*1e6;

		//Creation of data array for Plot d1, [Fe'] vs cell biomass
		var d1 = [];
		for(var i=0;i< resolution ;i++){
			d1.push([biomassmolest*1e6*i*stepsize, -(1e12/2)*(Krho*edta*kf - feedta*kdprime + i*stepsize*rhomax - Math.sqrt(Math.pow(Krho,2)*Math.pow(edta,2)*Math.pow(kf,2) + 2*Krho*edta*feedta*kdprime*kf + Math.pow(feedta,2)*Math.pow(kdprime,2) + Math.pow(i*stepsize,2)*Math.pow(rhomax,2) + 2*(Krho*edta*kf -feedta*kdprime)*i*stepsize*rhomax))/(edta*kf)])
		};
		//Creation of data array for Plot d2, [Fe'] vs cell density
		var d2 = [];
		for(var i=0;i< resolution ;i++){
			d2.push([ Math.exp(i*lnstepsize)/1e3,(-1e12/2)*(Krho*edta*kf - feedta*kdprime + Math.exp(i*lnstepsize)*rhomax - Math.sqrt(Math.pow(Krho,2)*Math.pow(edta,2)*Math.pow(kf,2) + 2*Krho*edta*feedta*kdprime*kf+ Math.pow(feedta,2)*Math.pow(kdprime,2) + Math.pow(Math.exp(i*lnstepsize),2)*Math.pow(rhomax,2)+ 2*(Krho*edta*kf -feedta*kdprime)*Math.exp(i*lnstepsize)*rhomax))/(edta*kf)
])
		};

			
			
		//warn if Fe' exceeds 400 nM, the solubility limit measuerd by Liu and Millero (2002)	
		if(abiotic > 4e-10){ $(function() {
        $( "#fewarning" ).dialog({
            height: 140,
            modal: true
        });
    });};
		//output the estimates of key values
		$('#feconc').text(abiotic.toPrecision(2));
		$('#density').text(failNmL.toPrecision(2));
		$('#biomass').text(failC.toPrecision(2));
		
		
		//plotting d1 and d2
		$.plot($("#plot1"),[d1], {
			hoverable:true,
			yaxis: {
				transform: function (v) {
					return Math.log(v)
					}, 
				ticks: [1,3, 10,30, 100,300, 1000]
				},	
			grid: {
			markings:[{xaxis: {from:  1, to:failC }, color: "#f8f8f8"},{xaxis: {from:  failC, to:250 }, color: "#999999"}  ]}});
		
		
		
		$.plot($("#plot2"),[d2], {
			yaxis:{min:0},
			xaxis: {
			transform: function (v) {
				return Math.log(v)},
				 min: 1,ticks: [[1e0,"1<sup> </sup>"],[1e1,"10<sup> </sup>"],[1e2,"10<sup>2</sup>"],[1e3,"10<sup>3</sup>"],[1e4,"10<sup>4</sup>"],[1e5,"10<sup>5</sup>"],[1e6,"10<sup>6</sup>"],[1e7,"10<sup>7</sup>"]]},
			grid: {
				markings:[{xaxis: {from: 0, to:failN/1000 }, color: "#F8F8F8"},{xaxis: {from: failN/1000, to:1e9 }, color: "#999999"}]}});
				;
	}

	