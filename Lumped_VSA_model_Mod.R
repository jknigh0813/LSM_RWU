
#Modified JoFlo
#Created: James Knighton
#Date: 10/21/2019
#JoFlo (Archibald et al. 2014) - modified to include:
#a. five vertically connected soil layers
#b. depth distribution of RWU (Bw, and RootDepth)
#c. RWU plasticity (EPCO)
#d. exponential parameter for gravity driven baseflow

Lumped_VSA_model_Mod <-
function(
	dateSeries, # Date series (format should be: 2003-02-36)
	P, # Rain + Snow melt (mm)
	Tmax, # Max daily temperature (C)
	Tmin, # Max daily temperature (C)
	Depth = 1000, # Average soil depth in watershed (mm) [don't need if AWC and SAT entered directly]
	SATper = NULL, # Porosity (fraction)
	AWCper = NULL, # Available Water Capacity (AWC) (fraction)
	percentImpervious = 0, # Percent of the watershed that is impervious
	no_wet_class = 10, # The number of wetness classes for saturated area designation
	Tp = 5, # Time to peak (hours)
	latitudeDegrees = 42.38, 
	albedo = 0.23, # Average albedo
	StartCond = avg, # Watershed conditions before first day of run (wet, dry, avg)
	PETin = NULL,# User has the option to enter PET values (mm/day)
	AWC = Depth*AWCper, # AWC as a depth (mm)
	SAT = Depth*SATper,# porosity as a depth (mm)
	SW1 = NULL, # Soil water on the first day (depth, mm)
	BF1 = 1, #  mm/day can use nearby watershed baseflow, ConvertFlowUnits(cfs,WA=W_Area)
	PETcap = 5, # Does not let PET get larger than this cut off (mm)
	rec_coef = 0.1, # based on a study by Weiler in NY state
  Bexp=2, #exponent
	Se_min = 78, # mm
	C1 = 3.1, # Coefficient relating soil water to Curve Number S
	Ia_coef = 0.05, # range ~ 0.05 - 0.2
	PreviousOutput = NULL, # Allows us to take previous model output to initiate the model
	runoff_breakdown = RunoffBreakdown(Tp, HrPrcDelay = (Tp/2-4)) ,
	Bw = 10, #plant rooting distribution shape parameter	
	EPCO = 1, #RWU plasticity index
	RootDepth = 1 #Maximum rooting depth (mm)
){
	#Soil Layer Available Water Capacity
	#Top 0.4 m spread between four upper layers
	#Bottom layer represents the remaining depth
	AWC_a = 100*AWCper
	AWC_b = 100*AWCper
	AWC_c = 100*AWCper
	AWC_d = 100*AWCper
	AWC_e = (Depth - 400)*AWCper 
	
	latitude<-latitudeDegrees*pi/180 ## latitude in radians

	### This is to allow us to use previous model output as an input - useful for calculating just
	### a few more days of modeled flow without having to run years of data ..
	if (!is.null(PreviousOutput)){# We will extend the input variables to include overlap with previous model input
		nr<- nrow(PreviousOutput)
		dateSeries<- c(PreviousOutput$Date[(nr-3):nr],dateSeries)
		P<- c(PreviousOutput$rain_snowmelt[(nr-3):nr], P)
		Tmax<- c(PreviousOutput$Tmax[(nr-3):nr], Tmax)
		Tmin<- c(PreviousOutput$Tmin[(nr-3):nr], Tmin)
		OutStart<- 5# The output that we report will not include previous output
	} else OutStart<- 1
	
	################################################
	Tav<-(Tmax+Tmin)/2
	day<-strptime(dateSeries,format = "%Y-%m-%d")$yday+1
	month<-strptime(dateSeries,format = "%Y-%m-%d")$mon+1

	## Potential Evapotranspiration 
	PET<-PET_fromTemp(Jday=day,Tmax_C=Tmax,Tmin_C=Tmin,AvgT=Tav,albedo=albedo,lat_radians=latitude)*1000## mm (Priestley-Taylor)
	PET[which(PET>PETcap)]<-PETcap#Sets a cap on PET estimates (usually ~ 5mm)
	ETo <- PET
	ET <- ETo    #  Initializing modeled actual ET, equal to ETo when deltaP > 0, but less than ETo on drier days
	
	#Initial Root Water Uptake (RWU) demand on each layer
	wz_a = (PET/(1 - exp(-1*Bw)))*(1 - exp(-1*Bw*100/RootDepth))
	wz_a[wz_a < 0] = 0
	if (RootDepth >= 200) {wz_b = (PET/(1 - exp(-1*Bw)))*(1 - exp(-1*Bw*200/RootDepth)) - wz_a}
	else {wz_b = rep(0,length(wz_a))}
	if (RootDepth >= 300) {wz_c = (PET/(1 - exp(-1*Bw)))*(1 - exp(-1*Bw*300/RootDepth)) - (wz_a + wz_b)}
	else {wz_c = rep(0,length(wz_a))}
	if (RootDepth >= 400) {wz_d = (PET/(1 - exp(-1*Bw)))*(1 - exp(-1*Bw*400/RootDepth)) - (wz_a + wz_b + wz_c)}
	else {wz_d = rep(0,length(wz_a))}
	if (RootDepth > 400) {wz_e = PET - (wz_a + wz_b + wz_c + wz_d)}
	else {wz_e = rep(0,length(wz_a))}
	wz_GW = a <- rep(0, length(wz_e))
	
	wz_a_prime = wz_a
	wz_b_prime = wz_b
	wz_c_prime = wz_c
	wz_d_prime = wz_d
	wz_e_prime = wz_e
	ET_a = wz_a	
	ET_b = wz_b	
	ET_c = wz_c	
	ET_d = wz_d	
	ET_e = wz_e	
	ET_GW = wz_GW
	
	#Initializing vectors for Water Budget Loop, and Runoff estimation
	SoilWater<-vector(length=length(P))##(mm)
	SoilWater_a<-vector(length=length(P))##(mm)
	excess_a<-vector(length=length(P))
	SoilWater_b<-vector(length=length(P))##(mm)
	excess_b<-vector(length=length(P))
	SoilWater_c<-vector(length=length(P))##(mm)
	excess_c<-vector(length=length(P))
	SoilWater_d<-vector(length=length(P))##(mm)
	excess_d<-vector(length=length(P))
	SoilWater_e<-vector(length=length(P))##(mm)
	excess_e<-vector(length=length(P))
	
	TM_S<-vector(length=length(P))##This is the daily time-step T-M storage, used to calc baseflow
	totQ<-vector(length=length(P))
	Se<-vector(length=length(P))
	Q2<-vector(length=length(P))
	baseflow<-vector(length=length(P))
	MaxWetClass<-vector(length=length(P))
	impervRunoff<-vector(length=length(P))
	OverlandFlow<-vector(length=length(totQ))
	ShallowInterflow<-vector(length=length(totQ))

	deltaP <- P - wz_a ## (mm) neg values indicate net evap
	impervIa<-Ia_coef * 5  ##Se = 5 mm for impervious areas (CN=98)
	impervPe<-deltaP - impervIa##  Effective precipitation on impervious areas
	impervPe[which(impervPe < 0)]<-0
	Pe<-vector(length=length(P))# Effective Precipitation (non-impervious areas)

	# Setting up initial values
	TM_S[1]<-BF1/rec_coef
	SoilWater[1] <- AWC
	SoilWater_a[1] <- AWC_a
	SoilWater_b[1] <- AWC_b
	SoilWater_c[1] <- AWC_c
	SoilWater_d[1] <- AWC_d
	SoilWater_e[1] <- AWC_e
	excess_a[1] <- 0# Assume no runoff from the day before
	excess_b[1] <- 0# Assume no runoff from the day before
	excess_c[1] <- 0# Assume no runoff from the day before
	excess_d[1] <- 0# Assume no runoff from the day before
	excess_e[1] <- 0# Assume no runoff from the day before
	Se[1]<- Se_min+C1*(AWC-SoilWater_a[1])
	totQ[1] <- Pe[1]*Pe[1]/(Pe[1]+Se[1])
	impervRunoff[1] <- impervPe[1]
	
	Ia<-Ia_coef*Se  # Initial abstraction = depth of precip before runoff starts

	## Thornthwaite-Mather Function and Runoff Generation 
	for(i in (2):length(deltaP)){  
		
		
		#2. Infiltration and runoff
		#2.1: Runoff
		if (deltaP[i]-Ia[i-1]>0){Pe[i] <- deltaP[i]-Ia[i-1]  
		#2.2: no runoff
		} else Pe[i]<-0
		totQ[i] <- Pe[i]*Pe[i]/(Pe[i]+Se[i-1]) ## Effective storage is from previous day
		
		####################################################################################
		#3. Soil Water Balance Layer A
		#3.1: Drying
		if(deltaP[i]<=0){  #DRYING CONDITION
		SoilWater_a[i]<- SoilWater_a[i-1]*exp(deltaP[i]/AWC_a)
		ET_a[i] <- max(0, SoilWater_a[i-1] - SoilWater_a[i])} 
		
		#3.2: Wetting below Theta_fc
		else if (deltaP[i] + SoilWater_a[i-1] <= AWC_a){
		SoilWater_a[i]<- deltaP[i] + SoilWater_a[i-1]}
		
		#3.3: Wetting above Theta_fc
		else { 
		SoilWater_a[i] <- AWC_a  ## So overall SW cannot exceed AWC
		excess_a[i] <- deltaP[i] + SoilWater_a[i-1] - AWC_a}
		
		####################################################################################
		#4. Soil Water Balance Layer B
		#4.1: Drying
		if (RootDepth >= 200) {wz_b_prime[i] = wz_b[i] + EPCO*max(0,wz_a[i] - ET_a[i])}
		else {wz_b_prime[i] = 0}
		
		if((excess_a[i] - wz_b_prime[i])<=0){  #DRYING CONDITION
		SoilWater_b[i]<- SoilWater_b[i-1]*exp((excess_a[i] - wz_b_prime[i])/AWC_b)
		ET_b[i] <- max(0, SoilWater_b[i-1] - SoilWater_b[i])}   ## amount that gets evaporated (mm)
		
		#4.2: Wetting below Theta_fc
		else if ((excess_a[i] - wz_b_prime[i]) + SoilWater_b[i-1] <= AWC_b){
		SoilWater_b[i]<- (excess_a[i] - wz_b_prime[i]) + SoilWater_b[i-1]
		
		#4.3: Wetting above Theta_fc
		} else { 
		SoilWater_b[i] <- AWC_b  ## So overall SW cannot exceed AWC
		excess_b[i] <- deltaP[i] + SoilWater_b[i-1] - AWC_b}
		
		####################################################################################
		#4. Soil Water Balance Layer C
		if (RootDepth >= 300) {wz_c_prime[i] = wz_c[i] + EPCO*max(0,(wz_a[i] + wz_b[i]) - (ET_a[i] + ET_b[i]))}
		else {wz_c_prime[i] = 0}

		#4.1: Drying
		if((excess_b[i] - wz_c_prime[i])<=0){  #DRYING CONDITION
		SoilWater_c[i]<- SoilWater_c[i-1]*exp((excess_b[i] - wz_c_prime[i])/AWC_c)
		ET_c[i] <- max(0, SoilWater_c[i-1] - SoilWater_c[i])}  ## amount that gets evaporated (mm)
		
		#4.2: Wetting below Theta_fc
		else if ((excess_b[i] - wz_c_prime[i]) + SoilWater_c[i-1] <= AWC_c){
		SoilWater_c[i]<- (excess_b[i] - wz_c_prime[i]) + SoilWater_c[i-1]
		
		#4.3: Wetting above Theta_fc
		} else { 
		SoilWater_c[i] <- AWC_c  ## So overall SW cannot exceed AWC
		excess_c[i] <- deltaP[i] + SoilWater_c[i-1] - AWC_c}		
		
		####################################################################################
		#4. Soil Water Balance Layer D
		if (RootDepth >= 400) {wz_d_prime[i] = wz_d[i] + EPCO*max(0,(wz_a[i] + wz_b[i] + wz_c[i]) - (ET_a[i] + ET_b[i] + ET_c[i]))}
		else {wz_d_prime[i] = 0}
		
		#4.1: Drying
		if((excess_c[i] - wz_d_prime[i])<=0){  #DRYING CONDITION
		SoilWater_d[i]<- SoilWater_d[i-1]*exp((excess_c[i] - wz_d_prime[i])/AWC_d)
		ET_d[i] <- max(0, SoilWater_d[i-1] - SoilWater_d[i])  ## amount that gets evaporated (mm)
		
		#4.2: Wetting below Theta_fc
		} else if ((excess_c[i] - wz_d_prime[i]) + SoilWater_d[i-1] <= AWC_d){
		SoilWater_d[i]<- (excess_c[i] - wz_d_prime[i]) + SoilWater_d[i-1]
		
		#4.3: Wetting above Theta_fc
		} else { 
		SoilWater_d[i] <- AWC_d  ## So overall SW cannot exceed AWC
		excess_d[i] <- deltaP[i] + SoilWater_d[i-1] - AWC_d}	

		####################################################################################
		#4. Soil Water Balance Layer E
		if (RootDepth > 400) {wz_e_prime[i] = wz_e[i] + EPCO*max(0,(wz_a[i] + wz_b[i] + wz_c[i] + wz_d[i]) - (ET_a[i] + ET_b[i] + ET_c[i] + ET_d[i]))}
		else {wz_e_prime[i] = 0}

		#4.1: Drying
		if((excess_d[i] - wz_e_prime[i])<=0){  #DRYING CONDITION
		SoilWater_e[i]<- (SoilWater_e[i-1])*exp((excess_d[i] - wz_e_prime[i])/(AWC_e))
		ET_e[i] <- max(0, SoilWater_e[i-1] - SoilWater_e[i])  ## amount that gets evaporated (mm)
		
		#4.2: Wetting below Theta_fc
		} else if ((excess_d[i] - wz_e_prime[i]) + SoilWater_e[i-1] <= AWC_e){
		SoilWater_e[i]<- (excess_d[i] - wz_e_prime[i]) + SoilWater_e[i-1]
		
		#4.3: Wetting above Theta_fc
		} else { 
		SoilWater_e[i] <- AWC_e  ## So overall SW cannot exceed AWC
		excess_e[i] <- deltaP[i] + SoilWater_e[i-1] - AWC_e}	

		####################################################################################
		#4. Groundwater ET
		wz_GW[i] = EPCO*max(0,(wz_a[i] + wz_b[i] + wz_c[i] + wz_d[i] + wz_e[i]) - (ET_a[i] + ET_b[i] + ET_c[i] + ET_d[i] + ET_d[i] + ET_e[i]))
		ET_GW[i] = max(0,min(TM_S[i-1] - max(0.1,Depth - RootDepth),wz_GW[i]))
		TM_S[i] = TM_S[i-1] - ET_GW[i]
		
		#Total Soil AET
		ET[i] = ET_a[i] + ET_b[i] + ET_c[i] + ET_d[i] + ET_e[i] + ET_GW[i]

		SoilWater[i] = SoilWater_a[i] + SoilWater_b[i] + SoilWater_c[i] + SoilWater_d[i] + SoilWater_e[i]
		
		#Initial abstraction conditioned on catchment wetness
		Se[i] <- Se_min + C1*(AWC-SoilWater[i])
		Ia[i]<-Ia_coef*Se[i]

		#Gravity driven baseflow
		baseflow[i]<-rec_coef*(TM_S[i]^Bexp)
		if ((excess_e[i])>=totQ[i]){# Ensure mass-balance, since curve number is empirical
		TM_S[i]<-max(0,TM_S[i] - baseflow[i] + excess_e[i] - totQ[i])  
		} else {
		TM_S[i]<-max(0,TM_S[i] - baseflow[i])
		}

		impervRunoff[i]<-impervPe[i] 
	}
	### End of Water Balance Loop ############################################################################

  #Not currently used (10/21/19 JK)
	# Use the coefficients generated from time of concentration
	runoff_breakdown[which(runoff_breakdown < 0.01)] <- 0
	if (OutStart==1){  # Initialize these values if not already taken from previous output
	OverlandFlow[1:4]<-totQ[1:4]*runoff_breakdown[1]
	ShallowInterflow[1]<-0# Assuming no runoff in previous days before this model run
	ShallowInterflow[2]<-totQ[1]*runoff_breakdown[2]
	ShallowInterflow[3]<-totQ[1]*runoff_breakdown[3] + totQ[2]*runoff_breakdown[2]
	ShallowInterflow[4]<-totQ[1]*runoff_breakdown[4] + totQ[2]*runoff_breakdown[3] + totQ[3]*runoff_breakdown[2]
	}

	#Write output files
	modeled_flow<-(totQ)*(1-0.01*percentImpervious) + impervRunoff*(0.01*percentImpervious) + baseflow#totQ+baseflow
	rain_snowmelt<-P
	Streamflow<-data.frame(Date=as.Date(as.character(dateSeries)),  Tmax, Tmin, rain_snowmelt, modeled_flow, baseflow, deltaP, Pe, Se, Ia, SoilWater_a, SoilWater_b, SoilWater_c, SoilWater_d, SoilWater_e, PET, wz_a_prime, wz_b_prime, wz_c_prime, wz_d_prime, wz_e_prime, wz_GW, ET, ET_a, ET_b, ET_c, ET_d, ET_e, ET_GW, TM_S, totQ)[OutStart:length(P),] 
	return(Streamflow)
}
