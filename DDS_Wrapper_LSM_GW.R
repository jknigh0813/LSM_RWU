
#Dynamically Dimensioned Search (DDS) hooked to JoFlo
#Created: James Knighton
#Date: 10/21/2019
#Function executing the Dynamically Dimensioned Search (DDS) calibration routine for JoFlo


DDS_Calibrate_LSM <- function(gage, numIter)  
{
  
  source("//nfs/jknighton-data/Lumped_VSA_model_Mod.R")
  
  # Returns numIter length list of entries to be peturbed
  probPeturb<-function(x, numIter){
    # Input is xBounds & numIter.  
    # Returns numIter entry list with the indices which will be peturbed
    xDims<-nrow(x)
    probabilityVector<-1-log(1:numIter)/log(numIter)
    peturbIdx<-apply(matrix(unlist(lapply(probabilityVector, function(x) as.logical(rbinom(xDims, 1, x)))), byrow=TRUE, ncol=xDims), 1, which)
    return(peturbIdx)
  }
  
  #Read in model settings file
  infile = "//nfs/jknighton-data/HCDN_Sites_Plasticity_Estimate.csv"
  USGS = read.csv(infile)
  USGS = USGS[USGS$Skip == 0,]
  
  #Read in met data      
  # MetData <- read.csv(paste("//nfs/jknighton-data/MetForcing_CPC/",USGS$STAID[gage],"_MetData.csv",sep=""),head=FALSE)
  MetData <- read.csv(paste("//nfs/jknighton-data/MetForcing_Livneh/",USGS$STAID[gage],"_MetData.csv",sep=""),head=FALSE)
  colnames(MetData) <- c("Precip_mm","Tmax_C","Tmin_C")
  MetData$Date = seq(as.Date("1999/01/01"), as.Date("2010/12/31"),"days")
  
  #Calculate average daily air temperature
  MetData$Tavg_C = (MetData$Tmax_C + MetData$Tmin_C)/2
  
  #Set catchment latitude
  latitudeDegrees = USGS$LAT_GAGE[gage] #decimal degrees
  latitudeRadians<-latitudeDegrees*pi/180 ## latitude in radians
  
  #Get USGS daily discharge data (2000 - 2010)  
  if (nchar(USGS$STAID[gage]) == 7){Obs = get_usgs_gage(paste("0",USGS$STAID[gage], sep=""), begin_date = "2000-01-01", end_date="2010-12-31")}
  if (nchar(USGS$STAID[gage]) == 8){Obs = get_usgs_gage(USGS$STAID[gage], begin_date = "2000-01-01", end_date="2010-12-31")}
  Obs$flowdata$flow = Obs$flowdata$flow/(1000*Obs$area)
  
  #Set calibration parameter boundaries
  xBounds.df = data.frame(matrix(ncol=2,nrow=10))
  colnames(xBounds.df)<-c("min", "max")
  
  #DDS parameters
  r= 0.2

  #rec_coef
  xBounds.df$min[1] = 0.000001
  xBounds.df$max[1] = 0.01
  
  #Bexp
  xBounds.df$min[2] = 1.5
  xBounds.df$max[2] = 2.5
  
  #Se_min
  xBounds.df$min[3] = 1
  xBounds.df$max[3] = 200
  
  #C1
  xBounds.df$min[4] = 0.01
  xBounds.df$max[4] = 15
  
  #Ia_Coeff
  xBounds.df$min[5] = 0.01
  xBounds.df$max[5] = 0.2
  
  #percentimpervios
  xBounds.df$min[6] = 0
  xBounds.df$max[6] = 5
  
  #light ext.
  xBounds.df$min[7] = 0
  xBounds.df$max[7] = 1
  
  #Bw
  xBounds.df$min[8] = -10
  xBounds.df$max[8] = 10
  
  #EPCO
  xBounds.df$min[9] = 0
  xBounds.df$max[9] = 1
  
  #Root Depth
  xBounds.df$min[10] = 100
  xBounds.df$max[10] = USGS$Thickness[gage]
  
  # Generate initial first guess
  x_init<-c(0.01, 2, 80, 0.02, 0.049, 3, 0.01, 0.5, 0.01, 300)
  x_best = data.frame(x_init)
  NSE_init = -9999
  MAE_best <- 10000
  NSE_best <-NSE_init
  Pbias_best <-NSE_init
  
  peturbIdx<-probPeturb(xBounds.df, numIter)
  # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if beyond boundaries
  sigma<-xBounds.df$max - xBounds.df$min
  
  for (i in 1:numIter){
    # Set up test parameter values as x_test
    x_test<-as.matrix(x_best)
    
    # Get entries we will peturb
    idx<-peturbIdx[[i]]
    if (sum(idx) == 0) {idx = round(runif(1,1,5))}
    
    # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
    peturbVec<-rep(0, length(x_test))
    # Generate the required number of random normal variables
    N<-rnorm(length(x_test), mean=0, sd=1)
    
    # Set up vector of peturbations
    peturbVec[idx]<-r*N[idx]*sigma[idx]
    
    # Temporary resulting x value if we peturbed it
    testPeturb<-x_test + peturbVec  
    # Find the values in testPeturb that have boundary violations.  Store the indices in boundaryViolationsIdx
    boundaryViolationIdx<-which(testPeturb<xBounds.df$min | testPeturb > xBounds.df$max)
    
    # Reset those violated indices to the opposite peturbation direction
    peturbVec[boundaryViolationIdx]<-(-1*r*N[boundaryViolationIdx]*sigma[boundaryViolationIdx])
    
    # Find values still at violations of min or max and set them to the minimum or maximum values
    x_test<-x_test + peturbVec
    minViolationIdx<-which(x_test<xBounds.df$min)
    maxViolationIdx<-which(x_test>xBounds.df$max)
    x_test[minViolationIdx]<-xBounds.df$min[minViolationIdx]
    x_test[maxViolationIdx]<-xBounds.df$max[maxViolationIdx]
    
    #Walter et al (2005) snow dynamics model
    snowmelt=SnowMelt(Date=MetData$Date, precip_mm=MetData$Precip_mm, Tmax_C=MetData$Tmax_C, Tmin_C=MetData$Tmin_C, lat_deg=latitudeDegrees, forest=x_test[7])
    
    #Precipitation = rainfall + snowmelt    
    MetData$precip.tot=snowmelt$SnowMelt_mm+snowmelt$Rain_mm
    
    #Archibald et al (2014), Archibald & Walter (2014), Knighton & Singh (2019)    
    Results <- Lumped_VSA_model_Mod(dateSeries = MetData$Date, P = snowmelt$SnowMelt_mm+snowmelt$Rain_mm, Tmax=MetData$Tmax_C, Tmin = MetData$Tmin_C,
                                    latitudeDegrees=latitudeDegrees, Depth = USGS$Thickness[gage], SATper = 0.5, AWCper = USGS$AWC[gage], 
                                    Tp = 5, albedo = 0.23, StartCond = "avg", BF1 = 1, PETcap = 100, rec_coef = x_test[1], Bexp = x_test[2], Se_min = x_test[3], 
                                    C1 = x_test[4], Ia_coef = x_test[5],percentImpervious=x_test[6],Bw = x_test[8], EPCO = x_test[9], RootDepth=x_test[10]) 
    
    #Merge simulated and observed discharge records
    Results$mdate = Results$Date
    AllData = merge(Results, Obs$flowdata, by="mdate")
    AllData$month =  as.numeric(format(AllData$mdate, "%m"))
    AllData$year =  as.numeric(format(AllData$mdate, "%Y"))
    
    AllData$modeled_flow[AllData$modeled_flow == 0] = 0.001
    AllData$flow[AllData$flow == 0] = 0.001
    
    #Compute objective functions
    #Daily NSE
    numer = 0
    denom = 0
    for (j in 1:nrow(AllData))
    {numer = numer + (log(AllData$modeled_flow[j]) - log(AllData$flow[j]))^2
    denom = denom + (log(AllData$flow[j]) - mean(log(AllData$flow)))^2}
    NSE_daily = 1 - numer/denom
    
    #Pbias
    Pbias = 100*(sum(AllData$modeled_flow) - sum(AllData$flow))/sum(AllData$flow)
    
    #MAE
    numer = 0
    denom = 0
    for (j in 1:nrow(AllData))
    {numer[j] = abs(log(AllData$modeled_flow[j]) - log(AllData$flow[j]))}
    MAE = mean(numer)
    
    #Check if this simulation is better
    if (MAE < MAE_best)
    {
      x_best = x_test
      MAE_best = MAE   
      NSE_best = NSE_daily
      Pbias_best = Pbias
      plotdata = AllData
    }
    print_str = paste("Eval:",i,"   MAE:",round(MAE,digits=5),"   MAE Best:",round(MAE_best,digits=5),"   NSE:", round(NSE_best,digits=4))
    print(print_str)
    
  }
  
  #Save results
  Output = matrix(nrow=1,ncol=14)
  Output[1] = USGS$STAID[gage]
  Output[2] = MAE_best
  Output[3] = NSE_best
  Output[4] = Pbias_best
  Output[5] = x_best[1]
  Output[6] = x_best[2]
  Output[7] = x_best[3]
  Output[8] = x_best[4]
  Output[9] = x_best[5]
  Output[10] = x_best[6]
  Output[11] = x_best[7]
  Output[12] = x_best[8]
  Output[13] = x_best[9]
  Output[14] = x_best[10]

  write.csv(Output,paste("//nfs/jknighton-data/Model_Val_Livneh_Cluster/ModelValidationResults_Livneh_",USGS$STAID[gage],".csv"),sep="")
  
  #return(Output)

}

