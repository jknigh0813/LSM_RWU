
#Run DDS -> JoFlo on SESYNC cluster
#Created: James Knighton
#Date:10/21/2019
#Sends the Dynamically Dimensioned Search calibration routine of JoFlo to the SESYNC cluster

library(EcoHydRology)
library(zoo)
library(rslurm)

source("/research-home/jknighton/LSM/DDS_Calibrate_LSM.R")
source("//nfs/jknighton-data/Lumped_VSA_model_Mod.R")

#Read in index
infile = "//nfs/jknighton-data/HCDN_Sites.csv"
USGS = read.csv(infile)
USGS = USGS[USGS$Skip == 0,]

#test = DDS_Calibrate_LSM(1,numIter=10)

pars = data.frame(gage = 1:264, numIter = 10000)

sjob <- slurm_apply(DDS_Calibrate_LSM, pars, jobname = 'test_apply',
                    nodes = 3, cpus_per_node = 8, submit = TRUE)