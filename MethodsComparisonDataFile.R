##########################
## Read in RPresence data files
## Save as RData file
##########################



###load data files
#all data combined
dta_all <- readPao('allmethods.pao')
#all camera trap data
dta_GC <- readPao('dta_CT.pao')
#canopy camera only data
dta_C <- readPao('dta_canopy.pao')
#ground camera only data
dta_G <- readPao('dta_ground.pao')
#line transect data
dta_L <- readPao('dta_LT.pao')
#line transect and canopy camera
dta_LC <- readPao('dta_LTCC.pao')
#line transect and ground camera 
dta_LG <- readPao('dta_LTGC.pao')
#ind canopy cameras
dta_ind_canopy <- readPao('dta_ind_canopy.pao')


save.image(file = "MethodsComparisonPublicationData.RData")
