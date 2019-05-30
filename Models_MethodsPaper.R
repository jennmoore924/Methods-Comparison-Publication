
rm(list=ls())

#load data files
load(file = "dta_all.RData") #all data combined by transect
load(file = "dta_LT.RData") #line transect data combined by transect
load(file = "dta_CT.RData") #all cameras combined by transect
load(file = "dta_CT_ground.RData") #ground cameras combined by transect
load(file = "dta_CT_canopy.RData") #canopy cameras combined by transect
load(file = "dta_LTCC.RData") #line transects and canopy cameras combined by transect
load(file = "dta_LTGC.RData") #line transects and ground cameras combined by transect
load(file = "dta_CT_ind_canopy.RData")  #individual canopy camera data

##species accumulation curves

#only load dta_LT and dta_CT files
library(reshape2)
library(vegan)
d <- rbind(lt5, ct4)
d2 <- aggregate(Count ~ Transect + Species, data = d, FUN = "sum")
d3 <- dcast(d2, Transect ~ Species, value.var = "Count")
d4 <- d3[ ,2:37]
alldata <- specaccum(d4, method = "random")

#dta_LT file
d2 <- aggregate(Count ~ Transect + Species, data = lt5, FUN = "sum")
d3 <- dcast(d2, Transect ~ Species, value.var = "Count")
d4 <- d3[ ,2:37]
linedata <- specaccum(d4, method = "random")

#dta_CT_ground file
d2 <- aggregate(Count ~ Transect + Species, data = ct5, FUN = "sum")
d3 <- dcast(d2, Transect ~ Species, value.var = "Count")
d4 <- d3[ ,2:37]
grounddata <- specaccum(d4, method = "random")

#dta_CT_canopy file
d2 <- aggregate(Count ~ Transect + Species, data = ct5, FUN = "sum")
d3 <- dcast(d2, Transect ~ Species, value.var = "Count")
d4 <- d3[ ,2:37]
canopydata <- specaccum(d4, method = "random")


#par(mfrow = c(2,2))
#plot(alldata, main = "All Data")
#plot(linedata, main = "Line Transect Data")
#plot(grounddata, main = "Ground Camera Data")
#plot(canopydata, main = "Canopy Camera Data")

par(mar = c(6, 6, 4, 2))
plot(alldata, xlab = "Number of Sites", ylab = "Number of Species", las = 1, lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, xlim = c(1,18), ylim = c(0,40))
plot(linedata, add = TRUE, col = "red", lty = 2, lwd = 2)
plot(grounddata, add = TRUE, col = "blue", lty = 3, lwd = 2)
plot(canopydata, add = TRUE, col = "limegreen", lty = 4, lwd = 2)
legend("topleft", bty = 'n', legend = c("all data", "ground cameras", "arboreal cameras", "line transects"), lty = c(1,3,4,2), col = c("black", "blue", "limegreen", "red"), cex=1.1)

######################
## Occupancy Models ##
######################
library(RMark)


#detection line transects
dtaLT <- dta_LT
indicesLT <- as.numeric(row.names(dtaLT[dtaLT$ch == "..",]))
dtaLT <- dtaLT[-indicesLT,]

#process data
pdLT=process.data(dtaLT,model="Occupancy", groups = c('transect2','obs2'))
#design data
ddLT=make.design.data(pdLT)
#initial model
mLT=mark(pdLT,ddLT,model="Occupancy",model.parameters=list(Psi=list(formula=~1),p=list(formula=~spef+tief)))
#detection 0.311
#occupancy 0.2218
#joint model detection 0.1404
mLTloc=mark(pdLT,ddLT,model="Occupancy",model.parameters=list(Psi=list(formula=~transect2),p=list(formula=~spef+tief)))


#canopy cameras
dtaCCT <- dta_CT_canopy
#process data
pdCCT=process.data(dtaCCT,model="Occupancy",groups = c('transect2','obs2'))
#design data
ddCCT=make.design.data(pdCCT)
#initial model
mCCT=mark(pdCCT,ddCCT,model="Occupancy",model.parameters=list(Psi=list(formula=~1),p=list(formula=~spef+tief)))
#detection 0.728
#occupancy 0.352
#joint model detection 0.7979
mCCTloc=mark(pdCCT,ddCCT,model="Occupancy",model.parameters=list(Psi=list(formula=~transect2),p=list(formula=~spef+tief)))


#ground cameras
dtaGCT <- dta_CT_ground
#process data
pdGCT=process.data(dtaGCT,model="Occupancy", groups = c('transect2','obs2'))
#design data
ddGCT=make.design.data(pdGCT)
#initial model
mGCT=mark(pdGCT,ddGCT,model="Occupancy",model.parameters=list(Psi=list(formula=~1),p=list(formula=~spef+tief)))
#detection 0.688
#occupancy 0.295
#joint model detection 0.7728
mGCTloc=mark(pdGCT,ddGCT,model="Occupancy",model.parameters=list(Psi=list(formula=~transect2),p=list(formula=~spef+tief)))


##detection 1
## overall detection for any species any method
dta1 <- rbind(dta_LT, dta_CT_canopy, dta_CT_ground)
dta1$obs2 <- as.factor(dta1$obs2)

#subset for method of interest
dta1 <- subset(dta1, dta1$method == "transect")
rownames(dta1) <- NULL

#scale effort covariates
dta1$stief1 <- scale(dta1$tief1)
dta1$stief2 <- scale(dta1$tief2)
dta1$sspef <- scale(dta1$spef)

#process data
#pd1=process.data(dta1,model="Occupancy", groups = c('transect2','obs2','method'))
pd1=process.data(dta1,model="Occupancy", groups = c('transect2','obs2'))
#design data
dd1=make.design.data(pd1)

#fix detection to 0 for animals that are not available for detection by each method
indices <- as.numeric(row.names(dta1[dta1$ch == ".." & dta1$transect2 == "1",]))
xfix <- dta1[indices,]

dd1$p$xfix <- 1
dd1$Psi$xfix <- 1

 for(i in 1:nrow(xfix)){   
   j=which(dd1$p$obs==xfix$obs2[i] & dd1$p$method==xfix$method[i]);  
   if (length(j)>0) dd1$p$xfix[j]=0
 }

for(i in 1:nrow(xfix)){   
  j=which(dd1$p$obs==xfix$obs2[i]);  
  if (length(j)>0) dd1$p$xfix[j]=0
}

 for(i in 1:nrow(xfix)){   
   j=which(dd1$Psi$obs==xfix$obs2[i] & dd1$Psi$method==xfix$method[i]);  
   if (length(j)>0) dd1$Psi$xfix[j]=0
 }

for(i in 1:nrow(xfix)){   
  j=which(dd1$Psi$obs==xfix$obs2[i]);  
  if (length(j)>0) dd1$Psi$xfix[j]=0
}

indices2p <- which(dd1$p$xfix == 0)
indices2psi <- which(dd1$Psi$xfix == 0)
valuesp <- rep(0, length(indices2p))
valuespsi <- rep(0, length(indices2psi))

psi.fixed=list(formula=~1, fixed=list(index=indices2psi, value=valuespsi))
p.fixed=list(formula=~1, fixed=list(index=indices2p,value=valuesp))

#fix detection to 0 for species that aren't available for detection
m1=mark(pd1,dd1,model="Occupancy",model.parameters=list(Psi=list(formula=~1), p=list(formula=~1)))
#detection transects: 0.08
#detection ground cameras: 0.79
#detection canopy cameras: 0.84
#2113.748
m1=mark(pd1,dd1,model="Occupancy",model.parameters=list(Psi=psi.fixed,p=p.fixed))

#remove species that aren't available for detection
#indices3 <- as.numeric(row.names(dta1[dta1$ch == "..",]))
#dta1a <- dta1[-indices3,]
#process data
#pd1a=process.data(dta1a,model="Occupancy", groups = c('transect2','obs2', 'method'))
#design data
#dd1a=make.design.data(pd1a)
#m1a=mark(pd1a,dd1a,model="Occupancy",model.parameters=list(Psi=list(formula=~method),p=list(formula=~method+spef+tief)))


################
## Figure 2 + 4
##species richness estimates depending on different combinations of data
##occupancy * number of available species for those methods = species richness
names <- c("Transect", "Ground", "Arboreal", "\nTransect \n+ Ground", "Transect \n+ Arboreal", "Ground \n+ Arboreal", "All")

par(mar=c(4,4.2,2,2)) #bottom, left, top, right
est <- c((8/21), (12/31), (12/23), (11/31), (12/30), (8/36), (0/36))
unique <- c((3/21), (13/31), (7/23), (18/31), (11/30), (25/36), (36/36))
data <- rbind(unique, est)
bp <- barplot(data, col = c("darkblue", "red"),
              beside = F, ylab = "Percent of Species Detected",
              ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2, las = 1)
mtext(text = names, side = 1, cex = 1.2, at = bp, line = 1.5)
freq <- c("n=21", "n=31", "n=23", "n=31", "n=30", "n=36", "n=36")
text(x = bp, y = est+unique, labels = freq, pos = 3, cex = 1.2)
mtext(text = "n=36", side = 3, at = 7.9, cex = 1.2, line = 0.2)


par(mfrow=c(2,1))
par(mar=c(1.5,5.5,2.5,2))
est <- c(4.33, 9.06, 8.11, 10.83, 9.27, 13.88, 14.74)
se <- c(1.05, 0.67, 0.60, 0.74, 0.66, 0.74, 0.76)

ylab1 <- c("Estimated \nSpecies Richness")
bp <- barplot(est, 
              cex.main = 1.5,beside = T, ylab = ylab1,
              ylim = c(0,16), cex.lab = 1.2, cex.axis = 1.2, las = 1)
means <- est
errors <- se
sub <- (means-errors)
plus <- means+errors
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "A.")

est <- c(0.34, 0.70, 0.73, 0.66, 0.72, 0.76, 0.74)
se <- c(0.08, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02)
par(mar=c(2.5,5.5,1.5,2))
ylab2 <- c("Estimated \nDetection Probaility")
bp <- barplot(est, 
              beside = T, ylab = ylab2,
              ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2, las = 1)
means <- est
errors <- se
sub <- (means-errors)
plus <- means+errors
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
mtext(text = names, side = 1, cex = 1.2, at = bp, line = 1.5)
title(adj = 0, "B.")
#text(seq(0.8,8,by=1.2), par("usr")[3]-0.075, 
#     srt = 15, adj= 1, xpd = TRUE,
#     labels = names, cex=1)




##########################
##detection 2
## diurnal species - transects vs. cameras
############################
dta2 <- rbind(dta_LT, dta_CT)
#subset for only diurnal species
sl <- read.csv("MammalList2017.csv", header = T)
diurnal <- sl[sl$Active == "Diurnal",]
dta2 <- subset(dta2, dta2$obs2 %in% diurnal$Species)

dta2$obs2 <- as.factor(dta2$obs2)
#scale effort covariates
dta2$stief1 <- scale(dta2$tief1)
dta2$stief2 <- scale(dta2$tief2)
dta2$sspef <- scale(dta2$spef)
#process data
pd2=process.data(dta2, model="Occupancy", groups = c('transect2','obs2', 'method'))
#design data
dd2=make.design.data(pd2)
#initial model
m2=mark(pd2,dd2,model="Occupancy",model.parameters=list(Psi=list(formula=~method),p=list(formula=~method+sspef+stief)))
#detection line transects:0.17
#detection cameras:0.87
#1145.031
#effort covariate not significant


#####################
##detection 3
##species both arboreal and terrestrial - ground cameras vs. canopy cameras
#####################

dta3 <- rbind(dta_CT_canopy, dta_CT_ground)
#subset for only species that are on the ground and in trees
sl <- read.csv("MammalList2017.csv", header = T)
both <- sl[sl$Location == "Both",]
dta3 <- subset(dta3, dta3$obs2 %in% both$Species)

dta3$obs2 <- as.factor(dta3$obs2)
#scale effort covariates
dta3$stief1 <- scale(dta3$tief1)
dta3$stief2 <- scale(dta3$tief2)
dta3$sspef <- scale(dta3$spef)
#process data
pd3=process.data(dta3, model="Occupancy", groups = c('transect2','obs2', 'method'))
#design data
dd3=make.design.data(pd3)
#initial model
m3=mark(pd3,dd3,model="Occupancy",model.parameters=list(Psi=list(formula=~method),p=list(formula=~method+sspef+stief)))
#detection ground cameras:0.63
#detection canopy cameras:0.69
#spatial effort significant



########################
##detection 4
## primates - transects vs. ground cameras vs. canopy cameras
########################

dta4 <- rbind(dta_LT, dta_CT_canopy, dta_CT_ground)
#subset for only diurnal species
sl <- read.csv("MammalList2017.csv", header = T)
all <- sl[sl$Active == "Diurnal" & sl$Location == "Both",]
all2 <- all[c(1:5,10,11),]
dta4 <- subset(dta4, dta4$obs2 %in% all2$Species)
dta4$obs2 <- as.factor(dta4$obs2)

#scale effort covariates
dta4$stief1 <- scale(dta4$tief1)
dta4$stief2 <- scale(dta4$tief2)
dta4$sspef <- scale(dta4$spef)
dta4$seffort1 <- scale(dta4$effort1)
dta4$seffort2 <- scale(dta4$effort2)

#process data
pd4=process.data(dta4, model="Occupancy", groups = c('transect2','obs2', 'method'))
#design data
dd4=make.design.data(pd4)
#initial model
m4null=mark(pd4,dd4,model="Occupancy",model.parameters=list(Psi=list(formula=~1),p=list(formula=~1)))
m4=mark(pd4,dd4,model="Occupancy",initial = m4null,model.parameters=list(Psi=list(formula=~method),p=list(formula=~method+sspef+stief)))
#detection line transects: 0.05
#detection ground cameras: 0.77
#detection canopy cameras: 0.88
#1238.741 AICc
#effort covariate not significant
m4sp=mark(pd4,dd4,model="Occupancy",model.parameters=list(Psi=list(formula=~method+obs2),p=list(formula=~method+sspef+stief)))



####################
##detection 5
######################
## just canopy cameras individually
dta5 <- dta_CT_ind_canopy

#subset for only arboreal/both species
sl <- read.csv("MammalList2017.csv", header = T)
all <- sl[sl$Location != "Terrestrial",]
dta5 <- subset(dta5, dta5$obs2 %in% all$Species)

dta5$obs2 <- as.factor(dta5$obs2)
dta5$transect2 <- as.factor(dta5$transect2)

dta5$stief1 <- scale(dta5$tief1)
dta5$stief2 <- scale(dta5$tief2)

#process data
pd=process.data(dta5, model="Occupancy", groups = c('transect2','obs2'))
#design data
dd=make.design.data(pd)
#initial model
m5=mark(pd,dd,model="Occupancy",model.parameters=list(Psi=list(formula=~1),p=list(formula=~obs2+height)))


#####
#combine figure 4-6 into one figure
#####
par(mfrow=c(2,3))
#library(rphylopic)

#duikerimg <- image_data('bb357997-4d4d-4e73-b853-d95a7a47ec68', size= '64')[[1]]
################
## Figure 4
#detection probability and species richness for diurnal species 
#comparing transects and cameras
#rnames2=rownames(m2$results$real)
#cindices2 <- which(grepl("camera", rnames2))
#lindices2 <- which(grepl("transect", rnames2))
#bp3mean <- c(m2$results$real$estimate[lindices2[1]], m2$results$real$estimate[cindices2[1]])
#bp3se <- c(m2$results$real$se[lindices2[1]], m2$results$real$se[cindices2[1]])
#occ3mean <- c(m2$results$real$estimate[lindices2[3]], m2$results$real$estimate[cindices2[3]])
#occ3se <- c(m2$results$real$se[lindices2[3]], m2$results$real$se[cindices2[3]])

#species richness
sr <- c(4.66, 8.37)
se <- c(1.22, 0.58)
par(mar=c(2,5,3,2)) #bottom, left, top, right
bp <- barplot(sr,
              beside = T, ylab = "Estimated Species Richness",
              ylim = c(0,10), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- sr
errors <- se
sub <- (means-errors)
plus <- means+errors
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "A.")
#add_phylopic_base(duikerimg, x = 0.1, y = 0.9, ysize = 0.15, alpha = 0.2, color = 'black')

#detection
par(mar=c(4.5,5,2.5,2)) #bottom, left, top, right
det <- c(0.17, 0.87)
detse <- c(0.17, 0.12)
bp <- barplot(det,
              beside = T, ylab = "Estimated Detection Probability",
              names.arg = c("Line Transect", "Camera Traps"),cex.names = 1.2,
              ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- det
errors <- detse
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "B.")
mtext(text = "Diurnal Species", side = 1, cex = 1, at = 1.3, line = 3)


################
## Figure 5
#detection probability and species richness for arb/terr species 
#comparing ground cameras and canopy cameras

#species richness
sr <- c(5.39, 5.65)
se <- c(0.58, 0.52)
par(mar=c(2,5,3,2)) #bottom, left, top, right
bp <- barplot(sr,
              beside = T, ylab = NA,
              ylim = c(0,10), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- sr
errors <- se
sub <- (means-errors)
plus <- means+errors
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "C.")

#detection
par(mar=c(4.5,5,2.5,2)) #bottom, left, top, right
det <- c(0.63, 0.69)
detse <- c(0.05, 0.04)
bp <- barplot(det,
              beside = T, ylab = NA,
              names.arg = c("Ground Camera", "Arboreal Camera"),cex.names =1.2,
              ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- det
errors <- detse
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "D.")
mtext(text = "Species using Terr. & Arb. Substrates", side = 1, cex = 1, at = 1.3, line = 3)


##########
## Figure 6
## all methods
## primate species (7)
## species richness and detection probability

#species richness
sr <- c(2.24, 2.74, 2.19)
se <- c(0.89, 0.42, 0.32)
par(mar=c(2,5,3,2)) #bottom, left, top, right
bp <- barplot(sr,
              cex.main = 1.5,beside = T, ylab = NA,
              ylim = c(0,10), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- sr
errors <- se
sub <- (means-errors)
plus <- means+errors
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "E.")

#detection
par(mar=c(4.5,5,2.5,2)) #bottom, left, top, right
det <- c(0.01, 0.89, 0.93)
detse <- c(0.02, 0.12, 0.08)
bp <- barplot(det,
              beside = T, ylab = NA,
              ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, las = 1)
means <- det
errors <- detse
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "F.")
mtext(text = c("Line \nTransect", "Ground \nCamera", "Arboreal \nCamera"), side = 1, cex = 0.8, at = bp, line = 1.5)
mtext(text = "Primates", side = 1, cex = 1, at = 2, line = 3)











################################################
#### Old plots - code
################################################

f6 <- read.csv("Figure6Data.csv", header = T)
f6$Species <- as.character(f6$Species)
est <- cbind(f6$est_transect, f6$est_ground, f6$est_canopy)
se <- cbind(f6$se_transect, f6$se_ground, f6$se_canopy)
par(mar = c(7, 4, 8, 2) + 0.2)
par(xpd=TRUE)
bp <- barplot(c(est[1,],est[2,],est[3,],est[4,],est[5,],est[6,],est[7,],est[8,],est[9,],est[10,],est[11,],est[12,],est[13,],est[14,]), 
              space = c(0,0,0,1,rep(c(0,0,1),12), 0,0),
              cex.main = 1.5,beside = T, 
              col = rep(c("blue", "red", "green"),14),
              ylim = c(0.0,1.0), cex.lab = 1.2, cex.axis = 1.2,
              las = 1, xlab = "", names.arg = NA)
names = c(f6$Species[1],NA,NA,f6$Species[2],NA,NA,f6$Species[3],NA,NA,f6$Species[4],NA,NA,f6$Species[5],NA,NA,f6$Species[6],NA,NA,f6$Species[7],NA,NA,f6$Species[8],NA,NA,f6$Species[9],NA,NA,f6$Species[10],NA,NA,f6$Species[11],NA,NA,f6$Species[12],NA,NA,f6$Species[13],NA,NA,f6$Species[14],NA,NA)
text(seq(1.6,57,by=1.33), par("usr")[3]-0.025, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = names, cex=0.70)

means <- c(est[1,],est[2,],est[3,],est[4,],est[5,],est[6,],est[7,],est[8,],est[9,],est[10,],est[11,],est[12,],est[13,],est[14,]) 
errors <-c(se[1,],se[2,],se[3,],se[4,],se[5,],se[6,],se[7,],se[8,],se[9,],se[10,],se[11,],se[12,],se[13,],se[14,])
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
legend(-2,1.25, c("Line Transects","Ground Cameras", "Canopy Cameras"), cex=1, fill=c("blue", "red", "green"), bty = 'n')



####################################
## Detection Probability Bar Plot ##
####################################

bp1mean <- c(mLT$results$real$estimate[1], mGCT$results$real$estimate[1], mCCT$results$real$estimate[1])
bp1se <- c(mLT$results$real$se[1], mGCT$results$real$se[1], mGCT$results$real$se[1])         

rnames=rownames(m1$results$real)
cindices <- which(grepl("canopy", rnames))
gindices <- which(grepl("ground", rnames))
lindices <- which(grepl("transect", rnames))
bp2mean <- c(m1$results$real$estimate[lindices[1]],m1$results$real$estimate[gindices[1]],m1$results$real$estimate[cindices[1]])
bp2se <- c(m1$results$real$se[lindices[1]],m1$results$real$se[gindices[1]],m1$results$real$se[cindices[1]])

rnames2=rownames(m2$results$real)
cindices2 <- which(grepl("camera", rnames2))
lindices2 <- which(grepl("transect", rnames2))
bp3mean <- c(m2$results$real$estimate[lindices2[1]], m2$results$real$estimate[cindices2[1]])
bp3se <- c(m2$results$real$se[lindices2[1]], m2$results$real$se[cindices2[1]])

rnames3=rownames(m3$results$real)
gindices3 <- which(grepl("ground", rnames3))
cindices3 <- which(grepl("canopy", rnames3))
bp4mean <- c(m3$results$real$estimate[gindices3[1]], m3$results$real$estimate[cindices3[1]])
bp4se <- c(m3$results$real$se[gindices3[1]], m3$results$real$se[cindices3[1]])

rnames4=rownames(m4$results$real)
lindices4 <- which(grepl("transect", rnames4))
cindices4 <- which(grepl("ground", rnames4))
gindices4 <- which(grepl("canopy", rnames4))
bp5mean <- c(m4$results$real$estimate[lindices4[1]], m4$results$real$estimate[gindices4[1]], m4$results$real$estimate[cindices4[1]])
bp5se <- c(m4$results$real$se[lindices4[1]], m4$results$real$se[gindices4[1]], m4$results$real$se[cindices4[1]])


bp <- barplot(c(bp1mean,bp2mean,bp3mean,bp4mean,bp5mean), space = c(0,0,0,1,0,0,1,0,1,0,1,0,0),
        main="Detection Probabilities", ylab="Detection Probability", 
        cex.main = 1.5,beside = T, 
        col = c("gray","orange", "blue", "gray", "orange", "blue", "gray", "darkred", "orange", "blue", "gray", "orange", "blue"),
        names.arg = c(NA,"All Species\nInd. Model", NA, NA, "All Species\nJoint Model", NA, "          Diurnal\n          Species", NA, "        Arb/Ter\n        Species", NA, NA, "Diurnal\nArb/Ter Species", NA),
        ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2)
means <- c(bp1mean, bp2mean, bp3mean, bp4mean, bp5mean)
errors <- c(bp1se, bp2se, bp3se, bp4se, bp5se)
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1


segments(bp, sub, bp, plus, lwd = 2)

arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)

legend("topleft", c("Line Transects","Ground Cameras","Canopy Cameras","All Cameras"), cex=1, 
       fill=c("gray", "orange", "blue", "darkred"))





