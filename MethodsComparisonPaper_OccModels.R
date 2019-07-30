####################################
## Methods Comparion Paper
## Line transects, ground cameras, arboreal cameras
## Moore et al. ()
## Updated April 30, 2019
## Occupancy Models
#####################################

rm(list=ls())
#load package
library(RPresence)

#read in data
load('MethodsComparisonPublicationData.RData')
#dta_all - all methods
#dta_L - line transect
#dta_G - ground cameras
#dta_C - canopy cameras
#dta_GC - ground and canopy cameras
#dta_LC - line transect + canopy camera
#dta_LG - line transect + ground camera 
#CeLh
#CeMi
#CeMo
#CoAn
#LoAl
#PaTr
#dta_ind - individual canopy cameras

######
#Figure 1 - study site map (created in ArcGIS)
######

######
#Figure 2 - raw species richness for each method/combination of methods
######
names <- c("Transect", "Ground", "Arboreal", "\nTransect \n+ Ground", "Transect \n+ Arboreal", "Ground \n+ Arboreal", "All")

par(mar=c(4,4.2,2.2,2)) #bottom, left, top, right
est <- c((8/20), (12/30), (11/22), (11/30), (12/29), (8/35), (0/35))
unique <- c((3/20), (13/30), (6/22), (18/30), (10/29), (24/35), (35/35))
data <- rbind(unique, est)
bp <- barplot(data, col = c("darkblue", "red"),
              beside = F, ylab = "Percent of Species Detected",
              ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2, las = 1)
mtext(text = names, side = 1, cex = 1.2, at = bp, line = 1.5)
freq <- c("n=20", "n=30", "n=22", NA, "n=29", "n=35", NA)
text(x = bp, y = est+unique, labels = freq, pos = 3, cex = 1.2)
mtext(text = "n=35", side = 3, at = 7.9, cex = 1.2, line = 0.2)
mtext(text = "n=30", side = 3, at = 4.3, cex = 1.2, line = 0.2)


#########################################
## All Species Analysis - Analysis #1  ##
#########################################
#comparison of all methods

#data from all methods
#run model based on different covariates to choose best model for these data
#then use best model but using data from different combination of methods
mods=list(); i=1

mods[[i]]=occMod(model=list(psi~Transect, p~mass),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~trail,    p~mass),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~access,   p~mass),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~minelev,  p~mass),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~maxelev,  p~mass),data=dta_all,type='so',randinit=10);i=i+1;

mods[[i]]=occMod(model=list(psi~Transect, p~mass+cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~trail,    p~mass+cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~access,   p~mass+cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~minelev,  p~mass+cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~maxelev,  p~mass+cat),data=dta_all,type='so',randinit=10);i=i+1;

mods[[i]]=occMod(model=list(psi~Transect, p~group),data=dta_all,type='so',randinit=10);i=i+1; 
mods[[i]]=occMod(model=list(psi~trail,    p~group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~access,   p~group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~minelev,  p~group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~maxelev,  p~group),data=dta_all,type='so',randinit=10);i=i+1;

mods[[i]]=occMod(model=list(psi~Transect, p~mass+group),data=dta_all,type='so',randinit=10);i=i+1; 
mods[[i]]=occMod(model=list(psi~trail,    p~mass+group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~access,   p~mass+group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~minelev,  p~mass+group),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~maxelev,  p~mass+group),data=dta_all,type='so',randinit=10);i=i+1;

mods[[i]]=occMod(model=list(psi~Transect, p~cat),data=dta_all,type='so',randinit=10);i=i+1; 
mods[[i]]=occMod(model=list(psi~trail,    p~cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~access,   p~cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~minelev,  p~cat),data=dta_all,type='so',randinit=10);i=i+1;
mods[[i]]=occMod(model=list(psi~maxelev,  p~cat),data=dta_all,type='so',randinit=10);

# create AIC table of model results and print
results<-createAicTable(mods); 
results$table$wgt=round(results$table$wgt,4)          # round off weights to 4 decimal places
results$table$modlike=round(results$table$modlike,4)  #  round off model likelihood values to 4 places
print(results$table)

#best model psi(trail)p(Species) - lowest AIC score
m_all=occMod(model=list(psi~maxelev, p~mass+group),data=dta_all,outfile = "allmethods",type='so',randinit=10)

#run this model for the rest of the sets of data combining different methods

#all cameras (ground and arboreal)
m_GC=occMod(model=list(psi~maxelev, p~mass+group),data=dta_GC,type='so',randinit=10)

#ground cameras + line transects
m_GL=occMod(model=list(psi~maxelev, p~mass+group),data=dta_LG,type='so',randinit=10)

#arboreal cameras + line transect
m_CL=occMod(model=list(psi~maxelev, p~mass+group),data=dta_LC,type='so',randinit=10)

#ground cameras only
m_G=occMod(model=list(psi~maxelev, p~mass+group),data=dta_G,type='so',randinit=10)

#arboreal cameras only
m_C=occMod(model=list(psi~maxelev, p~mass+group),data=dta_C,type='so',randinit=10)

#line transects only
m_L=occMod(model=list(psi~maxelev, p~mass+group),data=dta_L,type='so',randinit=10)


mean(m_all$real$psi$est[1:18]) #0.4343492 #0.1376938 #1, 3
mean(m_GC$real$psi$est[1:18]) #0.3978528 #0.1318527 #3, 2
mean(m_GL$real$psi$est[1:18]) #0.3810028 #0.1558962 #5, 4
mean(m_CL$real$psi$est[1:18]) #0.3958639 #0.1681748 #4, 6
mean(m_G$real$psi$est[1:18]) #0.2535904 #0.1173737 #7, 1
mean(m_C$real$psi$est[1:18]) #0.3477924 #0.1591870 #6, 5
mean(m_L$real$psi$est[1:18]) #0.4013121 #0.4081930 #2, 7
 
mean(unique(m_all$real$p$est)) #0.6814845 #0.2372511 #2, 1
mean(unique(m_GC$real$p$est)) #0.7060167 #0.2385052 #1, 2
mean(unique(m_GL$real$p$est)) #0.6145389 #0.3096234 #4, 4
mean(unique(m_CL$real$p$est)) #0.5400645 #0.2747944 #6, 3
mean(unique(m_G$real$p$est)) #0.6618811 #0.3450557 #3, 5
mean(unique(m_C$real$p$est)) #0.5776393 #0.4243497 #5, 7
mean(unique(m_L$real$p$est)) #0.2435697 #0.4075727 #7, 6

#######
#Figure 3
#Estimated species richness and detection probability by method or combination of methods
#######
names <- c("Transect", "Ground", "Arboreal", "\nTransect \n+ Ground", "Transect \n+ Arboreal", "Ground \n+ Arboreal", "All")
par(mfrow=c(2,1))
par(mar=c(1.5,5.5,2.5,2))

est_occ <- c(mean(m_L$real$psi$est[1:18]), mean(m_G$real$psi$est[1:18]), mean(m_C$real$psi$est[1:18]), 
         mean(m_GL$real$psi$est[1:18]), mean(m_CL$real$psi$est[1:18]), mean(m_GC$real$psi$est[1:18]), 
         mean(m_all$real$psi$est[1:18]))
se_occ <- c(sqrt(sum(m_L$real$psi$se[1:18]^2)), sqrt(sum(m_G$real$psi$se[1:18]^2)), sqrt(sum(m_C$real$psi$se[1:18]^2)), 
        sqrt(sum(m_GL$real$psi$se[1:18]^2)), sqrt(sum(m_CL$real$psi$se[1:18]^2)), sqrt(sum(m_GC$real$psi$se[1:18]^2)), 
        sqrt(sum(m_all$real$psi$se[1:18]^2)))

ylab1 <- c("Estimated \nOccupancy Probability")
bp <- barplot(est_occ, 
              cex.main = 1.5,beside = T, ylab = ylab1,
              ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2, las = 1)
means <- est_occ
errors <- se_occ
sub <- (means-errors)
for(i in 1:length(sub))if(sub[i]<0)sub[i]=0
plus <- means+errors
for(i in 1:length(plus))if(plus[i]>1)plus[i]=1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
title(adj = 0, "A.")

est_p <- c(mean(unique(m_L$real$p$est)), mean(unique(m_G$real$p$est)), mean(unique(m_C$real$p$est)), 
           mean(unique(m_GL$real$p$est)), mean(unique(m_CL$real$p$est)), mean(unique(m_GC$real$p$est)), 
           mean(unique(m_all$real$p$est)))
se_p <- c(sqrt(sum(unique(m_L$real$p$se)^2)), sqrt(sum(unique(m_G$real$p$se)^2)), sqrt(sum(unique(m_C$real$p$se)^2)), 
          sqrt(sum(unique(m_GL$real$p$se)^2)), sqrt(sum(unique(m_CL$real$p$se)^2)), sqrt(sum(unique(m_GC$real$p$se)^2)), 
          sqrt(sum(unique(m_all$real$p$se)^2)))

par(mar=c(2.5,5.5,1.5,2))
ylab2 <- c("Estimated \nDetection Probaility")
bp <- barplot(est_p, 
              beside = T, ylab = ylab2,
              ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2, las = 1)
means <- est_p
errors <- se_p
sub <- (means-errors)
for(i in 1:length(sub))if(sub[i]<0)sub[i]=0
plus <- means+errors
for(i in 1:length(plus))if(plus[i]>1)plus[i]=1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)
mtext(text = names, side = 1, cex = 1.2, at = bp, line = 1.5)
title(adj = 0, "B.")


######################################
## Primate Analysis - Analysis #2   ##
######################################
#files for 6 primate species
#CeLh,CeMi,CeMo,CoAn,LoAl,PaTr

#CeLh_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=CeLh,type='so.mm',randinit=10, outfile = "lhoest")
#CeMi_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=CeMi,type="so.mm",randinit=10, outfile = "blue")
#CeMo_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=CeMo,type='so.mm',randinit=10, outfile = "mona")
#CoAn_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=CoAn,type='so.mm',randinit=10, outfile = "lhoest")
#LoAl_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=LoAl,type='so.mm',randinit=10, outfile = "mangabey")
#PaTr_m=occMod(model=list(psi~1,theta~1,p~DEVICE),data=PaTr,type="so.mm",randinit=10, outfile = "chimp")

#CeLh: psi = 0.9465 (SE 0.0541)
#theta 1
#0.2641 (SE 0.0756) transect
#0.8218 (SE 0.0662) ground
#0.6457 (SE 0.0823) arboreal

#CeMi: psi = 0.9626 (SE 0.0563)
#theta 1
#0.3174 (SE 0.0799) transect
#0.3463 (SE 0.0818) ground
#0.6926 (SE 0.0822) arboreal

#mona: naive 0.0556
#only 1 detection 1 time by line transect

#colobus: psi = 0.6878 +\- 0.5665
#theta 1
#0.1615 (SE 0.1501) transect
#0.0404 (SE 0.0513) ground
#0.0404 (SE 0.0513) arboreal

#mangabey: naive 0.1111
#only 2 detections by line transect

#chimp: psi = 0.6944 +/- 0.3140
#theta 0.7406 (SE 0.6044)
#0.0540 (SE 0.0643) transect
#0.4861 (SE 0.3537) ground
#0.0540 (SE 0.0643) arboreal

#detection
par(mfrow=c(1,1))
par(mar=c(4.5,5,2.5,2)) #bottom, left, top, right
det <- c(0.2641, 0.8218, 0.6457, 
         0.3174, 0.3463, 0.6926, 
         0.1615, 0.0404, 0.0404, 
         0.0540, 0.4861, 0.0540)
detse <- c(0.0756, 0.0662, 0.0823,
           0.0799, 0.0818, 0.0822,
           0.1501, 0.0513, 0.0513, 
           0.0643, 0.3537, 0.0643)
bp <- barplot(det,
              beside = F, ylab = "Estimated Detection Probability",
              space = c(0,0,0,1,0,0,1,0,0,1,0,0),
              names.arg = c(NA, "C. lhoesti", NA, NA, "C. mitis", NA, NA, "C. angolensis", NA, NA, "P. troglodytes", NA),
              cex.names = 1.2,
              ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, las = 1, 
              col = c("red", "darkgreen", "blue"))
legend("topright", pch = 15, cex=1.2,
       legend = c("Line Transect", "Ground Camera", "Arboreal Camera"), 
       col = c("red", "darkgreen", "blue"))
means <- det
errors <- detse
sub <- (means-errors)
sub[sub < 0] <- 0
plus <- means+errors
plus[plus > 1] <- 1
segments(bp, sub, bp, plus, lwd = 2)
arrows(bp, sub, bp, plus,lwd = 2, angle = 90, code = 3, length = 0.05)


##################################################
## Arboreal Camera Height in Tree - Analysis #3 ##
##################################################
#dta_ind - individual arboreal camera file
dta_ind$unitcov$height <- as.numeric(dta_ind$unitcov$height)
ind_m=occMod(model=list(psi~maxelev + height, p~mass+group+height),data=dta_ind,type='so',randinit=10, outfile = "height")

