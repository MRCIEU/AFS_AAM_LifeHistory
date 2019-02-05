# This script is adapted from code provided by Jack Bowden, Suzi Gage and Amy Taylor
# Radial MR is from https://github.com/WSpiller/RadialMR and leave one out uses MRBase
# Note - have to clear R and load next dataset if want to use a different outcome (check right working directory as well)
# Note - be careful where code uses betaXG (not absolute) and BXG (absolute)

# Clear R
rm(list=ls(all=TRUE))

#################################################################################
# Pick Data
#################################################################################

##Change to your working directory (note forward slashes) 
setwd("M:/data/ukbiobank/_devs/UKBIOBANK_Phenotypes_App_6326/data/AAM123")  

# read in dataset want for exposure and whether 500k or 350k as outcome e.g.:
dat = read.csv ("dataset_AFB_500_AAM123.csv", header=T)

#################################################################################
# Prepare data for MR
#################################################################################

# if using AAM123
#dat$SNP <- as.vector(dat$SNP) 
#dat <- dat[!grepl("rs4946632", dat$SNP),] 
#dat <- dat[complete.cases(dat), ] 
#dat <- na.omit(dat) 

# If doing no BMI SNP sensitivity analysis
#dat <- dat[!grepl("rs10938397", dat$SNP),] 
#dat <- dat[!grepl("rs12446632", dat$SNP),]
#dat <- dat[!grepl("rs2947411", dat$SNP),]
#dat <- dat[!grepl("rs3101336", dat$SNP),]
#dat <- dat[!grepl( "rs543874", dat$SNP),]
#dat <- dat[!grepl("rs7103411", dat$SNP),]
#dat <- dat[!grepl("rs7138803", dat$SNP),]
#dat <- dat[!grepl("rs7514705", dat$SNP),]
#dat <- dat[!grepl("rs8050136", dat$SNP),]

# checks:
#library(data.table)
#setDT(dat, key='SNP')[.('rs10144321')] # a snp that should still be in the data
#setDT(dat, key='SNP')[.('rs8050136')] # snp that shouldnt be
#is.na(setDT(dat, key='SNP')[.('rs8050136')]) # any missing data
dat <- dat[complete.cases(dat), ] # remove rows with missing data
dat <- na.omit(dat) # remove any na's

BetaXG   = as.numeric(data$beta_x) # genetic association with exposure
BetaYG   = as.numeric(data$beta_y) # standard errors 
seBetaYG = as.numeric(data$SE_y) # genetic association with outcome
seBetaXG = as.numeric(data$SE_x) # standard errors

BYG             = BetaYG*sign(BetaXG) # Pre-processing steps to ensure all  
BXG             = abs(BetaXG)         # gene--exposure estimates are positive


#################################################################################
# STATISTICS
#################################################################################

# Mean F stat
F   = BXG^2/seBetaXG^2
mF  = mean(F)
mF 

# I-squared function
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

# I-sq results
Isq(BXG,seBetaXG) # for an unweighted MR-Egger estimate 
Isq(BXG/seBetaYG,seBetaXG/seBetaYG) # for a weighted MR-Egger estimate

#################################################################################
# IVW approach
#################################################################################

IVWfit      = summary(lm(BYG ~ -1+BXG,weights=1/seBetaYG^2)) 

DF      = length(BYG)-1
IVWBeta = IVWfit$coef[1,1]
SE      = IVWfit$coef[1,2]/IVWfit$sigma
IVW_p   = 2*(1-pt(abs(IVWBeta/SE),DF))
IVW_CI  = IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE

# Cochran's Heterogeneity statistic
# for IVW method

DF      = length(BYG)-1
phi_IVW = IVWfit$sigma^2
QIVW    = DF*phi_IVW
Qp      = 1-pchisq(QIVW,DF)

phi_IVW
QIVW
Qp

# IVWResults = (point estimate, corrected standard error, 
# 95% Confidence interval, t-statistic, p-value) 
IVWResults = c(IVWBeta,SE,IVW_CI,IVWBeta/SE,IVW_p)
IVWResults_test <- as.data.frame(t(IVWResults))
print(IVWResults)

#################################################################################
# MR-Egger regression (with MAF corrected weights)  
#################################################################################

MREggerFit      = summary(lm(BYG ~ BXG,weights=1/seBetaYG^2))

# Inference with correct standard errors

MREggerBeta0   = MREggerFit$coef[1,1]
MREggerBeta1   = MREggerFit$coef[2,1]
SE0            = MREggerFit$coef[1,2]/MREggerFit$sigma
SE1            = MREggerFit$coef[2,2]/MREggerFit$sigma
DF             = length(BYG)-2
MRBeta0_p      = 2*(1-pt(abs(MREggerBeta0/SE0),DF))
MRBeta1_p      = 2*(1-pt(abs(MREggerBeta1/SE1),DF))
MRBeta0_CI     = MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
MRBeta1_CI     = MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1

#################################################################################
# Weighted Median
#################################################################################

## Function
weighted.median <- function(betaIV.in, weights.in) {
  betaIV.order = betaIV.in[order(betaIV.in)]
  weights.order = weights.in[order(betaIV.in)]
  weights.sum = cumsum(weights.order)-0.5*weights.order
  weights.sum = weights.sum/sum(weights.order)
  below = max(which(weights.sum<0.5))
  weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(weighted.est) }

weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
  med = NULL
  for(i in 1:1000){
    betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
    betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
    betaIV.boot = betaYG.boot/betaXG.boot
    med[i] = weighted.median(betaIV.boot, weights.in)
  }
  return(sd(med)) 
}

## Analysis
betaIV   = BYG/BXG  
weights  = (seBetaYG/BXG)^-2 
betaWM   = weighted.median(betaIV, weights) 
sebetaWM = weighted.median.boot(BXG, BYG, seBetaXG, seBetaYG, weights) 
t     = betaWM/sebetaWM
p     = 2*(1-pt(abs(t),length(BYG)-1))
WMresults = data.frame(Estimate=betaWM,Std.Error=sebetaWM,t,p)
CI.WM     = betaWM + c(-1,1)*sebetaWM
WMresults_test <- WMresults
WMresults_test$LCI <- betaWM-(sebetaWM*1.96)
WMresults_test$UCI <- betaWM+(sebetaWM*1.96)

#################################################################################

plot(BXG,BYG,pch=19,cex=1.4,main="Scatter plot",xlab="",ylab="")
mtext(side=2,expression(paste("Gene-outcome assoc.")),line=2,cex=1.5)
mtext(side=1,expression(paste("Gene-exposure assoc.")),line=2.8,cex=1.5)

MREggerFit  = lm(BYG ~ BXG,weights=1/seBetaYG^2)
IVWFit      = lm(BYG ~ -1+BXG,weights=1/seBetaYG^2)

lines(BXG,MREggerFit$fitted.values,col="blue",lwd=2)
lines(BXG,IVWFit$fitted.values,col="red")
lines(BXG,betaWM*BXG,col="black",lwd=2)

legend("topleft", c("IVW slope","MR-Egger slope","Weighted Median Slope"),
       col=c("red","blue","black"),lwd=3,cex=1.1,bty="n")

dev.print(device=pdf,file=plotfile[1])

MREggerResults     = matrix(nrow = 2,ncol = 6)
MREggerResults[1,] = c(MREggerBeta0,SE0,MRBeta0_CI,MREggerBeta0/SE0,MRBeta0_p)
MREggerResults[2,] = c(MREggerBeta1,SE1,MRBeta1_CI,MREggerBeta1/SE1,MRBeta1_p)

BetaIV = BYG/BXG
Xrange = range(BetaIV,MRBeta1_CI)
Yrange = range(BXG/seBetaYG)

plot(BetaIV,BXG/seBetaYG,pch=19,cex=1.4,main="Funnel plot",xlab="",ylab="",xlim=Xrange,ylim=Yrange)
lines(rep(IVWBeta,2),c(0,Yrange[2]),col="red",lwd=2)
lines(rep(MREggerBeta1,2),c(0,Yrange[2]),col="blue",lwd=2)
lines(rep(betaWM,2),c(0,Yrange[2]),col="black",lwd=2)
mtext(side=1,expression(paste("Causal estimate")),line=3.5,cex=1.5)
mtext(side=2,expression(paste("MAF-corrected Instrument strength")),line=2,cex=1.5)

points(IVWBeta,10,cex=2,col="red",pch=15)
lines(IVW_CI,rep(10,2),lwd=2,col="red")
points(MREggerBeta1,6,cex=2,col="blue",pch=22)
lines(MRBeta1_CI,rep(6,2),lwd=2,col="blue")

points(betaWM,8,cex=2,col="black",pch=23)
lines(CI.WM,rep(8,2),lwd=2,col="black")


points(BetaIV,BXG/seBetaYG,pch=19,cex=1.4)

legend("topright",c("IVW","MR-Egger","Weighted Median"),lwd=2,
       col=c("red","blue","black"),cex=1,bty="n",pch=c(15,22,23))

dev.print(device=pdf,file=plotfile[2])

#################################################################################
# Summary of results
#################################################################################

# IVW approach
#summary(IVWFit)
# IVW with corrected standard errors - Beta,SE,CI,t,p
IVWResults

# MR-Egger approach
#summary(MREggerFit)
print(MREggerResults) # MREggerResults = (point estimate, corrected standard error, 95% Confidence interval, t-statistic, p-value) for intercept (row 1) and slope (row 2).

# Weighted Median - may be different each time as bootstrapped
WMresults

# Cochran's Q statistic for IVW (and p-value)
QIVW 
Qp 

# Mean F statistic for IVW
mF  

# I^2_GX for MR-Egger
Isq(BXG,seBetaXG) # unweighted
Isq(BXG/seBetaYG,seBetaXG/seBetaYG) # weighted 

# Write results to dataframe
Results <- data.frame(Beta=numeric(), SE=numeric(), LCI=numeric(), UCI=numeric(), t=numeric(), p=numeric()) 
colnames(IVWResults_test) <- c("Beta", "SE", "LCI", "UCI", "t", "p") 
colnames(MREggerResults) <- c("Beta", "SE", "LCI", "UCI", "t", "p")
colnames(WMresults_test) <- c("Beta", "SE", "t", "p","LCI", "UCI")
Results <- rbind(Results, IVWResults_test) 
Results <- rbind(Results, MREggerResults)
Results <- rbind(Results, WMresults_test)
rownames(Results) <- c("IVW", "MR-Egger Intercept", "MR-Egger Slope", "Weighted Median")
Results <- format(round(Results[,1:6],4)) 
Results <- Results[,c(1,3,4,5,6,2)] 
Results_word <- Results[,c("Beta", "LCI", "UCI", "p")] # table just of beta and CI's for word

#################################################################################
# Plots of results 
#################################################################################

# FUNNEL PLOT:
plot(BetaIV,BXG/seBetaYG,pch=19,cex=1.4,main="SchizNEB",xlab="",ylab="")
lines(rep(IVWBeta,2),c(0,Yrange[2]),col="red",lwd=2)
#lines(rep(MREggerBeta1,2),c(0,Yrange[2]),col="blue",lwd=2) (took out Mr-Egger as we did not include it in the paper)
legend("topright",c("IVW"),lwd=2,
       col=c("red"),cex=1.6,bty="n",pch=c(15,22))
mtext(side=1,expression(paste("Causal estimate  ", hat(beta)[j])),line=3.5,cex=1.5)
mtext(side=2,expression(paste("Instrument strength   ",hat(gamma)[j]^"C")),line=1.5,cex=1.5)

points(IVWBeta,0.2*Yrange[2],cex=2,col="red",pch=15)
lines(IVW_CI,rep(0.2*Yrange[2],2),lwd=2,col="red")
#points(MREggerBeta1,0.1*Yrange[2],cex=2,col="blue",pch=22)
#lines(MRBeta1_CI,rep(0.1*Yrange[2],2),lwd=2,col="blue")

# PLOT
plot(BXG,BYG,pch=19,cex=1.4,main="SchizNEB",xlab="",ylab="")
mtext(side=2,expression(paste("Gene-outcome coeff:  ", hat(Gamma)[j])),line=2,cex=1.5)
mtext(side=1,expression(paste("Gene-exposure coeff:   ",hat(gamma)[j])),line=2.8,cex=1.5)

MREggerFit      = lm(BYG ~ BXG,weights=1/seBetaYG^2)
lines(BXG,MREggerFit$fitted.values,col="blue")

IVWFit = lm(BYG ~ -1+BXG,weights=1/seBetaYG^2)
lines(BXG,IVWFit$fitted.values,col="red")

legend("topleft",c("IVW","MR-Egger"),lwd=3,
       col=c("red","blue"),cex=1.6,bty="n")


#################################################################################
# Simex correction MR-Egger for Age at first sex  
#################################################################################
#load package
#library(simex)

# MR-Egger regression (weighted) 
#Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
#Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
#mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
#mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
#mod1<-summary(mod.sim1)
#mod2<-summary(mod.sim2)

# Plot the simex extrapolation
# for weighted change to modsim1 and unweighted modsim2..
#l = mod.sim2$SIMEX.estimates[,1]+1 
#b = mod.sim2$SIMEX.estimates[,3] 
#plot(l[-1],b[-1],ylab="",xlab="",pch=19,ylim=range(b),xlim=range(l)) 
#mtext(side=2,"Causal estimate",line=2.5,cex=1.5) 
#mtext(side=1,expression(1+lambda),line=2.5,cex=1.5) 
#points(c(1,1),rep(Fit2$coef[2],2),cex=2,col="blue",pch=19) 
#points(c(0,0),rep((mod.sim2$coef[2]),2),cex=2,col="blue",pch=3) 
#legend("bottomleft",c("Naive MR-Egger","MR-Egger (SIMEX)"), pch = c(19,3),cex=1.5,bty="n",col=c("blue","blue")) 
#lsq = l^2; f = lm(b~l+lsq) 
#lines(l,f$fitted)

#################################################################################
# Meta analysis of outcome summary statistics for age at first sex
#################################################################################

#meta-analyse the gene-outcome effects
#library(meta)
#meta.results <- metagen(data$BYG, data$seBetaYG, comb.fixed=T, sm="Beta")
#res<-summary(meta.results)
#res

#################################################################################
# RADIAL MR for age at first sex
#################################################################################

#install.packages("devtools")
#install_github("WSpiller/RadialMR")
#library(devtools)
#library(RadialMR)

# format data for radial analysis
#RSID <- as.vector(data$SNP)
#radial_data <- format_radial(BXG,BYG,seBXG,seBYG,RSID)

# Ivw radial
#IVW <- ivw_radial(radial_data, alpha=0.01, summary=TRUE)
#IVW

# Egger radial
#egger <- egger_radial(radial_data, alpha=0.01, summary=TRUE)
#egger

# Plot radial
#plot_radial((c(IVW, egger)), radial_scale=TRUE, show_outliers=TRUE, scale_match=TRUE)

#################################################################################
# Leave one out analysis for age at first sex
#################################################################################

#format the exposure dataset like so...
# $ id.exposure           : chr "1009"
# $ SNP                   : chr "rs2075677"
# $ effect_allele.exposure: chr "A"
# $ other_allele.exposure : chr "G"
# $ eaf.exposure          : num 0.774
# $ beta.exposure         : num 0.021
# $ se.exposure           : num 0.004
# $ pval.exposure         : num 1.88e-08
# $ samplesize.exposure   : int 298420
# $ ncase.exposure        : int NA
# $ ncontrol.exposure     : int NA
# $ units.exposure        : chr "SD"
# $ exposure              : chr "Subjective well being || id:1009"
# $ pval_origin.exposure  : chr "reported"
# $ data_source.exposure  : chr "mrbase"
# $ mr_keep.exposure      : logi TRUE

# making vars i dont have exposure
effect_allele.exposure <- "A"
effect_allele.outcome <- "A"
other_allele.exposure <- "A"
eaf.exposure <- NA
beta.exposure <- BXG
se.exposure <- seBetaXG
exposure <- "exposure"
SNP <- as.vector(data$SNP)
id.exposure<-NA
samplesize.exposure<-NA 
ncase.exposure<- NA
ncontrol.exposure<-NA
units.exposure<-"UNIT"
pval_origin.exposure<-"reported"
data_source.exposure<-"mydata"
mr_keep.exposure<-TRUE
pval.exposure <-NA

# making vars i dont have outcome
beta.outcome <- BYG
se.outcome <- seBetaYG
effect_allele.outcome <- "A"
effect_allele.outcome <- "A"
other_allele.outcome <- "A"
eaf.outcome <- NA
outcome <- "outcome"
SNP <- as.vector(data$SNP)
id.outcome<-NA
samplesize.outcome<-NA 
ncase.outcome<- NA
ncontrol.outcome<-NA
units.outcome<-"UNIT"
pval_origin.outcome<-"reported"
data_source.outcome<-"mydata"
mr_keep.outcome<-TRUE
pval.outcome <-NA

# forming exposure data
exposure_dat <- data.frame(id.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, ncase.exposure, ncontrol.exposure, units.exposure, exposure, pval_origin.exposure, data_source.exposure, mr_keep.exposure)
exposure_dat$SNP <- as.vector(SNP)
exposure_dat$effect_allele.exposure <- as.character(effect_allele.exposure)
exposure_dat$other_allele.exposure <- as.character(other_allele.exposure)
exposure_dat$units.exposure <- as.character(units.exposure)
exposure_dat$exposure <- as.character(exposure)
exposure_dat$pval_origin.exposure <- as.character(pval_origin.exposure)
exposure_dat$data_source.exposure <- as.character(data_source.exposure)

# forming outcome data
outcome_dat <- data.frame(id.outcome, SNP, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome, ncase.outcome, ncontrol.outcome, units.outcome, outcome, pval_origin.outcome, data_source.outcome, mr_keep.outcome)
outcome_dat$SNP <- as.vector(SNP)
outcome_dat$effect_allele.outcome <- as.character(effect_allele.outcome)
outcome_dat$other_allele.outcome <- as.character(other_allele.outcome)
outcome_dat$units.outcome <- as.character(units.outcome)
outcome_dat$outcome <- as.character(outcome)
outcome_dat$pval_origin.outcome <- as.character(pval_origin.outcome)
outcome_dat$data_source.outcome <- as.character(data_source.outcome)


dat <- harmonise_data(exposure_dat, outcome_dat, action=1)

#Performs series of MR methods on dataset
#res<-mr(dat)

#Produces scatter plot showing estimates from multiple methods
#plot1<-mr_scatter_plot(res, dat)

#Produces forest plot of individual ratio estimates
#res_single <- mr_singlesnp(dat)
#plot2<- mr_forest_plot(res_single)
#plot2

#Produces plot showing leave-one-out estimates
res_loo <- mr_leaveoneout(dat)
plot3 <- mr_leaveoneout_plot(res_loo)
plot3

#Produces funnel plot for sensitivity analyses
plot4 <- mr_funnel_plot(res_single)
plot4

#Produces an HTML report showing results of MR analysis
mr_report(dat)



