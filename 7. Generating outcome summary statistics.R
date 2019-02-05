# Note: sometimes have to exit R inbetween loops
# Note: this is hard coded for column numbers - adjust read in data and column no.s for varying exposures

rm(list=ls())

#################################################################################
# Data prep
#################################################################################

# To use 500k or 350k
finaldataAll <- read.csv("finaldata_AAM389.csv", header=T) 
#finaldata350 <- subset(finaldataAll, finaldataAll$split==2) 
#finaldataAll <- finaldata350
dim(finaldataAll)

# Determining which columns to use for SNPs
str(finaldataAll) 
names(finaldataAll[48:352]) # check names of all are SNPs
numberofSNPs <- finaldataAll[48:352]
names(finaldataAll)[48:352] <- sub("_.*", "", names(finaldataAll[,48:352]))

library(gdata)
keep(finaldataAll)
keep(finaldataAll, sure=T)
ls()
attach(finaldataAll) 

#################################################################################
# Summary stats NEB 
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$num_kids ~ finaldataAll[,i] +  + birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}
SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) # create exposure file from mismatched and GWAS betas and se
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_NEB_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats AFB
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$age_fbirth_multi ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_AFB_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats ALB 
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$age_lbirth ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_ALB_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats birthrange
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$birthrange_m ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_birthrange_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats years of education
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$edu_years ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_edu_years_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats years of age left education
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$age_edu ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_age_edu_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats alcohol intake
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$alc_intake ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_alc_intake_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats num_partners
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- lm(finaldataAll$num_partners ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_num_partners_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats Childlessness 
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- glm(finaldataAll$childless ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial(link = logit), data = finaldataAll) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_childless_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")

#################################################################################
# Summary stats risk 
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- glm(finaldataAll$risk ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial(link = logit), data = finaldataAll) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_risk_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")


#################################################################################
# Summary stats ever_smok 
#################################################################################

beta_y <- c()
SE_y <- c()
df_y <- c()

for (i in 48:352) { 
  model <- glm(finaldataAll$ever_smok ~ finaldataAll[,i]  +birthyear + PC1 + PC2 + PC3 +PC4 +PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial(link = logit), data = finaldataAll) 
  beta <- model$coefficients[[2]] 
  se <- coef(summary(model))[2, 2] 
  df <- model$df.residual
  c(beta_y, beta)-> beta_y 
  c(SE_y, se) -> SE_y
  c(df_y, df) -> df_y
}

SNP <- names(finaldataAll)[48:352]

sum_stats <- data.frame(cbind(SNP, beta_y, SE_y, df_y))
sum_stats$SNP <- as.character(sum_stats$SNP)
sum_stats$beta_y <- as.character(sum_stats$beta_y)
sum_stats$SE_y <- as.character(sum_stats$SE_y)
sum_stats$df_y <- as.character(sum_stats$df_y)

exposure <- read.csv(file="exposure.csv", header=T) 
exposure$SNP <- as.character(exposure$SNP)

dataset <- merge(exposure, sum_stats, by.sum_stats=c("SNP"), sum_stats=T)
write.csv(file="dataset_ever_smok_350_AAM389_withdrawals.csv", dataset, quote=F, row.names=F)

quit("no")