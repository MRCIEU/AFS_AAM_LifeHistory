
#################################################################################
# Re-coding SNPs to use in analysis
#################################################################################
# "create_recoded_snps" script not always counting the effect allele
# I copied the column titles from the .raw file into exel and transposing, to compare the counted allele and effect allele using the IF command
# If I got a 2 for the if command - they did not match
# Note: this is very hard coded for column numbers

rm(list=ls())

# Convert the output of the "create_recoded_snps" shell script (.raw file) into txt and read in
dat <- read.table("snps_recoded.txt", header=T)

# Read in excel file where you have identified mismatches as 2 in column 6
mismatched <- read.csv("mismatched.csv", header=T)

# because I had separated the SNP name and its counted allele into separate cols for IF command, this puts them back together to match the .raw/.txt file again
mismatched$snpstomatch<- paste(mismatched[,1], "_", mismatched[,2], sep = "") 
mismatched2<-mismatched[which(mismatched[,6]==2),]

# Have a look at which snps were mismatches
head(data[,which(colnames(dat)%in%mismatched2[,9])]) 
matched<- 2-dat[,which(colnames(dat)%in%mismatched2[,9])] 

# check if the replication is now counting effect allele not the minor allele by comparing in excel
str(dat)
str(matched)

# replace replication data mismatched columns with new reverse coded ones
dat2 <- subset(dat2[,-which(colnames(dat)%in%mismatched2[,9])]) 
dat2 <- lapply(dat2[,1:length(dat2)], as.numeric) 

# stick matched and replication2 (the subsetted replicationdata without the mismatched SNPS) back together
genetic <- data.frame(dat2, matched)

# remove pallindromic SNPs
mismatched[nrow(mismatched)+1,] <- c(rep(NA,8), "IID") 
mismatched[nrow(mismatched)+1,] <- c(rep(NA,8), "FID") 
tail(mismatched)
genetic <- subset(genetic[,which(colnames(genetic)%in%mismatched[,9])])

#################################################################################
# Merging genetic and observational
#################################################################################

# read in observational from previous script
famfinalEX <- read.csv("/data_ID_excluded.csv", header=T)

# merge observational and genetic
finaldata <- merge(famfinalEX, genetic, by.famfinalEX=c("IID"), famfinalEX=T)

# write to file - pick which line to run depending exposure
write.csv(finaldata, file="finaldata_AFS.csv", row.names=F, quote=F)
write.csv(finaldata, file="finaldata_AAM123.csv", row.names=F, quote=F)
write.csv(finaldata, file="finaldata_AAM389.csv", row.names=F, quote=F)
