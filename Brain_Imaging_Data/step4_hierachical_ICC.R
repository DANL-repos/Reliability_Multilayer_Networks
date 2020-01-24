rm(list = ls())

library(nlme)

setwd('./example_data/hierarchical_ICC')

measure<-'flexibility'
scan<-'60min'
nROI<-200

# Input dynamic measure data matrix
inputFileName<-paste0(measure,'_', scan, '_gamma1.05_omega2.5_for_hierarchical_ICC.txt')

# Label corresponding to the data matrix
dataLabel<-read.csv('labels_hierarchical_ICC.csv')
dataLabel$session<-factor(dataLabel$session)

# Define output file name
outName_ICC_session<-paste0(measure,'_', scan, '_gamma1.05_omega2.5_ICC_session.txt')
outName_ICC_condition<-paste0(measure,'_', scan, '_gamma1.05_omega2.5_ICC_condition.txt')

#################################################################
dynamicData<-read.table(inputFileName, header=FALSE, row.names=NULL)
dataMatrix<-data.matrix(dynamicData)

icc_session_all_ROIs <- vector(mode="numeric", length=nROI)
icc_condition_all_ROIs<- vector(mode="numeric", length=nROI)
for (i in 1:nROI) {
  dataCombined<-cbind(dataLabel, dataMatrix[, i])
  colnames(dataCombined)[4] <- "y"
fm <- lme(y ~ 1, random=~1|subID/condition, data=dataCombined)
output <- summary(fm)
lme_var <- VarCorr(output)
sigma_sub <- as.numeric(lme_var[2, 'Variance'])
sigma_con <- as.numeric(lme_var[4, 'Variance'])
sigma_r <- as.numeric(lme_var['Residual','Variance'])
                                
sigma_all = sigma_sub + sigma_con + sigma_r
icc_session <- (sigma_sub + sigma_con)/sigma_all
icc_condition <- sigma_sub/(sigma_sub + sigma_con)

icc_session_all_ROIs[i]<-icc_session
icc_condition_all_ROIs[i]<-icc_condition

}

write.table(icc_condition_all_ROIs, file=outName_ICC_condition, col.names=FALSE, row.names=FALSE)
write.table(icc_session_all_ROIs, file=outName_ICC_session, col.names=FALSE, row.names=FALSE)

