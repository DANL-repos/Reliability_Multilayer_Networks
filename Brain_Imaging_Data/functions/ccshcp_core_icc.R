# Calcuate ICC for metrics.
# args[1]: The input matlab .mat file. Has: 'AllData','CovW','CovB','time','sID','nSubj','nSession','nNodes'
# args[2]: The name for the output matlab .mat file.
# Original written by Ting Xu (xutingxt@gmail.com) 2015.12.16, revised to save out statistics of covariances 2017.02.14
# Revised by Chao-Gan Yan (ycg.yan@gmail.com) to interface with y_ICC_Image_LMM_CallR.m 2016.11.15
# Revised by Xi-Nian Zuo (zuoxinian@gmail.com) to handle surface-based metrics 2017.01.08
# Key Laboratory of Behavioral Science and Magnetic Resonance Imaging Research Center
# Institute of Psychology, Chinese Academy of Sciences, Beijing, China
# 
# ## ----------------------------- Model Instruction ----------------------------- ##
# Intra-Class Correlation (ICC) Calculation with Linear Mixed Model (LME)
# Author: Ting Xu (xutingxt@gmail.com), Date: Dec 27, 2016
# Chinese Academy of Sciences, China; Child Mind Institute, NY USA
# ## ----------------------------------------------------------------------------- ##
# This function performs the Linear Mixed Model (LME) to calculate ICC, the within- and 
# between- subjects variability, and the corresponding group averaged map across all input 
# images. It's required to install R (https://www.r-project.org/). Alternatively you may 
# install the following R packages:
# 
# install.packages("nlme")
# install.packages("R.matlab")
# 
# More details about the LME in 'nlme' R package can be found at
# https://cran.r-project.org/web/packages/nlme/index.html
# 
# The basic linear mixed model are built to estimate ICC with ReML (Restricted Maximum 
# Likelihood) estimation. The Linear Mixed Model allows missing data included. In addition, 
# ReML method avoids the negative ICC estimation. Random effects are set for the intercept 
# in the model. 
# 
# The basic model without covariates is as follow:
# 
# Y(i,j) = mu + gamma(i) + e(i,j)
# 
# e(i,j) ~ i.i.d. N(0, delta-square)
# gamma(i) ~ i.i.d N(0, tao-square)
# 
#         
# Y(i,j) is the dependent variable (e.g., functional connectivity, ALFF, ReHo, etc.) for the 
# i-th session in j-th subject. The random effects (gamma) and the residual term (e) are
# assumed to be independent identical normally distributed (i.i.d) with mean zero and 
# variance delta-square and tao-square. ICC is estimated as the follow:
# 
#                    between-individual variation  
# ICC =  -------------------------------------------------------------
#        (between-individual variation + within-individual variation)
#        
#                 tao-square
#     =  -----------------------------
#         (tao-square + delta-square)
#         
# 
# In addition, the model can incorporate the within- and between- individual covariates. 
# When the within-individual covariates are included, the random effect are also set for 
# the within-individual covariates but constrains no correlations between any of the random 
# effects or residuals.
# 
# Input: .mat
# sID    : one vector, the subject ID for all the observations. Make sure the order of observations below (AllData, CovW, CovB) is the same as sID. 
# AllData: the number of observations X the number of Nodes.
# CovW   : the number of observations X the number of within individual covariances, e.g. age, sex (see the Notice below for more information)
# CovB   : the number of observations X the number of between indiviudal covariances, e.g. meanFD, global-metric  
# time   : one vector, the scan time label for the observations (e.g. 1 for the first session, 2 for the second session). 
#          This is not used in the model, but it is good to be written out to label the input data 
# nSubj  : the number of subjects
# nNodes : the number of Nodes (e.g., vertex on surface, voxel in volume)
# nSession: the number of sessions
#
# Output: .mat
# intercept_metric: nNodes X 1, The intercept of the mixed model.
#                   If the covariances are setup appropriately (or no covariances), this would be the estimation of average map across the samples.
# t_value_intercept, p_value_intercept: nNodes x 1, the t and p statistics of the intercept, test whether the intercept is different than zero.
# t_value_CovB, p_value_CovB: nNodes X the number of between individual covariances
#                             The t and p statistics of the beta for between individual covariances, test whether the beta is different than zeros.
#                             If more than two covariances are included in the model, the order is the same as the inputs.
# t_value_CovW, p_value_CovW: nNodes X the number of within individual covariances
#                             The t and p statistics of the beta for within individual covariances, test whether the beta is different than zeros.
#                             If more than two covariances are included in the model, the order is the same as the inputs.
#
# 
# Notice: 
# 1) The demean step should and is applied outside this function in 'y_ICC_Image_LMM_CalR.m'.
# Usually, the grand-mean centering is required for quantitative variables. 
# 
# 2) When including categorical variables as covariates (for instances, sex), make sure 
# dummy coding the variables (zeros and ones). 
# 
# 3) Be cautious with covariates, in most of the cases, the demean or dummy coding steps 
# does not affect the ICC estimation. However, the intercept (estimated average map) 
# and its interpretation is dependent on how you deal with covariates (demean, dummy coding).
# References:
# [1] Zuo XN, Xu T, Jiang L, Yang Z, Cao XY, He Y, Zang YF, Castellanos FX, Milham MP.
#     2013. Toward reliable characterization of functional homogeneity in the human brain: 
#     preprocessing, scan duration, imaging resolution and computational space. 
#     Neuroimage 65, 374-386.
# [2] Xu T, Yang Z, Jiang L, Xing XX, Zuo XN. 2015. A Connectome Computation System for 
#     discovery science of brain. Science Bulletin 60:86-95. 
# ==========================================================================================


library("nlme")
library("R.matlab")

args = commandArgs(trailingOnly=TRUE)
InputName <- args[1]
OutputName <- args[2]

mat <- readMat(InputName)

nsubj <- mat[["nSub"]]
nsession <- mat[["nSes"]]
nNodes <- mat[["nROI"]]
# load data
metric <- mat[["AllData"]]
# load nuisance
CovW <- mat[["CovW"]]
CovB <- mat[["CovB"]]
subid <- mat[['sID']]
xvisit <- mat[['time']]

# initial setup
icc = matrix(0, nNodes, 1)
t_value_intercept = matrix(0, nNodes, 1)
p_value_intercept= matrix(0, nNodes, 1)
intercept_metric = matrix(0, nNodes, 1)
varb = matrix(0, nNodes, 1)
varw = matrix(0, nNodes, 1)
varb_rate = matrix(0, nNodes, 1)
varw_rate = matrix(0, nNodes, 1)
t_value_CovB = matrix(0, nNodes, ncol(CovB))
t_value_CovW = matrix(0, nNodes, ncol(CovW))
p_value_CovB = matrix(0, nNodes, ncol(CovB))
p_value_CovW = matrix(0, nNodes, ncol(CovW))

for (n in 1:nNodes){
  if (n%%100 == 0) print(sprintf("Compute mixed model for %4.2f present ...", n/nNodes*100))
  y = metric[,n]
  dataframe = data.frame(y,subid)
  if (length(CovB)!=0) dataframe = data.frame(dataframe,CovB)
  if (length(CovW)!=0) dataframe = data.frame(dataframe,CovW)
  ColumnNames <-names(dataframe)
  CovB_names <- ColumnNames[3:(2+ncol(CovB))]
  CovW_names <- ColumnNames[(3+ncol(CovB)): length(ColumnNames)]
  
  lmeExpression <- "y ~ "
  if (length(ColumnNames)>3) {
    for (iColumn in 3:(length(ColumnNames)-1)){
      lmeExpression <- sprintf('%s%s+', lmeExpression,ColumnNames[iColumn])
    }
  }
  if (length(ColumnNames)>2) {
    lmeExpression <- sprintf('%s%s', lmeExpression,ColumnNames[length(ColumnNames)])
  }
  if (length(ColumnNames)==2) {
    lmeExpression <- sprintf('%s1', lmeExpression)
  }
  if (length(CovW)==0) {
    lmeExpression <- sprintf('fm <- lme(%s, random = ~ 1 |subid, data = dataframe)', lmeExpression)
  } else {
    WithinCovName <- sprintf('%s', ColumnNames[(length(ColumnNames)-ncol(CovW)+1)])
    if (ncol(CovW)>1) {
      for (iColumn in (length(ColumnNames)-ncol(CovW)+2):(length(ColumnNames))){
        WithinCovName <- sprintf('%s+%s', WithinCovName, ColumnNames[iColumn])
      }
    }
    lmeExpression <- sprintf('fm <- lme(%s, random = list(subid = pdDiag(~ %s)), data = dataframe)', lmeExpression, WithinCovName)
  } 
  
  try({eval(parse(text=lmeExpression))
  output <- summary(fm)
  sigma_r = output$sigma^2
  sigma_b = getVarCov(fm)[1,1]
  icc[n] = sigma_b/(sigma_r+sigma_b)
  t_value_intercept[n] = output$tTable['(Intercept)', 't-value']
  p_value_intercept[n] = output$tTable['(Intercept)', 'p-value']
  intercept_metric[n] = output$tTable['(Intercept)', 'Value']
  tmp_idx = rownames(output$tTable) %in% CovW_names
  t_value_CovW[n,] = output$tTable[tmp_idx, 't-value']
  p_value_CovW[n,] = output$tTable[tmp_idx, 'p-value']
  tmp_idx = rownames(output$tTable) %in% CovB_names
  t_value_CovB[n,] = output$tTable[tmp_idx, 't-value']
  p_value_CovB[n,] = output$tTable[tmp_idx, 'p-value']
  
  
  varb[n] = sigma_b
  varw[n] = sigma_r
  vary = var(y)
  varb_rate[n] = sigma_b / vary
  varw_rate[n] = sigma_r / vary
  
  }, silent = FALSE)
}

writeMat(OutputName, icc = icc, varb = varb, varw = varw, varb_rate = varb_rate, varw_rate = varw_rate, intercept_metric = intercept_metric, 
         t_value_intercept = t_value_intercept, t_value_CovB = t_value_CovB, t_value_CovW = t_value_CovW, 
         p_value_intercept = p_value_intercept, p_value_CovB = p_value_CovB, p_value_CovW = p_value_CovW)



