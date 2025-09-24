# install.packages("remotes")
# install.packages("glmnet")
# install.packages("gbm")
# install.packages("pec")
# remotes::install_github("dushoff/shellpipes", ref = "main", force = TRUE)
# remotes::install_github("CYGUBICKO/glmnetpostsurv", dependencies = T)
# remotes::install_github("CYGUBICKO/satpred", dependencies = TRUE)



require(dplyr)
require(ggplot2)
require(caTools)
require(glmnet)
library(survival)
library(caret)
library(caret)

library(gbm)
library(pec)
library(remotes)
library(reticulate)
# install_pycox(pip = TRUE, install_torch = TRUE)
# install_keras(pip = TRUE, install_tensorflow = TRUE)
library(shellpipes)
library(glmnetsurv)
library(survivalmodels)
library(tensorflow)
library(keras)
library(satpred)


satpredtheme()

# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("survcomp")

require(survcomp)


