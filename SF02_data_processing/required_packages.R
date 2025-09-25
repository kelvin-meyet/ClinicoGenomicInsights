if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("cBioPortalData", quietly = TRUE))
  install.packages("cBioPortalData")

BiocManager::install("cBioPortalData")
require(cBioPortalData)
require(BiocManager)
require(AnVIL)


install.packages("mice",'GGally','grpreg', dependencies = T)

library(readr)
library(janitor)
library(dplyr)
library(tidyr)
require(DataExplorer)
require(ggplot2)
require(ggthemes)
require(DataExplorer)
library(mice)
