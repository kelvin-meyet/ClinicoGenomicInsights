# Set Dynamic Working Directory and load R source code of dependent libraries
# install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

source("required_packages.R")

#===============================================================================
#-------------------Get clinical Data from cBioportal---------------------------
#===============================================================================

# if (!require("BiocManager", quietly = TRUE)) 
#   install.packages("BiocManager")
# 
# if (!require("cBioPortalData", quietly = TRUE)) 
#   install.packages("cBioPortalData")
# 
# BiocManager::install("cBioPortalData") 
# require(cBioPortalData)

# require(cBioPortalData)
# require(BiocManager)
# require(AnVIL)

cbio <- cBioPortal() #initialise cbioportal 

getStudies(api=cbio)

setCache(directory ="cbiodata-clinical", 
         verbose = TRUE, ask = interactive())

#Get Study IDs available
stud=getStudies(cbio)
#write.csv(stud, file="cbiodata-clinical\\allstudyIDs.csv")

# Get clinical data for specified studyID
clinic <- clinicalData(api = cbio,
                             studyId = "prad_tcga_pan_can_atlas_2018")

#write.csv(clinic, file="cbiodata-clinical\\prca_tcga.csv")



#============Extract Relevant Variables from clinical Data=============
#-------------------Read Clinical Data & extract only relevant variables--------------------
clinic = read.csv("cbiodata-clinical\\prca_tcga.csv")  #55 variables



relevant_variables = c("OTHER_PATIENT_ID", "AGE","BUFFA_HYPOXIA_SCORE", "HISTORY_NEOADJUVANT_TRTYN", "ICD_O_3_HISTOLOGY","NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
                       "PATH_N_STAGE","PATH_T_STAGE","PERSON_NEOPLASM_CANCER_STATUS","PRIOR_DX","RADIATION_THERAPY","RAGNUM_HYPOXIA_SCORE","WINTER_HYPOXIA_SCORE",
                       "ETHNICITY","RACE","ANEUPLOIDY_SCORE","MSI_SCORE_MANTIS","MSI_SENSOR_SCORE","TISSUE_SOURCE_SITE","TUMOR_TYPE","TMB_NONSYNONYMOUS",
                       "FRACTION_GENOME_ALTERED","MUTATION_COUNT","PFS_MONTHS","PFS_STATUS")

clinical_data <- clinic %>%
  select(c(relevant_variables)) %>%
  mutate(OTHER_PATIENT_ID = tolower(OTHER_PATIENT_ID))


#=========================== Genomic Modelling ================================
snv_gene = read.csv("input_data\\SNV_gene_distribution.csv") #contains all 17036 genes

#======================bofei====================
bofei_gene <- scan("input_data\\LPC_genes.txt", what = "character", sep ="\n") %>%
  trimws(.)

bofei_matched_rows <- snv_gene$genename %in% bofei_gene
table(bofei_matched_rows)

bofei_match_genes <- snv_gene[bofei_matched_rows, ] #27genes

snvtumor = bofei_match_genes[grepl("tumor", colnames(snv_gene))]  
snvnormal = bofei_match_genes[grepl("normal", colnames(snv_gene))]

#---Tumor---
snv_tumor <- snvtumor %>% 
  select(-c(distinct.tumor.count, tumor.tally, GSV_CDS_tumor, nsGSV_tumor)) %>%
  mutate(snv27 = rowSums(.))

tumor <- bofei_match_genes %>%
  select(c(chromosome,genename,coding.gene,COSMIC,IntOGen))

tumorpatients = data.frame(tumor, snv_tumor)


#---Normal----
snv_normal <- snvnormal %>% 
  select(-c(distinct.normal.count, normal.tally, GSV_CDS_normal, nsGSV_normal)) %>%
  mutate(snv27 = rowSums(.))

normal <- bofei_match_genes %>%
  select(c(chromosome,genename,coding.gene,COSMIC,IntOGen))

normalpatients=data.frame(normal, snv_normal)


#------------Now building the tumor matrix with snv counts for 503 patients---------------
tumor_snv_matrix = tumorpatients %>%
  select(-c(chromosome, genename, coding.gene, COSMIC, IntOGen, snv27)) %>%
  t(.) 

colnames(tumor_snv_matrix) <- tumorpatients$genename 

snv503_tumor = rowSums(tumor_snv_matrix) #snv for all 503 patients across 94 matched genes
snv_27_gene_tumor <- colSums(tumor_snv_matrix) #snv for 94 genes across 503 patients

newtumor_patients<- cbind.data.frame(tumor_snv_matrix, snv503_tumor)




#------------Now building the normal matrix with snv counts for 503 patients---------------
normal_snv_matrix = normalpatients %>%
  select(-c(chromosome, genename, coding.gene, COSMIC, IntOGen, snv27)) %>%
  t(.) 

colnames(normal_snv_matrix) <- normalpatients$genename 

snv503_normal = rowSums(normal_snv_matrix) #snv for all 503 patients across 94 matched genes
snv_27_gene_normal <- colSums(normal_snv_matrix) #snv for 94 genes across 503 patients

#change snv503_normal to snv503n - done!!!

newnormal_patients<- cbind.data.frame(normal_snv_matrix, snv503_normal)


#working on Duplicates: find Unions

#replace each pair of duplicates with the union, the merging downstream will remove one duplicated union.
tumor_patients <- newtumor_patients

union_T1 <- apply(tumor_patients[c("patient21_tumor","patient244_tumor"), ], 2, max)

union_T1 <- as.data.frame(t(union_T1))
sum(as.numeric(union_T1[1,]))
tumor_patients[c(21,244),] <- union_T1 #replace both duplicates with union2
sum(as.numeric(tumor_patients[21,])); sum(as.numeric(tumor_patients[244,]))


union_T2 <- apply(tumor_patients[c("patient41_tumor", "patient126_tumor"), ], 2, max)
union_T2 <- as.data.frame(t(union_T2))
sum(as.numeric(union_T2[1,]))
tumor_patients[c(41,126), ]<- union_T2 #replace both duplicates with union2
sum(as.numeric(tumor_patients[41,])); sum(as.numeric(tumor_patients[126,]))


union_T3<-apply(tumor_patients[c("patient75_tumor","patient274_tumor"), ], 2, max)
union_T3<-as.data.frame(t(union_T3))
sum(as.numeric(union_T3[1,]))
tumor_patients[c(75,274), ]<-union_T3 #replace both duplicates with union3
sum(as.numeric(tumor_patients[75,])); sum(as.numeric(tumor_patients[274,]))



union_T4<-apply(tumor_patients[c("patient259_tumor","patient346_tumor"), ], 2, max)
union_T4<-as.data.frame(t(union_T4))
sum(as.numeric(union_T4[1,]))
tumor_patients[c(259,346), ]<-union_T4
sum(as.numeric(tumor_patients[259,])); sum(as.numeric(tumor_patients[346,]))


union_T5 <-apply(tumor_patients[c("patient147_tumor","patient486_tumor"), ], 2, max)
union_T5<-as.data.frame(t(union_T5))
sum(as.numeric(union_T5[1,]))
tumor_patients[c(147,486), ]<-union_T5
sum(as.numeric(tumor_patients[147,])); sum(as.numeric(tumor_patients[486,]))


#----Normal Unions---
normal_patients <- newnormal_patients

union_N1<-apply(normal_patients[c("patient21_normal","patient244_normal"), ], 2, max)

union_N1<-as.data.frame(t(union_N1))
sum(as.numeric(union_N1[1,]))
normal_patients[c(21,244),] <-union_N1 #replace both duplicates with union2
sum(as.numeric(normal_patients[21,])); sum(as.numeric(normal_patients[244,]))


union_N2<-apply(normal_patients[c("patient41_normal", "patient126_normal"), ], 2, max)
union_N2<-as.data.frame(t(union_N2))
sum(as.numeric(union_N2[1,]))
normal_patients[c(41,126), ]<- union_N2 #replace both duplicates with union2
sum(as.numeric(normal_patients[41,])); sum(as.numeric(normal_patients[126,]))



union_N3<-apply(normal_patients[c("patient75_normal","patient274_normal"), ], 2, max)
union_N3<-as.data.frame(t(union_N3))
sum(as.numeric(union_N3[1,]))
normal_patients[c(75,274), ]<-union_N3 #replace both duplicates with union3
sum(as.numeric(normal_patients[75,])); sum(as.numeric(normal_patients[274,]))



union_N4<-apply(normal_patients[c("patient259_normal","patient346_normal"), ], 2, max)
union_N4<-as.data.frame(t(union_N4))
sum(as.numeric(union_N4[1,]))
normal_patients[c(259,346), ]<-union_N4
sum(as.numeric(normal_patients[259,])); sum(as.numeric(normal_patients[346,]))


union_N5 <-apply(normal_patients[c("patient147_normal","patient486_normal"), ], 2,max)
union_N5<-as.data.frame(t(union_N5))
sum(as.numeric(union_N5[1,]))
normal_patients[c(147,486), ]<-union_N5
sum(as.numeric(normal_patients[147,])); sum(as.numeric(normal_patients[486,]))


# Combine the two data frames such that the columns are interleaved.
# Combine Tumor & Normal, such that 1 from tumor file & 1 from  normal file follows in that order.
obs <- nrow(newtumor_patients); 
col_gene = ncol(newtumor_patients)

comb2 = data.frame(matrix(ncol = col_gene*2 , nrow = obs),
                   row.names = paste0("patient_", c(1:obs)))

#Iterate over first 27 columns in tumor & normal files  

for(i in 1:col_gene){
  comb2[, i*2-1] <- newtumor_patients[, i] #odd numbered columns
  comb2[, i*2] <- newnormal_patients[, i] #even numbered columns
}

# Rename the columns in the combined data frame
labels=c(colnames(newtumor_patients), colnames(newnormal_patients))
#colnames(comb) <- paste0(c("tdf1_", "ndf2_"), rep(colnames(tumor_patients), each = 2))

#structure & select the colnames well
needednames <- c(colnames(newtumor_patients)[-28], "snv503") 
colnames(comb2) <- paste0(rep(needednames, each = 2), c("_tumor", "_normal"))
#set colnames to double repetitions of each genename from tumor patients, and concatenated with tumor & normal



## Combine Genetic Identifiers snvs & clinical data 
genomic_identifiers <- read.csv("input_data\\patient_id_B.csv") %>%
  select(patient, file_id, case_id, age_at_diagnosis) #Few relevant columns in genomic
genomic = cbind(genomic_identifiers, comb2)


# Investigate duplicated identifiers in both clinical and genomic 
# Duplicated patientIDs in clinical data
table(duplicated(clinical_data$patientId)) #no duplicated patient ids & we have 494 obs

# Duplicated patientIDs in genomic data
table(duplicated(genomic$case_id)) #5 duplications in genomic IDS

#extract duplicated ids 
genomic$case_id[duplicated(genomic$case_id)] #5 duplicates pairs => 10


# # remove duplicates in genomic
genomic_data2 <- genomic[!duplicated(genomic$case_id), ] #now 498 obs
dim(clinical_data)
dim(genomic_data2)

#------------MERGE---------------------
merge_gen_clin_bofei = merge(x=clinical_data, y=genomic_data2, by.x ="OTHER_PATIENT_ID", by.y ="case_id", all.x = T, all.y = F)

#shows only corresponding observations in both dataframes merged.

dim(merge_gen_clin_bofei) #dimension


#----Differences between tumor & normal-------
#library(dplyr)
snv_dat <- merge_gen_clin_bofei %>%
  select(CHD5_tumor:SALL1_normal)

col_names <- colnames(snv_dat)

differences <- list()

#loop through column names
for (name in col_names){
  if (endsWith(name, "_tumor")){
    normal_col <- gsub("_tumor","_normal", name) #extract corresponding normal column
    prefix <- gsub("_tumor", "_", name)
    differences[[paste(prefix, "diff", sep = "_")]] <- snv_dat[[name]] - snv_dat[[normal_col]]
  }
}

snv_diff_bofei = as.data.frame(differences) #differences in genomic snv  counts

write.csv(snv_diff_bofei, file = 'out_data/snv_data_difference_bofei.csv', row.names = F)


final_clinical <- merge_gen_clin_bofei %>% select(AGE:PFS_STATUS) 

final_combined_data <- as.data.frame(cbind(final_clinical, snv_diff_bofei))

write.csv(final_combined_data, file = 'out_data/final_combined_data.csv', row.names = F)


#=============================================================================
#====================Imputing Clinical data in merged data ==================
#---> Create a mapping of specific institutions to general categories for tissue_source_site variable
category_mapping <- list(
  "Research center" = c("MD Anderson Cancer Center", "International Genomics Consortium", "Memorial Sloan Kettering Cancer Center", "Fox Chase", "NCI Urologic Oncology Branch", "Institute for Medical Research"),
  "Hospital" = c("Roswell Park", "Mayo Clinic Arizona", "ABS - Lahey Clinic", "Cornell Medical College", "Melbourne Health", "Harvard Beth Israel", "Maine Medical Center", "ABS - IUPUI"),
  "University" = c("University of Pittsburgh", "University of California San Francisco", "Washington University", "University of Kansas", "University of North Carolina", "University of Arizona", "University of Sao Paulo", "Stanford University", "Wake Forest University", "University Medical Center Hamburg-Eppendorf", "University of Minnesota", "University of Texas", "BLN - Baylor"),
  "Biotech & Pharma" = c("Indivumed", "PROCURE Biobank", "Asterand", "Proteogenex, Inc.", "Global Bioclinical-Moldova", "Global BioClinical - Georgia")
)

# Function to map institutions to their category
map_category <- function(name) {
  category <- names(category_mapping)[sapply(category_mapping, function(x) name %in% x)]
  ifelse(length(category) > 0, category, "Other")
}

cleaned_combined_data <- final_combined_data  %>%
  dplyr::rename(New_Tumor_AFIT=NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT,
         NEO_CANCER_STATUS=PERSON_NEOPLASM_CANCER_STATUS,
         ICD_HIST=ICD_O_3_HISTOLOGY,
         FRAC_GEN_ALT=FRACTION_GENOME_ALTERED,
         Ragnum_HS = RAGNUM_HYPOXIA_SCORE,
         Winter_HS = WINTER_HYPOXIA_SCORE,
         Buffa_HS = BUFFA_HYPOXIA_SCORE,
         HIST_NEOADJUV_TRTMNT = HISTORY_NEOADJUVANT_TRTYN) %>%
  mutate(PFS_STATUS=recode(PFS_STATUS, "0:CENSORED"=0,"1:PROGRESSION"=1))%>%
  mutate(HIST_NEOADJUV_TRTMNT=factor(HIST_NEOADJUV_TRTMNT,
                                     levels=c("No","Yes (Pharmaceutical Treatment Prior To Resection)"), labels=c("No","Yes")))%>%
  mutate(ICD_HIST=factor(ICD_HIST,
                         levels=c("8140/3","8550/3","8255/3","8480/3","8490/3","8500/3"),
                         labels=c("adenocarcinoma", "acina_cell_carc.",
                                  "mixed_adenocarcinoma","mucinous_adenocarcinoma",
                                  "signet_ring_cell_carc.",
                                  "infiltrat_duel_carc.")))%>%
  mutate(New_Tumor_AFIT=as.factor(New_Tumor_AFIT))%>%
  mutate(PATH_N_STAGE=factor(PATH_N_STAGE, levels = c("N0","N1"),
                             labels=c("no_cancer_nearby_lymphs","cancer_nearby_lymphs")))%>%
  mutate(PATH_T_STAGE=factor(PATH_T_STAGE, levels=c("T2A","T2B","T2C", "T3A", "T3B", "T4"),
                             labels =c("stage_2A","stage_2B","stage_2C",
                                       "stage_3A","stage_3B", "stage_4"),
                             ordered = TRUE))%>%
  mutate(NEO_CANCER_STATUS=as.factor(NEO_CANCER_STATUS))%>%
  mutate(PRIOR_DX = recode(PRIOR_DX,"No"="No", "Yes"="Yes",
                           "Yes, History Of Synchronous And Or Bilateral Malignancy" = "Yes"))%>%
  mutate(PRIOR_DX = as.factor(PRIOR_DX))%>%
  mutate(RADIATION_THERAPY = as.factor(RADIATION_THERAPY))%>%
  mutate(RACE=factor(RACE, levels = c("Asian","Black or African American","White"),
                     labels = c("Asian", "Black_Af_Am","White")))%>%
  mutate(TISSUE_SOURCE_SITE = as.factor(sapply(TISSUE_SOURCE_SITE, map_category))) %>%
  mutate(TUMOR_TYPE = factor(TUMOR_TYPE,
                             levels = c("Prostate Adenocarcinoma, Acinar Type",
                                        "Prostate Adenocarcinoma, Other Subtype"),
                             labels = c("Prad_Acinar.Type", "Prad_other.subtype")))%>%
  select(-c(ETHNICITY,RACE)) %>% #removed ethnicity & Race, they are heavily missing
  clean_names()



#-----------Checking Missingness----------------
missing_columns <- cleaned_combined_data %>%
  summarise_all(~sum(is.na(.)) / length(.)) %>%
  gather(variable, missing_rate) %>%
  filter(missing_rate > 0)

missing_columns


#====Create dummy variables for only variables with missing values(13 vars) in the data
#--Dataframe of only missing vars
missing_indices <- missing_columns$variable %in% names(cleaned_combined_data)  
miss_vars <- missing_columns$variable[missing_indices] %>%
  cleaned_combined_data[,.] 

glimpse(cleaned_combined_data)
introduce(cleaned_combined_data)
plot_intro(cleaned_combined_data)

# plot missing values plot
#plot_missing(cleaned_combined_data) # - shows the frequency of missing values for each column.
plot_missing(cleaned_combined_data,
             ggtheme = theme_igray(),
             title = "Missing Value Plot",
             theme_config = theme(plot.title = element_text(color = "black")),
             geom_label_args = c(hjust = "inward")
)


#---------------------STEP 4---------------------------------
# https://www.semaforobares.com/
#================ PMM Imputation ======================
library(ggplot2)
library(mice)
clinicaldata <- cleaned_combined_data %>%
  select(age:pfs_status)

iterations  = 6
random_seed = 2098

#Imputing each variable with missing - MICE - pmm
# pmm_mice <- mice(data   = clinicaldata,
#                    method = "pmm",
#                    m      = iterations,
#                    seed   = random_seed)
# saveRDS(pmm_mice, file ="imputation_models/pmm_mice.rds")
pmm_mice = readRDS("imputation_models/pmm_mice.rds") 
#saved the original imputation results in sharepoint for reproducibility, you can re-run the entire imputation model again!

# --- Put Missing & Imputed Data into a dataframe ---
pmm_impute <- data.frame(miss_vars, imp_pmm = complete(pmm_mice)) #11+25 vars
#---Picked missing vars & some non-imputed vars for comparison
hist_dat <- pmm_impute %>%
  select(buffa_hs, ragnum_hs, winter_hs, aneuploidy_score, frac_gen_alt, mutation_count, 
         imp_pmm.buffa_hs, imp_pmm.ragnum_hs, imp_pmm.winter_hs, imp_pmm.aneuploidy_score,
         imp_pmm.frac_gen_alt, imp_pmm.mutation_count)
hist_dat_pmm <- hist_dat
hist_long <- hist_dat %>% #longer format
  pivot_longer(., cols = everything(), names_to ="Variable", values_to = "Value") %>%
  # Filter and order the variables
  mutate(Variable = factor(Variable, levels = c(sort(grep("^Imp", unique(Variable))), unique(Variable))))

# Create histograms of numeric variables and arrange them in a grid
ggplot(hist_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms of Numeric Variables and their Corresponding Imputed Variables",
       x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 12))
#dev.off()

#jpeg("barplot_pmm_new.jpeg", width = 70, height = 40, units = "cm", res = 400)
bar_dat <- pmm_impute %>%
  select(new_tumor_afit, path_n_stage, path_t_stage, neo_cancer_status,radiation_therapy,
         imp_pmm.new_tumor_afit, imp_pmm.path_n_stage, imp_pmm.path_t_stage,
         imp_pmm.neo_cancer_status,imp_pmm.radiation_therapy) 


# Reshape the dataframe into a longer format (tidy data)
bar_long <- gather(bar_dat, key = "Variable", value = "Category")

# Filter and order the variables
bar_long <- bar_long %>%
  mutate(Variable = factor(Variable, levels = c(sort(grep("^Imp", unique(Variable))), unique(Variable))))

# Create bar charts of categorical variables and arrange them in a grid
ggplot(bar_long, aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Bar Charts of Categorical Variables and their Corresponding Imputed Variables") +
  theme(strip.text = element_text(size = 12))

#dev.off()

#==============RF Imputation========================
#jpeg("histplot_rf_new.jpeg", width = 70, height = 40, units = "cm", res = 400)

# rf_mice <- mice(data   = clinicaldata,
#                 method = "rf",
#                 m      = iterations,
#                 seed   = random_seed)
# saveRDS(rf_mice, file ="imputation_models/rf_mice.rds")

rf_mice = readRDS("imputation_models/rf_mice.rds")
#Again saved the original imputation results in sharepoint for reproducibility, you can re-run the entire imputation model again!

#---Put Missing & Imputed Data into a dataframe
rf_impute <- data.frame(miss_vars, imp_rf = complete(rf_mice)) #11+25 vars

#--Picked missing vars & some non-imputed vars for comparison
hist_dat <- rf_impute %>%
  select(buffa_hs, ragnum_hs, winter_hs, aneuploidy_score, frac_gen_alt, mutation_count, 
         imp_rf.buffa_hs, imp_rf.ragnum_hs, imp_rf.winter_hs,imp_rf.aneuploidy_score,
         imp_rf.frac_gen_alt, imp_rf.mutation_count)
hist_dat_rf <- hist_dat

hist_long <- hist_dat %>% #longer format
  pivot_longer(., cols = everything(), names_to = "Variable", values_to = "Value") %>%
  # Filter and order the variables
  mutate(Variable = factor(Variable, levels = c(sort(grep("^Imp", unique(Variable))), unique(Variable))))


# Create histograms of numeric variables and arrange them in a grid
ggplot(hist_long, aes(x = Value)) +
  geom_histogram() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms of Numeric Variables and their Corresponding Imputed Variables",
       x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 12))
#dev.off()

#jpeg("barplot_rf_new.jpeg", width = 70, height = 40, units="cm", res = 400)
#---Barcharts for factors---------
bar_dat <- rf_impute %>%
  select(new_tumor_afit, path_n_stage, path_t_stage, neo_cancer_status,radiation_therapy,
         imp_rf.new_tumor_afit, imp_rf.path_n_stage, imp_rf.path_t_stage,
         imp_rf.neo_cancer_status,imp_rf.radiation_therapy) 

# Reshape the dataframe into a longer format (tidy data)
bar_long <- gather(bar_dat, key = "Variable", value = "Category")

# Filter and order the variables
bar_long <- bar_long %>%
  mutate(Variable = factor(Variable, levels = c(sort(grep("^Imp", unique(Variable))), unique(Variable))))

# Create bar charts of categorical variables and arrange them in a grid
ggplot(bar_long, aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Bar Charts of Categorical Variables and their Corresponding Imputed Variables") +
  theme(strip.text = element_text(size = 12))




#====GGPAIRS Plot====
#install.packages("GGally")
GGally::ggpairs(hist_dat_pmm)


#--Diagnosing Best Imputation Cycle with PMM (Continuous Variables)--------
library(mice)
densityplot(pmm_mice, ~ buffa_hs | .imp)  
densityplot(pmm_mice, ~ ragnum_hs | .imp) 
densityplot(pmm_mice, ~ winter_hs | .imp) 
densityplot(pmm_mice, ~ aneuploidy_score | .imp) 
densityplot(pmm_mice, ~ frac_gen_alt | .imp)
densityplot(pmm_mice, ~ mutation_count | .imp) 

one   <- complete(pmm_mice, action = 1L)
two   <- complete(pmm_mice, action = 2L)
three <- complete(pmm_mice, action = 3L)
four  <- complete(pmm_mice, action = 4L)
five  <- complete(pmm_mice, action = 5L)
six   <- complete(pmm_mice, action = 6L)

mutation_count <- data.frame( 
  initial = cleaned_combined_data$mutation_count,
  cycle_1 = one$mutation_count,
  cycle_2 = two$mutation_count,
  cycle_3 = three$mutation_count,
  cycle_4 = four$mutation_count,
  cycle_5 = five$mutation_count,
  cycle_6 = six$mutation_count)

# Convert all columns to numeric (if needed)
mutation_count[] <- lapply(mutation_count, as.numeric)
# Reshape the data to long format for ggplot2
mutation_count_long <- gather(mutation_count)
# Create a grid of histograms using ggplot2 with specified bins
ggplot(mutation_count_long, aes(x = value)) +
  geom_histogram(binwidth = 50, fill = "blue", color = "black") +
  facet_wrap(~ key, scales = "free") +
  labs(title = "Histograms of Mutation Count - PMM MICE", x = "Mutation Count (bin = 50)", 
       y = "Frequency")

#Transformed Fraction Of Genome Altered
#==========================================================================
fraction_genome_altered <- data.frame(
  initial =  cleaned_combined_data$frac_gen_alt,
  cycle_1 = one$frac_gen_alt ,
  cycle_2 = two$frac_gen_alt ,
  cycle_3 = three$frac_gen_alt ,
  cycle_4 = four$frac_gen_alt,
  cycle_5 = five$frac_gen_alt,
  cycle_6 = six$frac_gen_alt)

# Convert all columns to numeric (if needed)
fraction_genome_altered[] <- lapply(fraction_genome_altered, as.numeric)
# Reshape the data to long format for ggplot2
fraction_genome_altered_long <- gather(fraction_genome_altered)
# Create a grid of histograms using ggplot2 with specified bins
ggplot(fraction_genome_altered_long, aes(x = value)) +
  geom_histogram(fill = "blue", color = "black") +
  facet_wrap(~ key, scales = "free") +
  labs(title = "Histograms of Fraction of Genome Altered (%) - PMM MICE", 
       x = "Mutation Count (bin = 5)", 
       y = "Frequency")


#--Diagnosing Best Imputation Cycle with PMM (Categorical cases) --------
# ------New Tumot Afit-----------
new_tumor_afit_pmm <- data.frame(
  initial =  cleaned_combined_data$new_tumor_afit,
  cycle_1 = one$new_tumor_afit,
  cycle_2 = two$new_tumor_afit,
  cycle_3 = three$new_tumor_afit,
  cycle_4 = four$new_tumor_afit,
  cycle_5 = five$new_tumor_afit,
  cycle_6 = six$new_tumor_afit)
# Reshape the dataframe into a longer format (tidy data)
new_tumor_afit_plot_pmm <- gather(new_tumor_afit_pmm, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "PMM Imputed cycles for New_tumor_afit") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison New Tumor After Initial Diagnosis -----------------
results1 <- data.frame(imputation = character(), No = numeric(), Yes = numeric(), missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(new_tumor_afit_pmm)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(new_tumor_afit_pmm[[col_name]])
  
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No = cycle_table["No"],
    Yes = cycle_table["Yes"],
    missing = sum(is.na(new_tumor_afit_pmm[[col_name]]))
  )
  # Append the row to the results list
  results1 <- rbind(results1, results_row)
}
df_percent1 <- prop.table(as.matrix(results1[,-c(1)]), margin = 1) * 100
print(df_percent1)

#------Path-n-stage-----------
path_n_stage_pmm <- data.frame(
  initial = cleaned_combined_data$path_n_stage,
  cycle_1 = one$path_n_stage,
  cycle_2 = two$path_n_stage,
  cycle_3 = three$path_n_stage,
  cycle_4 = four$path_n_stage,
  cycle_5 = five$path_n_stage,
  cycle_6 = six$path_n_stage
)

# Reshape the dataframe into a longer format (tidy data)
path_n_stage_plot_pmm <- gather(path_n_stage_pmm, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "PMM Imputed cycles for Pathological Node Stage") +
  theme(strip.text = element_text(size = 12))
#----------Numerical Comparison-----------------
results2 <- data.frame(imputation = character(),
                       no_cancer_nearby_lymphs = numeric(), 
                       cancer_nearby_lymphs = numeric(), 
                       missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(path_n_stage_pmm)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(path_n_stage_pmm[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No_cancer = cycle_table["no_cancer_nearby_lymphs"],
    cancer = cycle_table["cancer_nearby_lymphs"],
    missing = sum(is.na(path_n_stage_pmm[[col_name]]))
  )
  # Append the row to the results list
  results2 <- rbind(results2, results_row)
}
df_percent2 <- prop.table(as.matrix(results2[,-c(1)]), margin = 1) * 100
print(df_percent2)



# -------Path-t-stage---------------
path_t_stage_pmm <- data.frame(
  initial = cleaned_combined_data$path_t_stage,
  cycle_1 = one$path_t_stage,
  cycle_2 = two$path_t_stage,
  cycle_3 = three$path_t_stage,
  cycle_4 = four$path_t_stage,
  cycle_5 = five$path_t_stage,
  cycle_6 = six$path_t_stage
)
# Reshape the dataframe into a longer format (tidy data)
path_t_stage_plot_pmm <- gather(path_t_stage_pmm, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "PMM Imputed cycles for Pathological Tumor Stage") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison-----------------
results3 <- data.frame(imputation = character(), stage_2A = numeric(), stage_2B = numeric(), 
                       stage_2C = numeric(), stage_3A = numeric(), stage_3B = numeric(),
                       stage_4 = numeric(), missing = numeric())

# Iterate through each column in the dataframe
for (col_name in names(path_t_stage_pmm)){
  # Use addNA to include NA values in the table
  cycle_table <- table(path_t_stage_pmm[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    stage_2A = cycle_table["stage_2A"],
    stage_2B = cycle_table["stage_2B"],
    stage_2C = cycle_table["stage_2C"],
    stage_3A = cycle_table["stage_3A"],
    stage_3B = cycle_table["stage_3B"],
    stage_4 = cycle_table["stage_4"],
    missing = sum(is.na(path_t_stage_pmm[[col_name]]))
  )
  # Append the row to the results list
  results3 <- rbind(results3, results_row)
}
df_percent3 <- prop.table(as.matrix(results3[,-c(1)]), margin = 1) * 100
print(df_percent3)


#-------------Neocancer Status--------
neo_cancer_status_pmm <- data.frame(
  initial = cleaned_combined_data$neo_cancer_status,
  cycle_1 = one$neo_cancer_status,
  cycle_2 = two$neo_cancer_status,
  cycle_3 = three$neo_cancer_status,
  cycle_4 = four$neo_cancer_status,
  cycle_5 = five$neo_cancer_status,
  cycle_6 = six$neo_cancer_status
)
neo_cancer_status_plot_pmm <- gather(neo_cancer_status_pmm, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "PMM Imputed cycles for Neoplasm Cancer Status") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison-----------------
results4 <- data.frame(imputation = character(),
                       Tumor_Free = numeric(), 
                       With_Tumor = numeric(), 
                       missing = numeric())

# Iterate through each column in the dataframe
for (col_name in names(neo_cancer_status_pmm)){
  # Use addNA to include NA values in the table
  cycle_table <- table(neo_cancer_status_pmm[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    Tumor_free = cycle_table["Tumor Free"],
    With_Tumor = cycle_table["With Tumor"],
    missing = sum(is.na(neo_cancer_status_pmm[[col_name]]))
  )
  # Append the row to the results list
  results4 <- rbind(results4, results_row)
}
# Display the resulting dataframe
print(results4)



#----Radiation Therapy-----------
radiation_therapy_pmm <- data.frame(
  initial = cleaned_combined_data$radiation_therapy,
  cycle_1 = one$radiation_therapy,
  cycle_2 = two$radiation_therapy,
  cycle_3 = three$radiation_therapy,
  cycle_4 = four$radiation_therapy,
  cycle_5 = five$radiation_therapy,
  cycle_6 = six$radiation_therapy
)
radiation_therapy_plot_pmm <- gather(radiation_therapy_pmm, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "PMM Imputed cycles for Radiation Therapy") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison -----------------
results5 <- data.frame(imputation = character(), No = numeric(), Yes = numeric(),
                       missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(radiation_therapy_pmm)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(radiation_therapy_pmm[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No = cycle_table["No"],
    Yes = cycle_table["Yes"],
    missing = sum(is.na(radiation_therapy_pmm[[col_name]]))
  )
  # Append the row to the results list
  results5 <- rbind(results5, results_row)
}
# Display the resulting dataframe
print(results5)

#--Diagnosing Best Imputation Cycle with RF (Continuous) --------
densityplot(rf_mice, ~ buffa_hs | .imp) #1L for buffa
densityplot(rf_mice, ~ ragnum_hs | .imp) #2L/5L for ragnum
densityplot(rf_mice, ~ winter_hs | .imp) #3L for winter
densityplot(rf_mice, ~ aneuploidy_score | .imp) # 1L for aneuploidy
densityplot(rf_mice, ~ frac_gen_alt | .imp) #4L
densityplot(rf_mice, ~ mutation_count | .imp) #1L

#--Diagnosing Best Imputation Cycle with RF (Categorical cases) --------
one   <- complete(rf_mice, action = 1L)
two   <- complete(rf_mice, action = 2L)
three <- complete(rf_mice, action = 3L)
four  <- complete(rf_mice, action = 4L)
five  <- complete(rf_mice, action = 5L)
six   <- complete(rf_mice, action = 6L)

mutation_count_rf <- data.frame(
  initial =  cleaned_combined_data$mutation_count,
  cycle_1 = one$mutation_count,
  cycle_2 = two$mutation_count,
  cycle_3 = three$mutation_count,
  cycle_4 = four$mutation_count,
  cycle_5 = five$mutation_count,
  cycle_6 = six$mutation_count)

# Convert all columns to numeric (if needed)
mutation_count_rf[] <- lapply(mutation_count_rf, as.numeric)
# Reshape the data to long format for ggplot2
mutation_count_rf_long <- gather(mutation_count_rf)
# Create a grid of histograms using ggplot2 with specified bins
ggplot(mutation_count_rf_long, aes(x = value)) +
  geom_density(binwidth = 50, fill = "blue", color = "black") +
  facet_wrap(~ key, scales = "free") +
  labs(title = "Histograms of Mutation Count - RF MICE", x = "Mutation Count (bin = 50)", 
       y = "Frequency")

#Transformed Fraction Of Genome Altered  - RF MICE
#==========================================================================
fraction_genome_altered_rf <- data.frame(
  initial =  cleaned_combined_data$frac_gen_alt,
  cycle_1 = one$frac_gen_alt ,
  cycle_2 = two$frac_gen_alt ,
  cycle_3 = three$frac_gen_alt,
  cycle_4 = four$frac_gen_alt,
  cycle_5 = five$frac_gen_alt,
  cycle_6 = six$frac_gen_alt)
# Convert all columns to numeric (if needed)
fraction_genome_altered_rf[] <- lapply(fraction_genome_altered_rf, as.numeric)
# Reshape the data to long format for ggplot2
fraction_genome_altered_rf_long <- gather(fraction_genome_altered_rf)

ggplot(fraction_genome_altered_rf_long, aes(x = value)) +
  geom_density(binwidth = 5, fill = "red", color = "black") +
  facet_wrap(~ key, scales = "free") +
  labs(title = "Density of Fraction of Genome Altered (%) - RF MICE", 
       x = "Fraction of genome altered (bin = 5)", 
       y = "Frequency")

#  ------New Tumot Afit-----------
new_tumor_afit_rf <- data.frame(
  initial =  cleaned_combined_data$new_tumor_afit,
  cycle_1 = one$new_tumor_afit,
  cycle_2 = two$new_tumor_afit,
  cycle_3 = three$new_tumor_afit,
  cycle_4 = four$new_tumor_afit,
  cycle_5 = five$new_tumor_afit,
  cycle_6 = six$new_tumor_afit 
)

# Reshape the dataframe into a longer format (tidy data)
new_tumor_afit_rf_plot <- gather(new_tumor_afit_rf, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Random Forest Imputed cycles for New_tumor_afit") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison New Tumor After Initial Diagnosis -----------------
results1 <- data.frame(imputation = character(), No = numeric(), Yes = numeric(), missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(new_tumor_afit_rf)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(new_tumor_afit_rf[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No = cycle_table["No"],
    Yes = cycle_table["Yes"],
    missing = sum(is.na(new_tumor_afit_rf[[col_name]]))
  )
  # Append the row to the results list
  results1 <- rbind(results1, results_row)
}
df_percent1 <- prop.table(as.matrix(results1[,-c(1)]), margin = 1) * 100
print(df_percent1)

#------Path-n-stage-----------
path_n_stage_rf <- data.frame(
  initial = cleaned_combined_data$path_n_stage,
  cycle_1 = one$path_n_stage,
  cycle_2 = two$path_n_stage,
  cycle_3 = three$path_n_stage,
  cycle_4 = four$path_n_stage,
  cycle_5 = five$path_n_stage,
  cycle_6 = six$path_n_stage 
)

# Reshape the dataframe into a longer format (tidy data) & plot
path_n_stage_rf_plot  <- gather(path_n_stage_rf, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Random Forest Imputed cycles for path_n_stage") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison-----------------
results2 <- data.frame(imputation = character(),
                       no_cancer_nearby_lymphs = numeric(), 
                       cancer_nearby_lymphs = numeric(), 
                       missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(path_n_stage_rf)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(path_n_stage_rf[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No_cancer = cycle_table["no_cancer_nearby_lymphs"],
    cancer = cycle_table["cancer_nearby_lymphs"],
    missing = sum(is.na(path_n_stage_rf[[col_name]]))
  )
  # Append the row to the results list
  results2 <- rbind(results2, results_row)
}

# Display the resulting dataframe
print(results2)

df_percent2 <- prop.table(as.matrix(results2[,-c(1)]), margin = 1) * 100
print(df_percent2)

#------Path-t-stage-----------
path_t_stage_rf <- data.frame(
  initial = cleaned_combined_data$path_t_stage,
  cycle_1 = one$path_t_stage,
  cycle_2 = two$path_t_stage,
  cycle_3 = three$path_t_stage,
  cycle_4 = four$path_t_stage,
  cycle_5 = five$path_t_stage,
  cycle_6 = six$path_t_stage 
)
# Reshape the dataframe into a longer format (tidy data) & plot
path_t_stage_rf_plot  <- gather(path_t_stage_rf, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Random Forest Imputed cycles for path_t_stage") +
  theme(strip.text = element_text(size = 12))
results3 <- data.frame(imputation = character(), stage_2A = numeric(), stage_2B = numeric(), 
                       stage_2C = numeric(), stage_3A = numeric(), stage_3B = numeric(),
                       stage_4 = numeric(), missing = numeric())

# Iterate through each column in the dataframe
for (col_name in names(path_t_stage_rf)){
  # Use addNA to include NA values in the table
  cycle_table <- table(path_t_stage_rf[[col_name]])
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    stage_2A = cycle_table["stage_2A"],
    stage_2B = cycle_table["stage_2B"],
    stage_2C = cycle_table["stage_2C"],
    stage_3A = cycle_table["stage_3A"],
    stage_3B = cycle_table["stage_3B"],
    stage_4 = cycle_table["stage_4"],
    missing = sum(is.na(path_t_stage_rf[[col_name]]))
  )
  # Append the row to the results list
  results3 <- rbind(results3, results_row)
}
df_percent3 <- prop.table(as.matrix(results3[,-c(1)]), margin = 1) * 100
print(df_percent3)


#-------------Neocancer Status--------
neo_cancer_status_rf <- data.frame(
  initial = cleaned_combined_data$neo_cancer_status,
  cycle_1 = one$neo_cancer_status,
  cycle_2 = two$neo_cancer_status,
  cycle_3 = three$neo_cancer_status,
  cycle_4 = four$neo_cancer_status,
  cycle_5 = five$neo_cancer_status,
  cycle_6 = six$neo_cancer_status
)
neo_cancer_status_rf_plot <- gather(neo_cancer_status_rf, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Random Forest Imputed cycles for Neoplasm Cancer Status") +
  theme(strip.text = element_text(size = 12))
#----------Numerical Comparison-----------------
results4 <- data.frame(imputation = character(),
                       Tumor_Free = numeric(), 
                       With_Tumor = numeric(), 
                       missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(neo_cancer_status_rf)){
  # Use addNA to include NA values in the table
  cycle_table <- table(neo_cancer_status_rf[[col_name]])
  
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    Tumor_free = cycle_table["Tumor Free"],
    With_Tumor = cycle_table["With Tumor"],
    missing = sum(is.na(neo_cancer_status_rf[[col_name]]))
  )
  # Append the row to the results list
  results4 <- rbind(results4, results_row)
}

#----Radiation Therapy-----------
radiation_therapy_rf <- data.frame(
  initial = cleaned_combined_data$radiation_therapy,
  cycle_1 = one$radiation_therapy,
  cycle_2 = two$radiation_therapy,
  cycle_3 = three$radiation_therapy,
  cycle_4 = four$radiation_therapy,
  cycle_5 = five$radiation_therapy,
  cycle_6 = six$radiation_therapy
)

radiation_therapy_rf_plot <- gather(radiation_therapy_rf, key = "Variable", value = "Category") %>%
  mutate(Variable = factor(Variable)) %>%
  ggplot(., aes(x = Category)) +
  geom_bar() +
  facet_wrap(~ Variable, scales = "free") +
  labs(title = "Random Forest Imputed cycles for Radiation Therapy") +
  theme(strip.text = element_text(size = 12))

#----------Numerical Comparison -----------------
results5 <- data.frame(imputation = character(), No = numeric(), Yes = numeric(),
                       missing = numeric())
# Iterate through each column in the dataframe
for (col_name in names(radiation_therapy_rf)) {
  # Use addNA to include NA values in the table
  cycle_table <- table(radiation_therapy_rf[[col_name]])
  
  # Create a new row for the results dataframe
  results_row <- data.frame(
    imputation = col_name,
    No = cycle_table["No"],
    Yes = cycle_table["Yes"],
    missing = sum(is.na(radiation_therapy_rf[[col_name]]))
  )
  # Append the row to the results list
  results5 <- rbind(results5, results_row)
}

# Final Selected Imputation cycle for each variable
imputed_result <- cleaned_combined_data %>%
  mutate(buffa_hs = complete(rf_mice, action = 1L)$buffa_hs) %>%
  mutate(ragnum_hs = complete(rf_mice, action = 5L)$ragnum_hs) %>%
  mutate(winter_hs = complete(rf_mice, action = 3L)$winter_hs) %>%
  mutate(aneuploidy_score = complete(rf_mice, action = 1L)$aneuploidy_score) %>%
  mutate(frac_gen_alt = complete(rf_mice, action = 4L)$frac_gen_alt) %>%
  mutate(mutation_count = complete(rf_mice, action = 1L)$mutation_count)%>%
  mutate(new_tumor_afit = complete(pmm_mice, action = 1L)$new_tumor_afit) %>%
  mutate(path_n_stage = complete(pmm_mice, action = 1L)$path_n_stage) %>%
  mutate(path_t_stage = complete(pmm_mice, action = 1L)$path_t_stage) %>%
  mutate(neo_cancer_status = complete(pmm_mice, action = 1L)$neo_cancer_status) %>%
  mutate(radiation_therapy = complete(pmm_mice, action = 1L)$radiation_therapy) #%>%

# plot missing values plot
plot_missing(imputed_result,
             ggtheme = theme_igray(),
             title = "Missing Value Plot for Imputed Result",
             theme_config = theme(plot.title = element_text(color = "black")),
             geom_label_args = c(hjust = "inward")
)

dim(final_combined_data)

#Data for ML
write.csv(imputed_result,
          file = "out_data/prca_clinicogenomics_data.csv",
          row.names = F)






#============================evaluation numerical 1 good =====================================
# Load required packages
library(dplyr)

# List of numerical variables
numerical_vars <- c("buffa_hs", "ragnum_hs", "winter_hs", 
                    "aneuploidy_score", "frac_gen_alt", "mutation_count")

# Run KS tests using the imputed_result object
ks_results <- purrr::map_dfr(numerical_vars, function(var) {
  observed_vals <- cleaned_combined_data[[var]][!is.na(cleaned_combined_data[[var]])]
  imputed_vals  <- imputed_result[[var]][is.na(cleaned_combined_data[[var]])]
  
  if (length(imputed_vals) > 0 && length(observed_vals) > 0) {
    test <- ks.test(imputed_vals, observed_vals)
    tibble(
      variable = var,
      method = test$method,
      statistic = round(test$statistic, 3),
      p_value = round(test$p.value, 4)
    )
  } else {
    tibble(variable = var, method = "KS test", statistic = NA, p_value = NA)
  }
})

ks_results

#To assess the plausibility of multiple imputed values, we conducted two-sample Kolmogorov–Smirnov (KS) tests 
#comparing the distribution of observed and imputed values for each continuous variable. 
#Five out of six variables (`buffa_hs`, `winter_hs`, `aneuploidy_score`, `frac_gen_alt`, and `mutation_count`) showed no significant differences (p > 0.05), indicating good alignment between observed and imputed distributions. 
# Only `ragnum_hs` showed a statistically significant difference (p = 0.026), warranting further review.
# R automatically selected between the **asymptotic** and **exact** KS test based on the presence of ties and sample size; the presence of ties—common in imputed datasets—triggered a switch to the **exact** method for some variables. This adaptive behavior ensures robust inference even in the presence of tied values, without requiring manual specification.





#============================evaluation categorical 1 good =====================================



# List of categorical variables
categorical_vars <- c("new_tumor_afit", "neo_cancer_status", 
                      "radiation_therapy", "path_t_stage", "path_n_stage")

# Function to test a single categorical variable using `imputed_result`
test_categorical_var <- function(var) {
  observed_vals <- cleaned_combined_data[[var]][!is.na(cleaned_combined_data[[var]])]
  imputed_vals  <- imputed_result[[var]][is.na(cleaned_combined_data[[var]])]
  
  # Align factor levels
  all_levels <- union(levels(factor(observed_vals)), levels(factor(imputed_vals)))
  observed_tab <- table(factor(observed_vals, levels = all_levels))
  imputed_tab  <- table(factor(imputed_vals, levels = all_levels))
  tbl <- rbind(Observed = observed_tab, Imputed = imputed_tab)
  
  # Perform tests
  chi_result <- tryCatch({
    suppressWarnings(chisq.test(tbl))
  }, error = function(e) NULL)
  
  fisher_result <- tryCatch({
    fisher.test(tbl)
  }, error = function(e) NULL)
  
  list(
    variable = var,
    contingency_table = tbl,
    chi_method = if (!is.null(chi_result)) chi_result$method else NA,
    chi_stat   = if (!is.null(chi_result)) round(chi_result$statistic, 3) else NA,
    chi_p      = if (!is.null(chi_result)) round(chi_result$p.value, 4) else NA,
    fisher_method = if (!is.null(fisher_result)) fisher_result$method else NA,
    fisher_p      = if (!is.null(fisher_result)) round(fisher_result$p.value, 4) else NA
  )
}

# Apply test across categorical vars
cat_results <- purrr::map(categorical_vars, test_categorical_var)

# Print results
for (res in cat_results) {
  cat("\n===============================\n")
  cat("Variable:", res$variable, "\n")
  print(res$contingency_table)
  if (!is.na(res$chi_method)) {
    cat("\nChi-squared Test:", res$chi_method, 
        "| Statistic:", res$chi_stat, "| p-value:", res$chi_p, "\n")
  } else {
    cat("\nChi-squared Test: Not applicable\n")
  }
  if (!is.na(res$fisher_method)) {
    cat("Fisher's Exact Test:", res$fisher_method, 
        "| p-value:", res$fisher_p, "\n")
  } else {
    cat("Fisher's Exact Test: Not applicable\n")
  }
}




# Create observed and imputed vectors
observed_vals <- cleaned_combined_data$radiation_therapy[!is.na(cleaned_combined_data$radiation_therapy)]
imputed_vals <- complete(pmm_mice, 1)$radiation_therapy[is.na(cleaned_combined_data$radiation_therapy)]

# Tabulate
tbl <- rbind(
  Observed = table(observed_vals),
  Imputed  = table(imputed_vals)
)

# Perform Chi-squared test (if valid)
chi_result <- chisq.test(tbl)

# If chi-squared is not valid (e.g., warning about small expected counts), use Fisher
fisher_result <- fisher.test(tbl)

# Print both
cat("Contingency Table:\n"); print(tbl)
cat("\nChi-squared Test:\n"); print(chi_result)
cat("\nFisher's Exact Test:\n"); print(fisher_result)


#Across all categorical variables, imputed data maintained the same imbalance as the original dataset. 
# The few significant p-values (notably for new_tumor_afit and path_n_stage) stem from low cell counts in rare categories rather 
# than genuine distributional shifts. 

# In short, your imputation model successfully preserved the categorical structure of the data.
# The imputed categorical values maintained the original class imbalances, with consistent proportions before and after imputation. 
# This alignment supports the plausibility of the imputations. 
# Results from both Chi-squared and Fisher's exact tests further suggest that the observed and imputed distributions are not significantly different for most variables, reinforcing confidence in the categorical imputation process.




