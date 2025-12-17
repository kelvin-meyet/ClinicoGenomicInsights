#====================================================
#==========Exploratory and Univariate Analaysis =====
#====================================================
# Dynamic directory based on the location of the research folder
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

#---> Libraries <----
require(dplyr)
require(ggplot2)
require(caTools)
require(glmnet)
library(survival)
library(kernelshap)
library(shapviz)
library(survex)
library(survival)
library(survminer)

#----load clinical and genomic data 
#require(dplyr)
clinical_all_genomic <- read.csv("cleaned_imputed_data_for_ml-integrated.csv")
dat_new <- clinical_all_genomic

#---convert integer to numeric  and character to categorical--------
# -- Integer to numeric -- good
integer_cols <- sapply(dat_new, is.integer)
dat_new[integer_cols] <- lapply(dat_new[integer_cols], as.numeric) 
# Convert character columns to categorical
char_cols <- sapply(dat_new, is.character)
dat_new[char_cols] <- lapply(dat_new[char_cols], factor)
#------Arrange data-----------
endpoints <- dat_new %>% 
  dplyr::select(pfs_months, pfs_status)
all_vars <- dat_new %>% 
  dplyr::select(-c(pfs_months, pfs_status))
pfs_status <- endpoints$pfs_status
pfs_months <- endpoints$pfs_months

#------------------------------------------------------------------------------
#----------- Exploratory with clinical-all-genomic ----------------------------
# Calculate the number of censored records
data_viz <- dat_new
num_censored <- nrow(data_viz) - sum(data_viz$pfs_status)   #all patients - patients with progression free status
cat(sprintf("%.1f%% of records are censored\n", (num_censored / nrow(data_viz)) * 100))


png(filename = "Figure2.png", width = 1920, height = 1080, units = "px", res = 300)
ggplot(data_viz, aes(x = pfs_months,
                     fill = factor(pfs_status, 
                                   labels = c("Censored Events", "Observed Progression")))) +
  geom_histogram(alpha = 0.7, bins = 30, position = "stack", color = "white") +
  labs(x = "Time (months)", y = "Frequency",
       title = "Time Distribution for Censored and Disease Progression Patients") + 
  #scale_fill_brewer(palette = "Oranges", name = "Status", labels = c("Censored", "Progress_Free"))
  scale_fill_brewer(palette = "Set1", name = "Status", labels = c("Censored", "Progression")) +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(t=10, r=5, b=10, l=60)
    )
dev.off()


#=====================Fivenumber Summary====================
summary_stats <- data_viz %>%
  group_by(pfs_status) %>%
  summarise(
    n = n(),
    mean = mean(pfs_months, na.rm = TRUE),
    stdev = sd(pfs_months, na.rm = TRUE),
    min = min(pfs_months, na.rm = TRUE),
    q1 = quantile(pfs_months, 0.25, na.rm = TRUE),
    median = median(pfs_months, na.rm = TRUE),
    q3 = quantile(pfs_months, 0.75, na.rm = TRUE),
    max = round(max(pfs_months, na.rm = TRUE), 2)
  )
#relabel the pfs_status values for clarity
summary_stats$pfs_status <- factor(summary_stats$pfs_status, 
                                   levels = c(0,1), 
                                   labels = c("Censored", "Progression"))
print(summary_stats)



#----- Partition Data for train and test data exploration----------
#require(caTools)
set.seed(1) 
comb_dat <- cbind(all_vars, pfs_months, pfs_status)
sample = sample.split( Y= comb_dat$pfs_status, SplitRatio = 0.70)
train_viz <- train <- subset(comb_dat, sample == TRUE)
test_viz <- test  <- subset(comb_dat, sample == FALSE)

#Distribution of obs in Train and Test Set 
# Combine train and test data for easy comparison
train_viz$Set <- "Train"
test_viz$Set <- "Test"
combined_data_viz <- bind_rows(train_viz, test_viz)


##Plot 1 - good 
png(filename = "Figure A.png", width = 1920, height = 1080, units = "px", res = 300)

combined_data_viz$Set <- factor(combined_data_viz$Set,
                                levels = c("Train", "Test"))
ggplot(combined_data_viz, 
       aes(x = pfs_months,
           y = ..density.., 
           fill = Set)) +
  geom_histogram(fill = "#1f78b4", alpha = 0.75, bins = 30, color = "black") +
  facet_wrap(~Set, ncol = 1,
             labeller = labeller(Set = c(Train = "Training Data", Test = 'Test Data'))) +
  labs(
    title = "Probability Histograms of PFS Time by Dataset",
    x = "Time (months)",
    y = "Probability Density"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    #plot.margin = margin(t=10, r=5, b=10, l=60)
    #strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"   # removes legend
  )
dev.off()



# Plot 2: Proportion of `pfs_status`
pfs_status_summary <- combined_data_viz %>%
  group_by(Set, pfs_status) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count))

ggplot(pfs_status_summary, aes(x = Set, y = percentage, fill = as.factor(pfs_status))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  labs(
    title = "Distribution of Progression Free Survival Status",
    x = "Dataset",
    y = "Proportion of Events",
    fill = "PFS Status"
  ) +
  scale_fill_brewer(palette = "Set1", name = "Status", labels = c("Censored", "Progression")) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


#========================survival curve===================================
#---------Radiation Therapy ----------
png(filename = "Figure4.png", width = 1920, height = 1080, units = "px", res = 300)
comb_dat$pfs_status = factor(ifelse(comb_dat$pfs_status == 1, "progression", "censored"))
fit_rad_therapy <- survfit(Surv(comb_dat$pfs_months, comb_dat$pfs_status=="progression") ~ comb_dat$radiation_therapy,
                           data=comb_dat)
ggsurvplot(fit_rad_therapy,
           legend.title     = "Radiation Therapy (RT):",
           legend.labs      = c("No", "Yes"),
           conf.int         = TRUE,
           pval             = TRUE,
           pval.size        = 2.5,
           tables.height    = 0.25,
           tables.theme     = theme_cleantable() + 
             theme(axis.text.y = element_text(size=8)),
           risk.table       = TRUE, # Add risk table
           risk.table.col   = "strata", # Change risk table color by groups
           risk.table.fontsize = 2.5,  # smaller font size
           linetype         = "strata", # Change line type by groups
           surv.median.line = "hv",  # Specify median survival
           ggtheme          = theme_bw(), # Change ggplot2 theme
           palette          = c("#00AFBB", "#FC4E07")) #"#2D2D2D"  #E7B800 -- "#E7B800", "#2E9FDF"
dev.off()


#---------Neoplasm cancer status----------
png(filename = "Figure B.png", width = 1920, height = 1080, units = "px", res = 300)
fit_neoplasm_cancer <- survfit(Surv(comb_dat$pfs_months, comb_dat$pfs_status=="progression") ~ comb_dat$neo_cancer_status,
                               data = comb_dat)
ggsurvplot(fit_neoplasm_cancer,
           legend.title     = "Neoplasm Cancer Status:",
           legend.labs      = c("Tumor Free", "With Tumor"),
           conf.int         = TRUE,
           pval             = TRUE,
           pval.size        = 2.5,
           tables.height    = 0.2,
           tables.theme     = theme_cleantable() + 
             theme(axis.text.y = element_text(size=8)),
           risk.table       = TRUE, # Add risk table
           risk.table.col   = "strata", # Change risk table color by groups
           risk.table.fontsize = 2.5, #smaller font size
           linetype         = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme          = theme_bw(), # Change ggplot2 theme
           palette          = c("#00AFBB", "#FC4E07"))
dev.off()

#---------History of Neo-adjuvant treatment----------
png(filename = "Figure C.png", width = 1920, height = 1080, units = "px", res = 300)
fit_hist_neoadjuv_trtmnt <- survfit(Surv(comb_dat$pfs_months, comb_dat$pfs_status=="progression") ~ comb_dat$hist_neoadjuv_trtmnt,
                                    data = comb_dat)
ggsurvplot(fit_hist_neoadjuv_trtmnt,
           legend.title     = "Neoadjuvant Treatment History:",
           legend.labs      = c("No", "Yes"),
           conf.int         = TRUE,
           conf.int.style   = "step", 
           pval             = TRUE,
           pval.size        = 2.5,
           tables.height    = 0.2,
           xlab             = "Time (months)",
           tables.theme     = theme_cleantable() + theme(axis.text.y = element_text(size=8)),
           risk.table.fontsize = 2.5, 
           risk.table       = TRUE, 
           risk.table.col   = "strata", 
           linetype         = "strata", 
           surv.median.line = "hv", 
           ggtheme          = theme_bw(), 
           palette          = c("#00AFBB", "#FC4E07"))
dev.off()




#---------New Tumor After Initial Treatment----------
png(filename = "Figure3.png", width = 1920, height = 1080, units = "px", res = 300)
fit_new_tumor_afit <- survfit(Surv(comb_dat$pfs_months, comb_dat$pfs_status=="progression") ~ comb_dat$new_tumor_afit,
                              data = comb_dat)
ggsurvplot(fit_new_tumor_afit,
           legend.title     = "New Tumor After Initial Treatment:",
           legend.labs      = c("No", "Yes"),
           conf.int         = TRUE,
           pval             = TRUE,
           pval.size        = 2.5,
           tables.height    = 0.2,
           tables.theme     = theme_cleantable() + 
             theme(axis.text.y = element_text(size=8)),
           risk.table       = TRUE, # Add risk table
           risk.table.col   = "strata", # Change risk table color by groups
           risk.table.fontsize = 2.5, 
           linetype         = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme          = theme_bw(), # Change ggplot2 theme
           palette          = c("#00AFBB", "#FC4E07"))
dev.off()




