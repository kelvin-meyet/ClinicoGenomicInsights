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
clinical_all_genomic <- read.csv("prca_clinicogenomics_data.csv")
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
png(filename = "Figure S1.png", width = 1920, height = 1080, units = "px", res = 300)

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

#====================================================================
#======================== REVISED KM CURVE PLOTS =================
#====================================================================
comb_dat$pfs_status = factor(ifelse(comb_dat$pfs_status == 1, "progression", "censored"))
#=====REVISION - KM PLOTS FOR PUBLICATION ======
# library(survival)
# library(survminer)
# library(ggplot2)

#----------------------------
# 1) Fit KM - NTAIT
#-----------------------------
fit_new_tumor_afit <- survfit(
  Surv(pfs_months, pfs_status == "progression") ~ new_tumor_afit,
  data = comb_dat
)


# Log-rank p-value (italic p)
sd <- survdiff(
  Surv(pfs_months, pfs_status == "progression") ~ new_tumor_afit,
  data = comb_dat
)
p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

# Annotate ggplot 
p_label <- if (p < 0.0001) {
  "italic(p) < 0.0001"
} else {
  paste0("italic(p) == ", signif(p, 3))
}


risk_table_theme <- theme_cleantable(base_size = 11) +
  theme(
    axis.title.y = element_blank(),  # removes the left-side variable label
    axis.text.y  = element_text(size = 10),
    plot.margin  = margin(2, 8, 6, 8)
  )

# Build plot (thin CI band, clean)

g <- ggsurvplot(fit_new_tumor_afit, data = comb_dat,
  xlab = "Time (months)",
  ylab = "Survival probability",
  legend.title = "New Tumor After Initial Treatment (NTAIT):",
  legend.labs  = c("No", "Yes"),
  conf.int       = TRUE,
  conf.int.alpha = 0.10,     # thinner/lighter CI band
  size           = 0.85,     # line thickness (publication-friendly)
  censor.size    = 3.0,
  surv.median.line = "hv",
  risk.table          = TRUE,
  risk.table.height   = 0.22,
  risk.table.col      = "strata",
  risk.table.fontsize = 3.2,
  ggtheme      = theme_bw(), #pub_theme,
  tables.theme = risk_table_theme,
  palette = c("#00AFBB", "#FC4E07")
)

# 5. Add italicized p to the main panel
g$plot <- g$plot +
  annotate(
    "text",
    x = 10, y = 0.20,      # adjust to avoid overlap
    label = p_label,
    parse = TRUE,
    size = 3
  )
print(g)

# 6) Export to png
png("Figure3.png", width = 1920, height = 1200, res = 300)
print(g)
dev.off()




#-------------------------------------------------------
# 2) Fit KM - Radiation Therapy (RT)
#-------------------------------------------------------
fit_radiation_therapy <- survfit(
  Surv(pfs_months, pfs_status == "progression") ~ radiation_therapy,
  data = comb_dat
)

# Log-rank p-value (italic p)
sd <- survdiff(
  Surv(pfs_months, pfs_status == "progression") ~ radiation_therapy,
  data = comb_dat
)
p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

# Annotate p value 
p_label <- if (p < 0.0001) {
  "italic(p) < 0.0001"
} else {
  paste0("italic(p) == ", signif(p, 3))
}

# Risk Table
risk_table_theme <- theme_cleantable(base_size = 11) +
  theme(
    axis.title.y = element_blank(),  # removes the left-side variable label
    axis.text.y  = element_text(size = 10),
    plot.margin  = margin(2, 8, 6, 8)
  )


# Build plot (thin CI band, clean)
g2 <- ggsurvplot(fit_radiation_therapy, data = comb_dat,
                xlab = "Time (months)",
                ylab = "Survival probability",
                legend.title = "Radiation Therapy (RT):",
                legend.labs  = c("No", "Yes"),
                conf.int       = TRUE,
                conf.int.alpha = 0.10,     # thinner/lighter CI band
                size           = 0.85,     # line thickness (publication-friendly)
                censor.size    = 3.0,
                surv.median.line = "hv",
                risk.table          = TRUE,
                risk.table.height   = 0.22,
                risk.table.col      = "strata",
                risk.table.fontsize = 3.2,
                ggtheme      = theme_bw(), #pub_theme,
                tables.theme = risk_table_theme,
                palette = c("#00AFBB", "#FC4E07")
)

# Add italicized p to the main panel
g2$plot <- g2$plot +annotate("text",
                             x = 10, y = 0.20,      # adjust to avoid overlap
                             label = p_label,
                             parse = TRUE,
                             size = 3)
print(g2)

# Export to png
png("Figure4.png", width = 1920, height = 1200, res = 300)
print(g2)
dev.off()




#---------------------------------------------
# 3) Fit KM - neoplasm_cancer
#----------------------------------------------
fit_neoplasm_cancer <- survfit(
  Surv(pfs_months, pfs_status == "progression") ~ neo_cancer_status,
  data = comb_dat
)


# Log-rank p-value (italic p)
sd <- survdiff(
  Surv(pfs_months, pfs_status == "progression") ~ neo_cancer_status,
  data = comb_dat
)
p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

# Annotate p-value 
p_label <- if (p < 0.0001) {
  "italic(p) < 0.0001"
} else {
  paste0("italic(p) == ", signif(p, 3))
}

#Risk table 
risk_table_theme <- theme_cleantable(base_size = 11) +
  theme(
    axis.title.y = element_blank(),  # removes the left-side variable label
    axis.text.y  = element_text(size = 10),
    plot.margin  = margin(2, 8, 6, 8)
  )

# Build plot (thin CI band, clean)
g3 <- ggsurvplot(fit_neoplasm_cancer, data = comb_dat,
                 xlab = "Time (months)",
                 ylab = "Survival probability",
                 legend.title = "Neoplasm Cancer Status (NCS):",
                 legend.labs  = c("Tumor Free", "With Tumor"),
                 conf.int       = TRUE,
                 #conf.int.style = "step",
                 conf.int.alpha = 0.10,     # thinner/lighter CI band
                 size           = 0.85,     # line thickness (publication-friendly)
                 censor.size    = 3.0,
                 surv.median.line = "hv",
                 risk.table          = TRUE,
                 risk.table.height   = 0.22,
                 risk.table.col      = "strata",
                 risk.table.fontsize = 3.2,
                 ggtheme      = theme_bw(), #pub_theme,
                 tables.theme = risk_table_theme,
                 palette = c("#00AFBB", "#FC4E07")
)

# Add italicized p to the main panel
g3$plot <- g3$plot +annotate("text",
                             x = 10, y = 0.20,      # adjust to avoid overlap
                             label = p_label,
                             parse = TRUE,
                             size = 3)
print(g3)

# Export to png
png("FigureS2.png", width = 1920, height = 1200, res = 300)
print(g3)
dev.off()






#------------------------------------------------------
# 4) Fit KM - History of Neo-adjuvant treatment -- Fix later
#------------------------------------------------------
#install.packages("remotes")
# remotes::install_version("ggplot2", version = "3.4.4")
# library(ggplot2)

fit_hist_neoadjuv_trtmnt <- survfit(
  Surv(pfs_months, pfs_status == "progression") ~ hist_neoadjuv_trtmnt,
  data = comb_dat
)

#Log-rank p-value (italic p)
sd <- survdiff(
  Surv(pfs_months, pfs_status == "progression") ~ hist_neoadjuv_trtmnt,
  data = comb_dat
)
p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

# Annotate p value 
p_label <- if (p < 0.0001) {
  "italic(p) < 0.0001"
} else {
  paste0("italic(p) == ", signif(p, 3))
}


risk_table_theme <- theme_cleantable(base_size = 11) +
  theme(
    axis.title.y = element_blank(),  # removes the left-side variable label
    axis.text.y  = element_text(size = 10),
    plot.margin  = margin(2, 8, 6, 8)
  )


# Build plot (thin CI band, clean)
g4 <- ggsurvplot(fit_hist_neoadjuv_trtmnt, data = comb_dat,
                 xlab = "Time (months)",
                 ylab = "Survival probability",
                 legend.title = "Neoadjuvant Treatment History (NTH):",
                 legend.labs  = c("No", "Yes"),
                 conf.int       = TRUE,
                 conf.int.style = "step", # "ribbon",
                 linetype = 1,
                 conf.int.alpha = 0.10,     # thinner/lighter CI band
                 size           = 0.85,     # line thickness (publication-friendly)
                 censor.size    = 3.0,
                 surv.median.line = "hv",
                 risk.table          = TRUE,
                 risk.table.height   = 0.22,
                 risk.table.col      = "strata",
                 risk.table.fontsize = 3.2,
                 linetype         = "strata",
                 ggtheme      = theme_bw(), #pub_theme,
                 tables.theme = risk_table_theme,
                 palette = c("#00AFBB", "#FC4E07")
)

# # Force the CI ribbon to use constant aesthetics (no varying linetype/linewidth)
# g4$plot <- g4$plot +
#   guides(linetype = "none")  # remove any lingering linetype mapping

# g4$plot$layers <- g4$plot$layers[!sapply(g4$plot$layers, function(x)
#   "GeomRibbon" %in% class(x$geom)
# )]

# Add italicized p to the main panel
g4$plot <- g4$plot + annotate("text",
                              x = 10, y = 0.20,      # adjust to avoid overlap
                              label = p_label,
                              parse = TRUE,
                              size = 2)

while (dev.cur() > 1) dev.off()
print(g4)

#Export 
png("FigureS3.png", width = 1920, height = 1200, res = 300)
print(g4)
dev.off()














