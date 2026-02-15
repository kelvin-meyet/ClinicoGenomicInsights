# Set your Working Directory and load R source code of dependent libraries
# install.packages("rstudioapi")
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

source("pcm_libraries.R")
source("rsf_libraries.R")
source("deepsurv_libraries.R")
source("deepsurv.R")

#-----------------------> Load Data <----------------------------
#--- Read any of the data below to run the respective ML Experiment Pipeline--
# clinical <- read.csv("prca_clinicogenomics_data.csv")%>%
#   select(age:pfs_status)


clinical_all_genomic <- read.csv("prca_clinicogenomics_data.csv")


# clinical_novel_genomic <- read.csv("cleaned_imputed_data_for_ml-integrated.csv")%>%
#   select(-c(fat3_diff, atm_diff, kmt2d_diff, foxa1_diff, tp53_diff, spop_diff,
#             smad4_diff,lrp1b_diff,  idh1_diff, ctnnb1_diff, braf_diff, kmt2c_diff))


# clinical_known_genomic <- read.csv("cleaned_imputed_data_for_ml-integrated.csv")%>%
#   select(-c(nkx3_1_diff, csmd3_diff, trrap_diff, chd4_diff, vwf_diff, ephb1_diff, herc2_diff, mcm3_diff,
#             spta1_diff, sall1_diff,  herc1_diff, rybp_diff, ttn_diff, chd5_diff, myh6_diff))



dat_new <- clinical_all_genomic  #----> update this based on the data options loaded 

#----------------------->Data Modelling<-----------------
#---convert integer to numeric and character to categorical-
integer_cols <- sapply(dat_new, is.integer)
dat_new[integer_cols] <- lapply(dat_new[integer_cols], as.numeric) 
# Convert character columns to categorical
char_cols <- sapply(dat_new, is.character)
dat_new[char_cols] <- lapply(dat_new[char_cols], factor)
#---Arrange data---
endpoints <- dat_new %>%  #we need status be integer & months be  but numeric here
  dplyr::select(pfs_months, pfs_status)
vars <- dat_new %>% 
  dplyr::select(-c(pfs_months, pfs_status))
pfs_status <- endpoints$pfs_status
pfs_months <- endpoints$pfs_months

#---- Partition Data ----
set.seed(1)
comb_dat <- cbind(vars, pfs_months, pfs_status)
sample = sample.split( Y= comb_dat$pfs_status, SplitRatio = 0.70)
train <- subset(comb_dat, sample == TRUE)
test  <- subset(comb_dat, sample == FALSE)
#---Train set --
X_train <- train %>% select(-c(pfs_months, pfs_status))
Y_train <- train %>% select(pfs_months, pfs_status)

#-- Test set --
X_test <- test %>% select(-c(pfs_months, pfs_status))
Y_test <- test %>% select(pfs_months, pfs_status)

#---Process Train Data---
#One-Hot Encode Categorical Columns
catego_columns <- sapply(X_train, is.factor)
categorical_columns <- colnames(X_train)[catego_columns]
train_categorical <- train %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  model.matrix(~ ., data = .) %>%    
  as.data.frame() 
num_columns_train <- sapply(X_train, is.numeric)
numeric_columns_train <- X_train[num_columns_train]
#combine onehot & scaled
X_train_processed <- cbind(numeric_columns_train, train_categorical)

#--- Process Test Data ---
# -- one-hot encode categorical 
test_categorical <- test %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  model.matrix(~., data = ., xlev = attr(train_categorical, "contrasts")) %>%
  as.data.frame()
num_columns_test <- sapply(X_test, is.numeric)
numeric_columns_test <- X_test[num_columns_test]
#combine one-hot & scaled test
X_test_processed <- cbind(numeric_columns_test, test_categorical)

dim(X_train_processed)
dim(X_test_processed) 
dim(Y_train)
dim(Y_test)

processed_train_data <- data.frame(cbind(X_train_processed, Y_train))
processed_test_data <- data.frame(cbind(X_test_processed, Y_test))

X_train_processed <- as.matrix(X_train_processed)
Y_train <- as.matrix(Y_train)
X_test_processed <- as.matrix(X_test_processed)
Y_test <- as.matrix(Y_test)


# Predefined fold assignment 
foldid <- sample(1:5, size = nrow(processed_train_data), replace = TRUE)
 
#================= Penalized COx Model ================
pfs_ytrain <- cbind(time=Y_train[, "pfs_months"], status=Y_train[, "pfs_status"])
pfs_ytest <- cbind(time=Y_test[, "pfs_months"], status=Y_test[, "pfs_status"])

# Run / loop over all aphas
alphas = c(0, 0.25, 0.5, 0.75, 1)
cv_all_alpha <- data.frame(iter      = integer(),
                           c_indexes = numeric(),
                           lambdas   = numeric(),
                           non_zer0_coef = integer(),
                           stringsAsFactors = FALSE)
for(i in alphas){
  cv_0  <- cv.glmnet(x            = X_train_processed,
                     y            = pfs_ytrain,     
                     foldid       = foldid,      
                     alpha        = i,
                     type.measure = "C",
                     family       = "cox",
                     standardize  = TRUE, 
                     relax        = TRUE, 
                     trace.it     = 0)
  
  alpha_parameter = i 
  c_index <- as.numeric(cv_0$cvm[cv_0$index[1,1]])
  lambda = as.numeric(cv_0$lambda.min)
  pred = as.numeric(cv_0$nzero[cv_0$index[1,1]])
  cv_all_alpha = rbind(cv_all_alpha, 
                       data.frame(alpha = alpha_parameter , c_indexes = c_index, lambdas = lambda, non_zero_coef = pred) )
}

print(cv_all_alpha)

#-- Final Fitted Model on Train data with optimal parameters ---
best_cv_metrics=cv_all_alpha[which.max(cv_all_alpha$c_indexes), ]
best_cv_metrics
bestlambda = best_cv_metrics$lambdas
bestalpha = best_cv_metrics$alpha
pcm_clinical_fit <- glmnet(x           = X_train_processed, 
                           y           = pfs_ytrain,
                           family      = "cox",
                           lambda      = bestlambda,  
                           alpha       = bestalpha,
                           intercept   = TRUE,
                           standardize = T)
#---selected variables---
vars_ridge <- coef(pcm_clinical_fit, s = bestlambda); vars_ridge
vars_ridge <- as.matrix(vars_ridge)
ridge_selected <- names(vars_ridge[vars_ridge != 0,])
selected_variables <- as.character(ridge_selected)
cat(selected_variables, sep = "\n")

coefs_clinical_2 <- as.data.frame.matrix((coef(pcm_clinical_fit, s = bestlambda))) %>%
  filter(.[,1] != 0) %>%
  arrange(desc((.[, 1]))) %>%
  round(., 3)

# Replace spaces with periods in variable names
# Adjust variable names to avoid extra backticks
variable_names <- rownames(coefs_clinical_2)
variable_names <- sapply(variable_names, function(x) {
  x <- gsub(" ", ".", x) # Replace spaces with periods
  if (x %in% names(processed_train_data)) {
    x  # Use name directly if it exists in the dataset
  } else if (grepl("[^a-zA-Z0-9_]", x)) {
    paste0("`", x, "`")  # Add backticks only when necessary
  } else {
    x
  }
})

formula_final <- paste("Surv(pfs_months, pfs_status) ~", paste(variable_names, collapse = " + "))
surv_formula22 <- as.formula(formula_final)

# Print the updated formula
print(surv_formula22)
missing_vars <- setdiff(variable_names, names(processed_train_data))
print(missing_vars)
#---> Final Cox model
cox_model_final<-coxph(surv_formula22, data=processed_train_data, x=T, model = T) 
summary(cox_model_final)


base_surv <- survfit(cox_model_final)
summary(base_surv)
  
baseline_df <- data.frame(
  time = base_surv$time,
  n_risk = base_surv$n.risk,
  n_event = base_surv$n.event,
  surv = base_surv$surv,
  cumhaz = base_surv$cumhaz
)

base_surv$cumhaz # The estimated cumulative baseline hazard
base_surv$surv # The estimated baseline survival prob 


# Concordant - Index on Train Data
surv_pred_ridge <- predict(pcm_clinical_fit, newx=X_train_processed, type ="response")
predicted_results<-assess.glmnet(surv_pred_ridge, newy = pfs_ytrain, family = "cox"); predicted_results
k=assess.glmnet(pcm_clinical_fit, newx = X_train_processed, newy =  pfs_ytrain, family = "cox") 
pred_C <- k$C[1]
pred_C

#---Concordant - Index on Test Data ---
surv_pred_ridge <- predict(pcm_clinical_fit, newx=X_test_processed, type ="response")
predicted_results<-assess.glmnet(surv_pred_ridge, newy = pfs_ytest, family = "cox"); predicted_results
k=assess.glmnet(pcm_clinical_fit, newx = X_test_processed, newy =  pfs_ytest, family = "cox") ##test_surv_response
pred_C <- k$C[1]
pred_C


#----Do the prediction many times to find the average c_index------------
#------Perform this 10 times with the chosen gamma----------
# results1 <- data.frame(iter = numeric(), train_concordance = numeric(), predictors = numeric(), 
#                        se = numeric(), pred_concordance = nuemric() )

iterations <- 10
results_new <- data.frame()
for (iter in 1:iterations){
  
  #-----------Training period------------------
  cv_check <- cv.glmnet(x = X_train_processed,  #x_train,
                        y = pfs_ytrain, #train_surv_response,
                        #foldid = foldid, #gave constant preds
                        nfolds = 5,
                        alpha = bestalpha,
                        type.measure = "C",
                        family = "cox",
                        standardize = TRUE ) #default
  alpha = bestalpha
  BEST_LAMBDA <- round(cv_check$lambda.min, 4)
  variables <-as.numeric(cv_check$nzero[cv_check$index[1,1]])
  BEST_CONCORDANCE <- round(cv_check$cvm[as.numeric(cv_check$index[1,1])], 4)
  SE = round(cv_check$cvsd[as.numeric(cv_check$index[1,1])], 4)
  
  fit_r <- glmnet(x           = X_train_processed, 
                  y           = pfs_ytrain, 
                  family      = "cox",
                  lambda      = BEST_LAMBDA,
                  alpha       = bestalpha,
                  standardize = T)
  
  lambdas_new = round(fit_r$lambda, 4) #fitted
  
  #---------------------  Predict on test data -------------------------
  surv_pred_clinical <- predict(fit_r, newx = X_test_processed, type = "link") #predictions 
  k <- round(as.numeric(assess.glmnet(surv_pred_clinical, newy = pfs_ytest, family = "cox")$C), 4) 
  
  results_new <- rbind(results_new,  
                       data.frame(iteration = iter, penalty = alpha, sel_lambda = BEST_LAMBDA, 
                                  fitted_lambda = lambdas_new, sel_vars = variables, std.error = SE, 
                                  train_concordance = BEST_CONCORDANCE, Predicted_Concordance = k)) 
}

results_pred_clinical <- results_new
average_pred<-mean(results_pred_clinical$Predicted_Concordance)
cat("The average train concordance after ",iter," runs is: ",mean(results_pred_clinical$train_concordance))
cat("The sd train concordance after ",iter," runs is: ",sd(results_pred_clinical$train_concordance))
cat("The average predicted concordance after ",iter," runs is: ",mean(results_pred_clinical$Predicted_Concordance))
cat("The sd predicted concordance after ",iter," runs is: ",sd(results_pred_clinical$Predicted_Concordance))

#--Survival Predictions---predict at time = 72 months
library(survex) 
#works for selected models coxph, rsf but dont see neural network models
explainer_cox <- explain(cox_model_final)
#------survival predictions--------
req_variable_names <- sapply(variable_names, function(x) {
  x <- gsub(" ", ".", x) # Replace spaces with periods
  if (x %in% names(processed_test_data)) {
    x  # Use name directly if it exists in the dataset
  } else if (grepl("[^a-zA-Z0-9_]", x)) {
    paste0("`", x, "`")  # Add backticks only when necessary
  } else {
    x
  }
})

X_test_variables <- unique(req_variable_names)
X_test_pred <- as.data.frame(X_test_processed)

# Ensure column names are formatted consistently (replace spaces with dots)
colnames(X_test_pred) <- gsub(" ", ".", colnames(X_test_pred))

# Subset only the required columns dynamically
X_test_pred_req <- X_test_pred[, X_test_variables, drop = FALSE]

pred_cox_surv <- predict(explainer_cox, newdata = X_test_pred_req, output_type="survival", times = 72) #max(test_data$pfs_months)
pred_cox_surv[c(2,10),]
mean(pred_cox_surv)
median(pred_cox_surv)



#=====================Fivenumber Summary PCM====================
pred_df <- tibble(surv_prob = as.vector(pred_cox_surv))

# Five-number summary + extras
fivenum_pcm_surv <- pred_df %>%
  summarise(
    n = n(),
    mean = mean(surv_prob, na.rm = TRUE),
    stdev = sd(surv_prob, na.rm = TRUE),
    min = min(surv_prob, na.rm = TRUE),
    q1 = quantile(surv_prob, 0.25, na.rm = TRUE),
    median = median(surv_prob, na.rm = TRUE),
    q3 = quantile(surv_prob, 0.75, na.rm = TRUE),
    max = max(surv_prob, na.rm = TRUE)
  )

print(fivenum_pcm_surv)



#risk--> does not matter time
pred_cox_risk <- predict(explainer_cox, newdata = X_test_pred_req, output_type="risk")
pred_cox_risk[c(2,10)]
mean(pred_cox_risk)
median(pred_cox_risk)

predictions_risk_matrix_cox <- cbind(X_test_pred_req, pred_cox_surv, pred_cox_risk) 
predictions_risk_matrix_cox

# write.csv(predictions_risk_matrix_cox, 
#           file="C:/Users/dumbl/OneDrive/Desktop/ijerp-code/Analysis/pcm_prob_risk_predictions.csv")


#=====================Fivenumber Summary PCM====================
pred_df <- tibble(surv_prob = as.vector(pred_cox_risk))

# Five-number summary + extras
fivenum_pcm_risk <- pred_df %>%
  summarise(
    n = n(),
    mean = mean(surv_prob, na.rm = TRUE),
    stdev = sd(surv_prob, na.rm = TRUE),
    min = min(surv_prob, na.rm = TRUE),
    q1 = quantile(surv_prob, 0.25, na.rm = TRUE),
    median = median(surv_prob, na.rm = TRUE),
    q3 = quantile(surv_prob, 0.75, na.rm = TRUE),
    max = max(surv_prob, na.rm = TRUE)
  )

print(fivenum_pcm_risk)







#=======================================  Random Survival Forest =================================
#Prep Data for RSF ---> train & test are dataframes
#===Train set====
X_train <- train %>% select(-c(pfs_months, pfs_status))
Y_train <- train %>% select(pfs_months, pfs_status)

#===Test set===
X_test <- test %>% select(-c(pfs_months, pfs_status))
Y_test <- test %>% select(pfs_months, pfs_status)

#=====Process Train Data===
#One-Hot Encode Categorical Columns
catego_columns <- sapply(X_train, is.factor)
categorical_columns <- colnames(X_train)[catego_columns]

train_categorical <- train %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  model.matrix(~ . -1, data = .) %>%    #Remove intercept term
  as.data.frame() %>%
  select(-hist_neoadjuv_trtmntNo)  #Remove redundant term
num_columns_train <- sapply(X_train, is.numeric)
numeric_columns_train <- X_train[num_columns_train]
#combine onehot & scaled
X_train_processed <- cbind(numeric_columns_train, train_categorical)

#===Process Test Data===
#one-hot encode categorical 
test_categorical <- test %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  model.matrix(~.-1, data = ., xlev = attr(train_categorical, "contrasts")) %>%
  as.data.frame() %>%
  select(-hist_neoadjuv_trtmntNo) #remove redundnt level
num_columns_test <- sapply(X_test, is.numeric)
numeric_columns_test <- X_test[num_columns_test]
#combine one-hot & scaled test
X_test_processed <- cbind(numeric_columns_test, test_categorical)

dim(X_train_processed)
dim(X_test_processed) 
dim(Y_train)
dim(Y_test)
processed_train_data <- data.frame(cbind(X_train_processed, Y_train))
processed_test_data <- data.frame(cbind(X_test_processed, Y_test))
X_train_processed <- as.matrix(X_train_processed)
Y_train <- as.matrix(Y_train)
X_test_processed <- as.matrix(X_test_processed)
Y_test <- as.matrix(Y_test)


results <- data.frame(nodesize = numeric(),  mtry = numeric(), ntree = numeric(), error = numeric())
trees <- c(seq(10,100, by=10), 150, 200, 250, 300) 
for (tree_val in trees) {
  for (fold in unique(foldid)) {
    train_data <- processed_train_data[foldid != fold, ]    
    test_data <- processed_train_data[foldid == fold, ]
    # Tune a random survival forest model
    tune0 <- tune(formula     = Surv(pfs_months, pfs_status) ~ .,
                  data        = train_data,
                  mtryStart   = 1,
                  nodesizeTry = c(1:9, seq(10, 40, by = 5)),
                  ntreeTry    = tree_val
    )
    
    opt          <- tune0$results[which.min(tune0$results[,3]), ] #optimal info
    mtry_val     <- as.numeric(opt[2])   #fetch optimal mtry
    nodesize_val <- as.numeric(opt[1])   # fetch optimal nodesize
    ntree_val    <- as.numeric(tree_val) # fetch optimal tree
    error_val    <- as.numeric(opt[3])   # fetch optimal cv error
    
    results <- rbind(results, data.frame(mtry = mtry_val, nodesize = nodesize_val, ntree = tree_val, error = error_val))
  }
}
rsf_clin_best_cv_results <- results
best_parameters_clin <- results[which.min(rsf_clin_best_cv_results$error), ]
best_parameters_clin

#=====================================================================
# ===============BUILD BEST RSF MODEL===========================
#==================================================================
rsf_clinical_fit <- rfsrc(formula        = Surv(pfs_months, pfs_status) ~ .,
                          data           = processed_train_data,
                          perf.type      = "default",
                          mtry           = best_parameters_clin$mtry,    
                          nodesize       = best_parameters_clin$nodesize,   
                          ntree          = best_parameters_clin$ntree,   
                          importance     = TRUE,
                          forest         = TRUE
                          )
rsf_clinical_fit

#==Variable Importance Minimal Depth
rsf_vars<-var.select(rsf_clinical_fit, method = "md") #good 
rsf_vars$topvars


#--survex---Survival Predictions---
#predict at time = 72 months
library(survex) 
#works for selected models coxph, rsf but dont see neural network models
#Build explainer model with only selected variables
rsf_vars_req<- rsf_vars$topvars
X_train_pred <- as.data.frame(X_train_processed)
colnames(X_train_pred) <- gsub(" ", ".", colnames(X_train_pred))
X_train_pred_req <- X_train_pred[, rsf_vars_req] 

#Final rsf model with selected variables
rsf_final_data <- cbind(X_train_pred_req, Y_train)
rsf_clinical_fit_final <- rfsrc(formula        = Surv(pfs_months, pfs_status) ~ .,
                                data           = rsf_final_data,
                                perf.type      = "default",
                                mtry           = best_parameters_clin$mtry,    
                                nodesize       = best_parameters_clin$nodesize,   
                                ntree          = best_parameters_clin$ntree,   
                                importance     = TRUE,
                                forest         = TRUE
)

#==Explainer Created==
explainer_rsf <- explain(rsf_clinical_fit_final, data = X_train_pred_req)

#Build test data-matrix & Predict probabilities
X_test_pred_rsf <- as.data.frame(X_test_processed)
colnames(X_test_pred_rsf) <- gsub(" ", ".", colnames(X_test_pred_rsf))
X_test_pred_rsf_req <- X_test_pred_rsf[, rsf_vars_req]

pred_rsf <- predict(explainer_rsf, newdata = X_test_pred_rsf_req, output_type="survival", times = 72 ) #max(test_data$pfs_months)
pred_rsf[c(2,10),]
mean(pred_rsf)
median(pred_rsf)

#=====================Fivenumber Summary RSF-surv====================
pred_df <- tibble(surv_prob = as.vector(pred_rsf))

# Five-number summary + extras
fivenum_rsf_surv <- pred_df %>%
  summarise(
    n = n(),
    mean = mean(surv_prob, na.rm = TRUE),
    stdev = sd(surv_prob, na.rm = TRUE),
    min = min(surv_prob, na.rm = TRUE),
    q1 = quantile(surv_prob, 0.25, na.rm = TRUE),
    median = median(surv_prob, na.rm = TRUE),
    q3 = quantile(surv_prob, 0.75, na.rm = TRUE),
    max = max(surv_prob, na.rm = TRUE)
  )

print(fivenum_rsf_surv)




#risk--> does not matter time
pred_rsf_risk <- predict(explainer_rsf, newdata = X_test_pred_rsf_req, output_type="risk")
pred_rsf_risk[c(2,10)]
mean(pred_rsf_risk)
median(pred_rsf_risk)



predictions_risk_matrix_rsf <- cbind(X_test_pred_rsf_req, pred_rsf, pred_rsf_risk) 
predictions_risk_matrix_rsf
# 
# write.csv(predictions_risk_matrix_rsf, 
#           file="C:/Users/dumbl/OneDrive/Desktop/ijerp-code/Analysis/rsf_prob_risk_predictions.csv")

#=====================Fivenumber Summary RSF-risk====================
pred_df <- tibble(surv_prob = as.vector(pred_rsf_risk))

# Five-number summary + extras
fivenum_rsf_risk <- pred_df %>%
  summarise(
    n = n(),
    mean = mean(surv_prob, na.rm = TRUE),
    stdev = sd(surv_prob, na.rm = TRUE),
    min = min(surv_prob, na.rm = TRUE),
    q1 = quantile(surv_prob, 0.25, na.rm = TRUE),
    median = median(surv_prob, na.rm = TRUE),
    q3 = quantile(surv_prob, 0.75, na.rm = TRUE),
    max = max(surv_prob, na.rm = TRUE)
  )

print(fivenum_rsf_risk)





#=================================PREDICTIONS========================================
# randomForestSRC package  = can only predict the survival function at times which are present in the training dataset.---- rsf_clinical_fit$time.interest

#---------Prediction accuracy of best model on train data---------------
pred_rsf_clinical_fit_train <- predict(object = rsf_clinical_fit, newdata = processed_train_data)
pred_rsf_clinical_fit_train
# 1 - requested performance error = accuracy index
pred_rsf_clinical_fit_c_index_train <- 1 - get.cindex(time = pred_rsf_clinical_fit_train$yvar[,1], censoring =pred_rsf_clinical_fit_train$yvar[,2],
                                                      predicted = pred_rsf_clinical_fit_train$predicted)
pred_rsf_clinical_fit_c_index_train


#---------Prediction accuracy of best model on test data---------------
pred_rsf_clinical_fit <- predict(object = rsf_clinical_fit, newdata = processed_test_data)
pred_rsf_clinical_fit
# 1 - requested performance error = accuracy index
pred_rsf_clinical_fit_c_index2 <- 1 - get.cindex(time = pred_rsf_clinical_fit$yvar[,1], censoring =pred_rsf_clinical_fit$yvar[,2],
                                                 predicted = pred_rsf_clinical_fit$predicted)
pred_rsf_clinical_fit_c_index2


#-------------------------Predict on test set for 10 iterations------------------------------
# Initialize an empty data frame to store results
results_all <- data.frame(iteration = numeric(),
                          nodesize_best = numeric(),
                          mtry_best = numeric(),
                          ntree_best = numeric(),
                          train_concordance = numeric(),
                          pred_concordance = numeric())

iterations <- 10
for (iter in 1:iterations) {
  # Initialize results data frame
  par_results <- data.frame(nodesize = numeric(), mtry = numeric(), ntree = numeric(), train_error = numeric())
  trees <- c(seq(10,100, by=10), 150, 200) 
  for (tree_val in trees) {
    for (fold in unique(foldid)) {
      train_data <- processed_train_data[foldid != fold, ]
      test_data <- processed_train_data[foldid == fold, ]
      tune0 <- tune(formula = Surv(pfs_months, pfs_status) ~ .,
                    data = train_data, 
                    mtryStart = 1,
                    nodesizeTry = c(1:9, seq(10, 40, by = 5)),
                    ntreeTry = tree_val)
      opt <- tune0$results[which.min(tune0$results[, 3]), ]
      mtry_val <- as.numeric(opt[2])
      nodesize_val <- as.numeric(opt[1])
      ntree_val <- as.numeric(tree_val)
      error_val <- as.numeric(opt[3])
      par_results <- rbind(par_results, 
                           data.frame(mtry = mtry_val, nodesize = nodesize_val, ntree = tree_val, train_error = error_val)) 
    }
  }
  # Optimal parameters 
  best_par <- par_results[which.min(par_results$train_error), ]
  
  #---final model with tuned hyper-parameters------
  final_model <- rfsrc(formula = Surv(pfs_months, pfs_status) ~ ., 
                       data = processed_train_data,
                       mtry = best_par$mtry,
                       nodesize = best_par$nodesize,
                       ntree = best_par$ntree)
  # Predict on test data - model never seen & train data
  rsf_pred <- predict(object = final_model, newdata = processed_test_data)
  rsf_pred_train <- predict(object = final_model, newdata = processed_train_data)
  
  train_pred <- rsf_pred_train$predicted
  test_pred <- rsf_pred$predicted
  c_index <- get.cindex(time = processed_test_data$pfs_months, censoring = processed_test_data$pfs_status, 
                        predicted = rsf_pred$predicted)
  c_index_train <- get.cindex(time = processed_train_data$pfs_months, censoring = processed_train_data$pfs_status, 
                              predicted = rsf_pred_train$predicted)
  
  # Append results to results_all
  results_all <- rbind(results_all, data.frame(iteration = iter, 
                                               nodesize_best = best_par$nodesize, 
                                               mtry_best = best_par$mtry, 
                                               ntree_best = best_par$ntree,
                                               train_concordance = 1 - c_index_train,  #1-best_par$train_error,
                                               pred_concordance = 1 - c_index))
}
rsf_pred_results_clinical <- results_all
average_pred<-mean(rsf_pred_results_clinical$pred_concordance)
cat("The average train concordance after ",iter," runs is: ",  mean(rsf_pred_results_clinical$train_concordance))
cat("The sd of train concordance after ",iter," runs is: ",  sd(rsf_pred_results_clinical$train_concordance))
cat("The average test concordance after ",iter," runs is: ",  mean(rsf_pred_results_clinical$pred_concordance))
cat("The sd of test concordance after ",iter," runs is: ",  sd(rsf_pred_results_clinical$pred_concordance))

# this code implements a robust method for tuning the hyperparameters of a random survival forest model using cross-validation
# and evaluating its performance on a separate test set. The repeated iterations ensure that the model's performance is 
# adequately assessed and helps in obtaining stable estimates of the model's predictive performance
#==================================================================================================================================









#================DeepSurv ===================================================
#===Cross Validation==
source("deepsurv_libraries.R")
source("deepsurv.R")
tensorflow::set_random_seed(1)

#Prep Data for Deepsurv Jinli ---> train & test are dataframes
#==Train set ==
X_train <- train %>% select(-c(pfs_months, pfs_status))
Y_train <- train %>% select(pfs_months, pfs_status)
#==Test set ==
X_test <- test %>% select(-c(pfs_months, pfs_status))
Y_test <- test %>% select(pfs_months, pfs_status)

#======== Process Training data ===========
num_columns <- sapply(X_train, is.numeric)
numeric_columns <- colnames(X_train)[num_columns]
#Min-Max Scaling parameters of training data
scale_params <- train %>%
  select(all_of(numeric_columns)) %>%
  summarise(across(everything(), list(min=min, max=max)))
#Apply Min-Max scaling on training data
X_train[numeric_columns] <- train %>%
  select(all_of(numeric_columns)) %>%
  mutate(across(everything(), ~ (. - scale_params[[paste0(cur_column(), "_min")]]) / 
                  (scale_params[[paste0(cur_column(), "_max")]] - scale_params[[paste0(cur_column(),"_min")]])
  ))

#One-Hot Encode Categorical Columns
catego_columns <- sapply(X_train, is.factor)
categorical_columns <- colnames(X_train)[catego_columns]
#library(glmnet)
train_categorical <- train %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  { makeX(
    train = .,
    contrasts.arg = lapply(., function(x) contrasts(x, contrasts = FALSE)))
  } %>%as.data.frame()

#combine onehot & scaled
X_train_processed <- cbind(X_train[numeric_columns] , train_categorical)


#== Apply Transformation and process Test Data ===
#min-max scale on test data
X_test[numeric_columns] <- test %>%
  select(all_of(numeric_columns)) %>%
  mutate(across(everything(),
                ~ (. - scale_params[[paste0(cur_column(), "_min")]]) /
                  (scale_params[[paste0(cur_column(), "_max")]] - scale_params[[paste0(cur_column(), "_min")]])
  ))
#one-hot encode categorical 
test_categorical <- test %>%
  select(all_of(categorical_columns)) %>%
  mutate(across(everything(), as.factor)) %>%
  { makeX(
    train = .,
    contrasts.arg = lapply(., function(x) contrasts(x, contrasts = FALSE))
  )
  } %>%
  as.data.frame()

#combine one-hot & scaled test
X_test_processed <- cbind(test[numeric_columns], test_categorical)

dim(X_train_processed)
dim(X_test_processed) 
dim(Y_train)
dim(Y_test)
processed_train_data <- data.frame(cbind(X_train_processed, Y_train))
processed_test_data <- data.frame(cbind(X_test_processed, Y_test))
X_train_processed <- as.matrix(X_train_processed)
Y_train <- as.matrix(Y_train)
X_test_processed <- as.matrix(X_test_processed)
Y_test <- as.matrix(Y_test)

hyperparams <- expand.grid(
  num_layer = c(1),
  num_nodes = c(1),
  num_dropout = c(0.2, 0.3, 0.4),
  num_lr = c(0.0005, 0.001),
  batch_size = c(16),
  num_l2 = c(0.2,0.4,0.6)
)

# Perform K-Fold Cross-Validation for Grid Search
k_folds <- 5
folds <- sample(rep(1:k_folds, length.out = nrow(processed_train_data)))

# Grid Search with Cross-Validation
results <- data.frame(hyperparams)
results$Mean_C_Index <- NA
results$Std_C_Index <- NA

for (i in 1:nrow(hyperparams)) {
  params <- as.list(hyperparams[i, ])
  c_index_scores <- c()
  
  for (fold in 1:k_folds) {
    # Train/Validation Split
    train_indices <- which(folds != fold)
    val_indices <- which(folds == fold)
    
    X_train_fold <- X_train_processed[train_indices, ]
    Y_train_fold <- Y_train[train_indices, ]
    X_val_fold <- X_train_processed[val_indices, ]
    Y_val_fold <- Y_train[val_indices, ]
    
    # Build and Train Model
    deepsurv_model <- build_deepsurv(
      num_input = ncol(X_train_fold),
      num_layer = params$num_layer,
      num_nodes = params$num_nodes,
      string_activation = "selu",
      num_l2 = params$num_l2,
      num_dropout = params$num_dropout,
      num_lr = params$num_lr,
      lr_decay = 1e-2
    )
    
    # Add an early stopping callback
    early_stopping <- callback_early_stopping(
      monitor = "val_loss",       # Monitor validation loss
      patience = 20,              # Number of epochs with no improvement before stopping
      restore_best_weights = TRUE # Restore weights from the best epoch
    )
    
    history <- deepsurv_model %>% fit(
      x = X_train_fold,
      y = Y_train_fold,
      batch_size = params$batch_size,
      epochs = 100,
      validation_data = list(X_val_fold, Y_val_fold),
      verbose = 0,
      callbacks = list(early_stopping)
    )
    
    # Evaluate on Validation Set
    predictions <- deepsurv_model %>% predict(X_val_fold)
    c_index_scores <- c(c_index_scores, c_index(predictions, Y_val_fold[, 2], Y_val_fold[, 1]))
  }
  
  # Store Mean and Standard Deviation of C-Index Across Folds
  results$Mean_C_Index[i] <- mean(c_index_scores, na.rm = TRUE)
  results$Std_C_Index[i] <- sd(c_index_scores, na.rm = TRUE)
}

# Select Best Hyperparameters
best_params <- results[which.max(results$Mean_C_Index), ]
print(best_params)

# Train the final model with the best parameters
cat("Training final model with best hyperparameters...\n")
final_model <- build_deepsurv(
  num_input = ncol(X_train_processed),
  num_layer = best_params$num_layer,
  num_nodes = best_params$num_nodes,
  string_activation = "selu",
  num_l2 = best_params$num_l2,
  num_dropout = best_params$num_dropout,
  num_lr = best_params$num_lr,
  lr_decay = 1e-2)

print(final_model)


# Train the model using the entire training dataset
final_history <- final_model %>% fit(
  x = X_train_processed,
  y = Y_train,
  batch_size = best_params$batch_size,
  epochs = 100,
  validation_split = 0.25,  # Keep 25% of the training data for internal validation
  verbose = 1,
  callbacks =list(
    callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.5, patience = 5, min_lr = 1e-6),
    callback_early_stopping(monitor = "loss", patience = 10, restore_best_weights = TRUE)
  )
)


# Evaluate on the training set
cat("Evaluating final model on training set...\n")
final_train_predictions <- final_model %>% predict(X_train_processed)
final_train_c_index <- c_index(
  final_train_predictions, 
  Y_train[, "pfs_status"], 
  Y_train[, "pfs_months"]
)
cat("Final Train C-Index:", final_train_c_index, "\n")

# Evaluate on the test set
cat("Evaluating final model on test set...\n")
final_test_predictions <- final_model %>% predict(X_test_processed)
final_test_c_index <- c_index(
  final_test_predictions, 
  Y_test[, "pfs_status"], 
  Y_test[, "pfs_months"]
)
cat("Final Test C-Index:", final_test_c_index, "\n")

# Save the results
results_summary <- data.frame(
  Metric = c("Train C-Index", "Test C-Index"),
  Value = c(final_train_c_index, final_test_c_index)
)

valid_results <- subset(results, !is.na(Std_C_Index))
sorted_results <- valid_results[order(-valid_results$Mean_C_Index, valid_results$Std_C_Index), ]

# Fetch the top row
best_result <- sorted_results[1, ]

# Print the best hyperparameter set
cat("Best hyperparameter set with highest Mean_C_Index and lowest Std_C_Index:\n")
cat("Parameter set with low variability across folds:\n")
print(best_result)





#variable importance
variable_importance <- get_variable_importance(final_model, X_train_processed)
print(variable_importance)

##Plot Feature Attribution 
importance_df <- variable_importance[order(variable_importance$Importance, decreasing = T), ]
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.7) +
  coord_flip() +  # Flip axes to make feature names readable
  labs(title = "Clinical and All Genomic Features Attribution in DeepSurv Model",
       x = "Features",
       y = "Importance Score") +
  scale_fill_brewer(palette = "Set1") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) + theme_minimal()

#===Top feature from GBFA PLOT==
median_sensitivity <- round(median(importance_df$Importance), 4)
shortlisted_GBFA <- importance_df[importance_df$Importance > median_sensitivity,]
shortlisted_GBFA

# Plot histogram of shortlisted importance
ggplot(shortlisted_GBFA, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  coord_flip() +  # Flip axes to make feature names readable
  labs(title = "Clinical Features Attribution in DeepSurv Model above Median",
       x = "Features",
       y = "Importance Score") +
  scale_fill_brewer(palette = "Set1") +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) + theme_minimal()

#======Deepsurv Risk & Survival Predictions=======
median_sensitivity <- round(median(importance_df$Importance), 4)
shortlisted <- importance_df[importance_df$Importance > median_sensitivity,]
shortlisted_features <- shortlisted$Feature
shortlist_train <- X_train_processed[,shortlisted_features]
shortlist_test <- X_test_processed[, shortlisted_features]


final_model1 <- build_deepsurv(
  num_input = ncol(shortlist_train),
  num_layer = best_params$num_layer,
  num_nodes = best_params$num_nodes,
  string_activation = "selu",
  num_l2 = best_params$num_l2,
  num_dropout = best_params$num_dropout,
  num_lr = best_params$num_lr,
  lr_decay = 1e-2
)

predicted_log_risk_deepsurv <- final_model1 %>% predict(shortlist_test); predicted_log_risk_deepsurv
predicted_relative_risk_deepsurv <- exp(predicted_log_risk_deepsurv); predicted_relative_risk_deepsurv

mean(predicted_relative_risk_deepsurv)
median(predicted_relative_risk_deepsurv)
cat("Predicted Risk Score for Patient 2:", predicted_relative_risk_deepsurv[2], "\n")
cat("Predicted Risk Score for Patient 10:", predicted_relative_risk_deepsurv[10], "\n") 
deepsurv_risks <- data.frame(shortlist_test, predicted_log_risk_deepsurv, predicted_relative_risk_deepsurv) 
deepsurv_risks

# write.csv(deepsurv_risks, 
#           file="C:/Users/dumbl/OneDrive/Desktop/ijerp-code/Analysis/deepsurv_risks_predictions.csv")


#=====================Fivenumber Summary DeepSurv-risk====================
pred_df <- tibble(surv_risk = as.vector(predicted_relative_risk_deepsurv))

# Five-number summary + extras
fivenum_deepsurv_risk <- pred_df %>%
  summarise(
    n = n(),
    mean = mean(surv_risk, na.rm = TRUE),
    stdev = sd(surv_risk, na.rm = TRUE),
    min = min(surv_risk, na.rm = TRUE),
    q1 = quantile(surv_risk, 0.25, na.rm = TRUE),
    median = median(surv_risk, na.rm = TRUE),
    q3 = quantile(surv_risk, 0.75, na.rm = TRUE),
    max = max(surv_risk, na.rm = TRUE)
  )

print(fivenum_deepsurv_risk)



#======== Repeat 3 TIMES==============
# Initialize DataFrame to Store Results
results_summary <- data.frame(
  Run = integer(),
  Final_Train_C_Index = numeric(),
  Final_Test_C_Index = numeric()
)

# Number of repetitions
num_repeats <- 5

for (run in 1:num_repeats) {
  cat("\n===== RUN", run, "=====\n")
  
  # Set seed for reproducibility
  set.seed(123 + run )
  tensorflow::set_random_seed(123 + run)
  
  # Define the hyperparameter grid for this run
  hyperparams <- expand.grid(
    num_layer = c(1),
    num_nodes = c(1),
    num_dropout = c(0.2, 0.3, 0.4),
    num_lr = c(0.0005, 0.001),
    batch_size = c(16),
    num_l2 = c(0.2, 0.4, 0.6)
  )
  
  # Perform K-Fold Cross-Validation for Grid Search
  k_folds <- 5
  folds <- createFolds(Y_train[, "pfs_status"], k = k_folds, list = FALSE)
  
  # Grid Search with Cross-Validation
  results <- data.frame(hyperparams)
  results$Mean_C_Index <- NA
  results$Std_C_Index <- NA
  
  for (i in 1:nrow(hyperparams)) {
    params <- as.list(hyperparams[i, ])
    c_index_scores <- c()
    
    for (fold in 1:k_folds) {
      # Train/Validation Split
      train_indices <- which(folds != fold)
      val_indices <- which(folds == fold)
      
      X_train_fold <- X_train_processed[train_indices, ]
      Y_train_fold <- Y_train[train_indices, ]
      X_val_fold <- X_train_processed[val_indices, ]
      Y_val_fold <- Y_train[val_indices, ]
      
      # Build and Train Model
      deepsurv_model <- build_deepsurv(
        num_input = ncol(X_train_fold),
        num_layer = params$num_layer,
        num_nodes = params$num_nodes,
        string_activation = "selu",
        num_l2 = params$num_l2,
        num_dropout = params$num_dropout,
        num_lr = params$num_lr,
        lr_decay = 1e-2
      )
      history <- deepsurv_model %>% fit(
        x = X_train_fold,
        y = Y_train_fold,
        batch_size = params$batch_size,
        epochs = 100,
        validation_data = list(X_val_fold, Y_val_fold),
        verbose = 0
      )
      
      # Evaluate on Validation Set
      predictions <- deepsurv_model %>% predict(X_val_fold)
      c_index_scores <- c(c_index_scores, c_index(predictions, Y_val_fold[, 2], Y_val_fold[, 1]))
    }
    
    # Store Mean and Standard Deviation of C-Index Across Folds
    results$Mean_C_Index[i] <- mean(c_index_scores, na.rm = TRUE)
    results$Std_C_Index[i] <- sd(c_index_scores, na.rm = TRUE)
  }
  
  # Select Best Hyperparameters
  best_params <- results[which.max(results$Mean_C_Index), ]
  print(best_params)
  
  # Train the final model with the best parameters
  cat("Training final model with best hyperparameters...\n")
  final_model <- build_deepsurv(
    num_input = ncol(X_train_processed),
    num_layer = best_params$num_layer,
    num_nodes = best_params$num_nodes,
    string_activation = "selu",
    num_l2 = best_params$num_l2,
    num_dropout = best_params$num_dropout,
    num_lr = best_params$num_lr,
    lr_decay = 1e-2
  )
  
  # Train the model using the entire training dataset
  final_history <- final_model %>% fit(
    x = X_train_processed,
    y = Y_train,
    batch_size = best_params$batch_size,
    epochs = 100,
    validation_split = 0.35,  
    verbose = 1,
    callbacks = list(
      callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.5, patience = 5, min_lr = 1e-6),
      callback_early_stopping(monitor = "loss", patience = 10, restore_best_weights = TRUE)
    )
  )
  
  # Evaluate on the training set
  cat("Evaluating final model on training set...\n")
  final_train_predictions <- final_model %>% predict(X_train_processed)
  final_train_c_index <- c_index(
    final_train_predictions, 
    Y_train[, "pfs_status"], 
    Y_train[, "pfs_months"]
  )
  cat("Final Train C-Index:", final_train_c_index, "\n")
  
  # Evaluate on the test set
  cat("Evaluating final model on test set...\n")
  final_test_predictions <- final_model %>% predict(X_test_processed)
  final_test_c_index <- c_index(
    final_test_predictions, 
    Y_test[, "pfs_status"], 
    Y_test[, "pfs_months"]
  )
  cat("Final Test C-Index:", final_test_c_index, "\n")
  # Store results in the summary dataframe
  results_summary <- rbind(
    results_summary,
    data.frame(Run = run, 
               Final_Train_C_Index = final_train_c_index, 
               Final_Test_C_Index = final_test_c_index)
  )
}

print(results_summary)








#========================new===================

###############################################################################
# DeepSurv MAIN SCRIPT (Single-Run Grid Search + Final Fit + Evaluation)
# + Optional MULTI-RUN wrapper (if you later change your mind)
#
# REQUIREMENTS:
#   - You already have `train` and `test` data.frames in memory
#   - Both include: pfs_months (time), pfs_status (event: 1=event, 0=censored)
#   - You have saved the rewritten deepsurv.R in your working directory
#
# OUTPUTS:
#   SINGLE RUN:
#     - CV table across hyperparameter grid
#     - Best hyperparameters (by highest mean CV C-index; tie-break low SD)
#     - Final fitted model evaluated on TRAIN and TEST (C-index)
#     - Variable importance (GBFA) + plots
#     - Risk score predictions on TEST + five-number summary
#
#   MULTI-RUN (optional section at bottom):
#     - Repeats tuning+final evaluation over multiple seeds
#     - Summary mean±SD for Train/Test C-index
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(caret)
  library(ggplot2)
  library(tibble)
  library(keras)
  library(tensorflow)
})

#---------------------------#
# 0) Source your fixed DeepSurv functions
#---------------------------#
source("deepsurv.R")

#---------------------------#
# 1) Sanity checks
#---------------------------#
stopifnot(is.data.frame(train), is.data.frame(test))
stopifnot(all(c("pfs_months", "pfs_status") %in% names(train)))
stopifnot(all(c("pfs_months", "pfs_status") %in% names(test)))

train <- train %>% mutate(pfs_months = as.numeric(pfs_months),
                          pfs_status = as.integer(pfs_status))
test  <- test  %>% mutate(pfs_months = as.numeric(pfs_months),
                          pfs_status = as.integer(pfs_status))

if (!all(train$pfs_status %in% c(0, 1)) || !all(test$pfs_status %in% c(0, 1))) {
  stop("pfs_status must be coded 0/1 (0=censored, 1=event).")
}

#---------------------------#
# 2) Preprocessing: min-max numeric + one-hot categorical (train-anchored)
#---------------------------#
prep_design_matrices <- function(train_df, test_df,
                                 time_col = "pfs_months",
                                 event_col = "pfs_status") {
  
  X_train_raw <- train_df %>% select(-all_of(c(time_col, event_col)))
  X_test_raw  <- test_df  %>% select(-all_of(c(time_col, event_col)))
  
  # Coerce character -> factor
  X_train_raw <- X_train_raw %>% mutate(across(where(is.character), as.factor))
  X_test_raw  <- X_test_raw  %>% mutate(across(where(is.character), as.factor))
  
  is_num <- sapply(X_train_raw, is.numeric)
  num_cols <- names(X_train_raw)[is_num]
  cat_cols <- names(X_train_raw)[!is_num]
  
  # ---- Min-max scale numeric using TRAIN params only ----
  if (length(num_cols) > 0) {
    mins <- sapply(X_train_raw[num_cols], min, na.rm = TRUE)
    maxs <- sapply(X_train_raw[num_cols], max, na.rm = TRUE)
    ranges <- maxs - mins
    ranges[ranges == 0] <- 1
    
    X_train_num <- as.data.frame(Map(function(x, mn, rg) (x - mn) / rg,
                                     X_train_raw[num_cols], mins, ranges))
    X_test_num  <- as.data.frame(Map(function(x, mn, rg) (x - mn) / rg,
                                     X_test_raw[num_cols], mins, ranges))
    names(X_train_num) <- num_cols
    names(X_test_num)  <- num_cols
  } else {
    X_train_num <- data.frame()
    X_test_num  <- data.frame()
  }
  
  # ---- One-hot encode categorical using TRAIN levels ----
  if (length(cat_cols) > 0) {
    for (cc in cat_cols) {
      X_train_raw[[cc]] <- as.factor(X_train_raw[[cc]])
      X_test_raw[[cc]]  <- as.factor(X_test_raw[[cc]])
      levels(X_test_raw[[cc]]) <- levels(X_train_raw[[cc]])
    }
    
    fml <- as.formula(paste("~", paste(cat_cols, collapse = " + "), "- 1"))
    mm_train <- model.matrix(fml, data = X_train_raw)
    mm_test  <- model.matrix(fml, data = X_test_raw)
    
    X_train_cat <- as.data.frame(mm_train)
    X_test_cat  <- as.data.frame(mm_test)
  } else {
    X_train_cat <- data.frame()
    X_test_cat  <- data.frame()
  }
  
  X_train <- cbind(X_train_num, X_train_cat)
  X_test  <- cbind(X_test_num,  X_test_cat)
  
  # y is [time, event] to match fixed deepsurv.R
  y_train <- as.matrix(train_df %>% select(all_of(c(time_col, event_col))))
  y_test  <- as.matrix(test_df  %>% select(all_of(c(time_col, event_col))))
  colnames(y_train) <- c(time_col, event_col)
  colnames(y_test)  <- c(time_col, event_col)
  
  list(
    X_train = as.matrix(X_train),
    X_test  = as.matrix(X_test),
    y_train = y_train,
    y_test  = y_test
  )
}

prep <- prep_design_matrices(train, test, "pfs_months", "pfs_status")
X_train <- prep$X_train
X_test  <- prep$X_test
y_train <- prep$y_train
y_test  <- prep$y_test

cat("X_train dim:", paste(dim(X_train), collapse = " x "), "\n")
cat("X_test  dim:", paste(dim(X_test), collapse = " x "), "\n")

#---------------------------#
# 3) Hyperparameter grid (single-run tuning)
#    NOTE: Keep grid modest to reduce compute and instability.
#---------------------------#
# hyper_grid <- expand.grid(
#   num_layer   = c(1, 2),
#   num_nodes   = c(16, 32, 64),
#   num_dropout = c(0.2, 0.3, 0.4),
#   num_lr      = c(0.0005, 0.001),
#   batch_size  = c(16),
#   num_l2      = c(1e-4, 1e-3, 1e-2),
#   lr_decay    = c(1e-2)      # keep constant unless you truly want to tune it
# )

hyper_grid <- expand.grid(
  num_layer = c(1),
  num_nodes = c(1),
  num_dropout = c(0.2, 0.3, 0.4),
  num_lr = c(0.0005, 0.001),
  batch_size = c(16),
  num_l2 = c(0.2, 0.4, 0.6)
)

#---------------------------#
# 4) CV grid search (TRAIN only; stratified folds by event)
#---------------------------#
cv_grid_search <- function(X, y, hyper_grid,
                           k_folds = 5,
                           epochs = 100,
                           patience = 20,
                           seed = 123,
                           verbose_fit = 1) {
  
  set.seed(seed)
  tensorflow::set_random_seed(seed)
  
  folds <- caret::createFolds(y[, "pfs_status"], k = k_folds, list = FALSE)
  
  results <- hyper_grid
  results$Mean_C_Index <- NA_real_
  results$Std_C_Index  <- NA_real_
  
  for (i in seq_len(nrow(hyper_grid))) {
    params <- as.list(hyper_grid[i, ])
    c_scores <- numeric(0)
    
    for (fold in 1:k_folds) {
      tr_idx <- which(folds != fold)
      va_idx <- which(folds == fold)
      
      X_tr <- X[tr_idx, , drop = FALSE]
      y_tr <- y[tr_idx, , drop = FALSE]
      X_va <- X[va_idx, , drop = FALSE]
      y_va <- y[va_idx, , drop = FALSE]
      
      model <- build_deepsurv(
        num_input = ncol(X_tr),
        num_layer = params$num_layer,
        num_nodes = params$num_nodes,
        string_activation = "selu",
        num_l2 = params$num_l2,
        num_dropout = params$num_dropout,
        num_lr = params$num_lr,
        lr_decay = params$lr_decay
      )
      
      early_stopping <- callback_early_stopping(
        monitor = "val_loss",
        patience = patience,
        restore_best_weights = TRUE
      )
      
      model %>% fit(
        x = X_tr,
        y = y_tr,
        batch_size = params$batch_size,
        epochs = epochs,
        validation_data = list(X_va, y_va),
        verbose = verbose_fit,
        callbacks = list(early_stopping)
      )
      
      lp_va <- as.numeric(model %>% predict(X_va))
      c_scores <- c(c_scores,
                    c_index(lp_va,
                            time = y_va[, "pfs_months"],
                            event = y_va[, "pfs_status"]))
    }
    
    results$Mean_C_Index[i] <- mean(c_scores, na.rm = TRUE)
    results$Std_C_Index[i]  <- sd(c_scores, na.rm = TRUE)
    
    # Optional progress line
    cat("Grid", i, "of", nrow(hyper_grid),
        "Mean_C:", round(results$Mean_C_Index[i], 4),
        "SD:", round(results$Std_C_Index[i], 4), "\n")
  }
  
  # Best = highest mean; tie-break = lowest SD
  results <- results[order(-results$Mean_C_Index, results$Std_C_Index), , drop = FALSE]
  
  list(results = results, best_params = results[1, , drop = FALSE])
}

#---------------------------#
# 5) Fit final model (full TRAIN) + evaluate TRAIN and TEST
#---------------------------#
fit_and_evaluate <- function(X_train, y_train, X_test, y_test,
                             best_params,
                             epochs = 200,
                             validation_split = 0.25,
                             seed = 123,
                             verbose_fit = 0) {
  
  set.seed(seed)
  tensorflow::set_random_seed(seed)
  
  final_model <- build_deepsurv(
    num_input = ncol(X_train),
    num_layer = best_params$num_layer,
    num_nodes = best_params$num_nodes,
    string_activation = "selu",
    num_l2 = best_params$num_l2,
    num_dropout = best_params$num_dropout,
    num_lr = best_params$num_lr,
    lr_decay = best_params$lr_decay
  )
  
  callbacks <- list(
    callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.5, patience = 5, min_lr = 1e-6),
    callback_early_stopping(monitor = "loss", patience = 10, restore_best_weights = TRUE)
  )
  
  history <- final_model %>% fit(
    x = X_train,
    y = y_train,
    batch_size = best_params$batch_size,
    epochs = epochs,
    validation_split = validation_split,
    verbose = verbose_fit,
    callbacks = callbacks
  )
  
  lp_tr <- as.numeric(final_model %>% predict(X_train))
  lp_te <- as.numeric(final_model %>% predict(X_test))
  
  c_tr <- c_index(lp_tr, time = y_train[, "pfs_months"], event = y_train[, "pfs_status"])
  c_te <- c_index(lp_te, time = y_test[, "pfs_months"],  event = y_test[, "pfs_status"])
  
  list(model = final_model, history = history, c_train = c_tr, c_test = c_te)
}

#---------------------------#
# 6) SINGLE RUN: execute tuning + final evaluation
#---------------------------#
seed_single <- 123

cat("\n===== SINGLE RUN: GRID SEARCH CV =====\n")
gs <- cv_grid_search(
  X = X_train,
  y = y_train,
  hyper_grid = hyper_grid,
  k_folds = 5,
  epochs = 100,
  patience = 20,
  seed = seed_single,
  verbose_fit = 0
)

cv_results <- gs$results
best_params <- gs$best_params

cat("\n===== BEST HYPERPARAMETERS (SINGLE RUN) =====\n")
print(best_params)



cat("\n===== FIT FINAL MODEL + EVALUATE =====\n")
fe <- fit_and_evaluate(
  X_train = X_train, y_train = y_train,
  X_test  = X_test,  y_test  = y_test,
  best_params = best_params,
  epochs = 200,
  validation_split = 0.25,
  seed = seed_single,
  verbose_fit = 0
)

final_model <- fe$model
cat("Final Train C-index:", fe$c_train, "\n")
cat("Final Test  C-index:", fe$c_test, "\n")

# Optional: save CV table
write.csv(cv_results, "deepsurv_cv_results_single_run_clincal.csv", row.names = FALSE)

#---------------------------#
# 7) Variable importance (GBFA) + plots
#---------------------------#
cat("\n===== VARIABLE IMPORTANCE (GBFA) =====\n")
variable_importance <- get_variable_importance(final_model, X_train)
importance_df <- variable_importance %>% arrange(desc(Importance))
print(head(importance_df, 20))

ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "DeepSurv feature attribution (GBFA) — single run",
    x = "Features",
    y = "Importance score"
  ) +
  theme_minimal()

median_sensitivity <- median(importance_df$Importance, na.rm = TRUE)
shortlisted_GBFA <- importance_df %>% filter(Importance > median_sensitivity)
shortlisted_GBFA

ggplot(shortlisted_GBFA, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "DeepSurv features above median attribution — single run",
    x = "Features",
    y = "Importance score"
  ) +
  theme_minimal()



#---------------------------#
# 8) Risk predictions on TEST + five-number summary
#---------------------------#
cat("\n===== RISK PREDICTIONS (TEST) =====\n")
predicted_log_risk <- as.numeric(final_model %>% predict(X_test))
predicted_relative_risk <- exp(predicted_log_risk)

cat("Mean relative risk:", mean(predicted_relative_risk), "\n")
cat("Median relative risk:", median(predicted_relative_risk), "\n")

deepsurv_risks <- data.frame(
  predicted_log_risk = predicted_log_risk,
  predicted_relative_risk = predicted_relative_risk,
  pfs_months = y_test[, "pfs_months"],
  pfs_status = y_test[, "pfs_status"]
)

fivenum_deepsurv_risk <- tibble(surv_risk = predicted_relative_risk) %>%
  summarise(
    n = n(),
    mean = mean(surv_risk, na.rm = TRUE),
    stdev = sd(surv_risk, na.rm = TRUE),
    min = min(surv_risk, na.rm = TRUE),
    q1 = quantile(surv_risk, 0.25, na.rm = TRUE),
    median = median(surv_risk, na.rm = TRUE),
    q3 = quantile(surv_risk, 0.75, na.rm = TRUE),
    max = max(surv_risk, na.rm = TRUE)
  )

print(fivenum_deepsurv_risk)


# Optional: save predictions
# write.csv(deepsurv_risks, "deepsurv_risk_predictions_test_single_run.csv", row.names = FALSE)





###############################################################################
# OPTIONAL SECTION: MULTI-RUN REPEATS (If you later change your mind)
#   - Repeats CV tuning + final fit + test evaluation across multiple seeds
#   - Stores CV-best mean/SD and final train/test C-indices
###############################################################################

run_multiple_repeats <- function(X_train, y_train, X_test, y_test,
                                 hyper_grid,
                                 num_repeats = 5,
                                 base_seed = 123,
                                 k_folds = 5,
                                 cv_epochs = 100,
                                 final_epochs = 200,
                                 patience = 20,
                                 validation_split = 0.25,
                                 verbose_fit = 0) {
  
  results_summary <- data.frame(
    Run = integer(),
    Seed = integer(),
    Best_Mean_CV = numeric(),
    Best_SD_CV = numeric(),
    Final_Train_C_Index = numeric(),
    Final_Test_C_Index  = numeric()
  )
  
  models_by_run <- vector("list", length = num_repeats)
  best_params_by_run <- vector("list", length = num_repeats)
  cv_tables_by_run <- vector("list", length = num_repeats)
  
  for (r in 1:num_repeats) {
    seed_r <- base_seed + r
    cat("\n==============================\n")
    cat("MULTI-RUN", r, "Seed", seed_r, "\n")
    cat("==============================\n")
    
    gs_r <- cv_grid_search(
      X = X_train,
      y = y_train,
      hyper_grid = hyper_grid,
      k_folds = k_folds,
      epochs = cv_epochs,
      patience = patience,
      seed = seed_r,
      verbose_fit = verbose_fit
    )
    
    best_r <- gs_r$best_params
    fe_r <- fit_and_evaluate(
      X_train = X_train, y_train = y_train,
      X_test  = X_test,  y_test  = y_test,
      best_params = best_r,
      epochs = final_epochs,
      validation_split = validation_split,
      seed = seed_r,
      verbose_fit = verbose_fit
    )
    
    results_summary <- rbind(
      results_summary,
      data.frame(
        Run = r,
        Seed = seed_r,
        Best_Mean_CV = best_r$Mean_C_Index,
        Best_SD_CV   = best_r$Std_C_Index,
        Final_Train_C_Index = fe_r$c_train,
        Final_Test_C_Index  = fe_r$c_test
      )
    )
    
    models_by_run[[r]] <- fe_r$model
    best_params_by_run[[r]] <- best_r
    cv_tables_by_run[[r]] <- gs_r$results
    
    cat("Final Train C-index:", fe_r$c_train, "\n")
    cat("Final Test  C-index:", fe_r$c_test, "\n")
  }
  
  cat("\n===== MULTI-RUN SUMMARY =====\n")
  print(results_summary)
  
  cat("\nAggregate Train C-index mean:", mean(results_summary$Final_Train_C_Index), "\n")
  cat("Aggregate Train C-index sd  :", sd(results_summary$Final_Train_C_Index), "\n")
  cat("Aggregate Test  C-index mean:", mean(results_summary$Final_Test_C_Index), "\n")
  cat("Aggregate Test  C-index sd  :", sd(results_summary$Final_Test_C_Index), "\n")
  
  list(
    results_summary = results_summary,
    models_by_run = models_by_run,
    best_params_by_run = best_params_by_run,
    cv_tables_by_run = cv_tables_by_run
  )
}

# Example usage later:
# multi_out <- run_multiple_repeats(
#   X_train = X_train, y_train = y_train,
#   X_test = X_test, y_test = y_test,
#   hyper_grid = hyper_grid,
#   num_repeats = 5,
#   base_seed = 123,
#   k_folds = 5
# )





