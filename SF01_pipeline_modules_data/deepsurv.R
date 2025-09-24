library(tensorflow)
library(keras)
library(survival)
#BiocManager::install("survcomp")

neg_log_likelihd <- function(y_true, y_pred) {
  event <- y_true[, 1]
  time <- y_true[, 2]
  mask <- k_cast(time <= k_reshape(time, shape = c(-1, 1)), dtype = "float32")
  
  log_loss <- k_log(k_sum(mask * k_exp(y_pred), axis = 1))
  neg_log_loss <- -k_sum(event * (y_pred - log_loss))
  
  return(neg_log_loss / k_sum(event))
}


build_deepsurv <- function(num_input, 
                           num_layer, 
                           num_nodes, 
                           string_activation, 
                           num_l2, 
                           num_dropout, 
                           num_lr, 
                           lr_decay) {
  
  input_layer <- layer_input(shape = c(num_input),
                             dtype = "float32")
  
  hidden_layers <- input_layer
  
  for (i in 1:num_layer) {
    hidden_layers <- hidden_layers %>%
      layer_dense(units = num_nodes, 
                  activation = string_activation, 
                  # kernel_regularizer = regularizer_l2(num_l2),
                  kernel_initializer = "glorot_uniform") %>%
      layer_dropout(rate = num_dropout)
  }
  
  hidden_layers <- hidden_layers %>%
    layer_dense(units = 1, 
                activation = "linear", 
                kernel_regularizer = regularizer_l2(num_l2)) %>% 
    layer_activity_regularization(l2 = num_l2)
  
  model <- keras_model(input_layer, hidden_layers)
  
  model %>% compile(
    loss = neg_log_likelihd,
    optimizer = optimizer_nadam(learning_rate = num_lr, 
                                weight_decay = lr_decay, 
                                clipnorm=1.)
  )
}

c_index <- function(LP_pred, Y, E){
  c_index <- survcomp::concordance.index(LP_pred, Y, E)$c.index
  return(c_index)
}




#===========Define a function to get variable importance using gradients=====
#--Function to get the variable importance using gradients---
# # Load necessary libraries
# library(tensorflow)
# library(keras)

get_variable_importance <- function(model, input_data) {
  # Ensure input data is a TensorFlow constant tensor
  input_tensor <- tf$convert_to_tensor(as.matrix(input_data), dtype = tf$float32)
  
  # Use tf$GradientTape() to record operations for gradient calculation
  gradients_value <- NULL
  
  # Use tf_function() to ensure TensorFlow understands the context
  compute_gradients <- tf_function(function(input_tensor) {
    with(tf$GradientTape() %as% tape, {
      tape$watch(input_tensor)  # Start recording gradients on input tensor
      predictions <- model(input_tensor)  # Get model predictions
    })
    gradients <- tape$gradient(predictions, input_tensor)  # Calculate gradients
    return(gradients)
  })
  
  # Try to compute gradients using the wrapped function
  tryCatch({
    gradients_value <- compute_gradients(input_tensor)
    gradient_values <- keras::k_eval(gradients_value)  # Evaluate gradients
  }, error = function(e) {
    print("Error in computing gradients: ")
    print(e$message)
    return(NULL)
  })
  
  # Check if gradients were successfully computed
  if (is.null(gradient_values)) {
    stop("Failed to calculate gradients. Check TensorFlow setup and context.")
  }
  
  # Calculate mean absolute gradients for each feature across all samples
  feature_importance <- colMeans(abs(gradient_values))
  
  # Create a data frame with the actual feature names and their importance
  feature_names <- colnames(input_data)  # Extract feature names from input data
  
  if (is.null(feature_names)) {
    feature_names <- paste("Feature", 1:ncol(input_data), sep = "_")  # Default names if none found
  }
  
  importance_df <- data.frame(
    Feature = feature_names,  # Use actual feature names from input data
    Importance = feature_importance
  )
  
  # Sort by importance in descending order
  importance_df <- importance_df[order(-importance_df$Importance), ]
  
  return(importance_df)
}