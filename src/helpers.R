# ============================================================================
# Data Splitting: Training and Test Sets
# ============================================================================

# Simple function to split data into training and test sets
# Uses stratified sampling to ensure equal representation of treatment groups
split_data <- function(Y, P, X, 
                       tau_true = NULL, 
                       train_prop = 0.5, 
                       seed = 123) {
  # Input:
  #   Y: outcome variable (vector or matrix)
  #   P: treatment indicator (vector)
  #   X: covariates (matrix)
  #   tau_true: true treatment effects (optional, vector)
  #   train_prop: proportion of data for training (default 0.5)
  #   seed: random seed for reproducibility
  # Output:
  #   List containing training and test sets
  
  set.seed(seed)
  
  # Convert to matrices if needed
  Y <- as.matrix(Y)
  P <- as.matrix(P)
  X <- as.matrix(X)
  
  # Get sample size
  N <- nrow(X)
  
  # Convert P to vector for easier indexing
  P_vec <- as.vector(P)
  
  # Stratified sampling: split separately for each treatment group
  # Get indices for treated (P=1) and control (P=0) groups
  idx_treated <- which(P_vec == 1)
  idx_control <- which(P_vec == 0)
  
  # Sample from each group proportionally
  n_treated <- length(idx_treated)
  n_control <- length(idx_control)
  
  n_train_treated <- floor(n_treated * train_prop)
  n_train_control <- floor(n_control * train_prop)
  
  # Sample indices for training set from each group
  train_idx_treated <- sample(idx_treated, size = n_train_treated, replace = FALSE)
  train_idx_control <- sample(idx_control, size = n_train_control, replace = FALSE)
  
  # Combine training indices
  train_idx <- c(train_idx_treated, train_idx_control)
  test_idx <- setdiff(1:N, train_idx)
  
  # Split data - handle tau_true if provided
  if (!is.null(tau_true)) {
    # Convert tau_true to vector if needed
    tau_true <- as.vector(tau_true)
    
    S_tr <- list(
      Y = Y[train_idx, , drop = FALSE],
      P = P[train_idx, , drop = FALSE],
      X = X[train_idx, , drop = FALSE],
      tau_true = tau_true[train_idx]
    )
    
    S_te <- list(
      Y = Y[test_idx, , drop = FALSE],
      P = P[test_idx, , drop = FALSE],
      X = X[test_idx, , drop = FALSE],
      tau_true = tau_true[test_idx]
    )
  } else {
    # No tau_true provided
    S_tr <- list(
      Y = Y[train_idx, , drop = FALSE],
      P = P[train_idx, , drop = FALSE],
      X = X[train_idx, , drop = FALSE]
    )
    
    S_te <- list(
      Y = Y[test_idx, , drop = FALSE],
      P = P[test_idx, , drop = FALSE],
      X = X[test_idx, , drop = FALSE]
    )
  }
  
  return(list(
    S_tr = S_tr,
    S_te = S_te,
    train_idx = train_idx,
    test_idx = test_idx
  ))
}





PLS_fit_model <- function(Y_tr, 
                          X_tr, 
                          X_te,
                          X_eval, 
                          seed = 1231){
  
  set.seed(seed)
  
  
  
  ## Train the PLS model for the response (Y ~ X) and assignment (P ~ X) components
  pls_fit = plsr(Y_tr ~ X_tr, scale = TRUE, validation = "CV", seed = seed)
  
  # Extract RMSEP and find the minimum
  rmsep_cv <- RMSEP(pls_fit, estimate = "CV", seed = seed)
  rmsep_values <- as.vector(rmsep_cv$val[1, 1, ])
  
  # Find index of minimum (subtract 1 because index 1 corresponds to 0 components)
  optimal_ncomp <- which.min(rmsep_values) - 1
  
  if (optimal_ncomp > 0) {
    # Refit PLS with optimal number of components
    pls_fit_optimal <- plsr(Y_tr ~ X_tr,
                            ncomp = optimal_ncomp,
                            scale = TRUE,
                            validation = "none", 
                            seed = seed)
    
    # Compute scaling parameters from training data (what plsr does with scale=TRUE)
    center_values <- apply(X_tr, 2, mean)
    scale_values <- apply(X_tr, 2, sd)
    
    # Scale test data with train data mean and std. 
    X_te_scaled <- scale(X_te, center = center_values, scale = scale_values)
    X_tr_scaled <- scale(X_tr, center = center_values, scale = scale_values)
    X_eval_scaled <- scale(X_eval, center = center_values, scale = scale_values)
    
    
    
    # Step 2: Compute scores by multiplying scaled X by loadings from training model
    C_te <- X_te_scaled %*% loadings(pls_fit_optimal)[, 1:optimal_ncomp, drop = FALSE]
    C_tr <- X_tr_scaled %*% loadings(pls_fit_optimal)[, 1:optimal_ncomp, drop = FALSE]
    C_eval <- X_eval_scaled %*% loadings(pls_fit_optimal)[, 1:optimal_ncomp, drop = FALSE]
    
    
    # Add column names
    colnames(C_te) <- paste0("C_effect_", 1:optimal_ncomp)
    colnames(C_tr) <- paste0("C_effect_", 1:optimal_ncomp)
    colnames(C_eval) <- paste0("C_effect_", 1:optimal_ncomp)
    
    
    
    # Returns NULL if the optimal number of components is found to be zero
  } else {
    
    pls_fit_optimal <- NULL
    C_te <- NULL
    C_tr <- NULL
    C_eval <- NULL
    
  }
  
  return(list(C_te = C_te, 
              C_tr = C_tr, 
              C_eval = C_eval,
              optimal_ncomp = optimal_ncomp))
  
}



PLS_predict_components <- function(Y_tr, 
                                   P_tr,
                                   X_tr, 
                                   X_te,
                                   X_eval,
                                   seed = 231) {
  
  ## response/effect components
  pls_model_Y <- PLS_fit_model(Y_tr = Y_tr,
                               X_tr = X_tr, 
                               X_te = X_te,
                               X_eval = X_eval,
                               seed = seed)
  
  C_effect_tr <- pls_model_Y$C_tr
  C_effect_te <- pls_model_Y$C_te
  C_effect_eval <- pls_model_Y$C_eval

  pls_model_P <- PLS_fit_model(Y_tr = P_tr,
                               X_tr = X_tr, 
                               X_te = X_te,
                               X_eval = X_eval,
                               seed = seed)
  
  
  C_assignment_tr <- pls_model_P$C_tr
  C_assignment_te <- pls_model_P$C_te
  C_assignment_eval <- pls_model_P$C_eval

  
  return(list(C_effect_tr = C_effect_tr, 
              C_effect_te = C_effect_te,
              C_effect_eval = C_effect_eval,
              C_assignment_tr = C_assignment_tr, 
              C_assignment_te = C_assignment_te, 
              C_assignment_eval = C_assignment_eval))
  
}




fit_GRF <- function(X_matrix, 
                    Y_vector, 
                    P_vector, 
                    num.trees = 2000, 
                    min.node.size = 5,
                    honesty = TRUE,
                    seed = 16112){
  
  # Fit causal forest on the training data using components
  causal_forest_fit <- causal_forest(
    X = X_matrix,                # Features: effect components from training
    Y = Y_vector,                # Outcome (training)
    W = P_vector,                # Treatment (training)
    num.trees = num.trees,
    min.node.size = min.node.size,
    honesty = honesty, 
    seed = seed
  )
  
  
  return(causal_forest_fit)
  
}



forestPLS <- function(Y = Y, 
                      P = P, 
                      X = X, 
                      tau_true = tau_true, 
                      train_prop = 0.5, 
                      seed = 16112,
                      num.trees = 2000, 
                      min.node.size = 5,
                      internal_train_prop = 0.5,
                      honesty = TRUE, 
                      run_grf = TRUE) {
  
  # Split the data
  split_result <- split_data(Y = Y, P = P, X = X, 
                             tau_true = tau_true, 
                             train_prop = train_prop, 
                             seed = seed)
  
  ## Extract train and test variables
  Y_tr <- split_result$S_tr$Y
  Y_te <- split_result$S_te$Y
  X_tr <- split_result$S_tr$X
  X_te <- split_result$S_te$X
  P_tr <- split_result$S_tr$P
  P_te <- split_result$S_te$P
  
  
  ## Further split train data into two partitions where the first one uses PLS
  # to estimate components, the second one uses GRF with predicted components
  split_result_train <- split_data(Y = Y_tr, P = P_tr, X = X_tr, 
                             tau_true = NULL, 
                             train_prop = internal_train_prop, 
                             seed = seed)
  
  
  
  ## Extract train variables split into two
  Y_tr_1 <- split_result_train$S_tr$Y
  Y_tr_2 <- split_result_train$S_te$Y
  X_tr_1 <- split_result_train$S_tr$X
  X_tr_2 <- split_result_train$S_te$X
  P_tr_1 <- split_result_train$S_tr$P
  P_tr_2 <- split_result_train$S_te$P
  
  
                                     
  # Extract PLS components (here, we predict components for both tr_2 subsample that is used by the GRF
  # and evaluation set used for calculating conditional average policy effects)
  C_list <- PLS_predict_components(Y_tr = Y_tr_1,  P_tr = P_tr_1,
                                   X_tr = X_tr_1, X_te = X_tr_2, X_eval = X_te,
                                   seed = seed) 
  
  CATE_true_te <- split_result$S_te$tau_true
  
  
  
  # Combine components (handle NULL cases)
  C_tr_1 <- cbind(C_list$C_effect_tr, C_list$C_assignment_tr)
  C_tr_2 <- cbind(C_list$C_effect_te, C_list$C_assignment_te)
  C_te <- cbind(C_list$C_effect_eval, C_list$C_assignment_eval) 
  
  
  
  # Handle case where all components are NULL (ncomp = 0)
  if (is.null(C_tr_1) || is.null(C_tr_2)) {
    warning("No PLS components extracted (ncomp = 0). Using original covariates.")
    C_tr_1 <- X_tr_1
    C_tr_2 <- X_tr_2
    C_te <- X_te
  }
  
  # Fit causal forest with PLS components
  cf_model_C <- fit_GRF(X_matrix = C_tr_2, 
                        Y_vector = Y_tr_2, 
                        P_vector = P_tr_2,
                        num.trees = num.trees,
                        min.node.size = min.node.size,
                        honesty = honesty)
  
  # Predict treatment effects on test data
  tau_hat_comps <- predict(cf_model_C, newdata = C_te)$predictions
  rmse_C <- RMSE(CATE_true_te, tau_hat_comps)
  
  
  if (run_grf){
  ## Now fit for X 
  cf_model_X <- fit_GRF(X_matrix = X_tr_2, 
                        Y_vector = Y_tr_2, 
                        P_vector = P_tr_2, 
                        num.trees = num.trees, 
                        min.node.size = min.node.size,
                        honesty = honesty)
  
  # Predict treatment effects on test data
  tau_hat_X <- predict(cf_model_X, newdata = X_te)$predictions
  
  # Compute RMSE
  rmse_X <- RMSE(CATE_true_te, tau_hat_X)
  
  } else {
    
    tau_hat_X <- NULL 
    rmse_X <- NULL
    
  }
  

  # Return all requested outputs
  return(list(
    # Components (both types)
    C_effect_tr_1 = C_list$C_effect_tr,
    C_effect_tr_2 = C_list$C_effect_te,
    C_effect_te = C_list$C_effect_eval,
    C_assignment_tr_1 = C_list$C_assignment_tr,
    C_assignment_tr_2 = C_list$C_assignment_te,
    C_assignment_te = C_list$C_assignment_eval,
    C_tr_1 = C_tr_1,  # Combined components (train 1)
    C_tr_2 = C_tr_2,  # Combined components (train 2)
    C_te = C_te,     # Combined components (test/evaluation sample)
    
    # Train and test variables
    Y_tr = Y_tr,
    Y_te = Y_te,
    X_tr = X_tr,
    X_te = X_te,
    P_tr = P_tr,
    P_te = P_te,
    CATE_true_te = CATE_true_te,
    
    # Predicted effects
    tau_hat_comps = tau_hat_comps,  # Predictions using PLS components
    tau_hat_X = tau_hat_X,          # Predictions using full covariates
    
    # Test data RMSE
    rmse_X = rmse_X,                 # RMSE for full covariates
    rmse_C = rmse_C,                 # RMSE for PLS components
    
    # Fitted models
    cf_model_C = cf_model_C,         # Causal forest with components
    cf_model_X = cf_model_X,         # Causal forest with full covariates
    
    # Additional info
    split_result = split_result,    # Full split result
    C_list = C_list                 # Full component list
  ))
}







################ 
# -----------------------------
# KL divergence + Wasserstein helpers
# -----------------------------

safe_numeric <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

# KL(P||Q) using histogram on common bins (P=true, Q=pred)
kl_div_hist <- function(x_true, x_hat, nbins = 50L, eps = 1e-8) {
  x_true <- safe_numeric(x_true); x_hat <- safe_numeric(x_hat)
  if (length(x_true) < 5 || length(x_hat) < 5) return(NA_real_)
  
  pooled <- c(x_true, x_hat)
  
  # robust breaks via pooled quantiles (helps with heavy tails)
  probs <- seq(0, 1, length.out = nbins + 1L)
  brks <- unique(as.numeric(quantile(pooled, probs = probs, type = 8, names = FALSE)))
  
  # fallback if too many duplicates
  if (length(brks) < 10) {
    r <- range(pooled)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(NA_real_)
    brks <- seq(r[1], r[2], length.out = nbins + 1L)
  }
  
  hP <- hist(x_true, breaks = brks, plot = FALSE)
  hQ <- hist(x_hat,  breaks = brks, plot = FALSE)
  
  P <- hP$counts / sum(hP$counts)
  Q <- hQ$counts / sum(hQ$counts)
  
  P <- P + eps; Q <- Q + eps
  P <- P / sum(P); Q <- Q / sum(Q)
  
  sum(P * log(P / Q))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# KL + Wasserstein helpers
# -----------------------------
safe_numeric <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

# KL(P||Q) using histogram on common bins (P=true, Q=pred)
kl_div_hist <- function(x_true, x_hat, nbins = 50L, eps = 1e-8) {
  x_true <- safe_numeric(x_true); x_hat <- safe_numeric(x_hat)
  if (length(x_true) < 5 || length(x_hat) < 5) return(NA_real_)
  
  pooled <- c(x_true, x_hat)
  
  # robust breaks via pooled quantiles (helps with heavy tails)
  probs <- seq(0, 1, length.out = nbins + 1L)
  brks <- unique(as.numeric(quantile(pooled, probs = probs, type = 8, names = FALSE)))
  
  # fallback if too many duplicates
  if (length(brks) < 10) {
    r <- range(pooled)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(NA_real_)
    brks <- seq(r[1], r[2], length.out = nbins + 1L)
  }
  
  hP <- hist(x_true, breaks = brks, plot = FALSE)
  hQ <- hist(x_hat,  breaks = brks, plot = FALSE)
  
  P <- hP$counts / sum(hP$counts)
  Q <- hQ$counts / sum(hQ$counts)
  
  P <- P + eps; Q <- Q + eps
  P <- P / sum(P); Q <- Q / sum(Q)
  
  sum(P * log(P / Q))
}

# 1D Wasserstein-1 distance via quantile functions:
# W1(F,G) = ∫_0^1 |Q_F(u) - Q_G(u)| du  ≈ mean_u |Q_F(u)-Q_G(u)|
wasserstein_1d <- function(x_true, x_hat, M = 1000L) {
  x_true <- safe_numeric(x_true); x_hat <- safe_numeric(x_hat)
  if (length(x_true) < 2 || length(x_hat) < 2) return(NA_real_)
  
  u <- seq(0, 1, length.out = M)
  q1 <- as.numeric(quantile(x_true, probs = u, type = 8, names = FALSE))
  q2 <- as.numeric(quantile(x_hat,  probs = u, type = 8, names = FALSE))
  mean(abs(q1 - q2))
}
