##### SIMULATION ANALYSIS #####
library(pls)
library(ggplot2)
library(MASS)
library(reshape)
library(psych)
library(ggplot2)



##### High weight on features not 
##### influencing policy effect heterogeneity #####
simulate_data <- function(N){
  
  #set.seed(30)
  D = rbinom(N, 1, 0.5)
  'X <- cbind(runif(N, -1, 1),  runif(N, -1, 1), runif(N, -1, 1), 
                runif(N, -1, 1),  runif(N, -1, 1))'
  
  
  X <- cbind(rnorm(N, -1, 1),  rnorm(N, 1, 1), rnorm(N, 2, 1), 
             rnorm(N, 0, 1))
  #X <- as.matrix(X)
  Y <- 100*X[, 1] + 100*X[, 2] + D*X[, 3]*X[, 1] + rnorm(N, 0, 1)
  
  return(as.matrix(cbind(Y, D, X)))
}


##### The same with interaction terms in policy effect heterogeneity #####

simulate_data_group <- function(N){
  
  #set.seed(30)
  
  D = rbinom(N, 1, 0.5)
  'X <- cbind(runif(N, -1, 1),  runif(N, -1, 1), runif(N, -1, 1), 
                runif(N, -1, 1),  runif(N, -1, 1))'
  
  
  X <- cbind(rnorm(N, -1, 1),  rnorm(N, 1, 1), rnorm(N, 2, 1), 
             rnorm(N, 0, 1))
  #X <- as.matrix(X)
  Y <- 100*X[, 1] + 100*X[, 2] + D*(X[, 3] + 0.1*X[, 4] + 0.2*X[, 3]*X[, 4]) + rnorm(N, 0, 1)
  
  return(as.matrix(cbind(Y, D, X)))
}


##### Confounders present (various designs) #####

simulate_linear_back_door <- function(N, method = "lbd"){
  
  #set.seed(30)
  NX <- 5
  mean_X <- c(0, 1, -1, 2, 1)
  var_x <- diag(NX)
  X <- MASS::mvrnorm(N, c(mean_X), var_x)
  
  if (method == "lbd"){
  ## 1. Linear Back-Door
  D <- 0.1*(X[, 1])^2 - X[, 1] + rnorm(N, 0, 1)
  Y <- 0.5*D^2 - D*X[, 1] + rnorm(N, 0, 1)
  
  }
  else if (method == "nbd"){
    
    D <- X[, 1] - X[, 2] + 2*X[, 3] + rnorm(N, 0, 1)
    Y <- 100*X[, 1] + 100*X[, 2] - 100*X[, 3] + 3*D*X[, 1] + rnorm(N, 0, 1)
    
  }
  
  else if (method == "fd"){
    U <- rnorm(N, 0, 1)
    D <- U + rnorm(N, 2, 2)
    X[, 1] <- 2*D*X[, 2] + rnorm(N, 1, 2)
    Y <- X[, 1]^2 - X[, 1] + 100*U + 100*X[, 2] + rnorm(N, 0, 2)
    
  }
  
  else if (method == "niv"){ ## strong instrument, weak confounding
    
    Z_1 <- rnorm(N, -1, 1)
    Z_2 <- rnorm(N, 0, 1)
    U <- rnorm(N, 0, 1)
    D <- 3*Z_1 + 1.5*Z_2 + 0.5*U+ rnorm(N, 0, 1)
    Y <- 0.3*D^2 - 1.5*D*X[, 1] + U + rnorm(N, 0, 1)
    
  }
  
  else if (method == "liv"){ ## weak instrument, strong confounding
    
    Z_1 <- rnorm(N, -1, 1)
    Z_2 <- rnorm(N, 0, 1)
    U <- rnorm(N, 0, 1)
    D <- Z_1 - Z_2 + 0.5*U  + rnorm(N, 0, 1)
    Y <- 0.5*D*X[, 1] - U + rnorm(N, 0, 1)
    
  }
  
  return(cbind(Y, D, X))
  
}






#install.packages("glmnet")
library(glmnet)

##### Forest-PLS algorithm #####

forest_pls <- function(Y, D, X){
  
  
  sim_linear <- plsr(Y~X, scale=TRUE, validation="CV")
  ### Extract scores (components)
  X_scores <- as.data.frame(sim_linear$scores[, c(1, 2)])
  names(X_scores)[c(1, 2)] <- c("c1", "c2")
  X_scores <- as.matrix(X_scores)
  
  
  ### Run random forests for heterogeneity
  #install.packages("grf")
  library(grf)
  causal_model <- causal_forest(X_scores, Y, D, num.trees = 1000, seed = 0)
  oob_predictions <- predict(causal_model, estimate.variance = TRUE)
  head(oob_predictions)
  
  covariates_all <- as.data.frame(X)
  covariates_all$ate <- oob_predictions$predictions
  covariates_all$ate_var <- oob_predictions$variance.estimates
  
  covariates_all <- cbind(covariates_all, X_scores)

  
  return(covariates_all)
  
}


##### Causal forest #####

cf <- function(Y, D, X){
  
  ### causal_forest 
  causal_model_lbd <- causal_forest(X, Y, D, num.trees = 1000, seed = 0)
  oob_predictions_lbd_dir <- predict(causal_model_lbd, estimate.variance = TRUE)
  
  covariates_all <- as.data.frame(X)
  covariates_all$ate <- oob_predictions_lbd_dir$predictions
  covariates_all$ate_var <- oob_predictions_lbd_dir$variance.estimates
  
  return(covariates_all)
  
}


##### Plot all methods #####

#N = 100
plot_function_noconf <- function(N){
  
  list_ate_forest_pls <- list()
  list_ate_cf <- list()
  list_ate_true <- list()

  seed_numbers <- seq(50, 2500, length = 50)
  for (i in 1:50){
    set.seed(seed_numbers[i])
    data <- simulate_data(N = N)
    covariates_all <- forest_pls(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)] )
    covariates_all_cf <- cf(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)])
    
    
    
    list_ate_forest_pls[[i]] <- covariates_all$ate
    list_ate_cf[[i]] <- covariates_all_cf$ate
    list_ate_true[[i]] <- data[, 3]*data[, 5]
  }
  
  ## convert to a dataframe
  ate_forest_pls_data <- as.data.frame(do.call(cbind, list_ate_forest_pls))
  ate_forest_cf_data <- as.data.frame(do.call(cbind, list_ate_cf))
  ate_true <- as.data.frame(do.call(cbind, list_ate_true))
  
  
  means_ate_fpls <- rowMeans(ate_forest_pls_data)
  means_ate_cf <- rowMeans(ate_forest_cf_data)
  means_ate_true <- rowMeans(ate_true)
  covariates_all$ate_fpls <- means_ate_fpls
  covariates_all$ate_cf <- means_ate_cf
  covariates_all$ate_true <- means_ate_true
  data_all <- covariates_all
  

  if (N != 500){
  # Plotting multiple density plots
  
  plot <- ggplot(data_all) +
    geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
    geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
    geom_density(aes(x = ate_true, fill = "True"), alpha = 0.5) +
    
    labs(x = "Value", y = "Density") +
    scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                         "True" = "#967BB6")) +
    xlab("")  + xlim(-5, 1)+ ylim(0, 2) +
    theme_classic() +theme(legend.position = "none",
                                    axis.text=element_text(size=13),
                                    axis.title=element_text(size=13), 
                                    legend.text=element_text(size=13))
} else
    
    {plot <- ggplot(data_all) +
      geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
      geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
      geom_density(aes(x = ate_true, fill = "True"), alpha = 0.5) +
    
    labs(x = "Value", y = "Density") +
    scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                         "True" = "#967BB6")) +
    xlab("")  + xlim(-5, 1)+ylim(0, 2) +  theme_classic() +
      theme(legend.position = c(0.8, 0.7),
                                     axis.text=element_text(size=13),
                                     axis.title=element_text(size=13), 
                                     legend.text=element_text(size=13)) }
  
  return(plot)
}


plot_100 <- plot_function_noconf(N = 100)
plot_500 <- plot_function_noconf(N = 500)
plot_1000 <- plot_function_noconf(N = 1000)
plot_5000 <- plot_function_noconf(N = 5000)




##### GROUP HETEROGENEITY (with interactions) #####

## Here I just obtain data and run the methods similarly as in plots
## will need for the variable importance
obtain_data <- function(N){
  list_ate_forest_pls <- list()
  list_ate_cf <- list()
  list_ate_true <- list()
  seed_numbers <- seq(50, 2500, length = 50)
  for (i in 1:50){
    set.seed(seed_numbers[i])
    data <- simulate_data_group(N = N)
    covariates_all <- forest_pls(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)] )
    covariates_all_cf <- cf(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)])
    
    
    
    list_ate_forest_pls[[i]] <- covariates_all$ate
    list_ate_cf[[i]] <- covariates_all_cf$ate
    list_ate_true[[i]] <- data[, 5] + 0.5*data[, 6] + 0.2*data[, 5]*data[, 6]
  }
  
  ## convert to a dataframe
  ate_forest_pls_data <- as.data.frame(do.call(cbind, list_ate_forest_pls))
  ate_forest_cf_data <- as.data.frame(do.call(cbind, list_ate_cf))
  ate_true <- as.data.frame(do.call(cbind, list_ate_true))
  
  
  means_ate_fpls <- rowMeans(ate_forest_pls_data)
  means_ate_cf <- rowMeans(ate_forest_cf_data)
  means_ate_true <- rowMeans(ate_true)
  covariates_all$ate_fpls <- means_ate_fpls
  covariates_all$ate_cf <- means_ate_cf
  covariates_all$ate_true <- means_ate_true
  data_all <- covariates_all
  
  return(data_all)
}


##### GROUP HETEROGENEITY (with interaction and no confounders) #####
plot_function_noconf_group <- function(N){
  
  
  list_ate_forest_pls <- list()
  list_ate_cf <- list()
  list_ate_true <- list()
  seed_numbers <- seq(50, 2500, length = 50)
  for (i in 1:50){
    set.seed(seed_numbers[i])
    data <- simulate_data_group(N = N)
    covariates_all <- forest_pls(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)] )
    covariates_all_cf <- cf(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)])
    
    
    
    list_ate_forest_pls[[i]] <- covariates_all$ate
    list_ate_cf[[i]] <- covariates_all_cf$ate
    list_ate_true[[i]] <- data[, 5] + 0.5*data[, 6] + 0.2*data[, 5]*data[, 6]
  }
  
  ## convert to a dataframe
  ate_forest_pls_data <- as.data.frame(do.call(cbind, list_ate_forest_pls))
  ate_forest_cf_data <- as.data.frame(do.call(cbind, list_ate_cf))
  ate_true <- as.data.frame(do.call(cbind, list_ate_true))
  
  
  means_ate_fpls <- rowMeans(ate_forest_pls_data)
  means_ate_cf <- rowMeans(ate_forest_cf_data)
  means_ate_true <- rowMeans(ate_true)
  covariates_all$ate_fpls <- means_ate_fpls
  covariates_all$ate_cf <- means_ate_cf
  covariates_all$ate_true <- means_ate_true
  data_all <- covariates_all
  
  
  # Plotting multiple density plots
  if (N != 500){
  plot <- ggplot(data_all) +
    geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
    geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
    geom_density(aes(x = ate_true, fill = "True"), alpha = 0.5) +
    
    labs(x = "Value", y = "Density") +
    scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                         "True" = "#967BB6")) +
    xlab("")  + xlim(0, 5)+ ylim(0, 3)+
    theme_classic() + theme(legend.position = "none", 
                            axis.text=element_text(size=13),
                            axis.title=element_text(size=13), 
                            legend.text=element_text(size=13))
  } else 
    
  {
    plot <- ggplot(data_all) +
      geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
      geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
      geom_density(aes(x = ate_true, fill = "True"), alpha = 0.5) +
      
      labs(x = "Value", y = "Density") +
      scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                           "True" = "#967BB6")) +
      xlab("")   +xlim(0, 5)+ ylim(0, 3) +
      theme_classic() + theme(legend.position = c(0.8, 0.7),
                                axis.text=element_text(size=13),
                                axis.title=element_text(size=13), 
                                legend.text=element_text(size=13))
  }

  return(plot)
}


plot_70 <- plot_function_noconf_group(N = 70)
plot_500 <- plot_function_noconf_group(N = 500)
plot_1000 <- plot_function_noconf_group(N = 1000)
plot_5000 <- plot_function_noconf_group(N = 5000)





##### Variable Importance #####

varimp <- function(data){
  
  ## varimp of forest-pls when the outcome is ate
  model_ate_fpls <- regression_forest(data[, 1:4], data[, "ate_fpls"])
  varimp_fpls <- variable_importance(model_ate_fpls)
  
  ## varimp of causal-forest when the outcome is ate
  model_ate_cf <- regression_forest(data[, 1:4], data[, "ate_cf"])
  varimp_cf <- variable_importance(model_ate_cf)
  
  return(cbind(varimp_fpls, varimp_cf))
  
}


### Show what components consist of, and the treatment effects
data <- obtain_data(N = 100)
data_2 <- obtain_data(N = 1000)
data_3 <- obtain_data(N = 5000)


model_1 <- lm(c1 ~ V1 + V2 + V3 + V4, data = data)
model_2 <- lm(c1 ~ V1 + V2 + V3 + V4, data = data_2)
model_3 <- lm(c1 ~ V1 + V2 + V3 + V4, data = data_3)

model_1_c2 <- lm(c2 ~ V1 + V2 + V3 + V4, data = data)
model_2_c2 <- lm(c2 ~ V1 + V2 + V3 + V4, data = data_2)
model_3_c2 <- lm(c2 ~ V1 + V2 + V3 + V4, data = data_3)

stargazer::stargazer(model_1, model_1_c2, model_2, model_2_c2, model_3, model_3_c2)


summary(model_1)
summary(model_2)
summary(model_3)


summary(model_1_c2)

varimp_100 <- varimp(data = data)
varimp_1000 <- varimp(data = data_2)
varimp_5000 <- varimp(data = data_3)



stargazer::stargazer(varimp_100, varimp_1000, varimp_5000)
#ggpubr::ggarrange(plot_100, plot_500, nrow = 1)
###################






####  3rd CONFOUNDING  ####
## Linear back-door
plot_function_lbd <- function(N, method  = "lbd"){
  
  list_ate_forest_pls <- list()
  list_ate_cf <- list()
  list_ate_liv <- list()
  list_ate_lbd <- list()
  list_ate_nbd <- list()
  list_ate_fd <- list()
  
  seed_numbers <- seq(50, 2500, length = 50)
  for (i in 1:50){
    set.seed(seed_numbers[i])
    data <- simulate_linear_back_door(N = N, method = method)
    covariates_all <- forest_pls(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)] )
    covariates_all_cf <- cf(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)])
    
    
    
    list_ate_forest_pls[[i]] <- covariates_all$ate
    list_ate_cf[[i]] <- covariates_all_cf$ate
    list_ate_liv[[i]] <- 0.5*data[, 3]
    list_ate_lbd[[i]] <- 0.5 - data[, 5]
    list_ate_nbd[[i]] <- 3*data[, 3]
    list_ate_fd[[i]] <- 4*data[, 4]^2 - 2*data[, 4]
  }
  
  ## convert to a dataframe
  ate_forest_pls_data <- as.data.frame(do.call(cbind, list_ate_forest_pls))
  ate_forest_cf_data <- as.data.frame(do.call(cbind, list_ate_cf))
  ate_liv <- as.data.frame(do.call(cbind, list_ate_liv))
  ate_lbd <- as.data.frame(do.call(cbind, list_ate_lbd))
  ate_nbd <- as.data.frame(do.call(cbind, list_ate_nbd))
  ate_fd <- as.data.frame(do.call(cbind, list_ate_fd))
  
  
  
  
  
  means_ate_fpls <- rowMeans(ate_forest_pls_data)
  means_ate_cf <- rowMeans(ate_forest_cf_data)
  means_ate_liv <- rowMeans(ate_liv)
  means_ate_lbd <- rowMeans(ate_lbd)
  means_ate_nbd <- rowMeans(ate_nbd)
  means_ate_fd <- rowMeans(ate_fd)
  
  
  covariates_all$ate_fpls <- means_ate_fpls
  covariates_all$ate_cf <- means_ate_cf
  covariates_all$ate_liv <- means_ate_liv
  covariates_all$ate_lbd <- means_ate_lbd
  covariates_all$ate_nbd <- means_ate_nbd
  covariates_all$ate_fd <- means_ate_fd
  
  data_all <- covariates_all
  
  
  
  # Plotting multiple density plots

  if (method == "nbd"){
    
    
    if (N != 500)
    
    {plot <- ggplot(data_all) +
      geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
      geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
      geom_density(aes(x = ate_nbd, fill = "True"), alpha = 0.5) +
      
      labs(x = "Value", y = "Density") +
      scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                           "True" = "#967BB6")) +
      xlab("")  + xlim(-30, 20)+ ylim(0, 3) + 
      theme_classic() + theme(legend.position = "none",
                                axis.text=element_text(size=13),
                                axis.title=element_text(size=13), 
                                legend.text=element_text(size=13))
    
    
    } else 
        
      {plot <- ggplot(data_all) +
        geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
        geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
        geom_density(aes(x = ate_nbd, fill = "True"), alpha = 0.5) +
        
        labs(x = "Value", y = "Density") +
        scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                              "True" = "#967BB6")) +
        xlab("")  + xlim(-30, 20)+ylim(0, 3) +
    theme_classic() + theme(legend.position = c(0.8, 0.7),
                            axis.text=element_text(size=13),
                            axis.title=element_text(size=13), 
                            legend.text=element_text(size=13))}
 
  }
  
  else if (method == "fd"){
    
    if (N != 500)
    
    {plot <- ggplot(data_all) +
      geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
      geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
      geom_density(aes(x = ate_fd, fill = "True"), alpha = 0.5) +
      
      labs(x = "Value", y = "Density") +
      scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                           "True" = "#967BB6")) +
      xlab("")  + xlim(-20, 60)+ ylim(0, 1) + 
      theme_classic() + theme(legend.position = "none", 
                              axis.text=element_text(size=13),
                              axis.title=element_text(size=13), 
                              legend.text=element_text(size=13)) } 
    else
  
      {plot <- ggplot(data_all) +
        geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
        geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
        geom_density(aes(x = ate_fd, fill = "True"), alpha = 0.5) +
        
        labs(x = "Value", y = "Density") +
        scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                             "True" = "#967BB6")) +
        xlab("")  + xlim(-20, 60)+ ylim(0, 1) + 
        theme_classic() + theme(legend.position = c(0.8, 0.7), 
                                axis.text=element_text(size=13),
                                axis.title=element_text(size=13), 
                                legend.text=element_text(size=13)) }                           
         
  }
  
  else if (method == "liv"){
    
    
    if (N != 500){ 
    
    plot <- ggplot(data_all) +
      geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
      geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
      geom_density(aes(x = ate_liv, fill = "True"), alpha = 0.5) +
      
      labs(x = "Value", y = "Density") +
      scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                           "True" = "#967BB6")) +
      xlab("")  + xlim(-1, 1)+ ylim(0, 20) +
      theme_classic() + theme(legend.position = "none", 
                              axis.text=element_text(size=13),
                              axis.title=element_text(size=13), 
                              legend.text=element_text(size=13)) }
    else 
      
      plot <- ggplot(data_all) +
        geom_density(aes(x = ate_fpls, fill = "Forest-PLS"), alpha = 0.5) +
        geom_density(aes(x = ate_cf, fill = "Causal Forest"), alpha = 0.5) +
        geom_density(aes(x = ate_liv, fill = "True"), alpha = 0.5) +
        
        labs(x = "Value", y = "Density") +
        scale_fill_manual(name = "Policy Effect", values = c("Forest-PLS" = "#31906E", "Causal Forest" = "#F08080", 
                                                             "True" = "#967BB6")) +
        xlab("")  + xlim(-1, 1)+ ylim(0, 20) +
        theme_classic() + theme(legend.position = c(0.8, 0.7), 
                                axis.text=element_text(size=13),
                                axis.title=element_text(size=13), 
                                legend.text=element_text(size=13))
      
  }
  
  
  return(plot)
}


### We provide liv in the main part of the paper, and nbd and fd in appendix
plot_100 <- plot_function_lbd(N = 100, method = "nbd")
plot_500 <- plot_function_lbd(N = 500, method = "nbd")
plot_1000 <- plot_function_lbd(N = 1000, method = "nbd")
plot_5000 <- plot_function_lbd(N = 5000, method = "nbd")



plot_100 <- plot_function_lbd(N = 100, method = "fd")
plot_500 <- plot_function_lbd(N = 500, method = "fd")
plot_1000 <- plot_function_lbd(N = 1000, method = "fd")
plot_5000 <- plot_function_lbd(N = 5000, method = "fd")


## liv goes in the main part
plot_100 <- plot_function_lbd(N = 100, method = "liv")
plot_500 <- plot_function_lbd(N = 500, method = "liv")
plot_1000 <- plot_function_lbd(N = 1000, method = "liv")
plot_5000 <- plot_function_lbd(N = 5000, method = "liv")






####### Compare PLS and LASSO in the simulated experiment
## Lasso drops variables important for explaining policy effect heterogeneity
simulate_data_group_lasso <- function(N){
  
  set.seed(30)
  
  D = rbinom(N, 1, 0.5)
  'X <- cbind(runif(N, -1, 1),  runif(N, -1, 1), runif(N, -1, 1), 
                runif(N, -1, 1),  runif(N, -1, 1))'
  
  
  X <- cbind(rnorm(N, -1, 1),  rnorm(N, 1, 1), rnorm(N, 2, 1), 
             rnorm(N, 0, 1))
  #X <- as.matrix(X)
  Y <- 100*X[, 1] + 100*X[, 2] + D*(X[, 3] + 0.1*X[, 4] + 0.2*X[, 3]*X[, 4]) + rnorm(N, 0, 1)
  
  return(as.matrix(cbind(Y, D, X)))
}







# Read in data
data_for_lasso <- simulate_data_group_lasso(N = 1000)
cv_model <- cv.glmnet(data_for_lasso[, 3:ncol(data_for_lasso)], 
                      data_for_lasso[, 1], alpha = 1)

optimal_lambda <- cv_model$lambda.min


# Perform Lasso regression
lasso_model <- glmnet(data_for_lasso[, 3:ncol(data_for_lasso)], 
                      data_for_lasso[, 1], alpha = 1, lambda = optimal_lambda)
coef(lasso_model)

# Plot lambda versus coefficients
plot(lasso_model, xvar = "lambda", label = TRUE)




