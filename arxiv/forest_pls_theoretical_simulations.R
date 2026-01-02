############## SIMULATION ANALYSIS #########################
### No Confounding - High weight on other features


simulate_data_group <- function(N){
  
  set.seed(30)
  
  D = rbinom(N, 1, 0.5)
  'X <- cbind(runif(N, -1, 1),  runif(N, -1, 1), runif(N, -1, 1), 
                runif(N, -1, 1),  runif(N, -1, 1))'
  
  
  X <- cbind(rnorm(N, -1, 1),  rnorm(N, 1, 1), rnorm(N, 2, 1), 
             rnorm(N, 0, 1))
  #X <- as.matrix(X)
  Y <- 100*X[, 1] + 100*X[, 2] + D*(X[, 3] + 0.5*X[, 4] + 0.2*X[, 3]*X[, 4]) + rnorm(N, 0, 1)
  
  return(as.matrix(cbind(Y, D, X)))
}






#install.packages("glmnet")
library(glmnet)


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
  covariates_all$ate_std <- sqrt(oob_predictions$variance.estimates)
  
  
  return(covariates_all)
  
}




obtain_data_theory <- function(N){

data <- simulate_data_group(N = N)
covariates_all <- forest_pls(Y = data[, 1], D = data[, 2], X = data[, 3:ncol(data)] )
data_all <- covariates_all
data_all$ate_true <- data[, 5] + 0.5*data[, 6] + 0.2*data[, 5]*data[, 6]


return(data_all)

}



data_500 <- obtain_data_theory(N = 500)
data_500$ate_stand <- (data_500[, "ate"] - data_500[, "ate_true"])
data_500 <- as.data.frame(data_500)
mean(data_500$ate_stand)
sd(data_500$ate_stand)



data_1000 <- obtain_data_theory(N = 1000)
data_1000$ate_stand <- (data_1000[, "ate"] - data_1000[, "ate_true"])
data_1000 <- as.data.frame(data_1000)
mean(data_1000$ate_stand)
sd(data_1000$ate_stand)



data_10000 <- obtain_data_theory(N = 10000)
data_10000$ate_stand <- (data_10000[, "ate"] - data_10000[, "ate_true"])
data_10000 <- as.data.frame(data_10000)
mean(data_10000$ate_stand)
sd(data_10000$ate_stand)





data_100000<- obtain_data_theory(N = 100000)
data_100000$ate_stand <- (data_100000[, "ate"] - data_100000[, "ate_true"])
data_100000 <- as.data.frame(data_100000)
mean(data_100000$ate_stand)
sd(data_100000$ate_stand)




ggplot(data_500) +
  geom_density(aes(x= `ate_stand`, fill = "N = 500")) + xlim(-20, 20) + 
  theme_classic() +
  scale_fill_manual(name = "", values = c("N = 500" = "#D3D3D3")) + 
  theme(legend.position = "none", 
          axis.text=element_text(size=13),
          axis.title=element_text(size=13), 
          legend.text=element_text(size=13)) + xlab("") + ylab("Density")


ggplot(data_1000) +
  geom_density(aes(x= `ate_stand`, fill = "N = 1000")) + xlim(-20, 20) + 
  theme_classic() +
  scale_fill_manual(name = "", values = c("N = 1000" = "#D3D3D3"))  +
  theme(legend.position = "none", 
          axis.text=element_text(size=13),
          axis.title=element_text(size=13), 
          legend.text=element_text(size=13)) + xlab("") + ylab("")



ggplot(data_10000) +
  geom_density(aes(x= `ate_stand`, fill = "N = 10,000")) + xlim(-20, 20) + 
  theme_classic() +
  scale_fill_manual(name = "", values = c("N = 10,000" = "#D3D3D3"))  +
  theme(legend.position = "none", 
          axis.text=element_text(size=13),
          axis.title=element_text(size=13), 
          legend.text=element_text(size=13)) + xlab("")  + ylab("Density")
  


ggplot(data_100000) +
  geom_density(aes(x= `ate_stand`, fill = "N = 100,000")) + xlim(-20, 20) + 
  theme_classic() +
  scale_fill_manual(name = "", values = c("N = 100,000" = "#D3D3D3")) +
  theme(legend.position = "none", 
        axis.text=element_text(size=13),
        axis.title=element_text(size=13), 
        legend.text=element_text(size=13)) + xlab("") + ylab("")


