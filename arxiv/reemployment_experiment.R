#rm(list = ls())
library(ggplot2)
library(pls)
library(grf)
library(data.table)




setwd("C:/Users/Marie/OneDrive/Documents/forest_pls")
penn <- as.data.frame(read.table("penn_jae.dat", header=T ))
head(penn)

### Generate data
W <- as.matrix(as.numeric(penn$tg == 4))
Y <- as.matrix(log(penn$inuidur1))
X <- as.matrix(penn[, -c(2, 3, 4)])


##### RUN PLS and extract scores - use cross-validation for that
model_y <- plsr(Y~X, scale=TRUE, validation="CV")
plot(RMSEP(model_y, se = T)) ### 2 components
model <- RMSEP(model_y)
model$val


model_w <- plsr(W~X, scale=TRUE, validation="CV") ### since this is a randomized
## experiment, treatment is not really correlated with covariates, its balanced
# as shown also by chernozhukov. 

### Extract scores (components)
X_scores <- as.data.frame(model_y$scores[, c(1, 2)])
names(X_scores)[c(1, 2)] <- c("c1", "c2")
X_scores <- as.matrix(X_scores)



### Run random forests for heterogeneity
#install.packages("grf")
causal_model <- causal_forest(X_scores, Y, W, num.trees = 1000, seed = 0)
oob_predictions <- predict(causal_model, estimate.variance = TRUE)
head(oob_predictions)


### analyze treatment effect heterogeneity
covariates_all <- as.data.frame(cbind(penn[, -c(2, 3, 4)], X_scores))
covariates_all$ate <- oob_predictions$predictions
covariates_all$var_ate <- oob_predictions$variance.estimates
write.csv(covariates_all, "covariates_all.csv", row.names = F)


##rm(list = ls())
covariates_all <- fread("covariates_all.csv")
head(covariates_all)








### generate vigintiles
quantiles_c1 <- quantile(covariates_all$c1, seq(0, 1, 0.05))
quantiles_c2 <- quantile(covariates_all$c2, seq(0, 1, 0.05))

covariates_all$vigintiles_c1 <- 100
covariates_all$vigintiles_c2 <- 100


for (i in 2:length(quantiles_c1)){
  
  covariates_all$vigintiles_c1 <- ifelse(covariates_all$c1 >= quantiles_c1[i-1] &
                                           covariates_all$c1 <= quantiles_c1[i], 
                                         i-1, 
                                         covariates_all$vigintiles_c1) 
  
  covariates_all$vigintiles_c2 <- ifelse(covariates_all$c2 >= quantiles_c2[i-1] &
                                           covariates_all$c2 <= quantiles_c2[i], 
                                         i-1, 
                                         covariates_all$vigintiles_c2) 
}



summary(covariates_all$vigintiles_c1)
summary(covariates_all$vigintiles_c2)

library(data.table)
covariates_all <- as.data.table(covariates_all)
head(covariates_all)

#### Component 1
quantiles_data <- covariates_all[,  .(ate_med = quantile(ate, 0.5), 
                                      std_ate_med = quantile(sqrt(var_ate), 0.5), 
                                      ate_025 = quantile(ate, 0.025), 
                                      std_ate_025 = quantile(sqrt(var_ate), 0.025),
                                      ate_975 = quantile(ate, 1-0.025), 
                                      std_ate_975 = quantile(sqrt(var_ate), 1-0.025)
                                      ), .(vigintiles_c1)]


colors <- c(".5" = "#00827F", ".025" = "#967BB6", ".975" = "#797979")


ggplot(quantiles_data, aes( x = vigintiles_c1)) +
  geom_line(aes( y = ate_med, color = ".5"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_med - std_ate_med, 
                    ymax = ate_med + std_ate_med, 
                    color = ".5")) +
  geom_line(aes( y = ate_025, color = ".025"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_025 - std_ate_med, 
                    ymax = ate_025 + std_ate_med, 
                    color = ".025")) +
  geom_line(aes( y = ate_975, color = ".975"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_975 - std_ate_med, 
                    ymax = ate_975 + std_ate_med, 
                    color = ".975"))  +
  theme_classic() + scale_color_manual(values = colors, name = "Percentiles") +
  xlab("") + ylab("policy effect") +
  theme(axis.title = element_text(size = 13))






########## Component 2
quantiles_data <- covariates_all[,  .(ate_med = quantile(ate, 0.5), 
                                      std_ate_med = quantile(sqrt(var_ate), 0.5), 
                                      ate_025 = quantile(ate, 0.025), 
                                      std_ate_025 = quantile(sqrt(var_ate), 0.025),
                                      ate_975 = quantile(ate, 1-0.025), 
                                      std_ate_975 = quantile(sqrt(var_ate), 1-0.025)
                                      ), .(vigintiles_c2)]


colors <- c(".5" = "#00827F", ".025" = "#967BB6", ".975" = "#797979")


ggplot(quantiles_data, aes( x = vigintiles_c2)) +
  geom_line(aes( y = ate_med, color = ".5"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_med - std_ate_med, 
                    ymax = ate_med + std_ate_med, 
                    color = ".5")) +
  geom_line(aes( y = ate_025, color = ".025"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_025 - std_ate_med, 
                    ymax = ate_025 + std_ate_med, 
                    color = ".025")) +
  geom_line(aes( y = ate_975, color = ".975"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_975 - std_ate_med, 
                    ymax = ate_975 + std_ate_med, 
                    color = ".975"))  +
  theme_classic() + scale_color_manual(values = colors, name = "Percentiles") +
  xlab("") + ylab("policy effect") + 
  theme(axis.title = element_text(size = 13))



### Explain these components 
reg_c1 <- lm(X_scores[, 1] ~ X)
summary(reg_c1)

round(reg_c1$coefficients, 3)


reg_c2 <- lm(X_scores[, 2] ~ X)
summary(reg_c2)

round(reg_c2$coefficients, 3)


stargazer::stargazer(reg_c1, reg_c2, rownames = F, font.size = "small", 
                     no.space = T)


### describe scores
plot(RMSEP(model_y), main = "", legend = "topright", ylim = c(1.19, 1.225))


X_scores <- as.data.frame(X_scores)
ggplot(X_scores, aes(x = X_scores[, 1], y =  X_scores[, 2])) +
  geom_point(aes(color = X_scores[, 1]), alpha = 0.8, color = "#87AFC7", 
             size = 1.5) + theme_classic() +
  xlab("Component 1") + ylab("Component 2")  +
  theme(legend.position = "None")



########## RUN CAUSAL FOREST ##################################
causal_model_original <- causal_forest(X, Y, W, num.trees = 1000, seed = 0)
oob_predictions_orig <- predict(causal_model_original, estimate.variance = TRUE)
head(oob_predictions_orig)
covariates_all$ate_cf <- oob_predictions_orig$predictions
covariates_all$ate_var_cf <- oob_predictions_orig$variance.estimates


########## Component 1 #############
quantiles_data_cf <- covariates_all[,  .(ate_med = quantile(ate_cf, 0.5), 
                                      std_ate_med = quantile(sqrt(ate_var_cf), 0.5), 
                                      ate_025 = quantile(ate_cf, 0.025), 
                                      std_ate_025 = quantile(sqrt(ate_var_cf), 0.025),
                                      ate_975 = quantile(ate_cf, 1-0.025), 
                                      std_ate_975 = quantile(sqrt(ate_var_cf), 1-0.025)
), .(vigintiles_c1)]


colors <- c(".5" = "#00827F", ".025" = "#967BB6", ".975" = "#797979")


ggplot(quantiles_data_cf, aes( x = vigintiles_c1)) +
  geom_line(aes( y = ate_med, color = ".5"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_med - std_ate_med, 
                    ymax = ate_med + std_ate_med, 
                    color = ".5")) +
  geom_line(aes( y = ate_025, color = ".025"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_025 - std_ate_med, 
                    ymax = ate_025 + std_ate_med, 
                    color = ".025")) +
  geom_line(aes( y = ate_975, color = ".975"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_975 - std_ate_med, 
                    ymax = ate_975 + std_ate_med, 
                    color = ".975"))  +
  theme_classic() + scale_color_manual(values = colors, name = "Percentiles") +
  xlab("") + ylab("policy effect") +
  theme(axis.title = element_text(size = 13))




########## Component 2 #############
quantiles_data_orig_2 <- covariates_all[,  .(ate_med = quantile(ate_cf, 0.5), 
                                      std_ate_med = quantile(sqrt(ate_var_cf), 0.5), 
                                      ate_025 = quantile(ate_cf, 0.025), 
                                      std_ate_025 = quantile(sqrt(ate_var_cf), 0.025),
                                      ate_975 = quantile(ate_cf, 1-0.025), 
                                      std_ate_975 = quantile(sqrt(ate_var_cf), 1-0.025)
), .(vigintiles_c2)]


colors <- c(".5" = "#00827F", ".025" = "#967BB6", ".975" = "#797979")


ggplot(quantiles_data_orig_2, aes( x = vigintiles_c2)) +
  geom_line(aes( y = ate_med, color = ".5"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_med - std_ate_med, 
                    ymax = ate_med + std_ate_med, 
                    color = ".5")) +
  geom_line(aes( y = ate_025, color = ".025"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_025 - std_ate_med, 
                    ymax = ate_025 + std_ate_med, 
                    color = ".025")) +
  geom_line(aes( y = ate_975, color = ".975"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_975 - std_ate_med, 
                    ymax = ate_975 + std_ate_med, 
                    color = ".975"))  +
  theme_classic() + scale_color_manual(values = colors, name = "Percentiles") +
  xlab("") + ylab("policy effect") + 
  theme(axis.title = element_text(size = 13))


#### WHEN checking individual covariates, we find no heterogeneity ############

check_oroginal_covariates <- function(var){
data <- covariates_all[,  .(ate_med = quantile(ate_cf, 0.5), 
                                             std_ate_med = quantile(sqrt(ate_var_cf), 0.5), 
                                             ate_025 = quantile(ate_cf, 0.025), 
                                             std_ate_025 = quantile(sqrt(ate_var_cf), 0.025),
                                             ate_975 = quantile(ate_cf, 1-0.025), 
                                             std_ate_975 = quantile(sqrt(ate_var_cf), 1-0.025)
), .(var)]

return(ggplot(data, aes( x = var)) +
  geom_line(aes( y = ate_med, color = ".5"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_med - std_ate_med, 
                    ymax = ate_med + std_ate_med, 
                    color = ".5")) +
  geom_line(aes( y = ate_025, color = ".025"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_025 - std_ate_med, 
                    ymax = ate_025 + std_ate_med, 
                    color = ".025")) +
  geom_line(aes( y = ate_975, color = ".975"), size = 1.3) +
  geom_errorbar(aes(ymin = ate_975 - std_ate_med, 
                    ymax = ate_975 + std_ate_med, 
                    color = ".975"))  +
  theme_classic() + scale_color_manual(values = colors, name = "Percentiles") +
  xlab("") + ylab("policy effect") + 
  theme(axis.title = element_text(size = 13)))


}
check_oroginal_covariates(covariates_all$female)
check_oroginal_covariates(covariates_all$black)
check_oroginal_covariates(covariates_all$hispanic)
check_oroginal_covariates(covariates_all$agelt35)


