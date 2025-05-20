library(vars)
library(stringi)
library(dplyr)
library(MASS)
library(ggplot2)
library(mvtnorm)
df <- Canada


a <- VAR(y = df, p = 2, type = "const")

summary(a)

VARselect(y = df)

summary(a)

data(Canada) 
var.2c <- VAR(Canada, p = 2, type = "const") 
## Restrictions determined by thresh 
d <- restrict(var.2c, method = "ser") 

summary(var.2c)
coefficients(d)

fanchart(predict(d, n.ahead = 5))

fanchart(predict(var.2c, n.ahead = 10000))


## Restrictions set manually 
restrict <- matrix(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0), nrow=4, ncol=9, byrow=TRUE) 
a <- restrict(var.2c, method = "man", resmat = restrict)


coeff <- coefficients(a)
varcov_resid <- summary(a)$covres

names(c)

sapply(names(c), function(x) {
  x <- c[[x]][,1]
})

x2 <- bind_rows(c[["e"]][,1],c[["prod"]][,1], c[["rw"]][,1], c[["U"]][,1])

x2 <- rep(0, 4)

p <- 2
coeff_mat <- rep(0, 1 + p*length(names(coeff)))
names(coeff_mat) <- c("const", paste0(names(coeff),".l1"), paste0(names(coeff),".l2"))

coeff_mat <- bind_rows(coeff_mat, do.call('bind_rows', sapply(names(coeff), function(x){
  coeff[[x]][,1]
})))

coeff_mat[is.na(coeff_mat)] <- 0
coeff_mat <- data.frame(coeff_mat)
rownames(coeff_mat) <- names(coeff)

nsim <- 10000
burn_in_period <- 30
forecast_mat <- matrix(NA, nrow = length(names(coeff)), ncol = burn_in_period + p)
rownames(forecast_mat) <- names(coeff)

r <- nrow(df)
nv <- length(names(coeff))

forecast_mat[,1:p] <- t(df[(r - p + 1):r,])

forecast_mat <- as.matrix(forecast_mat)


#simulation

sim_list <- sapply(1:nsim, function(x) {
  
  #ith simulation
  for (j in (p+1):(burn_in_period + p)){
    forecast_mat[,j] <- coeff_mat[["const"]] + rmvnorm(1, sigma = varcov_resid)
    for (k in 1:p) {
      forecast_mat[,j] <- forecast_mat[,j] + 
        as.matrix(coeff_mat[-1][,((k-1)*nv + 1):(k*nv)]) %*% as.vector(forecast_mat[,(j - k)])
    }
  }
  forecast_mat 
}, simplify = F)

relative_lags <- rep(c(0,2,2,2), nv)
#select the simulations
sim_obs <- data.frame(do.call('rbind',sapply(1:nsim, function(x){
  sapply(1:nv, function(y){
    t(sim_list[[x]][,(p + burn_in_period - relative_lags[y])])[,y]
  })
}, simplify = F)))

#plot simulated vs actual density chart
library(reshape2)
df2 <- data.frame(df)
df2$data_type <- "Actual"
sim_obs$data_type <- "Simulated"
act_sim_obs <- rbind(df2, sim_obs)

act_sim_obs <- melt(act_sim_obs, id.vars = c("data_type"))

ggplot(data = act_sim_obs, aes(x = value, colour = data_type)) + 
  stat_density(geom = "line") + facet_wrap(variable ~ ., scales = "free") + 
  theme_bw() + labs(x = "", y = "Density", colour = "", title = "Density charts of Macro Economic Indiactors",
                    subtitle = "Actual vs Simulated") +
  theme(legend.position = "bottom", panel.grid = element_blank())

#plotting burn-in period charts