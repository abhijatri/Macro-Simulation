library(vars)
library(stringi)
library(dplyr)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(reshape2)


var_simulation <- function(df, p = 1, type = NULL, restrict_type = NULL, thresh = 2,
                            nsim = 100, burn_in_period = 30,
                            forecast_start_data = NULL, relative_lags = NULL) {
  
  r <- nrow(df)
  nv <- ncol(df)
  
  if (is.null(relative_lags) == TRUE) {
    relative_lags <- rep(0, ncol(df))
  } else {                                 
    relative_lags <- relative_lags
  }
  
  var.2c <- VAR(df, p = p, type = "const")
  
  if (is.null(type) == FALSE & is.null(restrict_type) == FALSE) {
    if (type == "restrict" & restrict_type == "ser") {
      a <- restrict(var.2c, method = restrict_type, thresh = thresh) 
    } else if (type == "restrict" & restrict_type == "man") {
      a <- restrict(var.2c, method = restrict_type, resmat = restrict_mat)  
    }} else if (is.null(type) == TRUE & is.null(restrict_type) == TRUE) {
      a <- var.2c
    } else {
      message("type and restrict_type combination not appropriately entered, will go ahead with the unrestricted VAR model")
      a <- var.2c
    }
  
  #defining default forecast start data - end points of data
  if(is.null(forecast_start_data) == TRUE) {
    forecast_start_data <- as.matrix(df[(r - p + 1):r,])
  }
  
  #collating the VAR model coefficients
  coeff <- coefficients(a)
  
  #storing the variance covariance matrix of the residuals to be used for drawing random samples from mvtnorm in the simulation
  varcov_resid <- summary(a)$covres
  
  coeff_mat <- rep(0, 1 + p*length(names(coeff)))
  
  nv <- length(names(coeff))
  names(coeff_mat)[1] <- "const"
  for (i in 1:p) {
    names(coeff_mat)[((i-1)*nv + 2): (i*nv + 1)] <- paste0(names(coeff),".l",i)
  }
  
  #create a generic coefficient matrix
  coeff_mat <- bind_rows(coeff_mat, do.call('bind_rows', sapply(names(coeff), function(x){
    j <- t(coeff[[x]][,1])
    j <- data.frame(j)
    colnames(j) <- rownames(coeff[[x]])
    j
  }, simplify = F)))
  
  coeff_mat[is.na(coeff_mat)] <- 0
  coeff_mat <- data.frame(coeff_mat)[-1,]
  rownames(coeff_mat) <- names(coeff)
  
  #define the forecast matrix
  forecast_mat <- matrix(NA, nrow = length(names(coeff)), ncol = burn_in_period + p)
  rownames(forecast_mat) <- names(coeff)
  
  forecast_mat[,1:p] <- t(forecast_start_data)
  
  forecast_mat <- as.matrix(forecast_mat)
  
  
  #simulation
  
  sim_list <- sapply(1:nsim, function(x) {
    message(paste("Simulation",x))
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
  
  #select the simulations
  sim_obs <- data.frame(do.call('rbind',sapply(1:nsim, function(x){
    sapply(1:nv, function(y){
      t(sim_list[[x]][,(p + burn_in_period - relative_lags[y])])[,y]
    })
  }, simplify = F)))
  
  
  #collating statistics between actuals and simulated
  
  
  #plot simulated vs actual density chart
  df2 <- data.frame(df)
  sim_obs2 <- data.frame(sim_obs)
  df2$data_type <- "Actual"
  sim_obs2$data_type <- "Simulated"
  act_sim_obs <- rbind(df2, sim_obs2)
  
  act_sim_obs <- melt(act_sim_obs, id.vars = c("data_type"))
  
  density_plot <- ggplot(data = act_sim_obs, aes(x = value, colour = data_type)) + 
    stat_density(geom = "line") + facet_wrap(variable ~ ., scales = "free") + 
    theme_bw() + labs(x = "", y = "Density", colour = "", title = "Density charts of Macro Economic Indiactors",
                      subtitle = paste0("Actual vs Simulated (",nsim,")")) +
    theme(legend.position = "bottom", panel.grid = element_blank())
  
  return_list <- list(model_summary = summary(a),
                      sim_obs = sim_obs,
                      sim_full = sim_list,
                      mean = list(actual = apply(df,2, mean), simulated = apply(sim_obs,2,mean)),
                      sd = list(actual = apply(df,2, sd), simulated = apply(sim_obs,2,sd)),
                      varcov = list(actual = var(df), simulated = var(sim_obs)),
                      corr =list(actual = cor(df), simulated = cor(sim_obs)),
                      density_plot = density_plot)
  return_list
}


p <- 2
type <- "restrict"
restrict_type <- "ser"
nsim <- 10000
burn_in_period <- 30
forecast_start_data <- as.matrix(Canada[(r - p + 1):r,])
                                 
#check what provides the optimum lag of endogeneous variables
varselect <- VARselect(y = df)
print(varselect)

l <- var_simulation(df = Canada, p = 2, nsim = 200, burn_in_period = 20, relative_lags = c(0,1,1,2),
                    type = "restrict", restrict_type = "ser", thresh = 5)                                 


#plotting burn-in period charts
                                 
                                 