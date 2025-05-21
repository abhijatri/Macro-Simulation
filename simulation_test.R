library(vars)
library(stringi)
library(dplyr)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(reshape2)


var_simulation <- function(df, p = 1, type = NULL, restrict_type = NULL, thresh = 2, restrict_mat = NULL,
                           nsim = 100, burn_in_period = 30,
                           forecast_start_data = NULL, relative_lags = NULL, varcov_mat = NULL) {
  
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
  if (is.null(varcov_mat) == TRUE) {
    varcov_resid <- summary(a)$covres
  } else {
    varcov_resid <- varcov_mat
  }
  
  
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
  
  #variance decomposition plot
  d <- fevd(a, n.ahead = burn_in_period)
  
  d <- do.call('rbind', sapply(names(coeff), function(x){
    f <- data.frame(d[[x]])
    f$var <- x
    f$t <- c(1:burn_in_period)
    f
  }, simplify = F))
  
  d <- melt(d, id.vars = c("var", "t"))
  
  var_decomp_plot <- ggplot(data = d, aes(x = t, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + facet_wrap(var ~ .) + 
    scale_fill_brewer(palette = 7, direction = 1)+ 
    theme_bw() + theme(legend.position = "bottom", panel.grid = element_blank()) + 
    labs(title = "Variance Decomposition (%)", x = "Forecast Length (t)", y = "% of Variance", fill = "")
  
  #plot burn-in-period analysis
  d <- do.call('rbind', sapply(1:nsim, function(x) {
    k <- data.frame(t(sim_list[[x]]))[-c(1:p),]
    k$sim_no <- x
    k$t <- as.numeric(row.names(k)) - p
    k
    
  }, simplify = F))
  
  d <- melt(d, id.vars = c("sim_no","t"))
  
  d_summ <- d %>%
    group_by(variable, t) %>%
    summarise(mean = mean(value, na.rm = T),
              sd = sd(value, na.rm = T)) %>%
    mutate(lcl = mean - sd,
           ucl = mean + sd) %>%
    select(variable,t, mean, lcl, ucl)
  
  d_summ <- melt(d_summ, id.vars = c("variable","t"))
  colnames(d_summ)[1] <- "var"
  d_summ$variable <- factor(d_summ$variable, levels = c("lcl", "mean", "ucl"))
  
  burn_in_plot <- ggplot(data = d_summ, aes(x = t, y = value, colour = variable)) + geom_line() +
    geom_point(size = 0.5, aes(colour = variable)) + 
    scale_colour_manual(values = c("red4", "green4", "red4")) + 
    facet_wrap(var ~. , scales = "free") + theme_bw() +
    labs(title = "Burn-in Period Analysis", x = "Forecast (t)", y = "Value", colour = "") + 
    theme(legend.position = "bottom")
  
  
  return_list <- list(model_object = a,
                      model_summary = summary(a),
                      sim_obs = sim_obs,
                      sim_full = sim_list,
                      mean = list(actual = apply(df,2, mean), simulated = apply(sim_obs,2,mean)),
                      sd = list(actual = apply(df,2, sd), simulated = apply(sim_obs,2,sd)),
                      varcov = list(actual = var(df), simulated = var(sim_obs)),
                      corr =list(actual = cor(df), simulated = cor(sim_obs)),
                      density_plot = density_plot,
                      var_decomp_plot = var_decomp_plot,
                      burn_in_plot = burn_in_plot)
  return_list
}

impulse_response_plot <- function(a, ci = 0.95, n.ahead = 10, runs = 100){
  coeff <- coefficients((a))
  cross <- expand.grid(impulse = names(coeff), response = names(coeff))
  ir <- do.call ('rbind', mapply(function(x, y) {
    d <- irf(x = a, impulse = x, response = y, ci = ci, n.ahead = n.ahead, runs = runs)
    d <- data.frame(cbind(unlist(d[["irf"]]), unlist(d[["Lower"]]), unlist(d[["Upper"]])))
    colnames(d) <- c("irf","lcl", "ucl")
    d$impulse <- x
    d$response <- y 
    d$t <- c(1:(n.ahead + 1))
    d
  },
  cross$impulse,  # names from first
  cross$response, SIMPLIFY = F))
  
  ir_long <- melt(ir, id.vars = c("impulse", "response", "t"))
  ir_long$variable <- factor(ir_long$variable, levels = c("lcl","irf", "ucl"))
  
  ggplot(data = ir_long, aes(x = t, y = value)) +
    geom_line(aes(linetype = variable, colour = variable)) + 
    geom_hline(yintercept = 0) + 
    scale_colour_manual(values = c("red4", "green4", "red4")) + 
    scale_linetype_manual(values = c(2,1,2)) + 
    facet_grid(impulse ~ response , scales = "free") + theme_bw() +
    labs(title = "Impulse Response", 
         subtitle = paste("CI =",ci),
         x = "Forecast (t)", y = "", linetype = "", colour = "") + 
    theme(legend.position = "bottom")
}


#check what provides the optimum lag of endogenous variables
varselect <- VARselect(y = au_me_data_4qtr_chg)
print(varselect)

restrict_mat <- matrix(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0), nrow=4, ncol=9, byrow=TRUE)

l <- var_simulation(df = au_me_data_log_8qtr_chg, p = 2, nsim = 10000, burn_in_period = 30,
                    type = "restrict", restrict_type = "ser", thresh = 2, restrict_mat = restrict_mat)                                 

impulse_response_plot(a = l$model_object, ci = 0.95, n.ahead = 30, runs = 100)
