library(ggplot2)
library(data.table)
library(dplyr)
library(readxl)
library(mgcv)
library(MASS)
library(plotly)
library(reshape2)

################################################################################
# Load Data
################################################################################

dat <- readxl::read_xlsx(file.path
                         (project_dir, "Data_PoF_sph_interaction_shared_v0.xlsx"),
                         col_names = TRUE,
                         sheet = "Col1_3D_JL",
                         .name_repair = "minimal")
setnames(dat,c("max(Rdot)/c_l","max(E_LKEnon_col)"),c("max_Rdot_c_l","max_E_LKEnon_col"))
# dat$PR <- NULL
dat$'dp/pbo' <- NULL
dat$'E_TEnon_col' <- NULL
dat <- reshape2::melt(data = dat, id.vars = c("D","RR","PR"))
setDT(dat)


###########################################################
## Linear Regression
###########################################################

# Initialize a data frame to store results
results <- data.frame(variable = character(), p1 = numeric(), p2 = numeric(), RMSE = numeric(), AIC = numeric(), BIC = numeric())

# List of variables (loop through each variable)
variables <- unique(dat$variable)

# List of p values (for optimization)
p_values <- seq(0.1, 4, by = 0.1) 
p_values <- c(p_values,-p_values)
dat_p <- expand.grid(p1 = p_values, p2 = p_values)
setDT(dat_p)

# Loop through each PR
for(PR_idx in unique(dat$PR)){
  # Loop through each variable
  for (var in variables) {
    
    cat("PR = ", PR_idx, ", var = ", var,"\n")
    # Filter data for each variable
    sub_dat_var <- dat[variable == var & PR == PR_idx, ]
    
    # Loop through each combination of p1 and p2
    for (i in 1:nrow(dat_p)) {
      # Fit the model with given p1 and p2
      p1 <- dat_p[i, p1]
      p2 <- dat_p[i, p2]
      lm_fit <- glm(value ~ I(RR^p1) * I(D^p2), data = sub_dat_var, family = "gaussian")
      
      # Calculate predicted values
      dat_grid <- expand.grid(D = unique(sub_dat_var$D), RR = unique(sub_dat_var$RR))
      dat_grid$value_pred <- predict(lm_fit, newdata = dat_grid)
      
      # Compare predicted values with actual values
      dat_predicted <- merge(sub_dat_var, dat_grid, by = c("D", "RR"))
      dat_predicted$error <- dat_predicted$value - dat_predicted$value_pred
      
      # Calculate RMSE
      RMSE <- sqrt(mean(dat_predicted$error^2))
      
      # Calculate AIC and BIC
      AIC_value <- AIC(lm_fit)
      BIC_value <- BIC(lm_fit)
      
      # Save results
      results <- rbind(results, data.frame(PR=PR_idx,variable = var, p1 = p1, p2 = p2, RMSE = RMSE, AIC = AIC_value, BIC = BIC_value))
    }
  }
}

opt_power <- results[,.SD[which.min(RMSE),.(p1,p2,RMSE)],by=c("PR","variable")]
opt_power

###########################################################
## Save the results
###########################################################

lm_coef_out_all <- data.table()
## Fit again with optimal parameters and create 3D plots
# Loop through each PR
for(PR_idx in unique(dat$PR)){
  # Loop through each variable to find best model based on RMSE
  for (var in variables) {
    cat("var = ", var,"\n")
    # Filter the results for the current variable
    sub_dat_var <- dat[variable == var & PR == PR_idx, ]
    var_results <- results[variable == var & PR==PR_idx, ]
    var_opt_power <- opt_power[variable == var & PR==PR_idx, ]
    
    # Refit the best model based on RMSE
    p1 <- var_opt_power$p1
    p2 <- var_opt_power$p2
    
    lm_fit <- glm(value ~ I(RR^{p1}) * I(D^{p2}), data = sub_dat_var, family = "gaussian")
    lm_coef <- data.frame(summary(lm_fit)$coefficients)
    
    # Create a data frame for each coefficient with the desired information
    lm_coef <- data.frame(summary(lm_fit)$coefficients)
    names(lm_coef) <- c("beta","se","t","p")
    lm_coef$coef <- rownames(lm_coef)
    setDT(lm_coef)
    lm_coef[,"LCI":=beta+qnorm(0.025)*se]
    lm_coef[,"UCI":=beta+qnorm(0.975)*se]
    lm_coef[, "label" := sprintf("%.3f [%.3f, %.3f] (p = %.3f)", beta, LCI, UCI, p)]
    
    lm_coef_out <- data.frame(t(lm_coef$label))
    names(lm_coef_out) <- lm_coef$coef
    lm_coef_out$PR <- PR_idx
    lm_coef_out$variable <- var
    lm_coef_out_all <- rbind(lm_coef_out_all,lm_coef_out)
    
    # Predict values using the fitted model
    dat_grid <- expand.grid(D = unique(sub_dat_var$D), RR = unique(sub_dat_var$RR))
    dat_grid$value_pred <- predict(lm_fit, newdata = dat_grid)
    
    # Merge actual and predicted values
    dat_predicted <- merge(sub_dat_var, dat_grid, by = c("D", "RR"))
    dat_predicted[ , error := value - value_pred]
    
    # Calculate RMSE for the refitted model
    RMSE <- sqrt(mean(dat_predicted$error^2))
    print(paste("Refitted RMSE: ", RMSE))
    
    # Create 3D plot for actual vs predicted values
    fig <- plot_ly() %>%
      add_trace(
        data = dat_predicted,
        x = ~D,
        y = ~RR,
        z = ~value,
        mode = 'markers',
        marker = list(color = 'blue', size = 5),
        name = 'Actual Values'
      ) %>%
      add_trace(
        data = dat_predicted,
        x = ~D,
        y = ~RR,
        z = ~value_pred,
        mode = 'markers',
        marker = list(color = 'red', size = 5),
        name = 'Predicted Values'
      ) %>%
      add_trace(
        x = dat_grid$D,
        y = dat_grid$RR,
        z = dat_grid$value_pred,
        type = 'mesh3d',
        opacity = 0.3,
        color = dat_grid$value_pred,
        colorscale = 'Viridis',
        showscale = TRUE,
        name = 'Fitted Surface'
      ) %>%
      layout(
        title = paste("3D Plot for", var),
        scene = list(
          xaxis = list(title = 'D'),
          yaxis = list(title = 'RR'),
          zaxis = list(title = 'Value')
        ),
        legend = list(
          x = 1.1,
          y = 0.5,
          traceorder = 'normal',
          font = list(size = 12),
          title = list(text = 'Legend')
        )
      )
    
    htmlwidgets::saveWidget(fig, file = file.path(project_dir,paste0(var,"_Plot_3D_PR_",PR_idx,".html")))
  }
}
