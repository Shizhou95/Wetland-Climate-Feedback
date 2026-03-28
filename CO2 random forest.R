# Wetland CO2 random forest model
# Code by: Shizhou Ma 
# Version Date: October. 2025
# This script fits a random forest model to wetland CO2 sequestration observations and generates:
#   1. Training vs testing prediction plots across repeated model runs
#   2. Variable importance summaries across repeated model runs
#   3. Partial dependence plots across repeated model runs
#Required package to run the code: tidyverse V2.0.0, tidymodels V1.2.0, ranger V0.17.0, 
#lubridate V1.9.4, vip V0.4.1, ggpubr V0.6.0, DALEXtra V2.3.0, doParallel V1.0.17, ggtexr V0.1.2, 
#patchwork V1.3.0, ggpmisc V0.6.1

# List of all packages needed
package_list <- c('tidyverse', 'tidymodels', 'ranger', 'lubridate', 'vip', 'ggpubr', 'DALEXtra','doParallel' , 'ggtext', 'patchwork', 'ggpmisc')
# Load all the packages
lapply(package_list, require, character.only = TRUE)

# SET YOUR PROJECT FOLDER HERE
project_dir <- "FILL-IN"


# Define input file in the project folder
input_file <- file.path(project_dir, "CO2input.csv")
# Set working directory
setwd(project_dir)
#read inputdata
inputdata <- read.csv(input_file)
inputdata <- inputdata %>%
  select(-Latitude, -Longitude)
inputdata <- inputdata %>% select(-BD) 

#filter out negative observations
dataformodel <- inputdata %>%
  filter(Cseq >= 0)%>%
  drop_na()


# log transform skewed dependent variable, shift a bit from 0 as well
vars_to_log <- c("Cseq")
#rename the log transformed variable
log_model <- dataformodel %>%
  mutate(across(.cols = all_of(vars_to_log), ~log(.x + .1))) %>%   
  rename_with( ~ str_c("Log_", vars_to_log), .cols = all_of(vars_to_log) ) 

#for future and Paleo
#log_model <- log_model %>% select(-NDVI, -IP)


# Define a function for ONE run with permutation importance ----
predict_RF_once <- function(df) {
  wetland_id_col <- df$Wetland.ID
  year_col       <- df$Year
  lat_col        <- df$Latitude
  long_col       <- df$Longitude
  group_col      <- df$group
  
  df_model <- df %>% select(-Wetland.ID, -Year, -group)
  
  split    <- initial_split(df_model, prop = 0.8)
  train_df <- training(split)
  test_df  <- testing(split)
  
  recipe_train <- recipe(Log_Cseq ~ ., data = train_df) %>% 
    step_center(all_predictors()) %>%
    step_scale(all_predictors())
  
  n.cores <- parallel::detectCores() - 1
  
  rf_mod <- rand_forest(
    mtry  = 4,
    trees = 150,
    min_n = 5
  ) %>%
    set_mode("regression") %>%
    set_engine("ranger", 
               importance   = "permutation",
               num.threads  = n.cores, 
               keep.inbag   = TRUE)
  
  wf <- workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(recipe_train)
  
  fit_wf <- fit(wf, data = train_df)
  
  # Permutation-based importance
  importance_vals <- fit_wf$fit$fit$fit$variable.importance
  
  # Generate predictions for the training set
  preds_train <- predict(fit_wf, new_data = train_df) %>%
    bind_cols(train_df) %>%
    mutate(dataset = "train")
  
  # Generate predictions for the testing set
  preds_test <- predict(fit_wf, new_data = test_df) %>%
    bind_cols(test_df) %>%
    mutate(dataset = "test")
  
  # Calculate performance metrics on training data
  r2_train <- 1 - sum((preds_train$Log_Cseq - preds_train$.pred)^2) /
    sum((preds_train$Log_Cseq - mean(preds_train$Log_Cseq))^2)
  lm_fit_train <- lm(Log_Cseq ~ .pred, data = preds_train)
  aic_train <- AIC(lm_fit_train)
  
  # Calculate performance metrics on testing data
  r2_test <- 1 - sum((preds_test$Log_Cseq - preds_test$.pred)^2) /
    sum((preds_test$Log_Cseq - mean(preds_test$Log_Cseq))^2)
  lm_fit_test <- lm(Log_Cseq ~ .pred, data = preds_test)
  aic_test <- AIC(lm_fit_test)
  
  # Return a list with separate outputs for training and testing
  list(
    importance = importance_vals,
    r2_train   = r2_train,
    aic_train  = aic_train,
    preds_train = preds_train,
    r2_test    = r2_test,
    aic_test   = aic_test,
    preds_test  = preds_test
  )
}

# Run the RF model 100 times and summarize variable importance ----
set.seed(123)
n_runs <- 100
rf_runs <- replicate(n_runs, predict_RF_once(log_model), simplify = FALSE)

# Extract predictions for training and testing from each run
all_preds_train <- lapply(rf_runs, function(run) run$preds_train)
all_preds_test  <- lapply(rf_runs, function(run) run$preds_test)

# Combine all predictions into separate data frames
all_preds_train_df <- bind_rows(all_preds_train, .id = "Run")  # 'Run' is a factor from "1" to "100"
all_preds_test_df  <- bind_rows(all_preds_test, .id = "Run")

# Calculate median R2 values across all runs
median_r2_train <- median(sapply(rf_runs, function(run) run$r2_train))
median_r2_test  <- median(sapply(rf_runs, function(run) run$r2_test))

# Create the plot for training data
plot_train <- ggplot(all_preds_train_df, aes(x = exp(.pred), y = exp(Log_Cseq))) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(group = Run),
    color = "blue",
    linewidth = 0.05,
    alpha = 0.05
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 0.3, 
    y = 1.5,
    label = sprintf("Training: R² median = %.2f", median_r2_train),
    hjust = 0,
    vjust = 1,
    color = "black",
    size = 4
  ) +
  scale_x_log10(
    breaks = c(0.5, 1, 2),
    labels = c("0.5", "1.0", "2.0"),
    limits = c(0.3, 2.5)  # Extended limits to ensure visibility
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2),
    labels = c("0.5", "1.0", "2.0"),
    limits = c(0.3, 2.5)
  ) +
  labs(
    title = "Training Data: Predictions vs Observations",
    x = "Cseq predictions (log scale)",
    y = "Cseq observations (log scale)"
  ) +
  theme_bw()

# Create the plot for testing data
plot_test <- ggplot(all_preds_test_df, aes(x = exp(.pred), y = exp(Log_Cseq))) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(group = Run),
    color = "red",
    linewidth = 0.05,
    alpha = 0.05
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 0.3, 
    y = 1.5,
    label = sprintf("Testing: R² median = %.2f", median_r2_test),
    hjust = 0,
    vjust = 1,
    color = "black",
    size = 4
  ) +
  # scale_x_log10(
  #   breaks = c(0.5, 1, 2),
  #   labels = c("0.5", "1.0", "2.0"),
  #   limits = c(0.3, 2.5)  # Extended limits to ensure visibility
  # ) +
  # scale_y_log10(
  #   breaks = c(0.5, 1, 2),
  #   labels = c("0.5", "1.0", "2.0"),
  #   limits = c(0.3, 2.5)
  # ) +
  scale_x_log10(labels = scales::number) +
  scale_y_log10(labels = scales::number) +
  labs(
    title = "Testing Data: Predictions vs Observations",
    x = "Cseq predictions (log scale)",
    y = "Cseq observations (log scale)"
  ) +
  theme_bw()

# Display the plots separately
plot_train
plot_test


#code for generating variable importance plot
# Extract variable importance from each run
all_importance <- map(rf_runs, ~ .x$importance)
importance_tibbles <- map(all_importance, ~ enframe(.x, "Variable", "Importance"))
importance_allruns <- bind_rows(importance_tibbles, .id = "Run")

# Summarize
importance_summary <- importance_allruns %>%
  group_by(Variable) %>%
  summarise(
    Mean_Importance = mean(Importance, na.rm = TRUE),
    SD_Importance   = sd(Importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Importance))

# Plot with error bars
ggplot(importance_summary, 
       aes(x = reorder(Variable, Mean_Importance), y = Mean_Importance)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean_Importance - SD_Importance,
                    ymax = Mean_Importance + SD_Importance),
                width = 0.4) +
  coord_flip() +
  labs(title = "",
       x = "Predictors", 
       y = "Variable Importance") +
  theme_minimal(base_size = 14)



# partial dependence plot----
get_pdp_once <- function(df) {
  
  #Keep track of non-predictor columns just to remove them from the model:
  wetland_id_col <- df$Wetland.ID
  year_col       <- df$Year
  lat_col        <- df$Latitude
  long_col       <- df$Longitude
  group_col      <- df$group
  
  #Subset strictly to the modeling columns
  df_model <- df %>% select(-Wetland.ID, -Year, -group)
  
  #Train/test split
  set.seed(sample.int(1e6, 1))  # new seed each time
  split    <- initial_split(df_model, prop = 0.8)
  train_df <- training(split)
  test_df  <- testing(split)
  
  #Recipe
  recipe_train <- recipe(Log_Cseq ~ ., data = train_df) %>% 
    step_center(all_predictors()) %>%
    step_scale(all_predictors())
  
  #Define the Random Forest Model
  n.cores <- parallel::detectCores() - 1
  rf_mod <- rand_forest(
    mtry  = 4, 
    trees = 150,
    min_n = 5
  ) %>%
    set_mode("regression") %>%
    set_engine("ranger", 
               importance   = "permutation", 
               num.threads  = n.cores, 
               keep.inbag   = TRUE)
  
  #Workflow & fit
  wf <- workflow() %>%
    add_recipe(recipe_train) %>%
    add_model(rf_mod)
  
  fit_wf <- fit(wf, data = train_df)
  
  #Create a DALEX explainer using training data
  # Because partial dependence is typically computed on the training (or a representative) set.
  explainer_rf <- explain_tidymodels(
    fit_wf,
    data = select(train_df, -Log_Cseq),
    y    = train_df$Log_Cseq,
    label = "RF"
  )
  
  #Compute partial dependence for all variables using `model_profile`
  #N = NULL => use all rows (or set N = some smaller number to speed up)
  pdp_rf <- model_profile(explainer_rf, N = NULL)
  
  # This object has multiple components; the partial dependence curves
  # are typically in `pdp_rf$agr_profiles` (averaged) or `pdp_rf$individual_profiles` (each row).
  pdp_data <- pdp_rf$agr_profiles %>%
    select(
      Variable = `_vname_`,
      x        = `_x_`,
      y_hat    = `_yhat_`
    )
  
  # return the aggregated partial dependence for all variables in a single tibble
  pdp_data
}

set.seed(123)  # For reproducibility across the whole experiment

n_runs <- 100

pdp_list <- map(1:n_runs, function(i) {
  pdp_data_i <- get_pdp_once(log_model) %>%
    mutate(iteration = i)
  pdp_data_i
})

# Combine them into one big data frame
pdp_allruns <- bind_rows(pdp_list)


pdp_summary <- pdp_allruns %>%
  group_by(Variable, x) %>%
  summarize(
    y_hat_mean = mean(y_hat, na.rm = TRUE),
    .groups = "drop"
  )


# Invert the log transformation for x (predictor) and y_hat (partial dependence)
pdp_allruns <- pdp_allruns %>%
  mutate(
    x_orig     = x,           # Original predictor scale
    y_hat_orig = exp(y_hat) - 0.1        # Original response scale (Cseq)
  )

pdp_summary <- pdp_summary %>%
  mutate(
    x_orig         = x,           # Original predictor scale
    y_hat_mean_orig = exp(y_hat_mean) - 0.1   # Mean partial dependence on original scale
  )

# Check for any negative values in the transformed data
summary(pdp_allruns$x_orig)
summary(pdp_allruns$y_hat_orig)

# Optionally, handle any negative values if they arise
# For example, set them to NA or a minimal positive value
pdp_allruns <- pdp_allruns %>%
  mutate(
    x_orig     = ifelse(x_orig < 0, NA, x_orig),
    y_hat_orig = ifelse(y_hat_orig < 0, NA, y_hat_orig)
  )

pdp_summary <- pdp_summary %>%
  mutate(
    x_orig         = ifelse(x_orig < 0, NA, x_orig),
    y_hat_mean_orig = ifelse(y_hat_mean_orig < 0, NA, y_hat_mean_orig)
  )

pdp_allruns <- pdp_allruns %>%
  group_by(Variable) %>%
  filter(x_orig >= quantile(x_orig, probs = 0.05, na.rm = TRUE),
         x_orig <= quantile(x_orig, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

pdp_summary <- pdp_summary %>%
  group_by(Variable) %>%
  filter(x_orig >= quantile(x_orig, probs = 0.05, na.rm = TRUE),
         x_orig <= quantile(x_orig, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

ggplot() +
  #All 100 PDP lines in faint blue
  geom_line(
    data = pdp_allruns,
    aes(x = x_orig, y = y_hat_orig, group = iteration),
    alpha = 0.2,
    color = "blue"
  ) +
  
  # # Mean Partial Dependence in solid red
  # geom_line(
  #   data = pdp_summary,
  #   aes(x = x_orig, y = y_hat_mean_orig),
  #   color = "red",
  #   size  = 1
  # ) +
  
  # Facet by Variable with free scales for both x and y axes
  facet_wrap(~ Variable, scales = "free", ncol=4) +
  
  # Styling
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme(legend.position = "none")







