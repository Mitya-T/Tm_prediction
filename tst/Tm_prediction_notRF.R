# =============================================================================
# G4 Tm PREDICTION PIPELINE (TOPOLOGY-AWARE VERSION)
# =============================================================================
# Features:
# - Canonical G4 motif parsing
# - Loop extraction
# - Loop permutation encoding
# - Tetrad estimation
# - Topology-aware features
# - Correlation filtering
# - Train/test split
# - Cross-validation
# - XGBoost regression
# - SHAP importance
# - Evaluation plots
# - Prediction output table
# =============================================================================
# setting the test directory:

setwd('/home/dimitri/Tm_prediction/tst/')
getwd()

# =============================================================================
# 1. LIBRARIES
# =============================================================================

packages <- c(
  "tidyverse",
  "stringr",
  "caret",
  "xgboost",
  "SHAPforxgboost",
  "pROC",
  "corrplot",
  "janitor",
  "readxl"
)

installed <- rownames(installed.packages())

# for(p in packages){
#   if(!p %in% installed){
#     install.packages(p)
#   }
# }


library(tidyverse)
library(stringr)
library(caret)
library(xgboost)
library(SHAPforxgboost)
library(corrplot)
library(janitor)
library(readxl)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

set.seed(42)

df <- read_excel("G4_All_updated.xlsx")

df <- df %>% clean_names()

df_sub <- df %>%
  select(
    type,
    name,
    sequence_5_3,
    conclusion,
    tm_c
  )

# =============================================================================
# 3. G4 FEATURE EXTRACTION
# =============================================================================

extract_g4_features <- function(sequence){
  
  # CHECK FOR NA FIRST
  if(is.na(sequence)){
    return(NULL)
  }
  
  seq <- toupper(sequence)
  
  seq_length <- nchar(seq)
  
  if(seq_length == 0){
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------
  # BASIC COMPOSITION
  # ---------------------------------------------------------------------------
  
  A_count <- str_count(seq, "A")
  T_count <- str_count(seq, "T")
  U_count <- str_count(seq, "U")
  G_count <- str_count(seq, "G")
  C_count <- str_count(seq, "C")
  
  gc_content <- (G_count + C_count) / seq_length
  
  g_richness <- G_count / seq_length
  
  # ---------------------------------------------------------------------------
  # G-RUN DETECTION
  # ---------------------------------------------------------------------------
  
  g_runs <- str_extract_all(seq, "G{2,}")[[1]]
  
  n_g_runs <- length(g_runs)
  
  g_run_lengths <- nchar(g_runs)
  
  max_g_run <- ifelse(length(g_run_lengths) > 0,
                      max(g_run_lengths), 0)
  
  mean_g_run <- ifelse(length(g_run_lengths) > 0,
                       mean(g_run_lengths), 0)
  
  tetrad_estimate <- ifelse(length(g_run_lengths) >= 4,
                            min(g_run_lengths[1:4]),
                            0)
  
  # ---------------------------------------------------------------------------
  # CANONICAL G4 PARSING
  # Pattern:
  # G3+N1-7G3+N1-7G3+N1-7G3+
  # ---------------------------------------------------------------------------
  
  canonical_pattern <- "(G{3,})(.{1,7}?)(G{3,})(.{1,7}?)(G{3,})(.{1,7}?)(G{3,})"
  
  match <- str_match(seq, canonical_pattern)
  
  if(!all(is.na(match))){
    
    loop1 <- nchar(match[,3])
    loop2 <- nchar(match[,5])
    loop3 <- nchar(match[,7])
    
    loop_vector <- c(loop1, loop2, loop3)
    
    loop_sum <- sum(loop_vector)
    
    loop_mean <- mean(loop_vector)
    
    loop_sd <- sd(loop_vector)
    
    loop_max <- max(loop_vector)
    
    loop_min <- min(loop_vector)
    
    loop_asymmetry <- loop_max - loop_min
    
    short_loops <- sum(loop_vector <= 2)
    
    # permutation encoding
    loop_pattern <- paste(loop_vector, collapse = "_")
    
  } else {
    
    loop1 <- NA
    loop2 <- NA
    loop3 <- NA
    
    loop_sum <- NA
    loop_mean <- NA
    loop_sd <- NA
    loop_max <- NA
    loop_min <- NA
    loop_asymmetry <- NA
    short_loops <- NA
    
    loop_pattern <- "unknown"
  }
  
  # ---------------------------------------------------------------------------
  # ADDITIONAL FEATURES
  # ---------------------------------------------------------------------------
  
  gg_freq <- str_count(seq, "GG") / max(seq_length - 1, 1)
  
  ggg_freq <- str_count(seq, "GGG") / max(seq_length - 2, 1)
  
  interruptions <- str_count(seq, "G[ATUC]G")
  
  pyrimidine_content <- (T_count + U_count + C_count) / seq_length
  
  purine_content <- (A_count + G_count) / seq_length
  
  return(data.frame(
    
    seq_length = seq_length,
    
    gc_content = gc_content,
    
    g_richness = g_richness,
    
    gg_freq = gg_freq,
    
    ggg_freq = ggg_freq,
    
    n_g_runs = n_g_runs,
    
    max_g_run = max_g_run,
    
    mean_g_run = mean_g_run,
    
    tetrad_estimate = tetrad_estimate,
    
    loop1 = loop1,
    loop2 = loop2,
    loop3 = loop3,
    
    loop_sum = loop_sum,
    loop_mean = loop_mean,
    loop_sd = loop_sd,
    loop_max = loop_max,
    loop_min = loop_min,
    
    loop_asymmetry = loop_asymmetry,
    
    short_loops = short_loops,
    
    interruptions = interruptions,
    
    purine_content = purine_content,
    
    pyrimidine_content = pyrimidine_content,
    
    loop_pattern = loop_pattern
    
  ))
  
}

# =============================================================================
# 4. FEATURE EXTRACTION FOR ALL SEQUENCES
# =============================================================================

cat("Extracting G4 topology-aware features...\n")

feature_list <- lapply(
  df_sub$sequence_5_3,
  extract_g4_features
)

features_df <- bind_rows(feature_list)

df_features <- bind_cols(df_sub, features_df)

which(is.na(df_sub$sequence_5_3))

# =============================================================================
# 5. FILTER TRAINING DATA
# =============================================================================

df_train <- df_features %>%
  filter(conclusion == "Forms a G4") %>%
  filter(!is.na(tm_c))

# =============================================================================
# 6. HANDLE CATEGORICAL VARIABLES
# =============================================================================

df_train$loop_pattern <- as.factor(df_train$loop_pattern)
df_train$type <- as.factor(df_train$type)

# =============================================================================
# 7. TRAIN / TEST SPLIT
# =============================================================================

set.seed(42)

train_index <- createDataPartition(
  df_train$tm_c,
  p = 0.8,
  list = FALSE
)

train_data <- df_train[train_index, ]
test_data  <- df_train[-train_index, ]

# =============================================================================
# 8. PREPARE MATRICES
# =============================================================================

x_train <- train_data %>%
  select(
    -tm_c,
    -name,
    -sequence_5_3,
    -conclusion
  )

x_test <- test_data %>%
  select(
    -tm_c,
    -name,
    -sequence_5_3,
    -conclusion
  )

# One-hot encoding
dummy_model <- dummyVars(~ ., data = x_train)

x_train_mat <- predict(dummy_model, x_train) %>% as.matrix()

x_test_mat <- predict(dummy_model, x_test) %>% as.matrix()

y_train <- train_data$tm_c
y_test  <- test_data$tm_c

# =============================================================================
# 9. REMOVE HIGHLY CORRELATED FEATURES
# =============================================================================

feature_variance <- apply(x_train_mat, 2, var)
zero_var_cols <- which(feature_variance == 0 | is.na(feature_variance))

x_train_mat <- x_train_mat[, -zero_var_cols]
x_test_mat  <- x_test_mat[, -zero_var_cols]

cor_matrix <- cor(x_train_mat)

high_cor <- findCorrelation(
  cor_matrix,
  cutoff = 0.90
)

if(length(high_cor) > 0){
  
  x_train_mat <- x_train_mat[, -high_cor]
  x_test_mat  <- x_test_mat[, -high_cor]
}

# =============================================================================
# 10. XGBOOST CROSS-VALIDATION
# =============================================================================

dtrain <- xgb.DMatrix(
  data = x_train_mat,
  label = y_train
)

params <- list(
  objective = "reg:squarederror",
  eta = 0.03,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8
)

cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 1000,
  nfold = 5,
  metrics = "rmse",
  early_stopping_rounds = 30,
  verbose = 1
)

# best_nrounds <- cv$evaluation_log[which.min(cv$evaluation_log$test_rmse_mean), "iter"]
best_nrounds <- cv$evaluation_log[which.min(cv$evaluation_log$test_rmse_mean)]$iter

best_nrounds

cat("\nBest number of rounds:", best_nrounds, "\n")

# =============================================================================
# 11. TRAIN FINAL MODEL
# =============================================================================

xgb_model <- xgboost(
  x = x_train_mat,
  y = y_train,
  nrounds = best_nrounds,
  objective = "reg:squarederror",
  learning_rate = 0.03,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8
)


# =============================================================================
# 12. PREDICTIONS
# =============================================================================

test_pred <- predict(
  xgb_model,
  x_test_mat
)

# =============================================================================
# 13. EVALUATION
# =============================================================================

rmse <- sqrt(mean((y_test - test_pred)^2))

r2 <- cor(y_test, test_pred)^2

mae <- mean(abs(y_test - test_pred))

cat("\n==============================\n")
cat("MODEL PERFORMANCE\n")
cat("==============================\n")

cat("RMSE:", round(rmse, 3), "°C\n")
cat("MAE :", round(mae, 3), "°C\n")
cat("R²  :", round(r2, 3), "\n")

# =============================================================================
# 14. EVALUATION PLOTS
# =============================================================================

eval_df <- data.frame(
  Actual = y_test,
  Predicted = test_pred
)

# Predicted vs actual
ggplot(eval_df, aes(Actual, Predicted)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "red"
  ) +
  
  theme_bw() +
  
  ggtitle("Predicted vs Actual Tm")

# Residuals
eval_df$residuals <- eval_df$Actual - eval_df$Predicted

ggplot(eval_df, aes(residuals)) +
  
  geom_histogram(
    bins = 30
  ) +
  
  theme_bw() +
  
  ggtitle("Residual Distribution")

# =============================================================================
# 15. FEATURE IMPORTANCE
# =============================================================================

importance_matrix <- xgb.importance(
  model = xgb_model
)

print(head(importance_matrix, 20))

xgb.plot.importance(
  importance_matrix[1:20]
)

# =============================================================================
# 16. SHAP VALUES
# =============================================================================

# Plot top 15 features by importance
print(xgb.plot.importance(
  importance_matrix[1:15]
))

# =============================================================================
# 17. PREDICT UNKNOWN Tm VALUES
# =============================================================================

# df_unknown <- df_features %>%
#   filter(conclusion == "Forms a G4") %>%
#   filter(is.na(tm_c))
# 
# if(nrow(df_unknown) > 0){
#   
#   # Convert to factor with SAME LEVELS as training data
#   df_unknown$loop_pattern <- factor(
#     df_unknown$loop_pattern, 
#     levels = levels(df_train$loop_pattern)
#   )
#   
#   # Also ensure type has same levels
#   df_unknown$type <- factor(
#     df_unknown$type,
#     levels = levels(df_train$type)
#   )
#   
#   x_unknown <- df_unknown %>%
#     select(
#       -tm_c,
#       -name,
#       -sequence_5_3,
#       -conclusion
#     )
#   
#   x_unknown_mat <- predict(
#     dummy_model,
#     x_unknown
#   ) %>%
#     as.matrix()
#   
#   # Match train columns
#   missing_cols <- setdiff(
#     colnames(x_train_mat),
#     colnames(x_unknown_mat)
#   )
#   
#   for(mc in missing_cols){
#     x_unknown_mat <- cbind(
#       x_unknown_mat,
#       0
#     )
#     colnames(x_unknown_mat)[ncol(x_unknown_mat)] <- mc
#   }
#   
#   x_unknown_mat <- x_unknown_mat[
#     ,
#     colnames(x_train_mat)
#   ]
#   
#   predicted_tm <- predict(
#     xgb_model,
#     x_unknown_mat
#   )
#   
#   df_unknown$predicted_tm <- predicted_tm
#   
#   cat("\n==============================\n")
#   cat("UNKNOWN Tm PREDICTIONS\n")
#   cat("==============================\n")
#   
#   print(
#     df_unknown %>%
#       select(
#         type,
#         name,
#         sequence_5_3,
#         predicted_tm
#       ) %>%
#       arrange(desc(predicted_tm))
#   )
# }

# Prepare features for prediction (all sequences)
df_features$loop_pattern <- factor(
  df_features$loop_pattern,
  levels = levels(df_train$loop_pattern)
)

df_features$type <- factor(
  df_features$type,
  levels = levels(df_train$type)
)

x_all <- df_features %>%
  select(
    -tm_c,
    -name,
    -sequence_5_3,
    -conclusion
  )

# One-hot encode using the training dummy model
x_all_mat <- predict(dummy_model, x_all) %>% as.matrix()

# Match train columns (add missing columns as 0)
missing_cols <- setdiff(
  colnames(x_train_mat),
  colnames(x_all_mat)
)

if(length(missing_cols) > 0){
  for(mc in missing_cols){
    x_all_mat <- cbind(x_all_mat, 0)
    colnames(x_all_mat)[ncol(x_all_mat)] <- mc
  }
}

# Reorder to match training columns
x_all_mat <- x_all_mat[, colnames(x_train_mat)]

# Predict Tm for all sequences
predicted_tm_all <- predict(xgb_model, x_all_mat)

# Add to original dataframe
df_sub$tm_predicted <- predicted_tm_all

cat("\n==============================\n")
cat("Tm PREDICTIONS ADDED\n")
cat("==============================\n")
cat("Added 'tm_predicted' column to df_sub\n")
cat("Sequences with measured Tm:\n")
print(nrow(df_sub %>% filter(!is.na(tm_c))))
cat("Sequences with predicted Tm (no measurement):\n")
print(nrow(df_sub %>% filter(is.na(tm_c))))


# =============================================================================
# 18. SAVE RESULTS
# =============================================================================

# write.csv(
#   df_unknown,
#   "G4_Tm_predictions.csv",
#   row.names = FALSE
# )

write.csv(
  df_sub,
  "G4_Tm_predictions.csv",
  row.names = FALSE
)


cat("\nResults saved:\n")
cat("  - G4_Tm_predictions.csv\n")

# =============================================================================
# END
# =============================================================================