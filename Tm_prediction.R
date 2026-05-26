setwd('/home/dimitri/Tm_prediction/')
getwd()

# install.packages("readxl")
library(readxl)

df <- read_excel("G4_All_updated.xlsx")

colnames(df)

install.packages("janitor")
library(janitor)

df <- df %>% clean_names()
colnames(df)

library(dplyr)

df_sub <- df %>% select(
  type,
  name,
  sequence_5_3,
  conclusion,
  tm_c
)

# install.packages("caret")
# install.packages("randomForest")
# install.packages("xgboost")


# Tm Prediction Model for G4-Forming RNA Sequences
# Train on sequences with known Tm values, predict for NA values

library(Biostrings)
library(dplyr)
library(tidyr)
library(caret)
library(randomForest)
library(xgboost)
library(tibble)

# ============================================================================
# 1. SEQUENCE FEATURE EXTRACTION
# ============================================================================

# Function to extract nucleotide composition features
extract_sequence_features <- function(sequence) {
  # Convert to uppercase
  seq <- toupper(sequence)
  
  # Basic composition
  A_count <- str_count(seq, "A")
  U_count <- str_count(seq, "U")
  G_count <- str_count(seq, "G")
  C_count <- str_count(seq, "C")
  
  seq_length <- nchar(seq)
  
  # If sequence is empty or invalid, return NA
  if (seq_length == 0) return(NULL)
  
  # GC content
  gc_content <- (G_count + C_count) / seq_length
  
  # AU content
  au_content <- (A_count + U_count) / seq_length
  
  # G-richness (important for G4)
  g_richness <- G_count / seq_length
  
  # Purine content (A + G)
  purine_content <- (A_count + G_count) / seq_length
  
  # Pyrimidine content (U + C)
  pyrimidine_content <- (U_count + C_count) / seq_length
  
  # Dinucleotide frequencies (selected G4-relevant ones)
  gg_freq <- str_count(seq, "GG") / (seq_length - 1)
  gc_freq <- str_count(seq, "GC") / (seq_length - 1)
  ga_freq <- str_count(seq, "GA") / (seq_length - 1)
  gu_freq <- str_count(seq, "GU") / (seq_length - 1)
  
  # G4 motif patterns (approximate: clusters of G's)
  # Count G-quartets (4 consecutive G's)
  gggg_count <- str_count(seq, "GGGG")
  
  # Runs of G's (using regex)
  g_runs <- gregexpr("G{3,}", seq)
  g_runs_count <- length(g_runs[[1]])
  
  return(list(
    seq_length = seq_length,
    gc_content = gc_content,
    au_content = au_content,
    g_richness = g_richness,
    purine_content = purine_content,
    pyrimidine_content = pyrimidine_content,
    gg_freq = gg_freq,
    gc_freq = gc_freq,
    ga_freq = ga_freq,
    gu_freq = gu_freq,
    gggg_count = gggg_count,
    g_runs = g_runs_count,
    a_count = A_count,
    u_count = U_count,
    g_count = G_count,
    c_count = C_count
  ))
}

# Function to extract features for all sequences
extract_all_features <- function(df) {
  features_list <- lapply(seq_along(df$sequence_5_3), function(i) {
    seq <- df$sequence_5_3[i]
    result <- tryCatch(
      extract_sequence_features(seq),
      error = function(e) NULL
    )
    
    # If extraction failed (NULL), create a list of NAs with correct structure
    if (is.null(result)) {
      result <- list(
        seq_length = NA,
        gc_content = NA,
        au_content = NA,
        g_richness = NA,
        purine_content = NA,
        pyrimidine_content = NA,
        gg_freq = NA,
        gc_freq = NA,
        ga_freq = NA,
        gu_freq = NA,
        gggg_count = NA,
        g_runs = NA,
        a_count = NA,
        u_count = NA,
        g_count = NA,
        c_count = NA
      )
    }
    return(result)
  })
  
  # Convert to dataframe - preserves all rows now
  features_df <- do.call(rbind, lapply(features_list, as.data.frame)) %>%
    as.data.frame()
  
  rownames(features_df) <- NULL  # Reset row names to be safe
  
  return(features_df)
}

# ============================================================================
# 2. DATA PREPARATION
# ============================================================================

# Extract features
cat("Extracting sequence features...\n")
sequence_features <- extract_all_features(df_sub)

# bad_row <- which(apply(sequence_features, 1, anyNA))
# bad_row

# Combine with original data
df_features <- df_sub %>%
  bind_cols(sequence_features)

# Filter for G4-forming sequences with known Tm values (training set)
df_train <- df_features %>%
  filter(conclusion == "Forms a G4") %>%
  filter(!is.na(tm_c)) %>%
  select(-type, -name, -sequence_5_3, -conclusion)  # Remove non-feature columns

# Prepare prediction data (G4-forming with unknown Tm)
df_predict <- df_features %>%
  filter(conclusion == "Forms a G4") %>%
  filter(is.na(tm_c)) %>%
  select(-type, -name, -sequence_5_3, -conclusion)

# Optional: Include "No G4" as negative examples if desired
df_no_g4 <- df_features %>%
  filter(conclusion == "No G4") %>%
  filter(!is.na(tm_c)) %>%
  select(-type, -name, -sequence_5_3, -conclusion)

cat("Training set size (Forms a G4 with known Tm):", nrow(df_train), "\n")
cat("Prediction set size (Forms a G4 with unknown Tm):", nrow(df_predict), "\n")
cat("Negative set size (No G4 with known Tm):", nrow(df_no_g4), "\n")

# xx <- df_sub[df_sub$conclusion == "No G4", ]


# ============================================================================
# 3. MODEL TRAINING
# ============================================================================

if (nrow(df_train) > 5) {
  cat("\nTraining Random Forest model...\n")
  
  # Remove rows with NaN or Inf
  df_train_clean <- df_train %>%
    filter(!if_any(everything(), ~is.nan(.) | is.infinite(.)))
  
  # Handle any remaining NA values
  df_train_clean <- df_train_clean %>% drop_na()
  
  # Separate Tm from features
  tm_values <- df_train_clean$tm_c
  X_train <- df_train_clean %>% select(-tm_c)
  
  # Train Random Forest
  set.seed(42)
  model_rf <- randomForest(
    x = X_train,
    y = tm_values,
    ntree = 100,
    mtry = sqrt(ncol(X_train)),
    importance = TRUE,
    na.action = na.omit
  )
  
  # Print feature importance
  cat("\nFeature Importance (Top 10):\n")
  importance_df <- data.frame(
    Feature = rownames(importance(model_rf)),
    Importance = importance(model_rf)[, "IncNodePurity"]
  ) %>%
    arrange(desc(Importance)) %>%
    head(10)
  print(importance_df)
  
  # Model performance on training set
  train_pred <- predict(model_rf, X_train)
  rmse_train <- sqrt(mean((tm_values - train_pred)^2))
  r2_train <- cor(tm_values, train_pred)^2
  
  cat("\nTraining Performance:\n")
  cat("RMSE:", round(rmse_train, 3), "\n")
  cat("R²:", round(r2_train, 3), "\n")
  
  # =========================================================================
  # 4. PREDICTION
  # =========================================================================
  
  if (nrow(df_predict) > 0) {
    cat("\nMaking predictions for unknown Tm values...\n")
    
    # Clean prediction data
    df_predict_clean <- df_predict %>%
      filter(!if_any(everything(), ~is.nan(.) | is.infinite(.))) %>%
      drop_na()
    
    if (nrow(df_predict_clean) > 0) {
      predicted_tms <- predict(model_rf, df_predict_clean)
      
      cat("Predictions made:", length(predicted_tms), "\n")
      cat("Predicted Tm range:", round(min(predicted_tms), 2), "-", 
          round(max(predicted_tms), 2), "°C\n")
    }
  }
  
} else {
  cat("Not enough training samples (n < 5). Cannot train model.\n")
}

# ============================================================================
# 5. ADD PREDICTIONS BACK TO ORIGINAL DATAFRAME
# ============================================================================

# Create output dataframe with predictions
df_output <- df_sub %>%
  mutate(tm_c_predicted = NA_real_)

# Add predicted values for sequences that had NA
if (nrow(df_train) > 5 && nrow(df_predict) > 0) {
  # Get indices of rows with NA tm_c that form G4
  na_indices <- which(df_sub$conclusion == "Forms a G4" & is.na(df_sub$tm_c))
  
  if (length(na_indices) > 0) {
    # Predict for all these rows
    features_for_pred <- df_features %>%
      slice(na_indices) %>%
      select(-type, -name, -sequence_5_3, -conclusion, -tm_c)
    
    features_for_pred_clean <- features_for_pred %>%
      filter(!if_any(everything(), ~is.nan(.) | is.infinite(.))) %>%
      drop_na()
    
    if (nrow(features_for_pred_clean) > 0) {
      predictions <- predict(model_rf, features_for_pred_clean)
      
      # Match back to original indices
      pred_indices <- na_indices[rownames(features_for_pred_clean) %>% 
                                   as.integer()]
      df_output$tm_c_predicted[pred_indices] <- predictions
    }
  }
}

# For sequences with known Tm, copy them to predicted column
df_output <- df_output %>%
  mutate(tm_c_predicted = coalesce(tm_c_predicted, tm_c))

cat("\n✓ Predictions complete!\n")
cat("Output dataframe 'df_output' contains 'tm_c_predicted' column\n")

# ============================================================================
# 6. OPTIONAL: MODEL DIAGNOSTICS
# ============================================================================

# Create a diagnostic plot (optional)
if (nrow(df_train) > 5) {
  pdf("tm_prediction_diagnostics.pdf", width = 12, height = 8)
  
  par(mfrow = c(2, 2))
  
  # Plot 1: Training predictions vs actual
  plot(tm_values, train_pred,
       main = "Training: Predicted vs Actual Tm",
       xlab = "Actual Tm (°C)", ylab = "Predicted Tm (°C)",
       pch = 16, col = rgb(0.2, 0.4, 0.8, 0.6))
  abline(0, 1, col = "red", lty = 2)
  
  # Plot 2: Residuals
  residuals <- tm_values - train_pred
  plot(train_pred, residuals,
       main = "Residuals vs Predicted",
       xlab = "Predicted Tm (°C)", ylab = "Residuals",
       pch = 16, col = rgb(0.8, 0.2, 0.2, 0.6))
  abline(h = 0, col = "red", lty = 2)
  
  # Plot 3: Feature importance
  importance_top <- importance(model_rf) %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    arrange(desc(IncNodePurity)) %>%
    head(10)
  
  barplot(importance_top$IncNodePurity,
          names.arg = importance_top$Feature,
          main = "Top 10 Feature Importance",
          ylab = "Importance (Node Purity)",
          las = 2, cex.names = 0.8)
  
  # Plot 4: Distribution of predictions
  hist(df_output$tm_c_predicted, breaks = 20,
       main = "Distribution of Predicted Tm Values",
       xlab = "Tm (°C)", ylab = "Frequency",
       col = rgb(0.2, 0.4, 0.8, 0.6))
  
  dev.off()
  cat("\n✓ Diagnostic plots saved to 'tm_prediction_diagnostics.pdf'\n")
}

# ============================================================================
# 7. SUMMARY STATISTICS
# ============================================================================

cat("\n" %+% strrep("=", 70) %+% "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 70) %+% "\n")
cat("Total sequences in df_sub:", nrow(df_sub), "\n")
cat("G4-forming sequences:", sum(df_sub$conclusion == "Forms a G4"), "\n")
cat("Non-G4 sequences:", sum(df_sub$conclusion == "No G4"), "\n")
cat("G4-forming with known Tm:", nrow(df_train), "\n")
cat("G4-forming with unknown Tm:", nrow(df_predict), "\n")

if (nrow(df_output) > 0) {
  cat("\nFinal Tm Statistics (all sequences):\n")
  cat("Mean Tm:", round(mean(df_output$tm_c_predicted, na.rm = TRUE), 2), "°C\n")
  cat("Median Tm:", round(median(df_output$tm_c_predicted, na.rm = TRUE), 2), "°C\n")
  cat("SD Tm:", round(sd(df_output$tm_c_predicted, na.rm = TRUE), 2), "°C\n")
  cat("Range:", round(min(df_output$tm_c_predicted, na.rm = TRUE), 2), "-",
      round(max(df_output$tm_c_predicted, na.rm = TRUE), 2), "°C\n")
}

cat("\n" %+% strrep("=", 70) %+% "\n")




