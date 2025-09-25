setwd("/Users/madisoncreach/Prediction_Project/Input_Datasets")

library(rrBLUP)

# Load phenotype and SNP data
phenotype_file <- "full_phenos_reduced_names_altered.csv"
full_phenos <- read.csv(phenotype_file)

expr <- read.csv("MI_NE_gene_expression.csv", row.names = 1, check.names = FALSE)
expr$genotype <- rownames(expr)
expr <- expr[, c(ncol(expr), 1:(ncol(expr)-1))]

# Load fold assignments and clean up
fold_file <- "fold_assignments_names_altered.csv"
fold_assignment <- read.csv(fold_file)
colnames(fold_assignment)[1] <- "genotype"
fold_assignment$fold <- fold_assignment$fold + 1

# Merge phenotype and fold data
phenos <- full_phenos
phenos <- phenos[ , -1]  # Drop GeneID column
colnames(phenos)[1] <- "genotype"
phenos <- merge(phenos, fold_assignment, by = "genotype")
phenos <- phenos[order(phenos$genotype),]

# Align exp to phenos
expr <- expr[expr$genotype %in% phenos$genotype, ]
expr <- expr[order(expr$genotype), ]
phenos <- phenos[phenos$genotype %in% expr$genotype, ]

# Clean markers
markers <- expr[, -1]  # remove genotype column
markers <- markers[, !(colnames(markers) %in% "genotype")]
markers <- markers[, colSums(is.na(markers)) == 0]

# Set parameters
ntaxa <- nrow(phenos)
kfolds <- 5
ind_train <- phenos$fold

# Prepare results dataframe
results_list <- list()

# Loop through all phenotype columns except genotype and fold
phenotype_cols <- setdiff(colnames(phenos), c("genotype", "fold"))

all_marker_effects <- list()

for (phenocol in phenotype_cols) {
  accuracy_rrBLUP <- numeric(kfolds)
  predicted_rrBLUP <- data.frame(predicted = numeric(), observed = numeric())
  
  for (i in 1:kfolds) {
    train_indices <- which(ind_train != i)
    test_indices <- which(ind_train == i)
    
    train_Pheno <- as.numeric(phenos[train_indices, phenocol])
    test_Pheno <- as.numeric(phenos[test_indices, phenocol])
    
    train_Markers <- as.matrix(markers[train_indices, ])
    test_Markers <- as.matrix(markers[test_indices, ])
    
    if (length(test_Pheno) == 0 || length(train_Pheno) == 0) {
      accuracy_rrBLUP[i] <- NA
      next
    }
    
    Pheno_answer <- mixed.solve(train_Pheno, Z = train_Markers, K = NULL, SE = FALSE, return.Hinv = FALSE)
    e <- as.matrix(Pheno_answer$u)
    marker_effects <- data.frame(
      Marker = colnames(train_Markers),
      Effect = as.vector(Pheno_answer$u)
    )
    all_marker_effects[[phenocol]] <- marker_effects
    
    pred_Pheno_valid <- test_Markers %*% e
    pred_Pheno <- pred_Pheno_valid + as.vector(Pheno_answer$beta)
    
    accuracy_rrBLUP[i] <- cor(pred_Pheno_valid, test_Pheno, use = "complete")
    
    
    predicted_rrBLUP <- rbind(predicted_rrBLUP, cbind(predicted = pred_Pheno_valid, observed = test_Pheno))
  }
  
  results_list[[phenocol]] <- list(
    fold_accuracies = accuracy_rrBLUP,
    mean_accuracy = mean(accuracy_rrBLUP, na.rm = TRUE)
  )
}

# Convert to summary dataframe
summary_df <- data.frame(
  phenotype = names(results_list),
  mean_accuracy = sapply(results_list, function(x) x$mean_accuracy)
)

print(summary_df)
write.csv(summary_df, "cv_fold_rrblup_exp_pcc_scores.csv", row.names = FALSE)

combined_effects_exp_df <- do.call(rbind, lapply(names(all_marker_effects), function(p) {
  df <- all_marker_effects[[p]]
  df$Phenotype <- p
  return(df)
}))

write.csv(combined_effects_exp_df, "/Users/madisoncreach/Prediction_Project/rrblup_feature_weights_exp.csv", row.names = FALSE)
