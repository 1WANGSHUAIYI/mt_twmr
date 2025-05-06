#################################################
Tissue1 <- read.csv("D:/GTEx/LD_clump/Adipose_Subcutaneous.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue2 <- read.csv("D:/GTEx/LD_clump/Brain_Frontal_Cortex_BA9.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue3 <- read.csv("D:/GTEx/LD_clump/Esophagus_Mucosa.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue4 <- read.csv("D:/GTEx/LD_clump/Heart_Left_Ventricle.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue5 <- read.csv("D:/GTEx/LD_clump/Kidney_Cortex.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue6 <- read.csv("D:/GTEx/LD_clump/Liver.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue7 <- read.csv("D:/GTEx/LD_clump/Lung.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue8 <- read.csv("D:/GTEx/LD_clump/Muscle_Skeletal.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue9 <- read.csv("D:/GTEx/LD_clump/Skin_Sun_Exposed_Lower_leg.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue10 <- read.csv("D:/GTEx/LD_clump/Thyroid.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue11 <- read.csv("D:/GTEx/LD_clump/Whole_Blood.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
data <- read.csv("~/wsy/GTEX/handoa.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
filtered_results <- list()
for (i in 1:11) {
  df <- get(paste0("Tissue", i))
  df$pval <- as.numeric(df$pval)  
  filtered_df <- df %>%
    filter(
      !grepl("^chrX", variant_id), 
      !grepl("^chrY", variant_id), 
      pval <= 0.005,              
      af >= 0.01 & af <= 0.99      
    )
  filtered_results[[paste0("Tissue", i, "_1")]] <- filtered_df
}
afiltered_results <- lapply(filtered_results, function(df) {
  df %>% select(Gene, rsid, slope, slope_se)
})                                        
for (i in seq_along(afiltered_results)) {
  prefix <- paste0("Tissue", i, "_")
  colnames(afiltered_results[[i]]) <- c("Gene", 
                                        "rsid", 
                                        paste0(prefix, "slope"), 
                                        paste0(prefix, "slope_se"))
}
a <- Reduce(function(x, y) merge(x, y, by = c("Gene", "rsid"), all = TRUE), afiltered_results)
b <- data %>% 
  dplyr::select(SNP, beta.exposure, se.exposure) %>%  
  left_join(a, by = c("SNP" = "rsid"))              
c <- b %>%
  filter(!is.na(Gene)) %>%          
  arrange(Gene)                    
slope_columns <- grep("^Tissue[1-9]_slope$|^Tissue10_slope$|^Tissue11_slope$", colnames(c), value = TRUE)

c$valid_rows <- apply(c[, slope_columns], 1, function(row) {
  non_na_values <- as.numeric(row[!is.na(row)])
  non_zero_values <- non_na_values[non_na_values != 0]
  
  if (length(non_zero_values) < 2) {
    return(FALSE)
  }
  return(TRUE)
})
D <- c
filtered_data <- D[D$valid_rows == TRUE, ]

filtered_data <- filtered_data[, -which(names(filtered_data) == "valid_rows")]



a1 <- filtered_data %>%
  dplyr::select(Gene, SNP, ends_with("_slope"))
a2 <- filtered_data %>%
  dplyr::select(Gene, SNP, ends_with("_slope_se"))  
a3 <- filtered_data %>%
  dplyr::select(Gene, SNP, beta.exposure)
a4 <- filtered_data %>%
  dplyr::select(Gene, SNP, se.exposure)

split_data <- split(a1, a1$Gene)
gene_matrices <- lapply(split_data, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices <- lapply(gene_matrices, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices <- lapply(gene_matrices, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data2 <- split(a2, a2$Gene)
gene_matrices2 <- lapply(split_data2, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  # 排除 rsid 和 Gene 列
  return(group_matrix)
})
gene_matrices2 <- lapply(gene_matrices2, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices2 <- lapply(gene_matrices2, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data3 <- split(a3, a3$Gene)
gene_matrices3 <- lapply(split_data3, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices3 <- lapply(gene_matrices3, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices3 <- lapply(gene_matrices3, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data4 <- split(a4, a4$Gene)
gene_matrices4 <- lapply(split_data4, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices4 <- lapply(gene_matrices4, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices4 <- lapply(gene_matrices4, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})




#####################################################################
K <- 11   
run_graph_lasso_with_l1 <- function(X, Y, W) {
  n <- nrow(X)
  k <- ncol(X)
  
  if (n == 1) {
    if (all(X == 0)) {
      return(list(
        estimated_beta = NA,
        se  = NA,
        pval = NA,
        method = "Wald Ratio",
        message = "SNP effect on expression is zero (X = 0)."
      ))
    } else {
      beta_wald <- as.numeric(Y) / as.numeric(X)
      se <- as.numeric(Yse) / abs(as.numeric(X))
      pval <- stats::pnorm(abs(beta_wald) / se, lower.tail=FALSE) * 2
      return(list(
        estimated_beta = beta_wald,
        se = se,
        pval = pval,
        method = "Wald Ratio",
        message = "Single SNP handled with Wald Ratio."
      ))
    }
  }
  
  initial_beta <- rep(0, k)  
  
  graph_lasso_penalty <- function(beta, W) {
    penalty <- 0
    for (i in 1:length(beta)) {
      for (j in 1:length(beta)) {
        if (i != j) {
          penalty <- penalty + W[i, j] * (beta[i] - beta[j])^2
        }
      }
    }
    return(penalty)
  }
  
  total_loss <- function(beta, X, Y, lambda1, lambda2, W) {
    residual <- Y - X %*% beta
    loss <- sum(residual^2)
    penalty_graph_lasso <- graph_lasso_penalty(beta, W)
    penalty_l1 <- sum(abs(beta))
    return(loss + lambda1 * penalty_graph_lasso + lambda2 * penalty_l1)
  }
  
  lambda1_vals <- c(0.003, 0.2)
  lambda2_vals <- c(0.008, 0.005)
  
  set.seed(42)
  folds <- sample(rep(1:5, length.out = n))
  best_loss <- Inf
  best_lambda1 <- NULL
  best_lambda2 <- NULL
  
  for (lambda1 in lambda1_vals) {
    for (lambda2 in lambda2_vals) {
      cv_loss <- 0
      for (fold in 1:5) {
        test_indices <- which(folds == fold)
        train_indices <- which(folds != fold)
        X_train <- X[train_indices, ]
        Y_train <- Y[train_indices]
        X_test <- X[test_indices, ]
        Y_test <- Y[test_indices]
        
        result_graph_l1 <- optim(par = initial_beta,
                                 fn = total_loss,
                                 X = X_train, Y = Y_train, lambda1 = lambda1, lambda2 = lambda2, W = W,
                                 method = "L-BFGS-B")
        
        beta_graph_l1 <- result_graph_l1$par
        predicted_Y <- X_test %*% beta_graph_l1
        cv_loss <- cv_loss + sum((Y_test - predicted_Y)^2)
      }
      if (cv_loss < best_loss) {
        best_loss <- cv_loss
        best_lambda1 <- lambda1
        best_lambda2 <- lambda2
      }
    }
  }
  
  result_graph_l1 <- optim(par = initial_beta,
                           fn = total_loss,
                           X = X, Y = Y, lambda1 = best_lambda1, lambda2 = best_lambda2, W = W,
                           method = "L-BFGS-B")
  
  beta_graph_l1 <- result_graph_l1$par
  fit_lasso <- glmnet(X, Y, alpha = 1, lambda = best_lambda2, start = beta_graph_l1)
  beta_final <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  result <- ifelse(beta_final == 0, beta_final, beta_graph_l1)
  
  return(list(estimated_beta = result, lambda1 = best_lambda1, lambda2 = best_lambda2))
}

results_list <- vector("list", 10)
for (i in 11:20) {
  X <- gene_matrices[[i]]
  Xse <- gene_matrices2[[i]]
  Y <- gene_matrices3[[i]]
  Yse <- gene_matrices4[[i]]
  p <- ncol(X)
  W <- matrix(1, nrow = K, ncol = K)
  # W <- W_list[[i]]
  result <- run_graph_lasso_with_l1(X, Y, W)
  results_list[[i]] <- result
}


##############################IVW#########################################
ivw_calculation_single_tissue <- function(beta_exposure, beta_outcome, se_exposure, se_outcome) {
  # Filter out NA values to ensure valid input    
  valid_indices <- which(!is.na(beta_exposure) & !is.na(beta_outcome) & !is.na(se_exposure) & !is.na(se_outcome))
  if (length(valid_indices) < 2) {    
    return(list(
      beta_ivw = NA,
      se_ivw = NA,
      p_ivw = NA,
      nsnp = NA,
      Q = NA,
      Q_df = NA,
      Q_pval = NA
    ))
  }
  
  beta_exposure <- beta_exposure[valid_indices]
  beta_outcome <- beta_outcome[valid_indices]
  se_outcome <- se_outcome[valid_indices]
  
  # Check for invalid data
  if (length(beta_exposure) < 2 || all(beta_exposure == 0) || all(beta_outcome == 0)) {
    return(list(
      beta_ivw = NA,
      se_ivw = NA,
      p_ivw = NA,
      nsnp = length(beta_exposure),
      Q = NA,
      Q_df = NA,
      Q_pval = NA
    ))
  }
  
  # Prepare weights and check validity
  weights <- 1 / se_outcome^2
  if (any(is.na(weights) | is.infinite(weights))) {
    return(list(
      beta_ivw = NA,
      se_ivw = NA,
      p_ivw = NA,
      nsnp = length(beta_exposure),
      Q = NA,
      Q_df = NA,
      Q_pval = NA
    ))
  }
  
  # Fit the model safely
  fit <- tryCatch(
    lm(beta_outcome ~ -1 + beta_exposure, weights = weights),
    error = function(e) return(NULL)
  )
  if (is.null(fit)) {
    return(list(
      beta_ivw = NA,
      se_ivw = NA,
      p_ivw = NA,
      nsnp = length(beta_exposure),
      Q = NA,
      Q_df = NA,
      Q_pval = NA
    ))
  }
  
  # Extract coefficients safely
  coefs <- coef(summary(fit))
  if (!"beta_exposure" %in% rownames(coefs)) {
    return(list(
      beta_ivw = NA,
      se_ivw = NA,
      p_ivw = NA,
      nsnp = length(beta_exposure),
      Q = NA,
      Q_df = NA,
      Q_pval = NA
    ))
  }
  b <- coefs["beta_exposure", "Estimate"]
  se <- coefs["beta_exposure", "Std. Error"] / min(1, summary(fit)$sigma)
  pval <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  
  # Calculate Q statistic
  Q_df <- length(beta_exposure) - 1
  Q <- summary(fit)$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  
  # Return results
  return(list(
    beta_ivw = b,
    se_ivw = se,
    p_ivw = pval,
    nsnp = length(beta_exposure),
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval
  ))
}

# Ensure all lists and variables are properly defined
result_list <- vector("list", length = 10531)  

for (g in 1:10531) {
  tissue_results <- vector("list") 
  # Check if gene_matrices dimensions are valid
  if (g <= length(gene_matrices) && g <= length(gene_matrices3)) {
    if (!is.null(dim(gene_matrices[[g]]))) {
      k <- ncol(gene_matrices[[g]])  # Use actual column number as k
      for (i in 1:k) { 
        if (i <= ncol(gene_matrices[[g]])) {  # Ensure i is within valid range
          X <- gene_matrices[[g]][, i] 
          Y <- gene_matrices3[[g]]  
          Xse <- gene_matrices2[[g]][, i]  
          Yse <- gene_matrices4[[g]]  
          # Ensure dimensions of outcome variables are valid
          if (length(Y) == nrow(gene_matrices[[g]]) && length(Yse) == nrow(gene_matrices2[[g]])) {
            result <- ivw_calculation_single_tissue(X, Y, Xse, Yse)
            tissue_results[[i]] <- result
          } else {
            cat("Warning: Outcome dimensions do not match for group g =", g, "and i =", i, "\n")
          }
        }
      }
      result_list[[g]] <- tissue_results
    } else {
      cat("Warning: Invalid matrix dimensions for group g =", g, "\n")
    }
  } else {
    cat("Warning: Index g =", g, "is out of bounds\n")
  }
}

##########################
p_ivw_matrix <- matrix(NA, nrow = length(result_list), ncol = 11)

for (i in seq_along(result_list)) {
  p_ivw_matrix[i, ] <- sapply(result_list[[i]], function(x) {
    if (is.list(x) && !is.null(x$p_ivw)) {
      x$p_ivw
    } else {
      NA 
    }
  })
}
positions <- which(p_ivw_matrix < 0.05, arr.ind = TRUE)
p_ivw_matrix <- cbind(p_ivw_matrix, output_matrix)


beta_ivw_matrix <- matrix(NA, nrow = length(result_list), ncol = 11)
for (i in seq_along(result_list)) {
  beta_ivw_matrix[i, ] <- sapply(result_list[[i]], function(x) {
    if (is.list(x) && !is.null(x$beta_ivw)) {
      x$beta_ivw
    } else {
      NA 
    }
  })
}
beta_ivw_matrix <- cbind(beta_ivw_matrix, output_matrix)

se_ivw_matrix <- matrix(NA, nrow = length(result_list), ncol = 11)

for (i in seq_along(result_list)) {
  se_ivw_matrix[i, ] <- sapply(result_list[[i]], function(x) {
    if (is.list(x) && !is.null(x$se_ivw)) {
      x$se_ivw
    } else {
      NA 
    }
  })
}
se_ivw_matrix <- cbind(se_ivw_matrix, output_matrix)


aa <- matrix(NA, nrow = 10531, ncol = 11)
for (i in 1:10531) {
  sublist <- results_list2[[i]]
  if (length(sublist) == 5) {
    pval <- results_list2[[i]][["pval"]]
    estimated_beta <- results_list2[[i]][["estimated_beta"]]
    pval[is.na(pval)] = 0
    if (any(pval < 0.05)) {
      beta_indices <- which(pval < 0.05)
      aa[i, beta_indices] <- estimated_beta[beta_indices]
    }
  } else if (length(sublist) == 3) {
    estimated_beta <- results_list2[[i]][["estimated_beta"]]
    aa[i, ] <- estimated_beta
  }
}
aa[is.infinite(aa)] <- 0
aa[is.na(aa)] <- 0

p <- apply(aa[, 1:11], 1, function(row) {
  if (all(row == 0)) 0 else 1
})

aa <- cbind(aa, p)


aaaaa <- matrix(nrow = length(results_list2), ncol = 2)
for (i in seq_along(results_list2)) {
  aaaaa[i, 1] <- paste0("list[", i, "]")
  aaaaa[i, 2] <- ifelse(length(results_list2[[i]]) == 5, 1, 2)
}
colnames(aaaaa) <- c("Name", "Category")


gene_names <- names(gene_matrices)
output_matrix <- matrix(nrow = length(gene_matrices), ncol = 2)
for (i in seq_along(gene_matrices)) {
  output_matrix[i, 1] <- gene_names[i]
  matrix_dims <- dim(gene_matrices[[i]])
  if (!is.null(matrix_dims) && matrix_dims[1] == 1 && matrix_dims[2] == 11) {
    output_matrix[i, 2] <- 1
  } else {
    output_matrix[i, 2] <- 2
  }
}
colnames(output_matrix) <- c("Name", "Category")
aaaaa<- cbind(aa, output_matrix)
m_df <- as.data.frame(aaaaa)
write.csv(m_df, file = "/home/user4/wsy/GTEX数据/LD_clump/jieguo0.00025.csv", row.names = FALSE)


##############IVWfdr#######
fdr_matrix <- t(apply(p_ivw_matrix, 1, function(row) {
  adjusted <- rep(NA, length(row))  
  non_na_indices <- !is.na(row)   
  adjusted[non_na_indices] <- p.adjust(row[non_na_indices], method = "fdr")
  return(adjusted)
}))
iivvww <- ifelse(!is.na(fdr_matrix) & fdr_matrix < 0.05, 1, 0)
new_columniivvww <- ifelse(rowSums(iivvww == 0) == ncol(iivvww), 0, 1)

iivvww <- cbind(iivvww, new_columniivvww)
iivvww <- cbind(iivvww, output_matrix)
m_df1 <- as.data.frame(iivvww)
write.csv(m_df1, file = "/home/user4/wsy/GTEX数据/LD_clump/jieguoiivvww.csv", row.names = FALSE)

m_df2 <-cbind(m_df1, m_df)
write.csv(m_df2, file = "/home/user4/wsy/GTEX数据/LD_clump/2.csv", row.names = FALSE)

#################################################################################################
library(data.table)
eQTLgen <- fread("/home/user4/wsy/GTEX数据/eQTLgen.txt", header = TRUE)
Tissue1 <- read.csv("/home/user4/wsy/GTEX数据/LD_clump/Brain_Frontal_Cortex_BA9.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue2 <- read.csv("/home/user4/wsy/GTEX数据/LD_clump/Lung.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue3 <- read.csv("/home/user4/wsy/GTEX数据/LD_clump/Muscle_Skeletal.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue4 <- read.csv("/home/user4/wsy/GTEX数据/LD_clump/Whole_Blood.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
data <- read.csv("~/wsy/GTEX数据/zdyyz.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
####bloodeqtl
setDT(eQTLgen)
eQTLgen_selected <- eQTLgen[, .(rsid = SNP, Gene, slope = beta, slope_se = se)]

filtered_results <- list()
for (i in 1:4) {
  df <- get(paste0("Tissue", i))
  df$pval <- as.numeric(df$pval)  
  filtered_df <- df %>%
    filter(
      !grepl("^chrX", variant_id), 
      !grepl("^chrY", variant_id), 
      pval <= 0.005,               
      af >= 0.01 & af <= 0.99      
    )
  filtered_results[[paste0("Tissue", i, "_1")]] <- filtered_df
}
afiltered_results <- lapply(filtered_results, function(df) {
  df %>% dplyr::select(Gene, rsid, slope, slope_se)
})                                        
for (i in seq_along(afiltered_results)) {
  prefix <- paste0("Tissue", i, "_")
  colnames(afiltered_results[[i]]) <- c("Gene", 
                                        "rsid", 
                                        paste0(prefix, "slope"), 
                                        paste0(prefix, "slope_se"))
}
#########
a <- Reduce(function(x, y) merge(x, y, by = c("Gene", "rsid"), all = TRUE), afiltered_results)
setDT(a)        
setDT(eQTLgen_selected)
a[, Gene_clean := sub("\\..*", "", Gene)]     
merged_data <- merge(a, eQTLgen_selected, by.x = c("rsid", "Gene_clean"), by.y = c("rsid", "Gene"), all.x = TRUE)
merged_data[, Gene_clean := NULL]
merged_data <- merged_data %>%
  mutate(Tissue4_slope = coalesce(as.numeric(Tissue4_slope), slope)) %>%
  select(-slope) 
merged_data <- merged_data %>%
  mutate(Tissue4_slope_se = coalesce(as.numeric(Tissue4_slope_se), slope_se)) %>%
  select(-slope_se) 
# merged_data <- merged_data %>%
#   select(-Tissue4_slope, -Tissue4_slope_se) %>%
#   rename(Tissue4_slope = slope, Tissue4_slope_se = slope_se)
a <- merged_data
b <- data %>% 
  dplyr::select(SNP, beta.exposure, se.exposure) %>%  
  left_join(a, by = c("SNP" = "rsid"))              
c <- b %>%
  filter(!is.na(Gene)) %>%         
  arrange(Gene)                    
slope_columns <- grep("^Tissue[1-4]_slope$|^Tissue10_slope$|^Tissue11_slope$", colnames(c), value = TRUE)

c$valid_rows <- apply(c[, slope_columns], 1, function(row) {
  non_na_values <- as.numeric(row[!is.na(row)])
  
  non_zero_values <- non_na_values[non_na_values != 0]
  
  if (length(non_zero_values) < 2) {
    return(FALSE)
  }
  sign_check <- all(sign(non_zero_values) == sign(non_zero_values[1]))
  
  return(sign_check)
  # return(TRUE)
})
D <- c
filtered_data <- D[D$valid_rows == TRUE, ]
filtered_data <- filtered_data[, -which(names(filtered_data) == "valid_rows")]



a1 <- filtered_data %>%
  dplyr::select(Gene, SNP, ends_with("_slope"))
a2 <- filtered_data %>%
  dplyr::select(Gene, SNP, ends_with("_slope_se"))  
a3 <- filtered_data %>%
  dplyr::select(Gene, SNP, beta.exposure)
a4 <- filtered_data %>%
  dplyr::select(Gene, SNP, se.exposure)

split_data <- split(a1, a1$Gene)
gene_matrices <- lapply(split_data, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices <- lapply(gene_matrices, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices <- lapply(gene_matrices, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data2 <- split(a2, a2$Gene)
gene_matrices2 <- lapply(split_data2, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices2 <- lapply(gene_matrices2, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices2 <- lapply(gene_matrices2, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data3 <- split(a3, a3$Gene)
gene_matrices3 <- lapply(split_data3, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices3 <- lapply(gene_matrices3, function(matrix) {
  matrix[] <- as.numeric(matrix) 
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices3 <- lapply(gene_matrices3, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})
split_data4 <- split(a4, a4$Gene)
gene_matrices4 <- lapply(split_data4, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("SNP", "Gene")])  
  return(group_matrix)
})
gene_matrices4 <- lapply(gene_matrices4, function(matrix) {
  matrix[] <- as.numeric(matrix)  
  matrix[is.na(matrix)] <- 0  
  return(matrix)
})
gene_matrices4 <- lapply(gene_matrices4, function(mat) {
  original_dim <- dim(mat)
  mat <- as.numeric(mat)
  if (!is.null(original_dim)) {
    mat <- matrix(mat, nrow = original_dim[1], ncol = original_dim[2])
  } else {
    mat <- matrix(mat, nrow = 1)
  }
  mat[is.na(mat)] <- 0
  return(mat)
})




#####################################################################
K <- 4   
run_graph_lasso_with_l1 <- function(X, Y, W) {
  n <- nrow(X)
  k <- ncol(X)
  
  if (n == 1) {
    if (all(X == 0)) {
      return(list(
        estimated_beta = NA,
        se  = NA,
        pval = NA,
        method = "Wald Ratio",
        message = "SNP effect on expression is zero (X = 0)."
      ))
    } else {
      beta_wald <- as.numeric(Y) / as.numeric(X)
      se <- as.numeric(Yse) / abs(as.numeric(X))
      pval <- stats::pnorm(abs(beta_wald) / se, lower.tail=FALSE) * 2
      return(list(
        estimated_beta = beta_wald,
        se = se,
        pval = pval,
        method = "Wald Ratio",
        message = "Single SNP handled with Wald Ratio."
      ))
    }
  }
  
  initial_beta <- rep(0, k)  
  
  graph_lasso_penalty <- function(beta, W) {
    penalty <- 0
    for (i in 1:length(beta)) {
      for (j in 1:length(beta)) {
        if (i != j) {
          penalty <- penalty + W[i, j] * (beta[i] - beta[j])^2
        }
      }
    }
    return(penalty)
  }
  
  total_loss <- function(beta, X, Y, lambda1, lambda2, W) {
    residual <- Y - X %*% beta
    loss <- sum(residual^2)
    penalty_graph_lasso <- graph_lasso_penalty(beta, W)
    penalty_l1 <- sum(abs(beta))
    return(loss + lambda1 * penalty_graph_lasso + lambda2 * penalty_l1)
  }
  
  lambda1_vals <- c(0.003, 0.2)
  lambda2_vals <- c(0.008, 0.005)
  
  set.seed(42)
  folds <- sample(rep(1:5, length.out = n))
  best_loss <- Inf
  best_lambda1 <- NULL
  best_lambda2 <- NULL
  
  for (lambda1 in lambda1_vals) {
    for (lambda2 in lambda2_vals) {
      cv_loss <- 0
      for (fold in 1:5) {
        test_indices <- which(folds == fold)
        train_indices <- which(folds != fold)
        X_train <- X[train_indices, ]
        Y_train <- Y[train_indices]
        X_test <- X[test_indices, ]
        Y_test <- Y[test_indices]
        
        result_graph_l1 <- optim(par = initial_beta,
                                 fn = total_loss,
                                 X = X_train, Y = Y_train, lambda1 = lambda1, lambda2 = lambda2, W = W,
                                 method = "L-BFGS-B")
        
        beta_graph_l1 <- result_graph_l1$par
        predicted_Y <- X_test %*% beta_graph_l1
        cv_loss <- cv_loss + sum((Y_test - predicted_Y)^2)
      }
      if (cv_loss < best_loss) {
        best_loss <- cv_loss
        best_lambda1 <- lambda1
        best_lambda2 <- lambda2
      }
    }
  }
  result_graph_l1 <- optim(par = initial_beta,
                           fn = total_loss,
                           X = X, Y = Y, lambda1 = best_lambda1, lambda2 = best_lambda2, W = W,
                           method = "L-BFGS-B")
  
  beta_graph_l1 <- result_graph_l1$par
  fit_lasso <- glmnet(X, Y, alpha = 1, lambda = best_lambda2, start = beta_graph_l1)
  beta_final <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
  result <- ifelse(beta_final == 0, beta_final, beta_graph_l1)
  
  return(list(estimated_beta = result, lambda1 = best_lambda1, lambda2 = best_lambda2))
}
results_list <- vector("list", 10)
for (i in 1:1242) {
  X <- gene_matrices[[i]]
  Xse <- gene_matrices2[[i]]
  Y <- gene_matrices3[[i]]
  Yse <- gene_matrices4[[i]]
  p <- ncol(X)
  W <- matrix(1, nrow = K, ncol = K)
  # W <- W_list[[i]]
  result <- run_graph_lasso_with_l1(X, Y, W)
  results_list[[i]] <- result
}




###########################
results_list2 <- results_list
aa <- matrix(NA, nrow = 1242, ncol = 4)
for (i in 1:1242) {
  sublist <- results_list2[[i]]
  if (length(sublist) == 5) {
    pval <- results_list2[[i]][["pval"]]
    estimated_beta <- results_list2[[i]][["estimated_beta"]]
    pval[is.na(pval)] = 0
    if (any(pval < 0.05)) {
      beta_indices <- which(pval < 0.05)
      aa[i, beta_indices] <- estimated_beta[beta_indices]
    }
  } else if (length(sublist) == 3) {
    estimated_beta <- results_list2[[i]][["estimated_beta"]]
    aa[i, ] <- estimated_beta
  }
}
aa[is.infinite(aa)] <- 0
aa[is.na(aa)] <- 0

p <- apply(aa[, 1:4], 1, function(row) {
  if (all(row == 0)) 0 else 1
})

aa <- cbind(aa, p)


aaaaa <- matrix(nrow = length(results_list2), ncol = 2)
for (i in seq_along(results_list2)) {
  aaaaa[i, 1] <- paste0("list[", i, "]")
  aaaaa[i, 2] <- ifelse(length(results_list2[[i]]) == 5, 1, 2)
}
colnames(aaaaa) <- c("Name", "Category")


gene_names <- names(gene_matrices)
output_matrix <- matrix(nrow = length(gene_matrices), ncol = 2)
for (i in seq_along(gene_matrices)) {
  output_matrix[i, 1] <- gene_names[i]
  matrix_dims <- dim(gene_matrices[[i]])
  if (!is.null(matrix_dims) && matrix_dims[1] == 1 && matrix_dims[2] == 4) {
    output_matrix[i, 2] <- 1
  } else {
    output_matrix[i, 2] <- 2
  }
}
colnames(output_matrix) <- c("Name", "Category")
aaaaa<- cbind(aa, output_matrix)
m_df <- as.data.frame(aaaaa)
write.csv(m_df, file = "/home/user4/wsy/GTEX数据/LD_clump/jieguozdyyz.csv", row.names = FALSE)







