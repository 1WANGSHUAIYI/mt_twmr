library("readxl")
data<- read_excel("D:/GTEx/GTEx_Analysis_v10_eQTL_updated/组织相似性.xlsx")
split_data <- split(data, data$Name)
new_list <- split_data[1:10]
gene_matrices <- lapply(new_list, function(group) {
  group_matrix <- as.matrix(group[, !names(group) %in% c("Name", "Description")])  # 排除 rsid 和 Gene 列
  return(group_matrix)
})
gene_matrices <- lapply(gene_matrices, function(matrix) {
  matrix[] <- as.numeric(matrix)  # 转换为数值型
  matrix[is.na(matrix)] <- 0  # 替换 NA 为 0
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
W_list <- list()
for (i in 1:5) {
  W_list[[i]] <- as.matrix(dist(t(gene_matrices[[i]]), method = "euclidean"))   
  W_list[[i]] <- 1 / (1 + W_list[[i]])
  diag(W_list[[i]]) <- 0
  W_list[[i]] <- W_list[[i]] / rowSums(W_list[[i]])
  W_list[[i]] <- (W_list[[i]]+ t(W_list[[i]])) / 2 
}
