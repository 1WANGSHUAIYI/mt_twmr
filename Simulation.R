N = c(50000)
Nexp1 = c(500)
Nexp2 = c(500)
K1 = 5
K2 = 5
pattern1 = 0
pattern2 = 0
G = 50
M = 15
M_gene = 12
prop_UHP = 1
h2_CHP = 0
h2_exp = 0.3
h2_meth = 0.3
h2_out_IV = c(0.05,0.1,0.15)
h2_asso = 0
h2_out = c(0.15)/10
h2_U=0.1
h2_U_exp=0.1
corr_gene = c(0.1)
corr_gene_cpg = 0.1
tissue_sharing = c(0.2,2)
seed = 1:20+20
rho = c(0)
threshold = c(1e-2)
cooccur = 0.5
density = 0.85
k_group1 = c(1)
k_group2 = c(1)
a_init = 0.1
b_init = 0.1
b_beta = 0.03
library(tidyverse)
library(pROC)
library(dplyr)
library(mr.raps)
library(MRPRESSO)
library(penalized)
library(GRAPPLE)
library(pbapply)
library(MendelianRandomization)
library(MRcML)
library(optparse)
library(data.table)
library(MASS)
library(glmnet)
library(quantreg)
library(robustbase)
library(mintMR)
library(mvtnorm)
library(MendelianRandomization)
par_data_generation_multiple_data <- function(args,ncl,run_MR=T,overlap=T) {
  token <- paste0(paste(names(args),args,sep = "_"),collapse = ".")
  
  N <- as.numeric(args["N"])
  Nexp1 <- as.numeric(args["Nexp1"])
  Nexp2 <- as.numeric(args["Nexp2"])
  
  K1 <- as.numeric(args["K1"])
  K2 <- as.numeric(args["K2"])
  
  pattern1 = as.numeric(args["pattern1"]) # high_prob
  pattern2 = as.numeric(args["pattern2"])
  
  G = as.numeric(args["G"]) # number of genes
  M = as.numeric(args["M"]) # number of IVs total
  M_gene = as.numeric(args["M_gene"]) #  number of IVs each gene
  
  prop_UHP = as.numeric(args["prop_UHP"]) # proportion of IV with UHP
  h2_CHP = as.numeric(args["h2_CHP"])
  
  h2_exp = as.numeric(args["h2_exp"])
  h2_meth = as.numeric(args["h2_meth"])
  
  h2_out_IV = as.numeric(args["h2_out_IV"])  # proportion of variation in outcome that can be explained by IVs
  h2_out = as.numeric(args["h2_out"])
  h2_asso = as.numeric(args["h2_asso"])
  h2_U = as.numeric(args["h2_U"])
  h2_U_exp = as.numeric(args["h2_U_exp"])
  
  corr_gene = as.numeric(args["corr_gene"])
  corr_gene_cpg = as.numeric(args["corr_gene_cpg"])
  
  tissue_sharing <- as.numeric(args["tissue_sharing"])
  
  cooccur <- as.numeric(args["cooccur"]) # change it to control the number of IVs in exposure
  seed = as.numeric(args["seed"])
  
  set.seed(seed)
  rho = as.numeric(args["rho"])
  threshold = as.numeric(args["threshold"])
  
  fixedBeta = as.numeric(args["fixedBeta"])
  fixedIVBeta = 0 # as.numeric(args["fixedIVBeta"])
  
  density <- as.numeric(args["density"]) # 0.95
  k_group1 <- as.numeric(args["k_group1"]) # 3
  k_group2 <- as.numeric(args["k_group2"]) # 3
  
  
  b_beta <- as.numeric(args["b_beta"])
  set.seed(seed)
  
  
  result <- true_beta <- matrix(NA,G,K1+K2)
  beta_Xlist <- beta_Ylist <- se_Xlist <- se_Ylist <- pval_list <- list()
  
  methods <- c("metaivw","ivw","egger","median","lasso","presso","grapple","conmix","raps","cml","cue","mix","clust",
               "mvivw","mvegger","mvmedian","mvlasso","mvcml","oracle","mintMR", "mvrobust")
  
  pval <- est <- vector(mode = "list", length = length(methods)) %>%
    lapply(function(x) matrix(NA, nrow = G, ncol = K1+K2))
  names(pval) <- names(est) <- methods
  
  # latent status across tissues are NOT correlated
  density <- rep(1-density,K1+K2)
  if(tissue_sharing <= 1e10) {
    disease_rel_exp <- K1 - 1:k_group1 + 1
    disease_rel_DNAm <- K1 + K2 - 1:k_group2 + 1
    density[-c(disease_rel_exp,disease_rel_DNAm)] <- 0.01 # used to be a fixed 0.01 in version v10 and before
    
    if(b_beta != 1) {# not the default value
      density[base::setdiff(1:K1,disease_rel_exp)] <- b_beta
    }
  }
  latent_status <- sapply(density,function(x) rbinom(G, 1, prob = x))
  
  
  
  
  if(K1 == 0) {
    colnames(latent_status) <- c(paste0("cpg",1:K2))
  }
  if(K2 == 0) {
    colnames(latent_status) <- c(paste0("gene",1:K1))
  }
  if(K1 != 0 & K2 != 0) {
    colnames(latent_status) <- c(paste0("gene",1:K1),paste0("cpg",1:K2))
  }
  
  rownames(latent_status) <- paste0("tissue",1:G)
  
  parres <- pblapply(1:G, function(repind) {
    set.seed(seed+repind)
    maf = runif(M,0.05,0.49) # minor allele frequency of IVs
    beta_causal_status <- latent_status[repind,]
    
    beta_causal <- numeric(K1+K2)
    
    adj_matrix <- matrix(0,nrow = M, ncol = K1+K2)
    adj_matrix <- apply(adj_matrix,2,function(x) {
      if(cooccur > 1){
        x[sample(1:cooccur,M_gene)] <- 1
      } else {
        x[sample(1:length(x),M_gene)] <- 1
      }
      return(x)
    })
    
    # step 0: generate Confounder
    U <- rnorm(N)
    
    # step 1: generate IVs
    IVs <- geno_generation(maf,M,rho = rho,n = N) %>% scale()
    LD <- cor(IVs)
    
    if(h2_CHP > 0) {
      beta <- rnorm(ncol(IVs),mean = 0,sd = sqrt(h2_CHP/ncol(IVs)))
      error <- rnorm(nrow(IVs),mean = 0,sd = sqrt(1 - h2_CHP))
      U <- IVs %*% beta + error
    }
    # step 2: generate gene expression
    
    Exposure <- generate_exposure_mult_consistency(N, Nexp1, Nexp2, IVs, h2_exp, h2_meth,M,
                                                   M_gene, K1, K2, adj_matrix,
                                                   diag_corr = corr_gene, off_diag_corr = corr_gene)
    
    Exposure <- add_association(Exps = Exposure[,1:K1], DNAm = Exposure[,setdiff(1:(K1+K2),1:K1)], K1 = K1, K2 = K2, h2_asso = h2_asso)
    Exposure <- add_correlation_structure(Exps = Exposure$Exps, DNAm = Exposure$DNAm,
                                          K1 = K1, K2 = K2,
                                          h2_exp = h2_exp, h2_meth = h2_meth,h2_U_exp = h2_U_exp,h2_asso = h2_asso,
                                          diag_corr = corr_gene_cpg, off_diag_corr = corr_gene_cpg)
    Exps <- Exposure$Exps; DNAm <- Exposure$DNAm
    Exposure <- cbind(Exps,DNAm)
    Exposure <- scale(Exposure)
    Exposure <- Exposure * sqrt(1 - h2_U_exp)
    diag(var(Exposure))
    beta_U_exp <- rnorm(ncol(Exposure),0,h2_U_exp)
    beta_U_exp <- rep(sqrt(0.1),ncol(Exposure))
    Exposure <- Exposure + lapply(beta_U_exp,function(x) x * U) %>% do.call(cbind,.)
    diag(var(Exposure))
    #################### Change 2024.1.5 #######################
    #################### Fix effect size #######################
    ############################################################
    corr_causal_eff <- 0
    beta_covariance <- diag(rep(sqrt(h2_out)),sum(beta_causal_status))
    correl_matrix <- matrix(1,sum(beta_causal_status),sum(beta_causal_status)) * corr_causal_eff + diag(rep(1,sum(beta_causal_status))) * (1-corr_causal_eff)
    beta_covariance <- beta_covariance %*% correl_matrix %*% beta_covariance
    ############################################################
    
    
    if(sum(beta_causal_status) > 0) {
      beta_causal[beta_causal_status == 1] <- mvtnorm::rmvnorm(1,mean = rep(0,sum(beta_causal_status)),sigma = beta_covariance)
      
      # every effect are from the Normal distribution with same variance
      beta_causal <- beta_causal/(sqrt(diag(var(Exposure%*%diag(beta_causal))))/sqrt(h2_out))
      beta_causal[is.na(beta_causal)] <- 0
      
      if(!is.na(fixedBeta)) {
        beta_causal[beta_causal_status == 1] <- fixedBeta
      }
    }
    
    if(tissue_sharing > 1) {
      random_index <- c(sample(1:K1),sample(1:K2+K1))
      beta_causal <- beta_causal[random_index]
      Exposure <- Exposure[,random_index]
    }
    
    true_beta[repind,] <- beta_causal
    
    error <- rnorm(N,0,1)
    
    if(is.na(fixedBeta)){
      error <- error / (sd(error) / sqrt(1-ifelse(sum(latent_status[repind,]) == 0,0,h2_out)-ifelse(prop_UHP==0,0,h2_out_IV)))
    }
    set.seed(seed+repind)
    ivIndex <- sample(1:ncol(IVs),ceiling(M*prop_UHP))
    beta_IV <- rnorm(length(ivIndex),fixedIVBeta,1)
    
    if(length(beta_IV) > 1) {
      coeff <- sd(IVs[,ivIndex] %*% beta_IV) / sqrt(h2_out_IV)
    } else {
      coeff <- sd(IVs[,ivIndex] * beta_IV) / sqrt(h2_out_IV)
    }
    
    
    if(length(ivIndex) > 1) {
      UHP <- IVs[,ivIndex] %*% beta_IV/coeff
    }
    if(length(ivIndex) == 1) {
      UHP <- IVs[,ivIndex] * beta_IV/coeff
    }
    if(length(ivIndex) == 0) {
      UHP <- 0
    }
    UHP <- UHP / sqrt(G)

    res <- list(IVs = IVs,
                adj_matrix = adj_matrix,
                ExposureEffect = Exposure %*% beta_causal,
                Exposure = Exposure,
                beta_causal = beta_causal,
                UHP = UHP,
                LD = LD,
                U = U)
    return(res)
  },cl=ncl)
  
  set.seed(seed)
  error <- rnorm(N,0,1)
  error <- error / (sd(error) / sqrt(1-h2_out*sum(latent_status)-h2_out_IV-h2_U))
  
  ###########################################################
  Out <- lapply(parres,function(x) x$ExposureEffect) %>% Reduce("+",.) +
    lapply(parres,function(x) x$UHP) %>% Reduce("+",.) +
    # U * sqrt(h2_U) +
    (lapply(parres,function(x) x$U) %>% Reduce("+",.))/sqrt(G) * sqrt(h2_U) +
    error
  
  
  true_beta <- lapply(parres,function(x) x$beta_causal) %>% do.call("rbind",.)
  
  parres <- pblapply(1:G, function(repind) {
    obj <- parres[[repind]]
    IVs <- obj$IVs
    Exposure <- obj$Exposure
    adj_matrix <- obj$adj_matrix
    LD <- obj$LD
    # step 5: perform linear reg to obtain summary statistics
    if(overlap) {
      summary_Y_IV <- t(apply(IVs,2,function(x) simplelm(y=Out,x=x))) %>% data.frame()
      sample_exp <- 1:Nexp1
      if(K1 == 0){
        summary_Exp_IV <- list()
      } else if (K1 == 1){
        summary_Exp_IV <- plyr::alply(matrix(Exposure[sample_exp,1:K1]),2,function(y)
          t(apply(IVs[sample_exp,],2,function(x) simplelm(y=y,x=x))))
      } else {
        summary_Exp_IV <- plyr::alply(Exposure[sample_exp,1:K1], 2,function(y)
          t(apply(IVs[sample_exp,],2,function(x) simplelm(y=y,x=x))))
      }
      
      
      sample_meth <- 1:Nexp2
      if(K2 == 0){
        summary_DNAm_IV <- list()
      } else if(K2 == 1){
        summary_DNAm_IV <- plyr::alply(matrix(Exposure[sample_meth,1:K2+K1]),2,function(y)
          t(apply(IVs[sample_meth,],2,function(x) simplelm(y=y,x=x))))
      } else {
        summary_DNAm_IV <- plyr::alply(Exposure[sample_meth,1:K2+K1],2,function(y)
          t(apply(IVs[sample_meth,],2,function(x) simplelm(y=y,x=x))))
      }
    } else {
      Nout = N - Nexp1 - Nexp2
      summary_Y_IV <- t(apply(IVs[1:Nout,],2,function(x) simplelm(y=Out[1:Nout,],x=x))) %>% data.frame()
      sample_exp <- Nout+1:Nexp1
      write.csv(Exposure[sample_exp, 1:K1], file = paste0("sample_exp_", repind, ".csv"), row.names = FALSE)
      if(K1 == 0){
        summary_Exp_IV <- list()
      } else if (K1 == 1){
        summary_Exp_IV <- plyr::alply(matrix(Exposure[sample_exp,1:K1]),2,function(y)
          t(apply(IVs[sample_exp,],2,function(x) simplelm(y=y,x=x))))
      } else {
        summary_Exp_IV <- plyr::alply(Exposure[sample_exp,1:K1], 2,function(y)
          t(apply(IVs[sample_exp,],2,function(x) simplelm(y=y,x=x))))
      }
      sample_meth <- Nout+Nexp1+1:Nexp2
      write.csv(Exposure[sample_meth, 1:K2 + K1], file = paste0("sample_meth_", repind, ".csv"), row.names = FALSE)
      if(K2 == 0){
        summary_DNAm_IV <- list()
      } else if(K2 == 1){
        summary_DNAm_IV <- plyr::alply(matrix(Exposure[sample_meth,1:K2+K1]),2,function(y)
          t(apply(IVs[sample_meth,],2,function(x) simplelm(y=y,x=x))))
      } else {
        summary_DNAm_IV <- plyr::alply(Exposure[sample_meth,1:K2+K1],2,function(y)
          t(apply(IVs[sample_meth,],2,function(x) simplelm(y=y,x=x))))
      }
      
    }
    
    summary_X_IV <- c(summary_Exp_IV,summary_DNAm_IV)
    beta_Y <- summary_Y_IV$Estimate
    beta_X <- lapply(summary_X_IV, function(x) data.frame(x)$Estimate) %>%
      do.call("cbind",.)
    se_Y <- summary_Y_IV$Std..Error
    se_X <- lapply(summary_X_IV, function(x) data.frame(x)$Std..Error) %>%
      do.call("cbind",.)
    
    pval_X <- 2*(1-pnorm(abs(beta_X/se_X)))
    adj_matrix <- 1*(pval_X < threshold)
    
    
    tissue_IV_num <- mean(colSums(adj_matrix))
    total_IV_num <- sum(rowSums(adj_matrix)>=1)
    # step 6: perform MVMR
    
    if(run_MR) {
      for (i in 1:(K1+K2)) {
        if (sum(adj_matrix[,i]) < 3) {
          next
        }
        print(i)
        bx <- apply(beta_X[adj_matrix[,i] == 1,] %>% data.matrix(),1,function(x) meta_iv(x,c(rep(Nexp1,K1),rep(Nexp2,K2))))
        sx <- apply(se_X[adj_matrix[,i] == 1,] %>% data.matrix(),1,function(x) meta_iv(x,c(rep(Nexp1,K1),rep(Nexp2,K2))))
        by <- beta_Y[adj_matrix[,i] == 1]
        sy <- se_Y[adj_matrix[,i] == 1]
        mr_obj <- mr_input(bx = bx,bxse = sx,by = by, byse = sy)
        # ========================= #
        # MR-IVW
        # ========================= #
        res_metaivw <- tryCatch(mr_ivw(mr_obj),silent = TRUE, error = function(x) return(NA))
        pval$metaivw[repind,i] <- tryCatch(res_metaivw$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$metaivw[repind,i] <- tryCatch(res_metaivw$Estimate,silent = TRUE, error = function(x) return(NA))
        
        bx <- beta_X[adj_matrix[,i] == 1,i]
        sx <- se_X[adj_matrix[,i] == 1,i]
        by <- beta_Y[adj_matrix[,i] == 1]
        sy <- se_Y[adj_matrix[,i] == 1]
        mr_obj <- mr_input(bx = bx,bxse = sx,by = by,byse = sy, correlation = LD[adj_matrix[,i] == 1,adj_matrix[,i] == 1])
        presso.df = data.frame(bx = bx, by = by,bxse = sx, byse = sy)
        grapple.df <- data.frame(gamma_exp = bx,gamma_out = by,se_exp = sx,se_out = sy,selection_pvals = 2*pnorm(-abs(bx/sx)))
        
        
        # common IVW
        res_ivw <- tryCatch(mr_ivw(mr_obj),silent = TRUE, error = function(x) return(NA))
        pval$ivw[repind,i] <- tryCatch(res_ivw$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$ivw[repind,i] <- tryCatch(res_ivw$Estimate,silent = TRUE, error = function(x) return(NA))
        
        
        # ========================= #
        # MR-Egger
        # ========================= #
        res_egger <- tryCatch(mr_egger(mr_obj),silent = TRUE, error = function(x) return(NA))
        pval$egger[repind,i] <- tryCatch(res_egger$Causal.pval,silent = TRUE, error = function(x) return(NA))
        est$egger[repind,i] <- tryCatch(res_egger$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # ========================= #
        # MR-Lasso
        # ========================= #
        res_lasso <- tryCatch(mr_lasso(mr_obj, lambda = 1),silent = TRUE, error = function(x) return(NA))
        pval$lasso[repind,i] <- tryCatch(res_lasso$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$lasso[repind,i] <- tryCatch(res_lasso$Estimate,silent = TRUE, error = function(x) return(NA))
      }
      
      exposure_index <- c(1:K1,1:K2+K1)
      
      bx <- beta_X[rowSums(adj_matrix[,exposure_index])>=1,exposure_index]
      sx <- se_X[rowSums(adj_matrix[,exposure_index])>=1,exposure_index]
      by <- beta_Y[rowSums(adj_matrix[,exposure_index])>=1]
      sy <- se_Y[rowSums(adj_matrix[,exposure_index])>=1]
      
      if(sum(rowSums(adj_matrix[,exposure_index])>=1) > length(exposure_index)) {
        mr_obj <- mr_mvinput(bx = bx,bxse = sx,
                             by = by,byse = sy,
                             correlation = LD[rowSums(adj_matrix[,exposure_index])>=1,
                                              rowSums(adj_matrix[,exposure_index])>=1])
        
        # mr_obj <- mr_mvinput(bx = bx,bxse = sx,by = by,byse = sy)
        # mvivw
        res_mvivw <- tryCatch(mr_mvivw(mr_obj),silent = TRUE, error = function(x) return(NA))
        pval$mvivw[repind,exposure_index] <- tryCatch(res_mvivw$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$mvivw[repind,exposure_index] <- tryCatch(res_mvivw$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # mvegger
        res_mvegger <- tryCatch(mr_mvegger(mr_obj),silent = TRUE, error = function(x) return(NA))
        pval$mvegger[repind,exposure_index] <- tryCatch(res_mvegger$Pvalue.Est,silent = TRUE, error = function(x) return(NA))
        est$mvegger[repind,exposure_index] <- tryCatch(res_mvegger$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # mvmedian
        res_mvmedian <- tryCatch(mr_mvmedian(mr_obj,iterations=100),silent = TRUE, error = function(x) return(NA))
        pval$mvmedian[repind,exposure_index] <- tryCatch(res_mvmedian$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$mvmedian[repind,exposure_index] <- tryCatch(res_mvmedian$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # mvlasso
        res_mvlasso <- tryCatch(mr_mvlasso(mr_obj, lambda = 1.8),silent = TRUE, error = function(x) return(NA))
        pval$mvlasso[repind,exposure_index] <- tryCatch(res_mvlasso$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$mvlasso[repind,exposure_index] <- tryCatch(res_mvlasso$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # mvrobust
        res_mvrobust <- tryCatch(mr_mvivw(mr_obj, robust = T, distribution = "t-dist"),silent = TRUE, error = function(x) return(NA))
        pval$mvrobust[repind,exposure_index] <- tryCatch(res_mvrobust$Pvalue,silent = TRUE, error = function(x) return(NA))
        est$mvrobust[repind,exposure_index] <- tryCatch(res_mvrobust$Estimate,silent = TRUE, error = function(x) return(NA))
        
        # # mvcML
        # ## rho_mat is pre-calculated by bivariate LDSC
        # res_mvcml <- tryCatch(mr_mvcML(mr_obj, n = min(N,Nexp1,Nexp2)),silent = TRUE, error = function(x) return(NA))
        # pval$mvcml[repind,exposure_index] <- tryCatch(res_mvcml$Pvalue,silent = TRUE, error = function(x) return(NA))
        # est$mvcml[repind,exposure_index] <- tryCatch(res_mvcml$Estimate,silent = TRUE, error = function(x) return(NA))
        
        
        # if(min(N,Nexp1,Nexp2) > 100) {
        #   rho_mat = diag(1,K1+K2+1)
        #   Sig_inv_l = invcov_mvmr(se_bx=sx,se_by=sy,rho_mat = rho_mat)
        #
        #   tryCatch({
        #     # MVcML_res = MVmr_cML_DP(b_exp=bx,
        #     #                         b_out=as.matrix(by),
        #     #                         se_bx=sx,
        #     #                         Sig_inv_l=Sig_inv_l,n = min(N,Nexp1,Nexp2),
        #     #                         num_pert = 50,
        #     #                         maxit = 10,
        #     #                         K_vec = seq(0,nrow(bx),by=4))
        #     MVcML_res = MVmr_cML_DP(b_exp=bx,
        #                             b_out=as.matrix(by),
        #                             se_bx=sx,
        #                             Sig_inv_l=Sig_inv_l,
        #                             n = min(N,Nexp1,Nexp2),
        #                             K_vec = 0:(ceiling(nrow(bx)/2)))
        #
        #     pval$mvcml[repind,] = (pnorm(-abs(MVcML_res$BIC_DP_theta/MVcML_res$BIC_DP_se))*2)[,1]
        #
        #     est$mvcml[repind,] <-
        #     pval$mvcml[repind,] = mr_mvcML(mr_obj, n = min(N,Nexp1,Nexp2))$Pvalue
        #   },silent = TRUE, error = function(x) return(NA))
        # }
        
        
        # oracle known status
        if(sum(obj$beta_causal != 0) > 0) {
          ind <- which(obj$beta_causal!=0)
          mr_obj <- mr_mvinput(bx = bx[,ind,drop=F],bxse = sx[,ind,drop=F],by = by,byse = sy)
          res_oracle <- tryCatch(mr_mvivw(mr_obj),silent = TRUE, error = function(x) return(NA))
          pval$oracle[repind,ind] <- tryCatch(res_oracle$Pvalue,silent = TRUE, error = function(x) return(NA))
          est$oracle[repind,ind] <- tryCatch(res_oracle$Estimate,silent = TRUE, error = function(x) return(NA))
        }
      }
    }
    
    beta_Xlist[[repind]] <- beta_X[rowSums(adj_matrix)>=1,]
    beta_Ylist[[repind]] <- beta_Y[rowSums(adj_matrix)>=1]
    se_Xlist[[repind]] <- se_X[rowSums(adj_matrix)>=1,]
    se_Ylist[[repind]] <- se_Y[rowSums(adj_matrix)>=1]
    LD <- LD[rowSums(adj_matrix)>=1,rowSums(adj_matrix)>=1]
    pval_list[[repind]] <- pval_X[rowSums(adj_matrix)>=1,]
    
    return(list(tissue_IV_num = tissue_IV_num,
                total_IV_num = total_IV_num,
                true_beta = true_beta[repind,],
                LD = LD,
                beta_Xlist = beta_Xlist[[repind]],
                beta_Ylist = beta_Ylist[[repind]],
                se_Xlist = se_Xlist[[repind]],
                se_Ylist = se_Ylist[[repind]],
                pval_list = pval_list[[repind]],
                pval = lapply(pval,function(x) x[repind,]),
                est = lapply(est, function(x) x[repind,])))
  },cl=ncl)
  
  
  
  tissue_IV_num <- sapply(parres,function(x) x$tissue_IV_num) %>% mean()
  total_IV_num <- sapply(parres,function(x) x$total_IV_num) %>% mean()
  # true_beta <- lapply(parres, function(x) x$true_beta) %>% do.call("rbind",.)
  LD <- lapply(parres, function(x) x$LD)
  beta_Xlist <- lapply(parres, function(x) x$beta_Xlist)
  beta_Ylist <- lapply(parres, function(x) x$beta_Ylist)
  se_Xlist <- lapply(parres, function(x) x$se_Xlist)
  se_Ylist <- lapply(parres, function(x) x$se_Ylist)
  p_list <- lapply(parres, function(x) x$pval_list)
  pval <- lapply(parres, function(x) x$pval) %>% reshape_list()
  est <- lapply(parres, function(x) x$est) %>% reshape_list()
  true_beta <- lapply(parres, function(x) x$true_beta) %>% do.call("rbind",.)
  latent_status <- 1 * (true_beta!=0)
  
  
  gamma_hat <- lapply(beta_Xlist,data.matrix)
  Gamma_hat <- lapply(beta_Ylist,matrix)
  se_g <- lapply(se_Xlist,data.matrix)
  se_G <- lapply(se_Ylist,matrix)
  
  names(pval) <- paste0("pval_",methods)
  names(est) <- paste0("est_",methods)
  
  result <- list(gamma_hat = gamma_hat,
                 Gamma_hat = Gamma_hat,
                 se_g = se_g,
                 se_G = se_G,
                 LD = LD,
                 beta_Xlist = beta_Xlist,
                 beta_Ylist = beta_Ylist,
                 se_Xlist = se_Xlist,
                 se_Ylist = se_Ylist,
                 p_list = p_list,
                 true_beta = true_beta,
                 latent_status = latent_status,
                 tissue_IV_num = tissue_IV_num,
                 total_IV_num = total_IV_num)
  result <- c(result,est,pval)
  return(result)
}
cor_mat_gen <- function(r1,r2,cut,K,diag = F) {
  cut <- mapply(function(x,k) c(x,k-x),cut,K,SIMPLIFY = F) %>% unlist()
  cor_mat <- matrix(0,nrow = sum(cut),ncol = sum(cut))
  cut <- c(0,cut)
  cut <- cumsum(cut)
  for (i in 1:(length(cut)-1)) {
    cor_mat[(cut[i]+1):cut[i+1],(cut[i]+1):cut[i+1]] <- r1[i]
  }
  for (i in 1:(length(cut)-3)) {
    if(diag) {
      diag(cor_mat[(cut[i]+1):cut[i+1],(cut[i+2]+1):cut[i+3]]) <- r2
    } else {
      cor_mat[(cut[i]+1):cut[i+1],(cut[i+2]+1):cut[i+3]] <- r2
    }
  }
  cor_mat[lower.tri(cor_mat)] <- 0
  cor_mat <- cor_mat + t(cor_mat)
  diag(cor_mat) <- 1
  return(cor_mat)
}

geno_generation <- function (maf, M, rho, n, type = "AR1") {
  SIGMA = matrix(nrow = M, ncol = M)
  for (i in 1:M) {
    for (j in 1:M) {
      if(type == "AR1") {
        SIGMA[i, j] = rho^(abs(i - j))
      }
      if(type == "exchangeable") {
        SIGMA[i, j] = rho^(ifelse(i - j == 0, 0, 1))
      }
    }
  }
  
  nsnp = M
  X = NULL
  
  AAprob = maf^2
  Aaprob = 2 * maf * (1 - maf)
  quanti = matrix(c(1 - Aaprob - AAprob, 1 - AAprob),
                  M, 2)
  Xt = mvtnorm::rmvnorm(n, mean = rep(0, M), sigma = SIGMA, method = "chol")
  Xt2 = matrix(0, n, M)
  for (j in 1:M) {
    cutoff = qnorm(quanti[j, ])
    Xt2[Xt[, j] < cutoff[1], j] = 0
    Xt2[Xt[, j] >= cutoff[1] & Xt[, j] < cutoff[2],
        j] = 1
    Xt2[Xt[, j] >= cutoff[2], j] = 2
  }
  X <- cbind(X, Xt2)
  
  return(X)
}

generate_exposure_mult_consistency <- function(N, Nexp1, Nexp2, IVs, h2_exp, h2_meth, M, M_gene, K1, K2, adj_matrix, diag_corr = 0, off_diag_corr = 0){
  beta_covariance <- diag(c(rep(sqrt(h2_exp/M),K1),rep(sqrt(h2_meth/M),K2)),nrow = K1+K2)
  correl_matrix <- (matrix(1,K1+K2,K1+K2) * diag_corr +
                      (1 - magic::adiag(matrix(1,K1,K1),matrix(1,K2,K2))) * (off_diag_corr - diag_corr) +
                      diag(rep(1,K1+K2)) * (1-diag_corr))
  beta_covariance <- beta_covariance %*% correl_matrix %*% beta_covariance
  eQTLbeta <- mvtnorm::rmvnorm(M,mean = rep(0,K1+K2),sigma = beta_covariance)
  Exposure <- lapply(1:(K1+K2), function(i) {
    beta <- eQTLbeta[,i]
    Exps = IVs %*% beta
    Exps[(N - Nexp1 - Nexp2) + 1:(Nexp1+Nexp2)] =
      IVs[(N - Nexp1 - Nexp2) + 1:(Nexp1+Nexp2),adj_matrix[,i]==1] %*%
      beta[adj_matrix[,i]==1]
    return(Exps)
  }) %>% do.call("cbind",.)
  return(Exposure)
}

add_correlation_structure <- function(Exps, DNAm, K1, K2,h2_exp, h2_meth, h2_U_exp, h2_asso,diag_corr = 0.8, off_diag_corr = 0.8){
  error_covariance <- diag(c(rep(sqrt(1 - h2_exp - h2_U_exp - h2_asso),K1),
                             rep(sqrt(1 - h2_meth - h2_U_exp - h2_asso),K2)),nrow = K1+K2)
  
  if(class(Exps)[1] == "numeric"){
    N <- length(Exps)
  } else {
    N <- nrow(Exps)
  }
  
  correl_matrix <- (matrix(1,K1+K2,K1+K2) * diag_corr +
                      (1 - magic::adiag(matrix(1,K1,K1),matrix(1,K2,K2))) * (off_diag_corr - diag_corr) +
                      diag(rep(1,K1+K2)) * (1-diag_corr))
  
  error_covariance <- error_covariance %*% correl_matrix %*% error_covariance
  
  if(K1+K2 > 1) {
    error <- mvtnorm::rmvnorm(N,mean = rep(0,K1+K2),sigma = error_covariance)
  } else {
    error <- rnorm(N, mean = 0, sd = sqrt(error_covariance)) %>% matrix()
  }
  Exps = Exps + error[,1:K1]
  DNAm = DNAm + error[,setdiff(1:(K2+K1),1:K1)]
  
  return(list(Exps = Exps, DNAm = DNAm))
}

add_association <- function(Exps, DNAm, K1, K2, h2_asso){
  beta_DNAm = matrix(rnorm(K1 * K2, mean = 0, sd = 1),nrow = K2)
  beta_Exps = matrix(rnorm(K1 * K2, mean = 0, sd = 1),nrow = K1)
  
  scale1 <- sqrt(diag(var(DNAm %*% beta_DNAm))) / sqrt(h2_asso/K1)
  scale2 <- sqrt(diag(var(Exps %*% beta_Exps))) / sqrt(h2_asso/K2)
  
  Exps_tmp <- Exps
  DNAm_tmp <- DNAm
  if(K2 > 0) {
    Exps = Exps_tmp + sweep(DNAm_tmp %*% beta_DNAm,2,scale1,"/")
  }
  if(K1 > 0) {
    DNAm = DNAm_tmp + sweep(Exps_tmp %*% beta_Exps,2,scale2,"/")
  }
  
  return(list(Exps = Exps, DNAm = DNAm))
}

meta_iv <- function(eff, w){
  sum(sqrt(w) * eff)/sqrt(sum(w))
}

reshape_list <- function(list_of_list) {
  att_names <- lapply(list_of_list,names) %>% unlist() %>% unique()
  reshaped <- lapply(att_names, function(x) lapply(list_of_list, function(y) y[[x]])) %>%
    lapply(function(x) do.call("rbind", x))
  names(reshaped) <- att_names
  return(reshaped)
}

simplelm <- function (x, y) {
  ## number of data
  n <- length(x)
  ## centring
  y0 <- sum(y) / length(y); yc <- y - y0
  x0 <- sum(x) / length(x); xc <- x - x0
  ## fitting an intercept-free model: yc ~ xc + 0
  xty <- c(crossprod(xc, yc))
  xtx <- c(crossprod(xc))
  slope <- xty / xtx
  rc <- yc - xc * slope
  ## Pearson estimate of residual standard error
  sigma2 <- c(crossprod(rc)) / (n - 2)
  ## standard error for slope
  slope_se <- sqrt(sigma2 / xtx)
  ## t-score and p-value for slope
  tscore <- slope / slope_se
  pvalue <- 2 * pt(abs(tscore), n - 2, lower.tail = FALSE)
  ## return estimation summary for slope
  c("Estimate" = slope, "Std. Error" = slope_se, "t value" = tscore, "Pr(>|t|)" = pvalue)
}



conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(IRanges::cor)
conflicted::conflicts_prefer(stats::var)
conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(stats::cov)
conflicted::conflicts_prefer(stats::cor)
conflicted::conflicts_prefer(mvtnorm::rmvnorm)
conflicted::conflicts_prefer(MVMRcML::MVmr_cML_DP)
conflicted::conflicts_prefer(base::rbind)
conflicted::conflicts_prefer(base::cbind)
conflicted::conflicts_prefer(dplyr::select)
setwd("D:/R/MINTMRCS/MVMR_ML")

option_list <- list(
  make_option(c("-N", "--N"), type="integer", default=50000, help="N", metavar="integer"),
  make_option("--Nexp1", type="integer", default=500, help="Nexp1", metavar="integer"),
  make_option("--Nexp2", type="integer", default=500, help="Nexp2", metavar="integer"),
  make_option(c("-k", "--K1"), type="integer", default=5, help="K1", metavar="integer"),
  make_option(c("-K", "--K2"), type="integer", default=5, help="K2", metavar="integer"),
  make_option("--pattern1", type="double", default=0, help="pattern1", metavar="double"),
  make_option("--pattern2", type="double", default=0, help="pattern2", metavar="double"),
  make_option(c("-G", "--G"), type="integer", default=50, help="G", metavar="integer"),
  make_option("--M", type="integer", default=15, help="M", metavar="integer"),
  make_option("--M_gene", type="integer", default=12, help="M_gene", metavar="integer"),
  make_option("--prop_UHP", type="double", default=1, help="prop_UHP", metavar="double"),
  make_option("--h2_CHP", type="double", default=0, help="h2_CHP", metavar="double"),
  make_option("--h2_exp", type="double", default=0.3, help="h2_exp", metavar="double"),
  make_option("--h2_meth", type="double", default=0.3, help="h2_meth", metavar="double"),
  make_option("--h2_out_IV", type="double", default=0.15, help="h2_out_IV", metavar="double"),
  make_option("--h2_asso", type="double", default=0, help="h2_asso", metavar="double"),
  make_option("--h2_out", type="double", default=0.015, help="h2_out", metavar="double"),
  make_option("--h2_U", type="double", default=0.1, help="h2_out", metavar="double"),
  make_option("--h2_U_exp", type="double", default=0.1, help="h2_out", metavar="double"),
  make_option("--corr_gene", type="double", default=0.1, help="corr_gene", metavar="double"),
  make_option("--corr_gene_cpg", type="double", default=0.1, help="corr_gene_cpg", metavar="double"),
  make_option("--tissue_sharing", type="integer", default=0.2, help="tissue_sharing", metavar="integer"),
  make_option("--seed", type="integer", default=21, help="seed", metavar="integer"),
  make_option("--rho", type="double", default=0, help="rho", metavar="double"),
  make_option("--threshold", type="double", default=1e-2, help="threshold", metavar="double"),
  make_option("--cooccur", type="double", default=0.5, help="cooccur", metavar="double"),
  make_option("--density", type="double", default=0.85, help="density", metavar="double"),
  make_option("--k_group1", type="integer", default=1, help="k_group1", metavar="integer"),
  make_option("--k_group2", type="integer", default=1, help="k_group2", metavar="integer"),
  make_option("--a_init", type="double", default=0.1, help="a_init", metavar="double"),
  make_option("--b_init", type="double", default=0.1, help="b_init", metavar="double"),
  make_option("--b_beta", type="double", default=0.03, help="b_beta", metavar="double")
)
args <- parse_args(OptionParser(option_list=option_list)) %>% unlist()

# option_list <- list(make_option(c("-n", "--number"), type="character", default="1"),
#                     make_option(c("-f", "--param_file"), type="character", default="params.txt"))
# args <- as.numeric(parse_args(OptionParser(option_list=option_list))$number)
# fn <- as.character(parse_args(OptionParser(option_list=option_list))$param_file)
# 
# params <- fread(paste0("script/",fn))
# args <- params[args,] %>% unlist()


unique_cache <- paste0("D:/R/MINTMRCS/MVMR_ML/script/", paste0(args,collapse = "_"))

if (!dir.exists(unique_cache)) {
  dir.create(unique_cache, recursive = TRUE)
}

dat <- par_data_generation_multiple_data(args,ncl=1,run_MR = T,overlap = F)
getwd()
save(dat, file = "dat3.RData")









dat[str_detect(names(dat),"pval")] <- lapply(dat[str_detect(names(dat),"pval")],function(x) x[,1:args["K1"]])
latent_status <- dat$latent_status[,1:args["K1"]]
#----------------------------------------------------#
# run mintMR
#----------------------------------------------------#
token = paste(args,collapse = "_")
set.seed(args["seed"])
print(args)
print(dat$total_IV_num)
print(dat$tissue_IV_num)
K1 = args["K1"]; K2 <- args["K2"]
K <- args["K2"]+args["K1"]
O <- 2; L = length(dat$gamma_hat); group <- list(1:args["K1"],1:args["K2"]+args["K1"])

############ derive parameter from data ################
b_beta_init <- mean(var(dat$est_mvivw[dat$latent_status==0])) * (K1+K2) / 2
b_beta_init
########################################################
opts = list(a_gamma = rep(0,L), b_gamma = rep(0,L),
            a_alpha = rep(0,L), b_alpha = rep(0,L),
            a_beta = rep(0,L), b_beta = rep(b_beta_init,L),
            a = args["a_init"], b = args["b_init"],
            maxIter = 4000, thin = 10, burnin = 1000)
set.seed(args["seed"])

CC <- floor((K1+K2)/4)
PC1 <- K1 - 2 - CC
PC2 <- K2 - 2 - CC
res_CCA_LD <- LDGibbsSampMultOmics(gammah = dat$gamma_hat,
                                   Gammah = dat$Gamma_hat,
                                   se1 = dat$se_g,
                                   se2 = dat$se_G,
                                   latent_status = dat$latent_status,
                                   corr_mat = dat$LD,
                                   group=group,
                                   opts = opts,
                                   low_rank = T,
                                   display_progress = T,
                                   CC=CC,PC1=PC1,PC2=PC2)
gibbs_result <- gibbs_output_process(res_CCA_LD)
gibbs_result_lsnolr <- gibbs_result

# single-gene, no latent status
res <- list()
for (j in 1:length(dat$gamma_hat)) {
  res_CCA_LD <- LDGibbsSingleGene(gammah = dat$gamma_hat[j],
                                  Gammah = dat$Gamma_hat[j],
                                  se1 = dat$se_g[j],
                                  se2 = dat$se_G[j],
                                  corr_mat = dat$LD[j],
                                  opts = opts,
                                  # low_rank = F,
                                  sgalprecompute = dat$tissue_IV_num,
                                  display_progress = F)
  nolr_gibbs_result <- gibbs_output_process(res_CCA_LD)
  res[[j]] <- nolr_gibbs_result$pval_mint
}
nolr_gibbs_result$pval_mint <- do.call("rbind",res)


res_CCA_LD <- LDGibbsSampMultOmics(gammah = dat$gamma_hat,
                                   Gammah = dat$Gamma_hat,
                                   se1 = dat$se_g,
                                   se2 = dat$se_G,
                                   latent_status = dat$latent_status,
                                   corr_mat = dat$LD,
                                   group=group,
                                   opts = opts,
                                   low_rank = F,
                                   display_progress = T,
                                   CC=2,PC1=2,PC2=2,oracle=T)
oracle_gibbs_result <- gibbs_output_process(res_CCA_LD)

## compare power
pow_comparison <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) sum(x[,1:K1][dat$latent_status[,1:K1] == 1] < 0.05,na.rm = T) /
           sum(dat$latent_status[,1:K1] == 1 & !is.na(x)))
pow_mint <- sum(gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 1] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 1)
pow_mint_nolr <- sum(nolr_gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 1] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 1)
pow_mint_oracle <- sum(oracle_gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 1] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 1)
pow_mint_lsnolr <- sum(gibbs_result_lsnolr$pval_mint[,1:K1][dat$latent_status[,1:K1] == 1] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 1)




pow_comparison$pval_mint <- pow_mint
pow_comparison$pval_mint_nolr <- pow_mint_nolr
pow_comparison$pval_mint_oracle <- pow_mint_oracle
pow_comparison$pval_mint_ldnolr <- pow_mint_lsnolr
pow_comparison <- data.frame(pow_comparison)
pow_data <- pow_comparison
## compare t1r
t1r_comparison <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) sum(x[,1:K1][dat$latent_status[,1:K1] == 0] < 0.05,na.rm = T) /
           sum(dat$latent_status[,1:K1] == 0 & !is.na(x)))
t1r_mint <- sum(gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 0] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 0)
t1r_mint_nolr <- sum(nolr_gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 0] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 0)
t1r_mint_oracle <- sum(oracle_gibbs_result$pval_mint[,1:K1][dat$latent_status[,1:K1] == 0] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 0)
t1r_mint_lsnolr <- sum(gibbs_result_lsnolr$pval_mint[,1:K1][dat$latent_status[,1:K1] == 0] < 0.05,na.rm = T)/
  sum(dat$latent_status[,1:K1] == 0)


t1r_comparison$pval_mint <- t1r_mint
t1r_comparison$pval_mint_nolr <- t1r_mint_nolr
t1r_comparison$pval_mint_oracle <- t1r_mint_oracle
t1r_comparison$pval_mint_lsnolr <- t1r_mint_lsnolr
t1r_comparison <- data.frame(t1r_comparison)
t1r_data <- t1r_comparison

## compare rmse
rmse_comparison <- dat[str_detect(names(dat),"est_")] %>%
  lapply(function(x) sqrt(mean((dat$true_beta[,1:K1] - x[,1:K1])^2,na.rm=T)))

rmse_mint <-sqrt(mean((dat$true_beta[,1:K1] - gibbs_result$beta_mint[,1:K1])^2,na.rm = T))
rmse_mint_nolr <- sqrt(mean((dat$true_beta[,1:K1] - nolr_gibbs_result$beta_mint[,1:K1])^2,na.rm = T))
rmse_mint_oracle <- sqrt(mean((dat$true_beta[,1:K1] - oracle_gibbs_result$beta_mint[,1:K1])^2,na.rm = T))
rmse_mint_lsnolr <- sqrt(mean((dat$true_beta[,1:K1] - gibbs_result_lsnolr$beta_mint[,1:K1])^2,na.rm = T))


rmse_comparison$est_mint <- rmse_mint
rmse_comparison$est_mint_nolr <- rmse_mint_nolr
rmse_comparison$est_mint_oracle <- rmse_mint_oracle
rmse_comparison$est_mint_lsnolr <- rmse_mint_lsnolr
rmse_comparison <- data.frame(rmse_comparison)
rmse_data <- rmse_comparison

args["tissue_IV_num"] <- dat$tissue_IV_num
args["total_IV_num"] <- dat$total_IV_num
args["nonzero_effects"] <- sum(dat$latent_status[,1:K1])
result <- list(args = args, pow = pow_data, t1r = t1r_data, rmse = rmse_data)
save(result, file = paste0("simulation/",paste0(args,collapse = "_"),".RData"))
system(paste0("rm -r ",unique_cache))




#########################################################################################################
load(file = "C:/Users/wang'shuai'yi/Desktop/mintMR数据生成/一/dat1.RData")

setwd("D:/R/MINTMRCS/MVMR_ML")
getwd()
library(readxl)
file_names <- sprintf("sample_exp_%d.csv", 1:50)
exposure_list <- lapply(file_names, function(file) {
  read.csv(file, header = TRUE, sep = ",")
})
exposure_list <- lapply(exposure_list, function(df) {
  as.matrix(sapply(df, as.numeric))
})
file_names2 <- sprintf("sample_meth_%d.csv", 1:50)
exposure_list2 <- lapply(file_names2, function(file) {
  read.csv(file, header = TRUE, sep = ",")
})
exposure_list2 <- lapply(exposure_list2, function(df) {
  as.matrix(sapply(df, as.numeric))
})
combined_exposure_list <- Map(function(mat1, mat2) {
  cbind(mat1, mat2)
}, exposure_list, exposure_list2)

W_list <- list()
for (i in 1:50) {
  W_list[[i]] <- as.matrix(dist(t(combined_exposure_list[[i]]), method = "euclidean"))   
  W_list[[i]] <- 1 / (1 + W_list[[i]])
  diag(W_list[[i]]) <- 0
  W_list[[i]] <- W_list[[i]] / rowSums(W_list[[i]])
  W_list[[i]] <- (W_list[[i]]+ t(W_list[[i]])) / 2
}


##################################################################################
library(mintMR)
library(glmnet)
load(file = "C:/Users/wang'shuai'yi/Desktop/mintMR数据生成/一/dat1.RData")
load(file = "D:/R/MINTMRCS/MVMR_ML/dat3_2.RData")
suppressMessages(library("tidyverse"))
library("ggplot2")
library("MendelianRandomization")
L <- length(dat$gamma_hat)
K1 <-  5
K2 <-  5
L  
group <- list(exposure1 = 1:K1,exposure2 = 1:K2+K1)
# column 1 to 5 in gamma_hat and se_g are from exposure 1
# column 6 to 10 in gamma_hat and se_g are from exposure 2
group
opts <- get_opts(L) 
names(opts)
# opts <- get_opts(L, maxIter = 1000)  
set.seed(1)
res_no_LD <- mintMR(gammah = dat$gamma_hat, 
                    Gammah = dat$Gamma_hat, 
                    se1 = dat$se_g, 
                    se2 = dat$se_G, 
                    group = group)




K <- K1 + K2
run_graph_lasso_with_l1 <- function(X, Y, W) {
  n <- nrow(X)
  k <- ncol(X)
  initial_beta <- rep(0, k)  
  
  graph_lasso_penalty <- function(beta, W) {
    penalty <- 0
    for (i in 1:length(beta)) {
      for (j in 1:length(beta)) {
        if (i!= j) {
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
  
  lambda1_vals <- c(0.003,0.04,0.2,0.0003)
  lambda2_vals <- c(0.008,0.01)
  
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
        train_indices <- which(folds!= fold)
        X_train <- X[train_indices,]
        Y_train <- Y[train_indices]
        X_test <- X[test_indices,]
        Y_test <- Y[test_indices]
        
        result_graph_l1 <- optim(par = initial_beta,
                                 fn = total_loss,
                                 X = X_train, Y = Y_train, lambda1 = lambda1, lambda2 = lambda2, W = W,
                                 method = "L-BFGS-B")
        
        beta_graph_l1 <- result_graph_l1$par
        # fit_lasso <- glmnet(X_train, Y_train, alpha = 1, lambda = lambda2, start = beta_graph_l1)
        # beta_final <- as.numeric(coef(fit_lasso, s = "lambda.min"))[-1]
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
for (i in 1:50) {
  X <- dat[["gamma_hat"]][[i]]
  Y <- dat[["Gamma_hat"]][[i]]
  p <- ncol(X)
  # W <- matrix(1, nrow = K, ncol = K)
  W <- W_list[[i]]
  result <- run_graph_lasso_with_l1(X, Y, W)
  results_list[[i]] <- result
}
estimated_betas <- lapply(results_list, function(x) x$estimated_beta)
estimated_betas <- do.call(rbind, estimated_betas)
formatted_matrix <- matrix(sapply(estimated_betas, function(x) if (x == 0) formatC(x, format = "f", digits = 0) else x), nrow = nrow(estimated_betas))
formatted_matrix <- matrix(sapply(estimated_betas, function(x) if (x == 0) formatC(x, format = "f", digits = 0) else formatC(x, format = "f", digits = 5)), nrow = nrow(estimated_betas))
num_rows <- nrow(formatted_matrix)
num_cols <- ncol(formatted_matrix)
row_names <- paste0("gene", 1:num_rows)
col_names <- paste0("tissue", 1:num_cols)
rownames(formatted_matrix) <- row_names
colnames(formatted_matrix) <- col_names






Z <- dat[["latent_status"]]
Z11 <- res_no_LD[["Pvalue"]]
Z11 <- ifelse(Z11 < 0.05, 1, 0)
formatted_matrix1 <- ifelse(formatted_matrix != 0, 1, 0)  
count_a_1 <- sum(formatted_matrix1 == 1 & Z == 1)
count_a_2 <- sum(Z11 == 1 & Z == 1)
count_a <- sum( Z == 1)
count_a_1 / count_a
count_a_2 / count_a     ##power
count_b_1 <- sum(formatted_matrix1 == 1 & Z == 0)
count_b_2 <- sum(Z11 == 1 & Z == 0)
count_b <- sum( Z == 0)
count_b_1 / count_b
count_b_2 / count_b      ##T1r
sum(formatted_matrix1 == 1 & Z == 1)/sum(formatted_matrix1 == 1)
sum(Z11 == 1 & Z == 1)/sum(Z11 == 1)            ##precision
(sum(formatted_matrix1 == 1 & Z == 1) + sum(formatted_matrix1 == 0 & Z == 0))/(sum(formatted_matrix1 == 1)+sum(formatted_matrix1 == 0))
(sum(Z11 == 1 & Z == 1) + sum(Z11 == 0 & Z == 0))/(sum(Z11 == 1)+sum(Z11 == 0))



# RMSE
# ZZ <- eta * beta_gk
ZZ <-dat[["true_beta"]]
if (!is.numeric(ZZ)) {
  ZZ <- as.numeric(ZZ)
}
formatted_matrix <- apply(res_no_LD[["Estimate"]], c(1,2), as.numeric)
error <- ZZ-formatted_matrix
mean(as.vector(error^2))
sqrt(mean(as.vector(error^2)))




pow_comparison <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) sum(x[,1:K][dat$latent_status[,1:K] == 1] < 0.05,na.rm = T) /
           sum(dat$latent_status[,1:K] == 1 & !is.na(x)))
t1r_comparison <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) sum(x[,1:K][dat$latent_status[,1:K] == 0] < 0.05,na.rm = T) /
           sum(dat$latent_status[,1:K] == 0 & !is.na(x)))
Precision <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) sum(x[,1:K][dat$latent_status[,1:K] == 1] < 0.05,na.rm = T) /
           sum(x < 0.05, na.rm = T))
Accurary <- dat[str_detect(names(dat),"pval")] %>%
  lapply(function(x) {
    num_positive_and_pval_less_than_005 <- sum(x[,1:K][dat$latent_status[,1:K] == 1] < 0.05, na.rm = T)
    num_negative_and_pval_greater_than_005 <- sum(x[,1:K][dat$latent_status[,1:K] == 0] > 0.05, na.rm = T)
    return((num_positive_and_pval_less_than_005 + num_negative_and_pval_greater_than_005) / sum(x < 5, na.rm = T))
  })
rmse_comparison <- dat[str_detect(names(dat),"est_")] %>%
  lapply(function(x) sqrt(mean((dat$true_beta[,1:K] - x[,1:K])^2,na.rm=T)))

###############################################################################

library(tidyverse)

data_power <- data.frame(
  power = c(0.87489393, 0.80457915, 0.60437274, 0.66247475, 0.72412764, 0.5898326, 0.62496361,
            0.82689799, 0.65495663, 0.48164995, 0.62572184, 0.6220358, 0.54237024, 0.5320039,
            0.79885017,0.55704571, 0.36343489,0.60845441, 0.54209588,0.4820748, 0.49459711),
  Group = factor(rep(c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW'), times = 3),
                 levels = c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW')), 
  Condition = factor(rep(c("UHP effects = 0.05", "UHP effects = 0.10", "UHP effects = 0.15"), each = 7)) 
)
data_T1r <- data.frame(
  T1r = c(0.012715417, 0.039073421,  0.060990802 , 0.14523068 , 0.091106975, 0.11809295, 0.119694022,
          0.022229626,0.025506715, 0.074086389, 0.20377299,0.11993957, 0.129952112, 0.128691122,
          0.034783431,0.022444925, 0.07660257, 0.26211672, 0.13245564, 0.1322775, 0.13500947),
  Group = factor(rep(c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW'), times = 3),
                 levels = c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW')), 
  Condition = factor(rep(c("UHP effects = 0.05", "UHP effects = 0.10", "UHP effects = 0.15"), each = 7)) 
)
data_rmse <- data.frame(
  rmse = c(0.014657113,0.025100952, 0.043818478, 0.048932418, 0.052204474, 0.060287823, 0.043818478,
           0.018202101, 0.028121551, 0.055955903, 0.064393324, 0.06587278,0.077753014, 0.056055903,
           0.021572689 , 0.029989789, 0.065760035, 0.076102458, 0.077074419, 0.092027223, 0.065760035),
  Group = factor(rep(c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW'), times = 3),
                 levels = c('MT-TWMR', 'MintMR', 'MVMR-Robust', 'MVMR-Lasso', 'MVMR-Median', 'MVMR-Egger', 'MVMR-IVW')), 
  Condition = factor(rep(c("UHP effects = 0.05", "UHP effects = 0.10", "UHP effects = 0.15"), each = 7)) 
)

ggplot(data_power, aes(x = Group, y = power, fill = Group)) + 
  geom_col() + 
  facet_grid(. ~ Condition) + 
  scale_fill_manual(values = c("#df9dc0", "#C69287", "#E79A90", "#EFBC91", "#FAE5B8", "#E4CD87", "#d6e1c2")) +
  labs(y = "Power", x = NULL) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey80", color = NA), 
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"), 
    panel.spacing = unit(2.5, "lines"),
    axis.text.x = element_text(angle = 30, hjust = 0.75)
  )

ggplot(data_T1r, aes(x = Group, y = T1r, fill = Group)) + 
  geom_col() + 
  facet_grid(. ~ Condition) + 
  scale_fill_manual(values = c("#df9dc0", "#C69287", "#E79A90", "#EFBC91", "#FAE5B8", "#E4CD87", "#d6e1c2")) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "yellow")+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red")+
  labs(y = "type-I error", x = NULL) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey80", color = NA), 
    strip.text = element_text(face = "bold"), 
    axis.title.y = element_text(size = 10, face = "bold"), 
    panel.spacing = unit(2.5, "lines"),
    axis.text.x = element_text(angle = 30, hjust = 0.75) 
  )

ggplot(data_rmse, aes(x = Group, y = rmse, fill = Group)) + 
  geom_col() + 
  facet_grid(. ~ Condition) + 
  scale_fill_manual(values = c("#df9dc0", "#C69287", "#E79A90", "#EFBC91", "#FAE5B8", "#E4CD87", "#d6e1c2")) +
  labs(y = "rmse", x = NULL) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey80", color = NA), 
    strip.text = element_text(face = "bold"), 
    axis.title.y = element_text(size = 10, face = "bold") , 
    panel.spacing = unit(2.5, "lines"),
    axis.text.x = element_text(angle = 30, hjust = 0.75) 
  )
