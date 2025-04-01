data <- read.csv("D:/GTEx/处理后/1Adipose_Subcutaneous.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
devtools::install_github("mrcieu/ieugwasr",force = TRUE)
library(ieugwasr)
devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
library(data.table)
library(R.utils)
#处理一下LD，以0.01开始
colnames(data)[2] <-"Gene"
colnames(data)[13] <-"rsid"
colnames(data)[7] <- "pval"
get_plink_exe()
genes <- unique(data$Gene)
i <-1
out <-data.frame()
for (i in 1:length(genes)){
  dat <- subset(data, Gene == genes[i])
  tryCatch(
    {clumtest6 <- ieugwasr::ld_clump_local(dat = dat,clump_r2=0.01,clump_kb=10000,clump_p=1,
                                           bfile = "D:/GTEx_Analysis_v8_eQTL/g1000_eur/g1000_eur", ## 欧洲的EUR
                                           plink_bin = "D:/R/R-4.3.1/library/plinkbinr/bin/plink_Windows.exe"
    )},
    warning=function(w){print("warning")},
    finally = {out <- rbind(out,clumtest6)}
  )
  if (i %% 10 == 0){print(i)
    print("############################又是10#############################################")
  }
  
}