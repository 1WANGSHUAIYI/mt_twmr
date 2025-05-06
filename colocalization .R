Tissue11 <- read.csv("D:/GTEx/LD_clump/Adipose_Subcutaneous.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue22 <- read.csv("D:/GTEx/LD_clump/Brain_Frontal_Cortex_BA9.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue33 <- read.csv("D:/GTEx/LD_clump/Esophagus_Mucosa.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue44 <- read.csv("D:/GTEx/LD_clump/Heart_Left_Ventricle.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue55 <- read.csv("D:/GTEx/LD_clump/Kidney_Cortex.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue66 <- read.csv("D:/GTEx/LD_clump/Liver.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue77 <- read.csv("D:/GTEx/LD_clump/Lung.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue88 <- read.csv("D:/GTEx/LD_clump/Muscle_Skeletal.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue99 <- read.csv("D:/GTEx/LD_clump/Skin_Sun_Exposed_Lower_leg.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue1010 <- read.csv("D:/GTEx/LD_clump/Thyroid.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue1111 <- read.csv("D:/GTEx/LD_clump/Whole_Blood.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
# data <- read.csv("D:/GTEx/handoa.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue1 <- read.csv("D:/GTEx/GTEx_v10/1Adipose_Subcutaneous.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue2 <- read.csv("D:/GTEx/GTEx_v10/1Brain_Frontal_Cortex_BA9.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue3 <- read.csv("D:/GTEx/GTEx_v10/1Esophagus_Mucosa.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue4 <- read.csv("D:/GTEx/GTEx_v10/1Heart_Left_Ventricle.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue5 <- read.csv("D:/GTEx/GTEx_v10/1Kidney_Cortex.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue6 <- read.csv("D:/GTEx/GTEx_v10/1Liver.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue7 <- read.csv("D:/GTEx/GTEx_v10/1Lung.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue8 <- read.csv("D:/GTEx/GTEx_v10/1Muscle_Skeletal.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue9 <- read.csv("D:/GTEx/GTEx_v10/1Skin_Sun_Exposed_Lower_leg.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue10 <- read.csv("D:/GTEx/GTEx_v10/1Thyroid.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
Tissue11 <- read.csv("D:/GTEx/GTEx_v10/1Whole_Blood.v10.eQTLs.signif_pairs.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)








library(data.table)
library(coloc)
library(TwoSampleMR)
library(locuscomparer)
library(ggplot2)
library(dplyr)


##
names(data)[names(data) == 'pval.exposure'] <- 'p_GWAS'



pqtl <- Tissue1
names(pqtl)[names(pqtl) == 'gene_id'] <- 'Gene'
# pqtl <- subset(pqtl ,Chr==Probe_Chr)
# ##UKB
# gene <- c("COMT","AGT","FN1","SERPINC1")
gene <- c("ENSG00000129158.11")
# 
pqtl <- pqtl[pqtl$Gene %in% gene,]


names(pqtl)[names(pqtl) == 'rs_id_dbSNP151_GRCh38p7'] <- 'rsid'
names(pqtl)[names(pqtl) == 'pval_nominal'] <- 'pval'

# IV <- data.frame()
# for (i in gene) {
#   gene1 <-pqtl[pqtl$Gene == i,]
#   clumtest6 <- ieugwasr::ld_clump_local(dat = gene1,clump_r2=0.1,clump_kb=1000,clump_p=1,
#                                         bfile="/home/user6/tyx/duibi/Data/1kg.v3/g1000_eur", ## 欧洲的EUR
#                                         plink_bin="/home/user6/R/x86_64-conda-linux-gnu-library/4.1/plinkbinr/bin/plink_Linux")
#   
#   IV <- rbind(IV,clumtest6)                                      
#   
# }




ngwas=463010
  case=54358
  
  #decode qtl
  neqtl=714
  
  
  
  
names(pqtl)[names(pqtl) == 'slope'] <- 'b_eqtl'
names(pqtl)[names(pqtl) == 'pval_nominal'] <- 'p_eqtl'
names(pqtl)[names(pqtl) == 'rsid'] <- 'SNP'

pqtl$MAF <- ifelse(pqtl$af>0.5,1-pqtl$af,pqtl$af)


IV_litt<- subset(Tissue11,Gene %in% "ENSG00000129158.11" )
pqtl_litt<- subset(pqtl,Gene %in% "ENSG00000129158.11" )


res <- data.frame()
for (i in 1:length(IV_litt$tss_distance)) {
  print(i)
  bp <- pqtl_litt$tss_distance[i]
  bp_up <- bp + 250000
  bp_low <- bp - 250000                                                                                                                                                          
  eqtl_merge <- subset(pqtl_litt, tss_distance > bp_low &  tss_distance< bp_up)
  gwas_merge <-  data[data$SNP %in% eqtl_merge$SNP,]
  eqtl_merge <- pqtl_litt[pqtl_litt$SNP %in% gwas_merge$SNP,]
  gwas_merge <- gwas_merge[!duplicated(gwas_merge$SNP),]
  eqtl_merge <- eqtl_merge[!duplicated(eqtl_merge$SNP),]
  input = merge(gwas_merge,eqtl_merge,by="SNP")
  input$varbeta <- input$se.exposure^2
  if (length(input$SNP) > 0){
    re= coloc.abf(dataset1 = list(pvalues=input$p_GWAS,snp=input$SNP,type="cc",N=ngwas,s=case/ngwas),
                  dataset2 = list(pvalues=input$pval,snp=input$SNP,beta=input$b_eqtl,varbeta=input$varbeta, type="quant",N=neqtl)
                  ,MAF = input$MAF
                  #,p12= 1e-4
    )}
  if (re$summary[6] > 0.5){
    re <- re$summary
    re <- c(re, i)
    res <- rbind(res,re)
    
  }
} 




