setwd("D:/brain structure")
library(readxl)
library(writexl)
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(stringr)
fs <- list.files("D:/brain structure/brain structure/",pattern = "ENIGMA3_mixed",full.names = TRUE) 
ids <- unlist(fs)
ids <- str_remove_all(ids,pattern="D:/brain structure/brain structure/")
ids <- str_remove_all(ids,pattern=".txt.gz")
exps <- list.files(paste0("exposure_new/"),pattern="csv",full.names = TRUE)
exps <- unlist(exps)
exps <- str_remove_all(exps,pattern="exposure_new/")
exps <- str_remove_all(exps,pattern=".csv")


for (exp in exps){if(!is.null(exp)){for (id in ids){
  exposure_dat <-fread(file=paste0("exposure_new/",exp,".csv"),header = T)
  out <- fread(file=paste0("D:/brain structure/brain structure/",id,".txt.gz"),header = T)
  outcome_dat<-merge(exposure_dat,out,by.x = "SNP",by.y = "SNP")#修改本地结局中对应SNP列的表头名字
  write.csv(outcome_dat,file = "q.csv")
  outcome_dat<-read_outcome_data(snps = exposure_dat$SNP,filename = "q.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "BETA1",
                                 se_col = "SE",
                                 eaf_col = "FREQ1",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 pval_col = "P",
                                 samplesize_col = "N")#修改结局中对应表头名
  outcome_dat$id.outcome <- id
  dat <- harmonise_data(exposure_dat,outcome_dat)
  res <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_ivw",
                                                              "mr_egger_regression",
                                                              "mr_weighted_median",
                                                              "mr_weighted_mode",
                                                              "mr_simple_mode",
                                                              "mr_wald_ratio")))
  
  write.csv(res,file = paste0("mendelian_new/","res",exp,"-",id,".csv"), row.names = FALSE)
  write.csv(dat,file = paste0("mendelian_new/","dat",exp,"-",id,".csv"), row.names = FALSE)
  print(paste(exp,"on",id,"is ok"))
}
}
}
