setwd("D:/brain structure/BNP/sensitivity_analysis")
library(TwoSampleMR)
library(data.table)
library(tidyverse)
dir.create(path = "BNP") 

exposure_dat <-fread(file=paste0("D:/brain structure/BNP/GCST90012082.csv"),header = T)
out <- fread(file=paste0("D:/brain structure/brain structure/ENIGMA3_mixed_se_wTHICK_Mean_supramarginal_thickavg_20200522.txt.gz"),header = T)
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
  outcome_dat$id.outcome <- 'supramarginal TH global weighted'
  dat <- harmonise_data(exposure_dat,outcome_dat)
  res <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_ivw",
                                                              "mr_egger_regression",
                                                              "mr_weighted_median")))
  write.csv(res,file = paste0("BNP/",'GCST90012082',"-res",".csv"), row.names = FALSE)
  write.csv(dat,file = paste0("BNP/",'GCST90012082',"-dat",".csv"), row.names = FALSE)
  
  #敏感性分析
  ##异质性检验
  het <- mr_heterogeneity(dat)
  het
  openxlsx::write.xlsx(het,file = paste0("BNP/",'GCST90012082',"-het.xlsx"), 
                       rowNames = FALSE)
  #run_mr_presso(dat,NbDistribution = 5000)#异质性检验Q_离群点
  res_presso <- TwoSampleMR::run_mr_presso(dat, NbDistribution = 5000)
  sink(paste0("BNP/",'GCST90012082',"_PRESSO.txt"),
       append=FALSE,split = FALSE)
  print(res_presso)
  sink()
  print(res_presso)
  #多效性检验
  pleio <- mr_pleiotropy_test(dat)
  pleio
  openxlsx::write.xlsx(pleio,file = paste0("BNP/",'GCST90012082',"-pleio.xlsx"), 
                       rowNames = FALSE)
  
  #森林图
  res_single <- mr_singlesnp(dat)
  forest <- mr_forest_plot(res_single)
  pdf(paste0("BNP/",'GCST90012082',"_forest.pdf"))
  print(forest[[1]])
  dev.off()
  #散点图
  mr_scatter <- mr_scatter_plot(res,dat)
  pdf(paste0("BNP/",'GCST90012082',"_scatter.pdf"))
  print(mr_scatter[[1]])
  dev.off()
  #漏斗图
  funnel <- mr_funnel_plot(res_single)
  pdf(paste0("BNP/",'GCST90012082',"_funnel.pdf"))
  print(funnel[[1]])
  dev.off()
  #逐个剔除检验
  single <- mr_leaveoneout(dat)
  leaveoneout <- mr_leaveoneout_plot(single)
  pdf(paste0("BNP/",'GCST90012082',"_leaveoneout.pdf"))
  print(mr_leaveoneout_plot(single))
  dev.off()

#汇总数据
library(readxl)
fs=list.files("BNP/",pattern = "-res",full.names = TRUE) 
df = map_dfr(fs,read.csv,header=T)
df_ivw <- subset(df,df$method=="Inverse variance weighted")
fs1 <- list.files("BNP/",pattern = "-het",full.names = TRUE)
df1 = map_dfr(fs1,read_xlsx)
df1_ivw <- subset(df1,df1$method=="Inverse variance weighted")
df1_mr_egger <- subset(df1,df1$method=="MR Egger")
fs2 <- list.files("BNP/",pattern = "-pleio",full.names = TRUE)
df2  <-  map_dfr(fs2,read_xlsx)

df3 <- 
  df %>% left_join(df1,by=c('id.outcome','id.exposure',
                            'outcome','exposure','method'))%>% 
  left_join(df2,by=c('id.outcome','id.exposure',
                     'outcome','exposure'))


df4 <- df3 %>% mutate(b=format(b,scientific = T,digits = 3)) %>% 
  mutate(pval=format(pval.x,scientific = T,digits = 3)) %>% 
  mutate(lo_ci=format(lo_ci,scientific = T,digits = 3)) %>% 
  mutate(up_ci=format(up_ci,scientific = T,digits = 3)) %>% 
  mutate(estimate= paste0(b, '(',lo_ci,
                          ',',up_ci,')'))
df5 <- df4 %>% select(exposure,id.outcome,nsnp,method,pval,estimate,Q,Q_pval,
                        egger_intercept,pval.y)
df5 <- df5 %>% mutate(Q=round(Q,2)) %>% 
  mutate(Q_pval=format(Q_pval,scientific = T,digits = 3)) %>% 
  mutate(egger_intercept=format(egger_intercept,scientific = T,digits = 3)) %>%
  mutate(pval.y=format(pval.y,scientific = T,digits = 3))

openxlsx::write.xlsx(df5,'sensitivity_analysis.xlsx',rowNames=F)