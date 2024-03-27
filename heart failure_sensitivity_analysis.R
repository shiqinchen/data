setwd("D:/brain structure/sensitivity_analysis")
library(TwoSampleMR)
library(data.table)
library(tidyverse)
dir.create(path = "heart_failure") 
ids <- c("ENIGMA3_mixed_se_wo_Mean_parsorbitalis_thickavg_20200522",
         "ENIGMA3_mixed_se_wo_Mean_lateralorbitofrontal_thickavg_20200522",
         "ENIGMA3_mixed_se_wTHICK_Mean_parsorbitalis_thickavg_20200522",
         "ENIGMA3_mixed_se_wo_Mean_temporalpole_thickavg_20200522",
         "ENIGMA3_mixed_se_wo_Mean_lingual_thickavg_20200522",
         "ENIGMA3_mixed_se_wo_Mean_precuneus_thickavg_20200522")
for (id in ids){
exposure_dat <-fread(file=paste0("D:/brain structure/exposure/","GCST90274223.csv"),header = T)
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
write.csv(res,file = paste0("heart_failure/",id,"-res",".csv"), row.names = FALSE)
write.csv(dat,file = paste0("heart_failure/",id,"-dat",".csv"), row.names = FALSE)

#敏感性分析
##异质性检验
het <- mr_heterogeneity(dat)
het
openxlsx::write.xlsx(het,file = paste0("heart_failure/",id,"-het.xlsx"), 
                     rowNames = FALSE)
#run_mr_presso(dat,NbDistribution = 5000)#异质性检验Q_离群点
 # res_presso <- TwoSampleMR::run_mr_presso(dat, NbDistribution = 5000)
 # sink(paste0("heart_failure/",id,"_PRESSO.txt"),
 #      append=FALSE,split = FALSE) 
 # print(res_presso)
 # sink()
 # print(res_presso)
#多效性检验
pleio <- mr_pleiotropy_test(dat)
pleio
openxlsx::write.xlsx(pleio,file = paste0("heart_failure/",id,"-pleio.xlsx"), 
                     rowNames = FALSE)

#森林图
res_single <- mr_singlesnp(dat)
forest <- mr_forest_plot(res_single)
pdf(paste0("heart_failure/",id,"_forest.pdf"))
print(forest[[1]])
dev.off()
#散点图
mr_scatter <- mr_scatter_plot(res,dat)
pdf(paste0("heart_failure/",id,"_scatter.pdf"))
print(mr_scatter[[1]])
dev.off()
#漏斗图
funnel <- mr_funnel_plot(res_single)
pdf(paste0("heart_failure/",id,"_funnel.pdf"))
print(funnel[[1]])
dev.off()
#逐个剔除检验
single <- mr_leaveoneout(dat)
leaveoneout <- mr_leaveoneout_plot(single)
pdf(paste0("heart_failure/",id,"_leaveoneout.pdf"))
print(mr_leaveoneout_plot(single))
dev.off()

print(paste(id,"is ok"))
}