#################
### CHI SQUARE TEST FOR ENRICHMENT, NEW VS. OLD - USED FOR TABLE 1
### R code for chi-square tests or enrichment
# Takes a table with the number of base pairs of non-B vs 
# non-B for two groups (in this case new and old)
# Written by Linnéa Smeds, September 2024

# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)

# Reading in the data and convert to tibble
statfile=paste("new_sequence/human_stat_table.tsv", sep="")
tib<-statfile %>% read.table(header=TRUE) %>% as_tibble() 

################################################################################
#  CHI-SQUARE TEST OF INDEPENDENCE

nonB_types=tib$nonB %>% unique()
chr_types=tib$Chr %>% unique()

chr<-c()
nb<-c()
pval<-c() 

# Go through the chromosome types and nonB types and perform the test for each
for (c in 1:length(chr_types)) {
  for (i in 1:length(nonB_types)) {
      X<-tib %>% filter(nonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(-Chr,-nonB,-Region)  %>% as.matrix() %>% as.table() %>% chisq.test()
      p<-X$p.value
      chr<-c(chr,chr_types[c])
      nb<-c(nb,nonB_types[i])
      pval<-c(pval,p)
  }
}
ntest=i*c
RES_TABLE <- tibble(chr,nb,pval) 

RES_TABLE %>% as_tibble() %>% mutate(sign=case_when(pval<0.05/24 ~ "yes",
                               TRUE ~ "no")) %>% print(n=24)
################################################################################
# SIMILAR FOR THE RESAMPLED SETS 
resamfile=paste("new_sequence/resample.human_stat_table.tsv", sep="")
rtib<-resamfile %>% read.table(header=TRUE) %>% as_tibble() 


resampling=rtib$Resample %>% unique()
nonB_types=rtib$nonB %>% unique()
chr_types=rtib$Chr %>% unique()

rsam<-c()
chr<-c()
nb<-c()
pval<-c() 

# Go through the resampling, the chromosome types and nonB types 
# and perform the test for each
for (r in 1:length(resampling)) {
  for (c in 1:length(chr_types)) {
    for (i in 1:length(nonB_types)) {
      X<-rtib %>% filter(nonB==nonB_types[i] & Chr==chr_types[c] & Resample==resampling[r])  %>% 
        select(-Resample,-Chr,-nonB,-Region)  %>% as.matrix() %>% as.table() %>% chisq.test()
      p<-X$p.value
      rsam<-c(rsam,resampling[r])
      chr<-c(chr,chr_types[c])
      nb<-c(nb,nonB_types[i])
      pval<-c(pval,p)
    }
  }
}
ntest=i*c*r
SUBSAMPLE <- tibble(rsam,chr,nb,pval) 
SUBSAMPLE %>% print(n=ntest)

SUBSAMPLE %>% as_tibble() %>% 
  mutate(sign=case_when(pval<0.001/ntest ~ "yes",
                        TRUE ~ "no")) %>% 
  group_by(chr,nb,sign) %>% 
  summarize() %>% print(n=ntest)


