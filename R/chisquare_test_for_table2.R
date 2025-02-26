#################
### CHI SQUARE TEST FOR ENRICHMENT, REPEATS VS. NON-REPEATS - USED FOR TABLE 2
### R code for chi-square tests or enrichment
# Takes a table with the number of base pairs of non-B vs 
# non-B for two groups (in this case repeats and non-repeats)
# Written by Linnéa Smeds, January 2025

# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)

# Reading in the data and convert to tibble
statfile=paste("repeats/7sp_repeat_vs_nonrepeat_contingency.txt", sep="")
tib<-statfile %>% read.table(header=TRUE) %>% as_tibble() 

################################################################################
#  CHI-SQUARE TEST OF INDEPENDENCE

nonB_types=tib$NonB %>% unique()
sp_types=tib$Species %>% unique()

sp<-c()
nb<-c()
pval<-c() 

# Go through the species and nonB types and perform the test for each
for (c in 1:length(sp_types)) {
  for (i in 1:length(nonB_types)) {
      X<-tib %>% filter(NonB==nonB_types[i] & Species==sp_types[c])  %>% 
        select(-Species,-NonB, -Type)  %>% as.matrix() %>% as.table() %>% chisq.test()
      p<-X$p.value
      sp<-c(sp,sp_types[c])
      nb<-c(nb,nonB_types[i])
      pval<-c(pval,p)
  }
}
ntest=i*c
RES_TABLE <- tibble(sp,nb,pval) 

RES_TABLE %>% as_tibble() %>% mutate(sign=case_when(pval<0.05/ntest ~ "yes",
                               TRUE ~ "no")) %>% print(n=ntest)
################################################################################
# SIMILAR FOR THE RESAMPLED SETS 
resamfile=paste("repeats/resample.10perc.stat_table.tsv", sep="")
rtib<-resamfile %>% read.table(header=TRUE) %>% as_tibble() 

resampling=rtib$Resample %>% unique()
nonB_types=rtib$nonB %>% unique()
sp_types=rtib$Species %>% unique()

rsam<-c()
sp<-c()
nb<-c()
pval<-c() 

# Go through the species, the resampling and nonB types 
# and perform the test for each
for (c in 1:length(sp_types)) {
  for (r in 1:length(resampling)) {
    for (i in 1:length(nonB_types)) {
      X<-rtib %>% filter(nonB==nonB_types[i] & Species==sp_types[c] & Resample==resampling[r])  %>% 
        select(-Resample,-Species,-nonB,-Region)  %>% as.matrix() %>% as.table() %>% chisq.test()
      p<-X$p.value
      rsam<-c(rsam,resampling[r])
      sp<-c(sp,sp_types[c])
      nb<-c(nb,nonB_types[i])
      pval<-c(pval,p)
    }
  }
}
ntest=i*c*r
SUBSAMPLE <- tibble(rsam,sp,nb,pval) 
SUBSAMPLE %>% print(n=ntest)

SUBSAMPLE %>% as_tibble() %>% 
  mutate(sign=case_when(pval<0.01/ntest ~ "yes",
                        TRUE ~ "no")) %>% 
  group_by(sp,nb,sign) %>% 
  summarize() %>% filter(sign=="no") %>% 
  print(n=ntest)


