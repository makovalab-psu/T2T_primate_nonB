################################################################################
### R code for plotting enrichment in functional regions, G4 motifs are 
### corrected for GC content. 
### Written by Linnéa Smeds Nov-2024

################################################################################
# Setting up, loading R libraries 
rm(list=ls())
require(tidyverse)
plotdir="plots/"

################################################################################
##### Reading in the data
sp="human"
file1<-paste("functional/",sp,"/enrichment_fullgenome.tsv", sep="")
file2<-paste("functional/",sp,"/enrichment_GQ_corr_fullgenome.tsv", sep="")

################################################################################
##### Convert to tibble 
tib_full <- file1 %>% read.table(header=TRUE) %>% as_tibble() %>%
    pivot_longer(c(-Region,-Class), names_to="NonB", values_to="Enrichment") %>% 
  mutate(Type="Corrected")

tib_GC <- file2 %>% read.table(header=TRUE) %>% as_tibble() %>%
  mutate(APR=0, DR=0, IR=0, MR=0, STR=0, Z=0) %>% 
  pivot_longer(c(-Region,-Class), names_to="NonB", values_to="Enrichment") %>%
  mutate(Type="Uncorrected")


comb_tib <- tib_full %>% bind_rows(tib_GC) %>% 
    mutate(PrintName=case_when(Class=="CDS_ProtCode" ~ "CDS - protein coding genes",  
                              Class=="CpG_Islands" ~ "CpG islands",    
                              Class=="Origin_of_Repl" ~ "Origins of replication",
                              Class=="Exons_ProtCode" ~ "Exons - protein coding genes",
                              Class=="Gene_NonProtCode"  ~"Non-protein coding genes",
                              Class=="Gene_ProtCode" ~ "Protein coding genes",
                              Class=="NGNR" ~ "Non-functional non-repetitive",
                              Class=="RNA_ProtCode" ~ "Transcript - protein coding genes",
                              Class=="3UTR" ~ "3'UTRs",
                              Class=="5UTR" ~ "5'UTRs",
                              TRUE ~ Class)) %>%  mutate(NonB=str_replace_all(NonB, c('GQ' = 'G4')))

# Gene_ProtCode, RNA_ProtCode, Exon_protCode are just summaries of introns+UTR+CDS
# We skip them for now 
comb_tib <- comb_tib %>% filter(Class!="Gene_ProtCode" & Class!="RNA_ProtCode"& Class!="Exons_ProtCode")
comb_tib$Type <- factor(comb_tib$Type, levels=c("Corrected", "Uncorrected"))
comb_tib$PrintName <- factor(comb_tib$PrintName, levels=c("Origins of replication", "Enhancers", "Promoters", "5'UTRs", "CDS - protein coding genes", "Introns", "3'UTRs", "Non-protein coding genes", "CpG islands", "Repeats", "Non-functional non-repetitive"))
comb_tib$NonB <- factor(comb_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))



nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45")
################################################################################
# PLOT FULL GENOME IN ONE FIGURE
set1<-comb_tib %>% filter(Type=="Corrected")
set2<-comb_tib %>% filter(Type=="Uncorrected")
p<-ggplot(set1, aes(x=PrintName, y=Enrichment, fill=NonB, color=NonB)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
  geom_bar(data=set2, aes(x=PrintName, y=Enrichment, fill=NonB), stat="identity", position=position_dodge(), color="black", alpha=0.8, show.legend = FALSE) +
  #facet_wrap(~Class, nrow=1, scale="free_x") +
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="top",
        legend.justification = "right",
        legend.text  = element_text(size = 8),
        legend.title=element_blank(),
        axis.text = element_text(margin = margin(r = 0, b=0)),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 55, hjust=1),
        plot.margin = margin(r = 30),
        strip.text = element_blank())+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + 
  ylab("Fold enrichment")+
  guides(fill = guide_legend(nrow = 1))
p

outfile=paste(plotdir,"FigS11_Enrichment_functional_G4corrected_",sp,".png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=4)
