################################################################################
### R code for plotting enrichment enrichment in different functional regions
### Written by Linnéa Smeds Oct-2024

################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
plotdir="plots/"
options(scipen = 999, decimals=2)

##### Reading in the data
sp="human"
file<-paste("functional/",sp,"/enrichment_fullgenome.tsv", sep="")

# Convert to tibble 
tib <- file %>% read.table(header=TRUE) %>% as_tibble() %>%
    pivot_longer(c(-Region,-Class), names_to="NonB", values_to="Enrichment") %>%
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
                              TRUE ~ Class)) %>%  
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))

# Gene_ProtCode, RNA_ProtCode, Exon_protCode are just summaries of introns+UTR+CDS
# We skip them for now 
DATA <- tib %>% filter(Class!="Gene_ProtCode" & Class!="RNA_ProtCode"& Class!="Exons_ProtCode")
DATA$PrintName <- factor(DATA$PrintName, levels=c("Origins of replication", "Enhancers", "Promoters", "5'UTRs", "CDS - protein coding genes", "Introns", "3'UTRs", "Non-protein coding genes", "CpG islands", "Repeats", "Non-functional non-repetitive"))
DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))

nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","#DB5829","#894B45" )
################################################################################
# PLOT ENRICHMENT (UNCORRECTED)
p<-ggplot(DATA, aes(x=PrintName, y=Enrichment, fill=NonB, color=NonB)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
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

outfile=paste(plotdir,"Fig5_Enrichment_functional_",sp, ".png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=4)
