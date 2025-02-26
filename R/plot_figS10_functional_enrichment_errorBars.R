################################################################################
### R code for plotting enrichment enrichment in different functional regions
### Written by Linnéa Smeds Oct-2024

################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
require(patchwork) 
plotdir="plots/"
options(scipen = 999, decimals=2)

##### Reading in the data
sp="human"
file<-paste("functional/",sp,"/enrichment_fullgenome.tsv", sep="")
fileVar10<-paste("functional/",sp,"/10perc.20rep.minmax.txt", sep="")
fileVar50<-paste("functional/",sp,"/50perc.20rep.minmax.txt", sep="")

vartib10<-fileVar10 %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))
vartib50<-fileVar50 %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))

# Convert to tibble 
# Gene_ProtCode, RNA_ProtCode, Exon_protCode are just summaries of introns+UTR+CDS
# We skip them for now 
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
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>% 
  filter(Class!="Gene_ProtCode" & Class!="RNA_ProtCode"& Class!="Exons_ProtCode")

tib$PrintName <- factor(tib$PrintName, levels=c("Origins of replication", "Enhancers", "Promoters", "5'UTRs", "CDS - protein coding genes", "Introns", "3'UTRs", "Non-protein coding genes", "CpG islands", "Repeats", "Non-functional non-repetitive"))
tib$NonB <- factor(tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))

# Combine with confidence intervals 
DATA10<-tib %>% inner_join(vartib10) %>%
  mutate(significance=if_else(Min<=1 & Max>=1, FALSE, TRUE))
DATA50<-tib %>% inner_join(vartib50) %>%
  mutate(significance=if_else(Min<=1 & Max>=1, FALSE, TRUE))


nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45")

################################################################################
# PLOT ENRICHMENT WITH ERROR BARS; 10% DATA, 20 REPL 
p1<-ggplot(DATA10, aes(x=PrintName, y=Enrichment, fill=NonB, color=NonB)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
  geom_errorbar(stat="identity", position=position_dodge(), aes(x=PrintName, ymin=Min, ymax=Max), colour="grey20", alpha=0.9)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="top",
        legend.justification = "right",
        legend.text  = element_text(size = 8),
        legend.title=element_blank(),
        plot.title = element_text(size=18, face="bold"),
        axis.text = element_text(margin = margin(r = 0, b=0)),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 55, hjust=1),
        plot.margin = margin(r = 30),
        strip.text = element_blank())+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14)) +
  xlab("") + 
  ylab("Fold enrichment")+ 
  ggtitle("A")+
  guides(fill = guide_legend(nrow = 1))
p1

# Code for adding asterisk to plot above (not used)
#  geom_text(aes(label = ifelse(significance, "*", ""), group = NonB), colour="grey20",
#            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt, show.legend = FALSE) +


################################################################################
# PLOT ENRICHMENT WITH ERROR BARS; 50% DATA, 20 REPL 

p2<-ggplot(DATA50, aes(x=PrintName, y=Enrichment, fill=NonB, color=NonB)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6, show.legend = FALSE) +
  geom_errorbar(stat="identity", position=position_dodge(), aes(x=PrintName, ymin=Min, ymax=Max), colour="grey20", alpha=0.9)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(size=18, face="bold"),
        axis.text = element_text(margin = margin(r = 0, b=0)),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 55, hjust=1),
        plot.margin = margin(r = 30),
        strip.text = element_blank())+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14)) +
  xlab("") + 
  ylab("Fold enrichment")+
  ggtitle("B")+
  guides(fill = guide_legend(nrow = 1))
p2

################################################################################
# COMBINE PLOTS
pcomb<-p1 + p2+ plot_layout(widths = c(1), heights= c(1.1,1)) 
pcomb
outfile=paste(plotdir,"FigS10_Enrichment_functional_DiffErrBars_",sp, ".png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=8)
