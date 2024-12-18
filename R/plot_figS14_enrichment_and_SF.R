################################################################################
### R code for plotting enrichment in centromeres belonging to different 
### suprachromosomal families.
### written by Linnéa Smeds 10-Dec-2024

#Note: Numbers per SF family were added manually.
################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork)  

plotdir="plots/"
options(scipen = 999, decimals=2)

# Input files 
pri_file="centromeres/pri_6sp_enrichment_merged.tsv"
alt_file="centromeres/alt_6sp_enrichment_merged.tsv"
priSF_file="centromeres/pri_6sp_SF_merged.tsv"
altSF_file="centromeres/alt_6sp_SF_merged.tsv"

# Reading in the data 
# Centromeres, fold enrichment
tib1<-pri_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
      pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
      mutate(Hap="pri")
tib2<-alt_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
      pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
      mutate(Hap="alt")
# SF annotation
tib3<-priSF_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
      mutate(Hap="pri")
tib4<-altSF_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Hap="alt")


DATA <- bind_rows(tib1,tib2) %>% inner_join(bind_rows(tib3,tib4)) %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))

# Check how many does have multiple SF annotated 
DATA %>% filter(str_detect(SF, ",")) %>% group_by(Species, LongChr) %>% summarize()
# Skip those 
DATA<-DATA %>% filter(!str_detect(SF, ","))


data_counts <- DATA %>% group_by(SF, Species, LongChr) %>% summarize() %>%
  ungroup() %>% group_by(SF) %>%summarise(count = n())

DATA<-DATA %>% left_join(data_counts, by = "SF")



nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","#DB5829","#894B45" )
################################################################################
# PLOT Centromere non-B enrichment comparison between SF types 

DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
DATA$SF <- factor(DATA$SF, levels=c("SF01","SF1", "SF2", "SF3", "SF4", "SF5"))

p<-ggplot(DATA, aes(x=SF, y=value, fill=NonB, color=NonB, width=count)) +
  geom_boxplot(show.legend=TRUE, varwidth = TRUE, alpha=0.6) +
  facet_wrap(~NonB, nrow=1)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  labs(x="Suprachromosomal family", y="Fold-enrichment") +
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=1)+
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'lightgray', colour="white"),
        strip.text = element_blank(),
        panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        legend.background = element_rect(color = NA),
        legend.position="top",
        legend.justification = "right",
        legend.title = element_blank(),
        plot.title = element_text(size=18),
       axis.text.x=element_text(angle = 55, vjust = 1, hjust=1, size=8),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))
p

outfile=paste(plotdir,"FigS14_Centromeres_per_SF.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=6)

