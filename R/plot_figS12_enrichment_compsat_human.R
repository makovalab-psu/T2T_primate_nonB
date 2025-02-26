################################################################################
### R code for plotting non-B DNA motifs enrichment in manually curated human 
### satellites and composite repeats.
### Written by Linnéa Smeds 30-Aug-2024

################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)

repeat_file <- "repeats/human_genome_compsat_enrichment.tsv"
lengths_file <-"repeats/human_genome_compsat_lengths.txt"
threshold <- 5000    # Example: 50000


################################################################################
# Processing files

# Read files and pivot table to put motifs in single column 
lentib=lengths_file %>% read.table(header=TRUE)
tib <- repeat_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
  left_join(lentib, by = c("Repeat" = "Repeat", "Type" = "Type")) %>% 
  pivot_longer(-c(Repeat,Len,Type), names_to="Motif") %>%
  mutate(PrintName=paste(Repeat," (",sprintf("%.2f",Len/1000000),"Mb)", sep="")) %>%
  drop_na() %>% mutate(head=case_when(Type=="newsat" ~ "A",
                                      Type=="compos" ~ "B"),
                       Motif=if_else(Motif=="GQ", "G4", Motif)) %>%
  mutate(PrintName = sub("SAT", "Sat", PrintName))

REPDATA <- tib %>% filter(Len>threshold) %>% mutate(logval=log2(value))
REPDATA$Motif <- factor(REPDATA$Motif, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))


################################################################################
# PLOTTING HORISONTAL 
  
scalestart=min(REPDATA %>% filter(value>0) %>%select(logval))
scaleend=max(REPDATA$logval)
scalebreak=-1*scalestart/(scaleend-scalestart)

# Removing repeats with less than threshold bp
p<-ggplot(REPDATA, aes(PrintName, fct_rev(Motif))) +
  geom_tile(aes(fill=logval), color = "white",lwd = 0.3,linetype = 1) +
  facet_wrap(head~., scales="free_x", ncol=1) +
  scale_fill_gradientn(name="log2(Enrichment)",
                       colors=c("#075AFF", "white","#FF0000"),
                       values=c(0, scalebreak, 1), na.value="#075AFF")+
  labs(x="", y="") +
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1, size=10),
        legend.text=element_text(size=10),
        legend.position="bottom",
        legend.title=element_text(size=10, vjust=1),
        legend.key.width  = unit(2, "lines"),
        legend.key.height = unit(0.7, "lines"),
        legend.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(colour = 'black', hjust=0, size=18, face="bold"),
        strip.background = element_rect(fill = 'white', colour="white"),
        panel.background = element_rect(fill = 'white', colour="white"))+
  scale_y_discrete(position="left")
p

outfile="plots/FigS12_enrichment_human_compsat.png"
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=8)
