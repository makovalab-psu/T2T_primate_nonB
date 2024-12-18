################################################################################
### R code for plotting enrichment in Centromeric and acrocentric repeats 
### Written by Linnéa Smeds Dec-2024
################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
require(ggtext)

################################################################################
# Read files 
repeat_file <- "repeats/7sp_genome_censat_enrichment.tsv"
lengths_file <-"repeats/7sp_genome_censat_lengths.txt"
threshold <- 50000    # Example: 50000


################################################################################
# Processing files

# Read files and pivot table to put motifs in single column 
lentib=lengths_file %>% read.table(header=TRUE)
tib <- repeat_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
      left_join(lentib, by = c("Repeat" = "Repeat", "Species" = "Species")) %>% 
        pivot_longer(-c(Repeat,Len,Species), names_to="Motif") %>%
     mutate(head=case_when(Species=="bonobo" ~ "<b>A.</b> Bonobo",
                                          Species=="chimp" ~ "<b>B.</b> Chimpanzee",
                                          Species=="human" ~ "<b>C.</b> Human",
                                          Species=="gorilla" ~ "<b>D.</b> Gorilla",
                                          Species=="borang" ~ "<b>E.</b> Bornean orangutan",
                                          Species=="sorang" ~ "<b>F.</b> Sumatran orangutan",
                                          Species=="siamang" ~ "<b>G.</b> Siamang"),
                           Motif=if_else(Motif=="GQ", "G4", Motif),
                           logval=log2(value)) %>% 
      mutate(FullName=case_when(Repeat=="active_asat" ~ paste("Active ","\U03B1","sat", sep=""),
      Repeat=="inactive_asat"  ~ paste("Inactive ","\U03B1","sat", sep=""),
      Repeat=="dhor" ~ paste("Divergent ","\U03B1","sat", sep=""),
      Repeat=="mon/hor" ~ "Monomeric/HOR",
      Repeat=="mixedsuperfamily" ~ paste("Mixed ","\U03B1","sat", sep=""),
      Repeat=="mon" ~ paste("Monomeric ","\U03B1","sat", sep=""),
      Repeat=="hsat1a" ~ "HSAT1A",
      Repeat=="hsat1b" ~ "HSAT1B",
      Repeat=="hsat2" ~ "HSAT2",
      Repeat=="hsat3" ~ "HSAT3",
      Repeat=="bsat" ~ "BSAT",
      Repeat=="gsat" ~ "GSAT",
      Repeat=="censat" ~ "Other cen sat",
      Repeat=="gap" ~ "GAP",
      Repeat=="subterm" ~ "Subterminal",
      Repeat=="rdna" ~ "rDNA",
      TRUE ~ NA)) %>% 
  mutate(PrintName=paste(FullName," (",sprintf("%.2f",Len/1000000),"Mb)", sep="")) %>%
  drop_na() %>% filter(Repeat!="gap") %>% filter(Repeat!="censat")

################################################################################
# PLOTTING HORISONTAL 
# Removing lines that are all NA and repeats with less than threshold bp
DATA <- tib %>% filter(Len>threshold) 
DATA$Motif <- factor(DATA$Motif, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
DATA$FullName <- factor(DATA$FullName, levels=c(paste("Active ","\U03B1","sat", sep=""), paste("Inactive ","\U03B1","sat", sep=""), paste("Divergent ","\U03B1","sat", sep=""), paste("Mixed ","\U03B1","sat", sep=""), "Monomeric/HOR", paste("Monomeric ","\U03B1","sat", sep=""), "HSAT1A", "HSAT1B", "HSAT2", "HSAT3", "BSAT", "GSAT", "Subterminal", "rDNA"))
DATA$PrintName = factor(DATA$PrintName, levels = unique(DATA$PrintName[order(DATA$FullName)]))
scalestart=min(DATA %>% filter(value>0) %>%select(logval))
scaleend=max(DATA$logval)
scalebreak=-1*scalestart/(scaleend-scalestart)

p<-ggplot(DATA, aes(PrintName, fct_rev(Motif))) +
  geom_tile(aes(fill=logval), color = "white",lwd = 0.3,linetype = 1) +
  facet_wrap(head~., scales="free_x", ncol=2) +
  scale_fill_gradientn(name="log2(Enrichment)",
                       colors=c("#075AFF", "white","#FF0000"),
                       values=c(0, scalebreak, 1), na.value = '#075AFF')+
  labs(x="", y="") +
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1, size=10),
        axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="right",
        legend.title=element_text(size=12, hjust=0, angle=90),
        legend.title.position="left",
        legend.key.width  = unit(1, "lines"),
        legend.key.height = unit(3, "lines"),
        legend.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
        strip.text = element_markdown(hjust=0, size=18),
        strip.background = element_rect(fill = 'white', colour="white"),
        panel.background = element_rect(fill = 'white', colour="white"),
        panel.spacing = unit(1, "lines"))+
   scale_y_discrete(position="left")
p
outfile="plots/FigS6_Enrichment_CenSat.png"
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=14)
