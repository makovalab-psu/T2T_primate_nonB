################################################################################
### R code for plotting GC content in centromeres vs background windows.
### written by Linnéa Smeds Oct-2024

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
setwd("~/Documents/T2T_nonB_paper/")
require(tidyverse)
plotdir="plots/"
#options(scipen = 999, decimals=2)

################################################################################
# Reading in the data, both centromeres and background data for testing significance 
file="centromeres/pri_6sp_enrichment_merged.tsv"
bgfile="centromeres/pri_6sp_background_enrichment_merged.tsv"
gcfile="centromeres/pri_6sp_GCcont_merged.tsv"
bggcfile="centromeres/pri_6sp_background_GCcont_merged.tsv"

################################################################################
# Make tibbles 
# Start with GC files, centro and background 
gctib<-gcfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centromere")
bggctib<-bggcfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="background")


# Main file with centromere densities, merge with GC directly
tib<-file %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(head=case_when(Species=="chimp" ~ "B",
                        Species=="bonobo" ~ "A",
                        Species=="human" ~ "C",
                        Species=="gorilla" ~ "D",
                        Species=="borang" ~ "E",
                        Species=="sorang" ~ "F"))
tib <- tib %>% mutate(textcol=if_else(value>0.8*max(tib$value) | value<0.3,"white","black")) %>% 
    mutate(Type="centromere") %>% inner_join(gctib)

# Background file, merge with GC directly
bgtib<-bgfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species,-Window), names_to="NonB", values_to="value") %>%
  mutate(head=case_when(Species=="chimp" ~ "B",
                        Species=="bonobo" ~ "A",
                        Species=="human" ~ "C",
                        Species=="gorilla" ~ "D",
                        Species=="borang" ~ "E",
                        Species=="sorang" ~ "F")) %>%
  mutate(Type="background") %>% inner_join(bggctib)

################################################################################
# PLOT GC CONTENT IN BACKGROUND VS CENTROMERES 

gc1<- tib %>% select(Species, Chr,Type,GCcont,head) %>% distinct()
gc2<-bgtib %>% select(Species,Chr,Type,GCcont,head) %>% distinct()
gc_comb<- bind_rows(gc1,gc2)

p<-ggplot(gc_comb, aes(x=GCcont, fill=Type, color=Type)) +
  geom_density(alpha=0.6) +
  facet_wrap(~head, ) +
  scale_fill_manual(values=c("darkgrey", "red"))+
  scale_color_manual(values=c("darkgrey", "red"))+
 # geom_vline(data=Line_tib, aes(xintercept = GCcont), color = "#DB5829")+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(colour = 'black', hjust=0, size=16, face="bold"),
        legend.title=element_blank(),
        legend.position = "bottom")+
  xlab("GC content") + 
  ylab("Window density") 
p

outfile=paste(plotdir,"FigS3_GC_Centromere_vs_Background_pri.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=7)

