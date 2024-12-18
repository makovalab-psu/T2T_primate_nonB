################################################################################
### R code for plotting OVERLAP OF NON-B DNA TYPES used for Figure 2
### Written by Linnéa Smeds 2024/11/12

# Note 1: This analysis is only performed on the primary haplotypes
# Note 2: Almost identical code for producing the non-human plots are found in
# plot_figS2_upset_and_tile.R

# Setting up, loading R libraries 
rm(list=ls())
require(tidyverse)
library(ggupset)
plotdir="plots/"
options(scipen = 999, decimals=2)
species="human"

# Reading in the data, create tibbles 
fileA=paste("overlap/",species,".summary.autosomes.txt", sep="")
fileX=paste("overlap/",species,".summary.chrX.txt", sep="")
fileY=paste("overlap/",species,".summary.chrY.txt", sep="")
filepw=paste("overlap/",species,".pairwise.summary.txt", sep="")

tibA<-fileA %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
  mutate(Types=strsplit(Comb, "-"), Chr="Autosomes", Head="A")
tibX<-fileX %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
  mutate(Types=strsplit(Comb, "-"), Chr="Chromosome X", Head="B")
tibY<-fileY %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
  mutate(Types=strsplit(Comb, "-"), Chr="Chromosome Y", Head="C")
tibP<-filepw %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(textcol=if_else(Frac>0.5,"white","black")) %>% 
  mutate(NB1 = str_replace_all(NB1, c('GQ' = 'G4')), 
         NB2 = str_replace_all(NB2, c('GQ' = 'G4'))) %>%
  mutate(Chr=case_when(Chr=="autosomes" ~ "Autosomes",
                       Chr=="chrX" ~ "Chromosome X",
                       Chr=="chrY" ~ "Chromsome Y"))
tibP$NB1=factor(tibP$NB1, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
tibP$NB2=factor(tibP$NB2, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))

################################################################################
# PART1 upset plot
# Only plot groups with more than 10kb intersect 
DATA <- tibA %>% bind_rows(tibX) %>% bind_rows(tibY) %>% 
  filter(Bp>10000) %>% arrange(desc(Bp))

pup<-ggplot(DATA, aes(x=Types, fill=Chr)) +
  geom_col(aes(y=Bp/1000000), color="darkred", fill="darkred", alpha=0.8, show.legend = FALSE) +
  facet_wrap(~Chr, nrow=3, scales="free_y")+
#  scale_fill_manual(values=c("#1964B0","#7BB0DF", "#B6DBFF"))+
  scale_x_upset(order_by="degree")+
  labs(x="") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "darkred",
                   combmatrix.panel.line.color = "darkred")+
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        panel.grid.major = element_line(colour = 'beige'),
        strip.background = element_rect(fill = 'white', colour="white"),
        plot.title = element_text(size=18),
        strip.text = element_markdown(hjust=0, size=13, face="bold"),
        panel.spacing = unit(1, "lines"),
        axis.ticks.x =element_blank())+
  scale_y_continuous(name='Mb')+
  labs(title = bold('A') ~ '')
pup

################################################################################
# PART2 heatmap with pairwise overlap

ppw<-ggplot(tibP, aes(x=NB2, y=fct_rev(NB1), fill=Frac)) +
  geom_tile(color="white", show.legend = FALSE) +
  geom_text(aes(label=round(Frac*100, digits=2), color=textcol), size=3, show.legend=FALSE) +
  facet_wrap(~Chr, nrow=3, scales='free') +
  labs(x="", y="") +
  scale_fill_gradient(low = "#fff4f0", high = "darkred", name="Overlap", na.value = 'white') +
  scale_color_manual(values=c("black","white")) +
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.grid.minor = element_line(colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(0, "lines"),
        legend.position="bottom",
        plot.title = element_text(size=18),
        legend.title=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(hjust=0, size=13, face="bold"))+
  labs(title = bold('B') ~ '')
ppw


################################################################################
# COMBINE Upsetplot with Heatmap

pcomb <- pup + ppw + plot_layout(widths = c(2.5, 1)) 
pcomb

outfile=paste(plotdir,"Figure2_upset_pairwise_",species,".png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=10)
