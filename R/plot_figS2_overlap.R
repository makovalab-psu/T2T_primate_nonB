################################################################################
### R code for plotting OVERLAP OF NON-B DNA TYPES used for Figure S2
### Written by Linnéa Smeds 2024/11/18

# Note 1: This analysis is only performed on the primary haplotypes
# Note 2: The (almost) identical code for producing the human figure is found in
# plot_fig3_overlap.R
# Note 3: The number of overlapping types in the bottom of the upset plot was 
# added manually

# Setting up, loading R libraries 
rm(list=ls())
require(tidyverse)
require(ggupset)
require(ggtext)
require(patchwork)
plotdir="plots/"
options(scipen = 999, decimals=2)

spvect=c("bonobo", "chimp", "gorilla", "borang", "sorang", "siamang")
# Make a hash with the species (needs R >= v4.2.0)
h1 <- hashtab()
sethash(h1, "bonobo", 1)
sethash(h1, "chimp", 2)
sethash(h1, "gorilla", 3)
sethash(h1, "borang", 4)
sethash(h1, "sorang", 5)
sethash(h1, "siamang", 6)
# Assign header letters 
head1_vec=c("A", "C", "E", "G", "I", "K")
head2_vec=c("B", "D", "F", "H", "J", "L")
suff_vec=c("Bonobo", "Chimpanzee", "Gorilla", "B. orangutan", "S. orangutan", "Siamang")

for (i in (1:length(spvect))) {
  
  species=spvect[i]

  # Reading in the data, create tibbles 
  fileA=paste("overlap/",species,".summary.autosomes.txt", sep="")
  fileX=paste("overlap/",species,".summary.chrX.txt", sep="")
  fileY=paste("overlap/",species,".summary.chrY.txt", sep="")
  filepw=paste("overlap/",species,".pairwise.summary.txt", sep="")
  tibA<-fileA %>% read.table(header=FALSE) %>% as_tibble() %>% 
    rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
    mutate(Types=strsplit(Comb, "-"), Chr="Autosomes", Sp=species)
  tibX<-fileX %>% read.table(header=FALSE) %>% as_tibble() %>% 
    rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
    mutate(Types=strsplit(Comb, "-"), Chr="Chromosome X", Sp=species)
  tibY<-fileY %>% read.table(header=FALSE) %>% as_tibble() %>% 
    rename(Comb=V1, Bp=V2) %>% mutate(Comb = str_replace_all(Comb, c('GQ' = 'G4'))) %>%
    mutate(Types=strsplit(Comb, "-"), Chr="Chromosome Y", Sp=species)
  tibP<-filepw %>% read.table(header=TRUE) %>% as_tibble() %>% 
    mutate(textcol=if_else(Frac>0.5,"white","black")) %>% 
    mutate(NB1 = str_replace_all(NB1, c('GQ' = 'G4')), 
           NB2 = str_replace_all(NB2, c('GQ' = 'G4'))) %>%
    mutate(Chr=case_when(Chr=="autosomes" ~ "Autosomes",
                         Chr=="chrX" ~ "Chromosome X",
                         Chr=="chrY" ~ "Chromsome Y"), Sp=species)
  tibP$NB1=factor(tibP$NB1, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
  tibP$NB2=factor(tibP$NB2, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
  
  ################################################################################
  # PART1 upset plot
  # Only plot groups with more than 10kb intersect 
  DATA1 <- tibA %>% bind_rows(tibX) %>% bind_rows(tibY) %>% 
    filter(Bp>10000) %>% arrange(desc(Bp)) %>% mutate(PrintStrip=paste(Chr,suff_vec[h1[[species]]], sep=" - "))
  
  pup<-ggplot(DATA1, aes(x=Types, fill=Chr)) +
    geom_col(aes(y=Bp/1000000), color="darkred", fill="darkred", alpha=0.8, show.legend = FALSE) +
    facet_wrap(~PrintStrip, nrow=3, scales="free_y")+
    scale_x_upset(order_by="degree")+
    labs(x="") +
    theme_combmatrix(combmatrix.panel.point.color.fill = "darkred",
                     combmatrix.panel.line.color = "darkred")+
    theme(panel.background = element_rect(fill = 'white', colour="white"),
          panel.grid.major = element_line(colour = 'beige'),
          strip.background = element_rect(fill = 'white', colour="white"),
          plot.title = element_text(size=18, face="bold"),
          strip.text = element_markdown(hjust=0, size=11, face="bold"),
          panel.spacing = unit(1, "lines"),
          axis.ticks.x =element_blank())+
    scale_y_continuous(name='Mb')+
    labs(title = head1_vec[h1[[species]]])
  pup
  
  ################################################################################
  # PART2 heatmap with pairwise overlap
  DATA2 <- tibP %>% mutate(PrintStrip=paste(Chr,suff_vec[h1[[species]]], sep=" - "))
  
  ppw<-ggplot(DATA2, aes(x=NB2, y=fct_rev(NB1), fill=Frac)) +
    geom_tile(color="white", show.legend = FALSE) +
    geom_text(aes(label=round(Frac*100, digits=2), color=textcol), size=3, show.legend=FALSE) +
    facet_wrap(~PrintStrip, nrow=3, scales='free') +
    labs(x="", y="") +
    scale_fill_gradient(low = "#fff4f0", high = "darkred", name="Overlap", na.value = 'white') +
    scale_color_manual(values=c("black","white")) +
    theme(panel.grid.major = element_line(colour = 'white'),
          panel.grid.minor = element_line(colour = 'white'),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(0, "lines"),
          legend.position="bottom",
          plot.title = element_text(size=18, face="bold"),
          legend.title=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.ticks.y=element_blank(),
          strip.background = element_rect(fill='white'),
          strip.text = element_text(hjust=0, size=11, face="bold"))+
    labs(title = head2_vec[h1[[species]]])
  ppw
  
  
  ################################################################################
  # COMBINE Upsetplot with Heatmap
  pcomb <- pup + ppw + plot_layout(widths = c(2.5, 1)) 
  pcomb
  
  outfile=paste(plotdir,"FigS2_overlap_",species,".png", sep="")
  ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=10,height=10)
  
} 
