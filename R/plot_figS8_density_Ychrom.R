################################################################################
### R code for plotting non-B motif density in Y chromosomes chromosomes
### written by Linnéa Smeds 11-Feb-2025

#Note: Arrow in legend was added manually
################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
require(ggtext)

# Input files (files with all species combined)
density_file_pri <- "densities/pri_6sp_merged.txt"
lengths_file_pri <- "densities/pri_6sp_lengths.txt"
centro_file1_pri <-"centromeres/pri_6sp_centromeres_merged.tsv"
centro_file2_pri <-"centromeres/pri_6sp_centromeres_unmerged.tsv"
sat_file_pri<-"centromeres/pri_6sp_cenSat_merged.tsv"
density_file_alt <- "densities/alt_6sp_merged.txt"
lengths_file_alt <- "densities/alt_6sp_lengths.txt"
centro_file1_alt <-"centromeres/alt_6sp_centromeres_merged.tsv"
centro_file2_alt <-"centromeres/alt_6sp_centromeres_unmerged.tsv"
sat_file_alt<-"centromeres/alt_6sp_cenSat_merged.tsv"
#acro_file<-"T2T_primate_nonB/helpfiles/Acrocentric.txt"

################################################################################
# DENFINE COLORS 
nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45")
cencol=c("red","orange", "darkred", "#ca6850", "#FFCC99", "#00DE60","#1B998B", "#0080FA", "#335189","#FA99FF", "#AC33C7", "#00CCCC", "grey", "grey20")

################################################################################
# Processing files and merge data into a single tibble

# Only PRI (sex chromosomes are placed in the primary haplotype)
len_tib_pri <- lengths_file_pri %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Species=V1, ChrLong=V2, Len=V3, Chr=V4) 
dens_tib_pri <- density_file_pri  %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, NonB=V6, ChrLong=V2, Start=V3, End=V4, Count=V5)
cen_tib_pri <- centro_file1_pri %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, ChrLong=V2, Cen_start=V3, Cen_end=V4)

pri<-dens_tib_pri %>% left_join(len_tib_pri) %>% left_join(cen_tib_pri)

# Combine density data
comb_tib <- pri %>% filter(NonB!="all" & Chr=="chrY") %>% 
  group_by(NonB, Species) %>% mutate(max=max(Count, na.rm=TRUE))  %>%  
  mutate(norm=Count/max) %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>%
  mutate(color=case_when(NonB=='APR' ~ nonBcol[1],
                         NonB=='DR' ~ nonBcol[2],
                         NonB=='STR' ~ nonBcol[3],
                         NonB=='IR' ~ nonBcol[4],
                         NonB=='MR' ~ nonBcol[5],
                         NonB=='G4' ~ nonBcol[6],
                         NonB=='Z' ~ nonBcol[7],
                         TRUE ~ 'black')) %>%
  mutate(PrintName=case_when(Species=="bonobo" ~  "Bonobo",
                             Species=="chimp" ~ "Chimp.", 
                             Species=="human" ~ "Human", 
                             Species=="gorilla" ~ "Gorilla", 
                             Species=="borang" ~ "B. orang.", 
                             Species=="sorang" ~ "S. orang.",
                            TRUE ~ NA))

comb_tib$Species <- factor(comb_tib$Species, levels=c("chimp", "bonobo", "human", "gorilla", "borang", "sorang"))
comb_tib$NonB <- factor(comb_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))

#Satellite data
sat_pri <- sat_file_pri %>% read.table(header=TRUE) %>% as_tibble() %>%
  left_join(len_tib_pri) %>% left_join(cen_tib_pri, relationship = "many-to-many") 

sat_tib <- sat_pri %>% filter(Chr=="chrY") %>%
  mutate(FullName=case_when(Sat=="active_asat" ~ paste("Active ","\U03B1","Sat", sep=""),
                            Sat=="inactive_asat"  ~ paste("Inactive ","\U03B1","Sat", sep=""),
                            Sat=="dhor" ~ paste("Divergent ","\U03B1","sat", sep=""),
                            Sat=="mon/hor" | Sat=="mixedsuperfamily" ~ paste("Mixed ","\U03B1","Sat", sep=""),
                            Sat=="mon" ~ paste("Monomeric ","\U03B1","Sat", sep=""),
                            Sat=="hsat1a" ~ "HSat1A",
                            Sat=="hsat1b" ~ "HSat1B",
                            Sat=="hsat2" ~ "HSat2",
                            Sat=="hsat3" ~ "HSat3",
                            Sat=="bsat" ~ "BSat",
                            Sat=="gsat" ~ "GSat",
                            Sat=="censat" ~ "Other cen sat",
                            Sat=="gap" ~ "GAP",
                            Sat=="subterm" ~ "Subterminal",
                            Sat=="rdna" ~ "rDNA",
                            TRUE ~ NA)) %>% filter(Sat!="ct") %>% 
  mutate(PrintName=case_when(Species=="bonobo" ~  "Bonobo",
                             Species=="chimp" ~ "Chimp.", 
                             Species=="human" ~ "Human", 
                             Species=="gorilla" ~ "Gorilla", 
                             Species=="borang" ~ "B. orang.", 
                             Species=="sorang" ~ "S. orang.",
                             TRUE ~ NA))



################################################################################
# PLOT ACROCENTRIC P-ARMS

sp<-c("Bonobo", "Chimp.", "Human", "Gorilla", "B. orang.", "S. orang.")
DATA<-comb_tib %>%
  select(Species,NonB,norm,Start,End,Chr,Cen_end,PrintName) %>% distinct()
DATA$PrintName<-factor(DATA$PrintName, levels=sp)
sat_tib$PrintName<-factor(sat_tib$PrintName, levels=sp)

p<-ggplot(DATA, aes(Start/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  geom_tile(aes(alpha=norm), color="transparent") +
  scale_alpha_continuous(range = c(0.05, 1), name="Norm. non-B density")+
  facet_grid(PrintName~., scales="free",  switch='y')+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol, name="")+
  labs(x="Mb", y="") +
  geom_line(data=sat_tib %>% filter(Sat!="rdna"), aes(x=Sat_start/1000000, y=10, color=FullName), linewidth=3, inherit.aes = FALSE) + 
  geom_rect(aes(xmin=0, ymin=10, ymax=10, xmax=Cen_end/1000000), color="white", fill="white", linewidth=5, show.legend = FALSE) +
  geom_rect(data=sat_tib %>% filter(Sat!="rdna"), aes(xmin=Sat_start/1000000, ymin=10, ymax=10, xmax=Sat_stop/1000000, color=FullName), fill="white", linewidth=3, inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_manual(breaks=c(paste("Active ","\U03B1","Sat", sep=""), paste("Inactive ","\U03B1","Sat", sep=""), paste("Divergent ","\U03B1","Sat", sep=""), paste("Mixed ","\U03B1","Sat", sep=""), paste("Monomeric ","\U03B1","Sat", sep=""),  "HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSat", "Other cen sat", "GAP", "Subterminal"), values=cencol, name="")+
  geom_segment(data=sat_tib %>% filter(Sat=="rdna"), aes(x=Sat_start/1000000, xend=Sat_stop/1000000, y=9.8, yend=9.8), color="black", linewidth=0.4, inherit.aes = FALSE, arrow=arrow(type="closed", length = unit(0.07, "cm"))) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=11),
        plot.title = element_markdown(hjust=0, size=16),
        legend.position="bottom",
        legend.box="vertical", 
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(1, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  #labs(title = PrintName)+
  scale_x_continuous(expand = c(0,0))+
  guides(fill=guide_legend(nrow=1, size=0.5),
         color=guide_legend(nrow=2))+
  scale_linetype_manual()
p

output<- paste("plots/FigS8_chrY_density.png", sep="")
ggsave(output, plot = p,dpi = 600,limitsize = TRUE, width=9,height=11)
