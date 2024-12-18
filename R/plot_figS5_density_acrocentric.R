################################################################################
### R code for plotting non-B motif density in p-arms or acrocentric chromosomes
### written by Linnéa Smeds 28-Nov-2024

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
acro_file<-"T2T_primate_nonB/helpfiles/Acrocentric.txt"

################################################################################
# DENFINE COLORS 
nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","#DB5829","#894B45")
cencol=c("red","orange", "darkred", "#ca6850", "#FFCC99", "#00DE60","#1B998B", "#0080FA", "#335189","#FA99FF", "#AC33C7", "#00CCCC", "grey", "grey20")

################################################################################
# Processing files and merge data into a single tibble

# PRI
len_tib_pri <- lengths_file_pri %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Species=V1, ChrLong=V2, Len=V3, Chr=V4) %>% mutate(Hap="pri")
dens_tib_pri <- density_file_pri  %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, NonB=V6, ChrLong=V2, Start=V3, End=V4, Count=V5) %>%
  mutate(Hap="pri")
cen_tib_pri <- centro_file1_pri %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, ChrLong=V2, Cen_start=V3, Cen_end=V4) %>% mutate(Hap="pri")
# ALT
len_tib_alt <- lengths_file_alt %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Species=V1, ChrLong=V2, Len=V3, Chr=V4) %>% mutate(Hap="alt")
dens_tib_alt <- density_file_alt  %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, NonB=V6, ChrLong=V2, Start=V3, End=V4, Count=V5) %>%
  mutate(Hap="alt")
cen_tib_alt <- centro_file1_alt %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Species=V1, ChrLong=V2, Cen_start=V3, Cen_end=V4) %>% mutate(Hap="alt")

pri<-dens_tib_pri %>% left_join(len_tib_pri) %>% left_join(cen_tib_pri)
alt<-dens_tib_alt %>% left_join(len_tib_alt) %>% left_join(cen_tib_alt)

# Acrocentric information
acro_tib <- acro_file %>% read.table(header=FALSE) %>% as_tibble() %>%
    rename(Species=V1, Hap=V2, Chr=V3) %>% mutate(Acro="Yes")

# Combine density data
comb_tib <- pri %>% bind_rows(alt) %>% 
  filter(NonB!="all" & Chr!="chr1522") %>% 
#  mutate(Cen_mid=(Cen_end-Cen_start)/2+Cen_start) %>% 
  filter(End<Cen_end+100000) %>%
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
  inner_join(acro_tib) %>% mutate(PrintName=paste(Chr,Hap, sep=" "))
comb_tib$Species <- factor(comb_tib$Species, levels=c("chimp", "bonobo", "human", "gorilla", "borang", "sorang"))
comb_tib$NonB <- factor(comb_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))

#Satellite data
sat_pri <- sat_file_pri %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Hap="pri") %>% 
  left_join(len_tib_pri) %>% left_join(cen_tib_pri, relationship = "many-to-many") 
sat_alt <- sat_file_alt %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Hap="alt") %>% 
  left_join(len_tib_alt) %>% left_join(cen_tib_alt, relationship = "many-to-many") 

sat_tib <- sat_pri %>% bind_rows(sat_alt) %>%
  mutate(FullName=case_when(Sat=="active_asat" ~ paste("Active ","\U03B1","sat", sep=""),
                            Sat=="inactive_asat"  ~ paste("Inactive ","\U03B1","sat", sep=""),
                            Sat=="dhor" ~ paste("Divergent ","\U03B1","sat", sep=""),
                            Sat=="mon/hor" | Sat=="mixedsuperfamily" ~ paste("Mixed ","\U03B1","sat", sep=""),
                            Sat=="mon" ~ paste("Monomeric ","\U03B1","sat", sep=""),
                            Sat=="hsat1a" ~ "HSAT1A",
                            Sat=="hsat1b" ~ "HSAT1B",
                            Sat=="hsat2" ~ "HSAT2",
                            Sat=="hsat3" ~ "HSAT3",
                            Sat=="bsat" ~ "BSAT",
                            Sat=="gsat" ~ "GSAT",
                            Sat=="censat" ~ "Other cen sat",
                            Sat=="gap" ~ "GAP",
                            Sat=="subterm" ~ "Subterminal",
                            Sat=="rdna" ~ "rDNA",
                            TRUE ~ NA)) %>%
  mutate(Sat_stop=if_else(Sat_stop>Cen_end & Sat_start<Cen_end, Cen_end, Sat_stop)) %>%
  filter(Sat_stop<=Cen_end) %>% mutate(PrintName=paste(Chr,Hap, sep=" ")) %>% inner_join(acro_tib)


################################################################################
# PLOT ACROCENTRIC P-ARMS

sp<-c("bonobo", "chimp", "human", "gorilla", "borang", "sorang")
plotname<-c("<b>A.</b> Bonobo", "<b>B.</b> Chimpanzee", "<b>C.</b> Human", "<b>D.</b> Gorilla", "<b>E.</b> Bornean orangutan", "<b>F.</b> Sumatran orangutan")

# Plot all as single pages
for (i in 1:6) {
  sub_tib <- comb_tib %>% filter(Species==sp[i]) %>% 
      select(Species,NonB,norm,Start,End,Chr,Cen_end,PrintName) %>% distinct()
  chr_vect <- sub_tib$Chr %>% sort() %>% unique()
  sub_tib$Chr <- factor(sub_tib$Chr, levels=chr_vect)
  sub_tib <- sub_tib %>% arrange(Chr)
  print_vect <- sub_tib$PrintName %>%  unique()
  sub_tib$PrintName <- factor(sub_tib$PrintName, levels=print_vect)
  sub_tib$NonB <- factor(sub_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
  sub_sat <- sat_tib %>% filter(Species==sp[i]) %>% filter(Sat!="ct")
  sub_sat$PrintName <- factor(sub_sat$PrintName, levels=print_vect)
  
  pacr<-ggplot(sub_tib, aes(Start/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
    geom_tile(aes(alpha=norm), color="transparent") +
    scale_alpha_continuous(range = c(0.05, 1), name="Norm. non-B density")+
    facet_grid(PrintName~., scales="free",  switch='y')+
    scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values =nonBcol, name="")+
    labs(x="Mb", y="") +
    geom_line(data=sub_sat %>% filter(Sat!="rdna"), aes(x=Sat_start/1000000, y=9, color=FullName), linewidth=3, inherit.aes = FALSE) + 
    geom_rect(aes(xmin=0, ymin=9, ymax=9, xmax=Cen_end/1000000), color="white", fill="white", linewidth=5, show.legend = FALSE) +
    geom_rect(data=sub_sat %>% filter(Sat!="rdna"), aes(xmin=Sat_start/1000000, ymin=9, ymax=9, xmax=Sat_stop/1000000, color=FullName), fill="white", linewidth=3, inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(breaks=c(paste("Active ","\U03B1","sat", sep=""), paste("Inactive ","\U03B1","sat", sep=""), paste("Divergent ","\U03B1","sat", sep=""), paste("Mixed ","\U03B1","sat", sep=""), paste("Monomeric ","\U03B1","sat", sep=""),  "HSAT1A", "HSAT1B", "HSAT2", "HSAT3", "BSAT", "GSAT", "Other cen sat", "GAP", "Subterminal"), values=cencol, name="")+
    geom_segment(data=sub_sat %>% filter(Sat=="rdna"), aes(x=Sat_start/1000000, xend=Sat_stop/1000000, y=8.8, yend=8.8), color="black", linewidth=0.4, inherit.aes = FALSE, arrow=arrow(type="closed", length = unit(0.07, "cm"))) +
    theme(panel.background = element_rect(fill = 'white', colour="white"),
          strip.background = element_rect(fill = 'white', colour="white"),
          strip.text = element_text(size=11),
          plot.title = element_markdown(hjust=0, size=16),
          legend.position="bottom",
          legend.box="vertical", 
          legend.text = element_text(size = 8),
          panel.grid.major = element_line(colour = 'white'),
          panel.spacing.y=unit(0.3, "lines"),
          axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
          axis.ticks.y =element_blank(),
          axis.text.y =element_blank())  +
    labs(title = plotname[i])+
    scale_x_continuous(expand = c(0,0))+
    guides(fill=guide_legend(nrow=1, size=0.5),
           color=guide_legend(nrow=2),
    )+
    scale_linetype_manual()
  pacr
  output<- paste("plots/FigS5_acrocentric_p-arm_",sp[i],".png", sep="")
  if(i==3){
    ggsave(output, plot = pacr,dpi = 600,limitsize = TRUE, width=9,height=8)
  }
  else{
    ggsave(output, plot = pacr,dpi = 600,limitsize = TRUE, width=9,height=12)
  }
  
}  
