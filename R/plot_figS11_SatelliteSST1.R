################################################################################
### R code for plotting non-B motifs occuring in the satellite SST1
### Written by Linnéa Smeds 12-Dec-2024

################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)

# Input files 
sf1_repeat_file <- "repeats/human/SST1/chr13.exactRep.bed"
Zfile <- "repeats/human/SST1/chr13.Z.bed"
DRfile <- "repeats/human/SST1/chr13.DR.bed"
MRfile <- "repeats/human/SST1/chr13.MR.bed"
TRIfile <- "repeats/human/SST1/chr13.TRI.bed"
G4file <- "repeats/human/SST1/chr13.GQ.bed"
STRfile <- "repeats/human/SST1/chr13.STR.bed"

################################################################################
# Processing files and merge data into a single tibble

sf1_tib <- sf1_repeat_file %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="monomer")

Z_tib <- Zfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="Z",Ypos=1)
MR_tib <- MRfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="MR",Ypos=4)
STR_tib <- STRfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="STR",Ypos=5)
G4_tib <- G4file %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="G4",Ypos=2)
DR_tib <- DRfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="DR",Ypos=6)
TRI_tib <- TRIfile %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Chr=V1, Start=V2, Stop=V3) %>% mutate(Type="nonB") %>% 
  mutate(NonB="TRI",Ypos=3)

NONB<-bind_rows(Z_tib,MR_tib,STR_tib,G4_tib,DR_tib,TRI_tib)
NONB$NonB <- factor(NONB$NonB, levels=c("DR", "STR", "MR", "TRI", "G4", "Z"))

nonBcol=c("#7BB0DF","#008A69","#F4A637", "darkorange","#DB5829","#894B45" )
#nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","#DB5829","#894B45" )
################################################################################
# PLOT CHR13, 6 repeat units

p<-ggplot(NONB, aes(color=NonB, fill=NonB)) +
  geom_hline(yintercept=0.25, linetype="solid", color = "black", linewidth=1)+
  geom_segment(aes(x=Start/1000000, xend=Stop/1000000, y=Ypos, yend=Ypos), linewidth=5) +
  scale_fill_manual(breaks=c("DR", "STR", "MR", "TRI", "G4", "Z"), values=nonBcol, name="")+
  scale_color_manual(breaks=c("DR", "STR", "MR", "TRI", "G4", "Z"), values=nonBcol, name="")+
  labs(x="Pos Chr13 (Mb)", y="") +
  geom_rect(data=sf1_tib, aes(xmin=Start/1000000, xmax=Stop/1000000, ymin=0, ymax=0.5), color="black", fill="white", linewidth=0.4, inherit.aes=FALSE) +
   theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.minor = element_line(colour = 'beige'),
        legend.position="right",
        legend.box="vertical", 
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(colour = 'beige'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  scale_x_continuous(expand = c(0,0), limits=c(12.302300,12.310800))+
  scale_y_continuous(limits=c(0,6), breaks=c(1,2,3,4,5,6), labels=c("Z", "G4", "TRI", "MR", "STR", "DR"))
p

output="plots/FigS11_Human_SST1_chr13_zoom.png"
ggsave(output, plot = p,dpi = 600,limitsize = TRUE, width=12,height=2)
