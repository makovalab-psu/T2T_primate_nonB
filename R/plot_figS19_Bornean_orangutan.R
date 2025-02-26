################################################################################
### R code for plotting alignment between Bornean orangutan chr10 haplotypes. 
### written by Linnéa Smeds 28-Oct-2024

rm(list=ls())
################################################################################
# Setting up, loading R libraries
require(tidyverse)
require(patchwork)  

# Input files 
density_pri <- "densities/borang_pri_centro_100kb.bed"
density_alt <- "densities/borang_alt_centro_100kb.bed"
align_file <- "centromeres/borang_pri/chr10_to_alt.aligned.centromere.txt"
centro_file <-"densities/alt_6sp_centromeres.txt"
windsize=100000

################################################################################
# Processing files and merge data into a single tibble

pri_tib <- density_pri %>% read.table(header=FALSE) %>% as_tibble() %>% 
            rename(NonB=V5, ChrLong=V1, Start=V2, End=V3, Count=V4) %>%
            mutate(hap="Pri")
alt_tib <- density_alt %>% read.table(header=FALSE) %>% as_tibble() %>% 
            rename(NonB=V5, ChrLong=V1, Start=V2, End=V3, Count=V4) %>%
            mutate(hap="Alt")

align_tib <- align_file %>% read.table(header=FALSE) %>% as_tibble() %>% 
            select(-V1,-V4) %>% rename(X4=V2, X3=V3, X1=V5, X2=V6) %>% 
            mutate(ID=row_number()) %>%
            pivot_longer(c(X1, X2, X3, X4), names_to="Xcoord", values_to="X")%>%
            mutate(Y=rep(seq(3,1,-2), each=2, 18))


cen_tib <- centro_file  %>% read.table(header=FALSE) %>% as_tibble() %>%
  rename(Species=V1, ChrLong=V2, Cen_start=V3, Cen_end=V4) %>% 
  filter(Species=="borang")

comb_tib <- pri_tib %>% bind_rows(alt_tib) %>% left_join(cen_tib) %>% 
  filter(NonB!="all") %>% group_by(NonB)  %>% mutate(NonB=if_else(NonB=="GQ", "G4", NonB)) %>%
  mutate(max=max(Count, na.rm=TRUE)) %>% mutate(norm=Count/max) 
  
comb_tib$NonB<-factor(comb_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))

nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45" )
################################################################################
# Plotting
# MAKE ONE PLOT FOR EACH HAPLOTYPE, THEN COMBINE 

### ~~~~~~~~~~~~~~~~~~ Alternative haplotype (from 29 to 38.7 Mb) 
sub_tib1 <- comb_tib %>% filter(hap=="Alt"& End<38800000)

p1<-ggplot(sub_tib1, aes(Start/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  geom_rect(data=sub_tib1 %>% filter(NonB=="APR"), aes(xmin=Cen_start/1000000, ymin=8.55, ymax=8.6, xmax=Cen_end/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_rect(data=sub_tib1 %>% filter(NonB=="APR"), aes(xmin=Cen_start/1000000, ymin=0.4, ymax=0.45, xmax=Cen_end/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent", height=.9, show.legend = TRUE) +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values =nonBcol, name="")+
  labs(x="", y="Alternative chr10") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        panel.grid.major = element_line(colour = 'white'),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.position="top",
        legend.box="horisontal", 
        legend.margin=margin())  +
  guides(fill=guide_legend(nrow=1))+
  scale_x_continuous(expand = c(0,0), limits=c(29,39), position = "top")+
  scale_y_discrete(expand = c(0,0))
p1


### ~~~~~~~~~~~~~~~~~~ Primary haplotype, from 29 to 35.5
sub_tib2 <- comb_tib %>% filter(hap=="Pri" & End<35600000)

p2<-ggplot(sub_tib2, aes(Start/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
   geom_tile(aes(alpha=norm), color="transparent", height=.9, show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1))+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values =nonBcol, name="")+
  labs(x="Mb", y="Primary chr10") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        panel.grid.major = element_line(colour = 'white'),
        plot.margin = margin(t = 15, r = 0, b = 0, l = 0),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())+
  guides(fill=guide_legend(nrow=1))+
  scale_x_continuous(expand = c(0,0), limits=c(29,39), position = "bottom")+
  scale_y_discrete(expand = c(0.1,0))
p2

### ~~~~~~~~~~~~~~~~~~~ Plot alignment 
p3<-ggplot(align_tib, aes(X/1000000, Y)) +
  geom_polygon(aes(group = ID), fill="lightgray", alpha=0.5)+
  labs(x="", y="") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        plot.margin = margin(t = -40, r = 0, b = 0, l = 0),
        axis.ticks =element_blank(),
        axis.text =element_blank())+
  guides(fill=guide_legend(nrow=1))+
  scale_x_continuous(expand = c(0,0), limits=c(29,39), position = "top")+
  scale_y_continuous(expand = c(0,0), limits=c(1,3))
p3



### ~~~~~~~~~~~~~~~~~~~
# JOIN PLOTS 
pcomb <- p1 + p3 + p2 + plot_layout(widths = c(1), heights= c(5,4,5)) 
pcomb
# Save to figure 
outfile<-"plots/FigS19_Bornean_orangutan_chr10.png"
ggsave(outfile, plot = pcomb,dpi = 600,limitsize = TRUE, width=6,height=8)
