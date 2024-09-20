#!/usr/bin/env Rscript
#################
### written by Linnéa Smeds 29-Aug-2024

rm(list=ls())
################################################################################
# Setting up, loading R libraries
require(tidyverse)
require(patchwork)  
#install.packages("ggh4x")
#require(ggh4x)

# Input files 
density_file <- "densities/pri_6sp_merged.txt"
lengths_file <- "densities/pri_6sp_lengths.txt"
centro_file <-"densities/pri_6sp_centromeres.txt"
xlimit <- 250000000 #(works well for human, but not for the others..)
outfile <- "plots/Figure3_primary_densities.jpg"
windsize=100000

################################################################################
# Processing files and merge data into a single tibble

len_tib <- lengths_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Species=V1, ChrLong=V2, Len=V3, Chr=V4)
dens_tib <- density_file  %>% read.table(header=FALSE) %>% as_tibble() %>% 
            rename(Species=V1, NonB=V6, ChrLong=V2, Start=V3, End=V4, Count=V5) %>%
            group_by(NonB, Species) %>% mutate(max=max(Count, na.rm=TRUE))
cen_tib <- centro_file  %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Species=1, ChrLong=V2, Cen_start=V3, Cen_end=V4)
comb_tib <- dens_tib %>% left_join(len_tib) %>% left_join(cen_tib) %>% 
  filter(NonB!="all" & Chr!="chr1522") %>%  mutate(norm=Count/max) %>% 
  mutate(maxlen=case_when(Species=="human" ~ 250000000,
                          Species=="chimp" ~ 260000000,
                          Species=="borang" ~ 260000000,
                          Species=="gorilla" ~ 300000000,
                          TRUE ~ 280000000))
comb_tib$Species <- factor(comb_tib$Species, levels=c("chimp", "bonobo", "human", "gorilla", "borang", "sorang"))

merged_tib<-comb_tib %>% group_by(Chr, Species) %>% mutate(maxstart=Start[which.max(Start)]) %>% ungroup() %>%
  mutate(Newstart=case_when((Chr=="chr24" | Chr=="chr23" | Chr=="chr22" | Chr=="chr21"| Chr=="chr20" | Chr=="chr19" | 
                               Chr=="chr18" | Chr=="chr17" | Chr=="chr16" | 
                               Chr=="chr15" | Chr=="chrY") ~ maxlen-maxstart+Start,
                            TRUE ~ Start), 
         Newend=case_when((Chr=="chr24" | Chr=="chr23" | Chr=="chr22" | Chr=="chr21"| Chr=="chr20" | Chr=="chr19" | 
                             Chr=="chr18" | Chr=="chr17" | Chr=="chr16" | 
                             Chr=="chr15" | Chr=="chrY") ~ maxlen-maxstart+End,
                          TRUE ~ End),
         Newcenstart=case_when((Chr=="chr24" | Chr=="chr23" | Chr=="chr22" | Chr=="chr21"| Chr=="chr20" | Chr=="chr19" | 
                                  Chr=="chr18" | Chr=="chr17" | Chr=="chr16" | 
                                  Chr=="chr15" | Chr=="chrY") ~ maxlen-maxstart+Cen_start,
                               TRUE ~ Cen_start),
         Newcenend=case_when((Chr=="chr24" | Chr=="chr23" | Chr=="chr22" | Chr=="chr21"| Chr=="chr20" | Chr=="chr19" | 
                                Chr=="chr18" | Chr=="chr17" | Chr=="chr16" | 
                                Chr=="chr15" | Chr=="chrY") ~ maxlen-maxstart+Cen_end,
                             TRUE ~ Cen_end),
         Righttext=case_when(((Chr=="chr24" | Chr=="chr23" | Chr=="chr22" | Chr=="chr21"| Chr=="chr20" | Chr=="chr19" | 
                                 Chr=="chr18" | Chr=="chr17" | Chr=="chr16" | 
                                 Chr=="chr15" | Chr=="chrY") & Newstart==maxlen-maxstart & NonB=="IR") ~ as.character(Chr),
                             TRUE ~ NA)) %>%
  mutate(Chr=case_when(Chr=="chr24" ~ 'chr5',
                       Chr=="chr23" ~ 'chr6',
                       Chr=="chr22" ~ 'chr7', 
                       Chr=="chr21" ~ 'chr8', 
                       Chr=="chr20" ~ 'chr9', 
                       Chr=="chr19" ~ 'chr10', 
                       Chr=="chr18" ~ 'chr11', 
                       Chr=="chr17" ~ 'chr12', 
                       Chr=="chr16" ~ 'chr13', 
                       Chr=="chr15" ~ 'chr14', 
                       Chr=="chrY" ~ 'chrX', 
                       TRUE ~ Chr)) 


################################################################################
# Plotting
# MAKE ONE PLOT FOR EACH SP, THEN COMBINE 

### ~~~~~~~~~~~~~~~~~~ Chimp
sub_tib1 <- merged_tib %>% filter(Species=="chimp")
chr_vect <- sub_tib1$Chr %>% unique()
sub_tib1$Chr <- factor(sub_tib1$Chr, levels=chr_vect)
maxl=sub_tib1$maxlen[1]

p1<-ggplot(sub_tib1, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  geom_rect(data=sub_tib1 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent", show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1))+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"))+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  ggtitle("A")+
  scale_x_continuous(expand = c(0,1))
p1

### ~~~~~~~~~~~~~~~~~~ Bonobo
sub_tib2 <- merged_tib %>% filter(Species=="bonobo")
chr_vect <- sub_tib2$Chr %>% unique()
sub_tib2$Chr <- factor(sub_tib2$Chr, levels=chr_vect)
maxl=sub_tib2$maxlen[1]

p2<-ggplot(sub_tib2, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  # geom_segment(data=merged_tib, aes(x=Newcenstart/1000000, yend=forcats::fct_rev(NonB), xend=Newcenend/1000000), color="red", alpha=0.5, linewidth=15) +
  geom_rect(data=sub_tib2 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent",  show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"), name="")+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  ggtitle("B")+
  scale_x_continuous(expand = c(0,1))
p2

### ~~~~~~~~~~~~~~~~~~ Human
sub_tib3 <- merged_tib %>% filter(Species=="human")
chr_vect <- sub_tib3$Chr %>% unique()
sub_tib3$Chr <- factor(sub_tib3$Chr, levels=chr_vect)
maxl=sub_tib3$maxlen[1]

p3<-ggplot(sub_tib3, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  # geom_segment(data=merged_tib, aes(x=Newcenstart/1000000, yend=forcats::fct_rev(NonB), xend=Newcenend/1000000), color="red", alpha=0.5, linewidth=15) +
  geom_rect(data=sub_tib3 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent",  show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"), name="")+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  ggtitle("C")+
  scale_x_continuous(expand = c(0,1))
p3

### ~~~~~~~~~~~~~~~~~~ Gorilla
sub_tib4 <- merged_tib %>% filter(Species=="gorilla")
chr_vect <- sub_tib4$Chr %>% unique()
sub_tib4$Chr <- factor(sub_tib4$Chr, levels=chr_vect)
maxl=sub_tib4$maxlen[1]

p4<-ggplot(sub_tib4, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  geom_rect(data=sub_tib4 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent",  show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"), name="")+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  ggtitle("D")+
  scale_x_continuous(expand = c(0,1))
p4

### ~~~~~~~~~~~~~~~~~~ Borang
sub_tib5 <- merged_tib %>% filter(Species=="borang")
chr_vect <- sub_tib5$Chr %>% unique()
sub_tib5$Chr <- factor(sub_tib5$Chr, levels=chr_vect)
maxl=sub_tib5$maxlen[1]

p5<-ggplot(sub_tib5, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB))) +
  geom_rect(data=sub_tib5 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent") +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"), name="")+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        legend.position="bottom",
        legend.box="vertical", 
        legend.margin=margin(),
         axis.text.y =element_blank())  +
  ggtitle("E")+
  scale_x_continuous(expand = c(0,1))+
  guides(fill=guide_legend(nrow=1))
p5


### ~~~~~~~~~~~~~~~~~~ Sorang
sub_tib6 <- merged_tib %>% filter(Species=="sorang")
chr_vect <- sub_tib6$Chr %>% unique()
sub_tib6$Chr <- factor(sub_tib6$Chr, levels=chr_vect)
maxl=sub_tib6$maxlen[1]

p6<-ggplot(sub_tib6, aes(Newstart/1000000, fct_rev(NonB), fill=fct_rev(NonB)), ) +
  geom_rect(data=sub_tib6 %>% filter(NonB=="APR"), aes(xmin=Newcenstart/1000000, ymin=8, ymax=8, xmax=Newcenend/1000000), color="red", linewidth=3, inherit.aes = FALSE) +
  geom_tile(aes(alpha=norm), color="transparent",  show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.1, 1), name="Norm. non-B density")+
  facet_grid(Chr~., scales="free",  switch='y')+
  geom_text(aes(label=Righttext), y="IR", size=4, angle=90, hjust=0.5, vjust=-1)+
  scale_fill_manual(breaks=c("APR", "DR", "GQ", "IR", "MR", "STR", "Z"), values =c("#DB5829","#7BB0DF", "#894B45", "#008A69", "#882D71", "#1964B0", "#386350"), name="")+
  labs(x="", y="") +
  xlim(0, maxl/1000000) +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'white', colour="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = 'white'),
        panel.spacing.y=unit(0.3, "lines"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.ticks.y =element_blank(),
        axis.text.y =element_blank())  +
  ggtitle("F")+
  scale_x_continuous(expand = c(0,1))
p6



### ~~~~~~~~~~~~~~~~~~ Combine all 6 plots 

pcomb <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(widths = c(1, 1, 1), heights= c(1,1)) 
pcomb
# Save to figure 
ggsave(outfile, plot = pcomb,dpi = 600,limitsize = TRUE, width=15,height=20)


