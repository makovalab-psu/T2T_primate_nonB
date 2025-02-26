################################################################################
### R code for plotting density of methylation scores and Quadron scores for 
# different G4s in LSAU
### written by Linnéa Smeds 21-Nov-2024

################################################################################
# Setting up, loading R libraries 
rm(list=ls())
require(tidyverse)
require(patchwork)
plotdir="plots/"
options(scipen = 999, decimals=1)

##### Reading in the data
sp="human"

fwdfile=paste("methylation/",sp,"/cluster/LSAU.sepchr.fwd.G4s.meth.txt", sep="")
revfile=paste("methylation/",sp,"/cluster/LSAU.sepchr.rev.G4s.meth.txt", sep="")

# Convert to tibble 
fwd_tib<-fwdfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  select(Type, Cluster, HG002_mean, CHM13_mean, G4_score) %>%
  mutate(Strand='forward')
rev_tib<-revfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  select(Type, Cluster, HG002_mean, CHM13_mean, G4_score) %>% 
  mutate(Strand='reverse')

# Merge forward and reverse, and add markings to chosen types
comb_tib <-fwd_tib %>% bind_rows(rev_tib) %>%
  rename(HG002=HG002_mean, CHM13=CHM13_mean) %>%
  pivot_longer(c(HG002,CHM13), names_to="Celltype", values_to="Mean") %>%
  mutate(Type = str_replace_all(Type, c('excl4-10' = 'excl4+10'))) %>%
  mutate(Group=paste(Strand,Type,Cluster, sep="_"), ) %>%
  mutate(Chosen=case_when(Group=="reverse_chr10_1" ~ "LS 1",
                          Group=="reverse_chr10_2" | Group=="reverse_chr4_3" ~ "LS 2",
                          Group=="reverse_excl4+10_2" ~ "LS 3",
                          TRUE ~ NA))

# Check numbers of G4 per cluster (no point plotting clusters with few sequences)
selected<-comb_tib %>% group_by(Group, Celltype) %>% summarize(n=n()) %>%
  filter(n>20) %>%select(Group) %>% unique()

################################################################################
# PLOT (combine with facet)

####### METHYLATION DISTRIBUTIONS 
DATA<-comb_tib %>% inner_join(selected)
TEXT <- DATA %>% select(Group, Chosen, Celltype) %>% 
  filter(Celltype=="HG002") %>% unique()


p1<-ggplot(DATA, aes(x=Mean, fill=Celltype, color=Celltype)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("#882D71", "#008A69"))+
  scale_color_manual(values=c("#882D71", "#008A69"))+
  geom_text(data=TEXT, aes(x=Inf,y=Inf,hjust=1,vjust=1.2,label=Chosen), color="red", size=8)+
  facet_wrap(~Group, scales="free_y") +
  xlab("Mean methylation score")+
  ylab("Density")+
  ggtitle("A")+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.background = element_rect(fill='lightgray'),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title=element_text(size=11),
        plot.title = element_text(size=18, face="bold"),
        strip.text = element_text(colour = 'black', size=11, face="bold"))
p1

####### QUADRON SCORE DISTRIBUTION 
# Extract Quadron score (use just one celltype, otherwise the scores are duplicated)
DATA2<-comb_tib %>% inner_join(selected) %>% filter(Celltype=="HG002") %>%
  select(Group, G4_score)

p2<-ggplot(DATA2, aes(x=G4_score)) + 
  geom_density(alpha=0.5, fill="black", color="black") +
  facet_wrap(~Group, scales="free_y") +
  geom_text(data=TEXT, aes(x=-Inf,y=Inf,hjust=0,vjust=1.2,label=Chosen), color="red", size=8)+
  xlab("Quadron score")+
  ylab("Density")+
  ggtitle("B")+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.background = element_rect(fill='lightgray'),
        axis.title=element_text(size=11),
        plot.title = element_text(size=18, face="bold"),
        strip.text = element_text(colour = 'black', size=11, face="bold"))
p2


####### COMBINE PLOTS 
pcomb <- p1 + p2 + plot_layout(heights = c(1, 1)) 
pcomb
outfile1=paste(plotdir,"FigS2_LSAU_G4_methylation.png", sep="")
ggsave(outfile1,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=12)

