################################################################################
### R code for plotting summary of non-B content for all species and chr types 
### written by Linnéa Smeds Nov-2024

################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
require(ggh4x)
plotdir="plots/"
options(scipen = 999, decimals=2)

################################################################################
##### Reading in the data
file<-paste("nonB_annotation/7sp_summary.txt", sep="")

# Convert to tibble 
DATA <- file %>% read.table(header=TRUE) %>% as_tibble() %>%
    mutate(MR=if_else(Type=="Mb", MR-TRI, MR)) %>%
    pivot_longer(c(-Species,-Chr,-Type), names_to="NonB", values_to="value") %>%
    pivot_wider(names_from=Type, values_from=value) %>%  
    mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>% 
    mutate(Head=case_when(Species=="chimp" ~ "Chimpanzee",
                        Species=="bonobo" ~ "Bonobo",
                        Species=="human" ~ "Human",
                        Species=="gorilla" ~ "Gorilla",
                        Species=="borang" ~ "Bornean orangutan",
                        Species=="sorang" ~ "Sumatran orangutan",
                        Species=="siamang" ~ "Siamang")) %>%
  group_by(Species,Chr) %>% mutate(tot=frac[NonB=="all"], place=sum(Mb[NonB!="all"])) %>%
  mutate(PrintNum=paste(sprintf("%.1f",tot*100),"%",sep="")) 

DATA$Head <- factor(DATA$Head, levels=c("Bonobo","Chimpanzee", "Human", "Gorilla", "Bornean orangutan", "Sumatran orangutan", "Siamang"))
DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))

nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","darkorange","#DB5829","#894B45")
################################################################################
# PLOT NON-B CONTENT AS HORISONTAL BARPLOT
SUB<-DATA %>% filter(NonB!="all") %>% mutate(PrintNum=if_else(NonB=="APR", PrintNum, NA)) 

p<-ggplot(SUB) + 
  geom_col(aes(Mb, fct_rev(Head), fill=fct_rev(NonB), color=fct_rev(NonB)), alpha=0.6, width=0.6)+
  geom_text(aes(label=PrintNum, x=place, y=fct_rev(Head)), hjust=-0.1)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"), values=nonBcol)+
  facet_grid(~Chr, scales="free") +
  force_panelsizes(cols = c(4, 1,1))+
  coord_cartesian(clip = "off") + #add this to print the full numbers
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.text=element_text(hjust=0,size=12, face='bold'),
        panel.spacing = unit(2.5, "lines"),
        legend.position="top",
        legend.justification = "right",
        legend.text  = element_text(size = 10),
        legend.title=element_blank(),
        plot.margin = margin(t = 20, r = 50, b = 40, l = 10),
        axis.text = element_text(margin = margin(r = 0, b=0)))+
  guides(fill = guide_legend(nrow = 1))+
  xlab("Mb") + 
  ylab("")+
  scale_x_continuous(expand = c(0, 0))
p

outfile=paste(plotdir,"Fig2.tiff", sep="")
ggsave(outfile,plot = p,scale = 1, device = "tiff", dpi = 600,limitsize = TRUE,width=10,height=6, compression = "lzw")
