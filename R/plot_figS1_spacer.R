################################################################################
### R code for plotting SPACER DISTRIBUTION

# Setting up, loading R libraries 
rm(list=ls())
require(tidyverse)
plotdir="plots/"
options(scipen = 999)

# Reading in the data
IRfile=paste("spacer/IR.spacers.txt", sep="")
IRtib<-IRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="A", Length=as.numeric(Length)) 
MRfile=paste("spacer/MR.spacers.txt", sep="")
MRtib<-MRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="B", Length=as.numeric(Length))

# Combine IR and MR
combtib<-IRtib %>% full_join(MRtib)


################################################################################
# PLOT (combine with facet)

outfile=paste(plotdir,"FigureS1_Spacer_human.jpg", sep="")

p<-ggplot(combtib, aes(x=Length)) +
  geom_histogram(binwidth=1,boundary=0, color='#1964B0', fill='#1964B0') +
  facet_wrap(~Type, scales="free") +
  labs(x="Spacer length")+
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.grid.minor = element_line(colour = 'white'),
        panel.background = element_rect(fill = '#edede6', colour = '#edede6'),
        panel.spacing = unit(2, "lines"),
        legend.position="bottom",
        legend.title=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(colour = 'black', hjust=0, size=12))
p

ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=6,height=3)
