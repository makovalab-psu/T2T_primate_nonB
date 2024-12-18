#################
### R code for plotting SPACER DISTRIBUTION
### Written by Linnéa Smeds August 2024


# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
plotdir="plots/"
options(scipen = 999)

# Reading in the data, convert to tibble
DRfile=paste("spacer/DR.spacers.txt", sep="")
DRtib<-DRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="A", Length=as.numeric(Length)) 
IRfile=paste("spacer/IR.spacers.txt", sep="")
IRtib<-IRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="B", Length=as.numeric(Length)) 
MRfile=paste("spacer/MR.spacers.txt", sep="")
MRtib<-MRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="C", Length=as.numeric(Length))

combtib<-DRtib%>% bind_rows(IRtib,MRtib)

# Check % of spacers has 0 length
IRfrac<-length(IRtib$Length[IRtib==0])/length(IRtib$Length)
DRfrac<-length(DRtib$Length[DRtib==0])/length(DRtib$Length)
MRfrac<-length(MRtib$Length[MRtib==0])/length(MRtib$Length)

################################################################################
# PLOT (combine with facet)

p<-ggplot(combtib, aes(x=Length)) +
  geom_histogram(binwidth=1, position = "identity", color='#1964B0', fill='#1964B0') +
  facet_wrap(~Type, scales="free") +
  labs(x="Spacer length", y="Count")+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(2, "lines"),
        legend.position="bottom",
        legend.title=element_blank(),
        #axis.ticks.x=element_blank(), 
        #axis.ticks.y=element_blank(),
        axis.text.y=element_text(angle=55),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(colour = 'black', hjust=0, size=14, face="bold"))
p

outfile=paste(plotdir,"FigS1_Spacer_human.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=3)
