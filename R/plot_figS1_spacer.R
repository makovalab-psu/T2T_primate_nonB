#################
### R code for plotting SPACER DISTRIBUTION
### Written by Linnéa Smeds August 2024


# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork) 
plotdir="plots/"
options(scipen = 999)

# Reading in the data, convert to tibble
DRfile=paste("spacer/DR.spacers.txt", sep="")
DRtib<-DRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="A", NonB="DR", Length=as.numeric(Length)) 
IRfile=paste("spacer/IR.spacers.txt", sep="")
IRtib<-IRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="B", NonB="IR", Length=as.numeric(Length)) 
MRfile=paste("spacer/MR.spacers.txt", sep="")
MRtib<-MRfile %>% read.table(header=FALSE) %>% as_tibble() %>% 
  rename(Length=V1) %>% mutate(Type="C", NonB="MR", Length=as.numeric(Length))

combtib<-DRtib%>% bind_rows(IRtib,MRtib)

# Check % of spacers has 0 length
IRfrac<-length(IRtib$Length[IRtib==0])/length(IRtib$Length)
DRfrac<-length(DRtib$Length[DRtib==0])/length(DRtib$Length)
MRfrac<-length(MRtib$Length[MRtib==0])/length(MRtib$Length)

# Cumulative values 
cumtib <- combtib %>% group_by(Type, Length) %>% summarize(n=n()) %>% 
  mutate(cum=cumsum(n)) %>% mutate(frac=cum/sum(n)) %>% mutate(max=max(n))

# Max for plotting 
maxtib<-cumtib %>% select(max,Type) %>% group_by(max,Type) %>% summarize()

nonBcol=c("#7BB0DF","#E9DC6D","#F4A637")

################################################################################
# PLOT (combine with facet)

p1<-ggplot(combtib, aes(x=Length)) +
  geom_histogram(binwidth=1, position = "identity", color='#1964B0', fill='#1964B0') +
  facet_wrap(~Type, scales="free") +
  labs(x="", y="Count")+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
     #   panel.spacing = unit(2, "lines"),
        plot.margin=margin(t=1,b=0,l=1,r=1, unit="lines"),
        legend.position="bottom",
        legend.title=element_blank(),
        axis.text.y=element_text(angle=55),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(colour = 'black', hjust=0, size=14, face="bold"))
p1

  
# MAKE A CUMULATIVE LINE PLOT
p2<-ggplot(cumtib, aes(x=Length, y=frac)) +
  geom_line(color='red') +
  facet_wrap(~Type, scales="free") +
  labs(x="Spacer length", y="Cumulative fraction")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
      #  panel.spacing = unit(1.80, "lines"),
        plot.margin=margin(t=0,b=1,l=1,r=1, unit="lines"),
        legend.position="bottom",
        legend.title=element_blank(),
        #axis.text.y=element_text(angle=55),
        strip.background = element_rect(fill='white'),
        strip.text = element_blank())
p2

pcomb<-p1 +  plot_spacer() + p2 + plot_layout(widths = c(1), heights= c(5,-0.2,5)) 
pcomb



h <- hist(x, plot = FALSE)
plot(h, col = "red", ylim = c(0, sum(h$counts)))
lines(h$mids, cumsum(h$counts))

outfile=paste(plotdir,"FigS1_Spacer_human.png", sep="")
ggsave(outfile,plot = p,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=3)

################################################################################
# CODE FOR ADDING BOTH HISTOGRAM AND LINE TO THE SAME PLOT USING LOOP
# This works for each part individually, but when combining them, the second 
# axes change so they are all transformed after the last nonB-type (MR). 

loop=c("DR", "IR", "MR")
name=c("A", "B", "C")
maxvect=c(maxtib$max[maxtib$Type=="A"],maxtib$max[maxtib$Type=="B"],maxtib$max[maxtib$Type=="C"])
plots <- list()

for(i in 1:length(loop)){
  DATA <- combtib %>% filter(NonB==loop[i])
  CUMDATA <- DATA %>% group_by(Type, Length) %>% summarize(n=n()) %>% mutate(cum=cumsum(n)) %>% mutate(frac=cum/sum(n))
  p<-ggplot(DATA, aes(x=Length)) +
    geom_histogram(binwidth=1, position = "identity", color='#1964B0', fill='#1964B0', show.legend = legvect[i]) +
    labs(x="Spacer length")+
    geom_line(data=CUMDATA, aes(x=Length,y=frac*maxvect[i]), color="red", linewidth=1)+
    scale_y_continuous(name = "Count", limits = c(0,maxvect[i]),
                       sec.axis = sec_axis(transform=~./maxvect[i], name="Cumulative fraction")) +
    theme(panel.grid.major = element_line(colour = 'beige'),
          panel.grid.minor = element_line(colour = 'beige'),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(2, "lines"),
          legend.position="bottom",
          plot.title=element_text(hjust=0, size=14, face="bold"),
          #axis.ticks.x=element_blank(), 
          #axis.ticks.y=element_blank(),
          axis.text.y.left=element_text(angle=55))+
    #       strip.background = element_rect(fill='white'),
    #      strip.text = element_text(colour = 'black', hjust=0, size=14, face="bold"))+
    ggtitle(name[i])
  p
  plots[[i]] <- p
}

pcomb1<-plots[[1]] + plots[[2]] + plots[[3]] + plot_layout(widths = c(1,1,1)) 
pcomb1

################################################################################
# CODE FOR ADDING BOTH HISTOGRAM AND LINE TO THE SAME PLOT, HARDCODED 
# THIS WORKS, BUT ONLY WHEN NAMING THE MAX VARIABLES DIFFERENTLY FOR EACH PLOT
# (when they are all named max, p1 and p2 changes after p3 is drawn)

# ==== THIS PART IS WHAT IS PRESENTLY USED FOR THE SUPPLEMENT FIGURE ====

loop=c("DR", "IR", "MR")
name=c("A", "B", "C")
maxvect=c(maxtib$max[maxtib$Type=="A"],maxtib$max[maxtib$Type=="B"],maxtib$max[maxtib$Type=="C"])

# DR
DATA <- combtib %>% filter(NonB==loop[1])
CUMDATA <- DATA %>% group_by(Type, Length) %>% summarize(n=n()) %>% mutate(cum=cumsum(n)) %>% mutate(frac=cum/sum(n))
max1=maxvect[1]
p1<-ggplot(DATA, aes(x=Length)) +
  geom_histogram(binwidth=1, position = "identity", color='#7BB0DF', fill='#7BB0DF') +
  labs(x="Spacer length")+
  geom_line(data=CUMDATA, aes(x=Length,y=frac*max1), color="red", linewidth=0.5, linetype=2)+
  scale_y_continuous(name = "Count", limits = c(-1,max1),
                     sec.axis = sec_axis(transform=~./max1, name="")) +
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title=element_text(hjust=0, size=14, face="bold"),
        plot.margin=margin(t=0,b=0,l=0,r=0, unit="lines"),
        axis.title.y.right=element_blank(),
        axis.text.y.left=element_text(angle=55))+
  ggtitle(name[1])
p1

# IR
DATA <- combtib %>% filter(NonB==loop[2])
CUMDATA <- DATA %>% group_by(Type, Length) %>% summarize(n=n()) %>% mutate(cum=cumsum(n)) %>% mutate(frac=cum/sum(n))
max2=maxvect[2]
p2<-ggplot(DATA, aes(x=Length)) +
  geom_histogram(binwidth=1, position = "identity", color='#7BB0DF', fill='#7BB0DF') +
  labs(x="Spacer length")+
  geom_line(data=CUMDATA, aes(x=Length,y=frac*max2), color="red", linewidth=0.5, linetype=2)+
  scale_y_continuous(name = "", limits = c(-1,max2),
                     sec.axis = sec_axis(transform=~./max2, name="")) +
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title=element_text(hjust=0, size=14, face="bold"),
        plot.margin=margin(t=0,b=0,l=0,r=0, unit="lines"),
        axis.title.y=element_blank(),
        axis.text.y.left=element_text(angle=55))+
  ggtitle(name[2])
p2

# MR
DATA <- combtib %>% filter(NonB==loop[3])
CUMDATA <- DATA %>% group_by(Type, Length) %>% summarize(n=n()) %>% mutate(cum=cumsum(n)) %>% mutate(frac=cum/sum(n))
max3=maxvect[3]
p3<-ggplot(DATA, aes(x=Length)) +
  geom_histogram(binwidth=1, position = "identity", color='#7BB0DF', fill='#7BB0DF') +
  labs(x="Spacer length")+
  geom_line(data=CUMDATA, aes(x=Length,y=frac*max3), color="red", linewidth=0.5, linetype=2)+
  scale_y_continuous(name = "", limits = c(-1,max3),
                     sec.axis = sec_axis(transform=~./max3, name="Cumulative fraction")) +
  theme(panel.grid.major = element_line(colour = 'beige'),
        panel.grid.minor = element_line(colour = 'beige'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title=element_text(hjust=0, size=14, face="bold"),
        plot.margin=margin(t=0,b=0,l=0,r=0, unit="lines"),
        axis.title.y.left=element_blank(),
        axis.text.y.left=element_text(angle=55))+
  ggtitle(name[3])
p3

pcomb<-p1 + plot_spacer() + p2 + plot_spacer() + p3 + plot_layout(widths = c(1,-0.05,1,-0.05,1)) 
pcomb

outfile=paste(plotdir,"FigS1_Spacer_human.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=7,height=3)



