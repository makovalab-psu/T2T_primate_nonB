################################################################################
### R code for plotting enrichment in Repeats for non-human apes  
### written by Linnéa Smeds 30-Aug-2024
################################################################################
# Setting up, loading R libraries
rm(list=ls())
require(tidyverse)
require(patchwork) 
require(ggtext)

repeat_file <- "repeats/6sp_genome_repeat_enrichment_renamed.tsv"
lengths_file <-"repeats/6sp_genome_repeat_lengths_renamed.txt"
outprefix <-  "plots/FigS13_6sp_enrichment_repeats"
threshold <- 50000    # Example: 50000


################################################################################
# Processing files

# Read files and pivot table to put motifs in single column 
lentib=lengths_file %>% read.table(header=TRUE)
tib <- repeat_file %>% read.table(header=TRUE) %>% as_tibble() %>% 
      left_join(lentib, by = c("Repeat" = "Repeat", "Species" = "Species")) %>% 
        pivot_longer(-c(Repeat,Len,Species), names_to="Motif") %>%
      mutate(PrintName=paste(Repeat," (",sprintf("%.2f",Len/1000000),"Mb)", sep="")) %>%
      mutate(PrintName = sub("SAT", "Sat", PrintName)) %>%
      drop_na() %>% mutate(head=case_when(Species=="bonobo" ~ "<b>A.</b> Bonobo",
                                          Species=="chimp" ~ "<b>B.</b> Chimpanzee",
                                          Species=="gorilla" ~ "<b>C.</b> Gorilla",
                                          Species=="borang" ~ "<b>D.</b> Bornean orangutan",
                                          Species=="sorang" ~ "<b>E.</b> Sumatran orangutan",
                                          Species=="siamang" ~ "<b>F.</b> Siamang"),
                           Motif=if_else(Motif=="GQ", "G4", Motif),
                           logval=log2(value), 
                           Class=case_when(str_detect(Repeat, "DNA") ~ "TEs",
                                          str_detect(Repeat, "LINE") ~ "TEs",
                                          str_detect(Repeat, "SINE") ~ "TEs",
                                          str_detect(Repeat, "LTR") ~ "TEs",
                                          str_detect(Repeat, "Helitron") ~ "TEs",
                                          str_detect(Repeat, "Retroposon") ~ "TEs",
                                          str_detect(Repeat, "Satellite") ~ "Satellites",
                                          str_detect(Repeat, "RNA") ~ "RNA",
                                          TRUE ~ "Other"))
tib$Motif <- factor(tib$Motif, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
tib$Class <- factor(tib$Class, levels=c("TEs", "Satellites", "RNA", "Other"))

 
################################################################################
# PLOTTING HORISONTAL AS TWO PLOTS TO GO ON TWO PAGES 
# Removing lines that are all NA and repeats with less than threshold bp
######  PART 1
species=c("bonobo", "chimp", "gorilla", "borang", "sorang", "siamang")
legvect<-c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
plots <- list()

scalestart=min(tib %>% filter(value>0) %>%select(logval))
scaleend=max(tib$logval)
scalebreak=-1*scalestart/(scaleend-scalestart)


for(i in 1:length(species)){
  DATA <- tib %>% filter(Len>threshold) %>% filter(Species==species[i])
  DATA$Motif <- factor(DATA$Motif, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))
  DATA$Class <- factor(DATA$Class, levels=c("TEs", "Satellites", "RNA", "Other"))
 
  p<-ggplot(DATA, aes(PrintName, fct_rev(Motif))) +
    geom_tile(aes(fill=logval), color = "white",lwd = 0.3,linetype = 1, show.legend = legvect[i]) +
    facet_grid(~Class, space="free", scales="free_x") +
    scale_fill_gradientn(name="log2(Enrichment)",
                         colors=c("#075AFF", "white","#FF0000"),
                         values=c(0, scalebreak, 1), na.value="#075AFF")+
    labs(x="", y="") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1, size=10),
          axis.text.y = element_text(size=12),
          plot.title=element_markdown(hjust=0, size=18),
          legend.text=element_text(size=12),
          legend.position="bottom",
          legend.title=element_text(size=12, vjust=1),
          legend.key.width  = unit(2, "lines"),
          legend.key.height = unit(0.7, "lines"),
          legend.margin = margin(t=0,r=0,b=0,l=0, unit="pt"),
          strip.text = element_text(colour = 'black', size=11, face="bold"),
         # strip.background = element_rect(fill = 'white', colour="white"),
          panel.background = element_rect(fill = 'white', colour="white"))+
    ggtitle(unique(DATA$head))+
    scale_y_discrete(position="left")
  p
  plots[[i]] <- p
}

pcomb1<-plots[[1]] + plots[[2]] + plots[[3]] + plot_layout(widths = c(1), heights= c(1,1,1)) 
pcomb2<-plots[[4]] + plots[[5]] + plots[[6]] + plot_layout(widths = c(1), heights= c(1,1,1)) 

outfile1=paste(outprefix,"_partABC.png", sep="")
ggsave(outfile1,plot = pcomb1,scale = 1,dpi = 600,limitsize = TRUE,width=12,height=14)
outfile2=paste(outprefix,"_partDEF.png", sep="")
ggsave(outfile2,plot = pcomb2,scale = 1,dpi = 600,limitsize = TRUE,width=12,height=14)


