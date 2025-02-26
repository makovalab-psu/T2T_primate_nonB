################################################################################
### R code for plotting enrichment in centromeres for alternative haplotypes
### written by Linnéa Smeds Aug-2024

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
setwd("~/Documents/T2T_nonB_paper/")
require(tidyverse)
require(cowplot)
require(ggtext)
plotdir="plots/"
options(scipen = 999, decimals=2)

################################################################################
# Reading in the data, both centromeres and background data for testing significance 
file="centromeres/alt_6sp_enrichment_merged.tsv"
bgfile="centromeres/alt_6sp_background_enrichment_merged.tsv"
gcfile="centromeres/alt_6sp_GCcont_merged.tsv"
bggcfile="centromeres/alt_6sp_background_GCcont_merged.tsv"

# Start with GC files, centro and background 
gctib<-gcfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="centro")
bggctib<-bggcfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Type="background")


# Main file with centromere densities, merge with GC directly
tib<-file %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(head=case_when(Species=="bonobo" ~ "<b>A.</b> Bonobo",
                        Species=="chimp" ~ "<b>B.</b> Chimpanzee",
                        Species=="human" ~ "<b>C.</b> Human",
                        Species=="gorilla" ~ "<b>D.</b> Gorilla",
                        Species=="borang" ~ "<b>E.</b> Bornean orangutan",
                        Species=="sorang" ~ "<b>F.</b> Sumatran orangutan"))%>%
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))
tib <- tib %>% mutate(textcol=if_else(value>0.8*max(tib$value) | value<0.3,"white","black")) %>% 
    mutate(Type="centro") %>% inner_join(gctib)

# Background file, merge with GC directly
bgtib<-bgfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species,-Window), names_to="NonB", values_to="value") %>% 
  mutate(Type="background") %>% inner_join(bggctib) %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4')))


# Use an ordered file to get the right number and order of chromosomes 
len_file="centromeres/pri_6sp_order.tsv"
len_tib <- len_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Species=V1, Chr=V2, Len=V3)

chr_vect <- len_tib$Chr %>% unique()
tib$Chr <- factor(tib$Chr, levels=chr_vect)
bgtib$Chr <- factor(bgtib$Chr, levels=chr_vect)


################################################################################
# CHECK IF CENTROMERE ENRICHMENT IS SIGNIFICANTLY DIFFERENT FROM BACKGROUND 

species_types=tib$Species %>% unique()
nonB_types=tib$NonB %>% unique()
chr_types=tib$Chr %>% unique()

sign_matrix<- matrix(nrow = 0, ncol = 4)
colnames(sign_matrix) <- c("Species", "Chr", "NonB", "Significance")


# Go through each cell and check significance 
for (s in 1:length(species_types)) {
  for (c in 1:length(chr_types)) {
    for (i in 1:length(nonB_types)) {
      set.seed(1234+i)
      background <- bgtib %>% filter(Species==species_types[s] & NonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(value)  %>% pull()
      observed <- tib %>% filter(Species==species_types[s] & NonB==nonB_types[i] & Chr==chr_types[c])  %>% 
        select(value)  %>% pull()
      if(length(background>0) & length(observed>0)) {
          lower_bound <- quantile(background, 0.025)
          upper_bound <- quantile(background, 0.975)
          # Check if the observed value is outside the bounds
          is_significant <- observed < lower_bound | observed > upper_bound
          if(is_significant==TRUE){
            sign=TRUE
          }
          else {
            sign=FALSE
          }
      }
      else {
        sign=NA
      }
      # Add a new row to the matrix
      sign_matrix <- rbind(sign_matrix, c(species_types[s], as.character(chr_types[c]), nonB_types[i], sign) )
    }
  } 
}
sign_tib <- sign_matrix %>% as_tibble()


################################################################################
# ADD SIGNIFICANCE TO THE CENTROMERE TIBBLE
# and add print column with "*" for significant values 

cen_tib<-tib %>% inner_join(sign_tib) %>% 
  mutate(PrintNum=if_else(Significance==TRUE, paste(sprintf("%.2f",value), "*", sep=""), paste(sprintf("%.2f",value)))) %>%
  mutate(PrintType=if_else(Significance==TRUE, "bold", "plain"))
cen_tib$Chr <- factor(cen_tib$Chr, levels=chr_vect)
cen_tib$NonB <- factor(cen_tib$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "TRI", "G4", "Z"))


################################################################################
# PLOT WITH FACET 

p<-ggplot(cen_tib, aes(x=NonB, y=forcats::fct_rev(Chr), fill=value))+
  geom_tile(color="white", show.legend=FALSE) +
  geom_text(aes(label=PrintNum, color=textcol, fontface=PrintType), size=2.5, show.legend=FALSE) +
  facet_wrap(~head, scales="free") +
  labs(x="", y="") +
  scale_fill_gradientn(name = "Enrichment",
                       colors=c("#075AFF", "white","#FF0000"),
                       values=c(0, 1/max(tib$value), 1),
                       breaks=c(-10 ,0, 1),
                       labels = c("0","1","10"),
                         na.value="#075AFF") +
  scale_color_manual(values=c("black","white")) +
  theme(panel.grid.major = element_line(colour = 'white'),
        panel.grid.minor = element_line(colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(0.5, "lines"),
        legend.position="bottom",
        legend.title=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        #legend.position="bottom",
       # legend.box="vertical", 
        legend.margin=margin(c(0,0,0,0)),
        strip.background = element_rect(fill='white'),
        strip.text = element_markdown(hjust=0, size=12))
p

#Fake plot to make legend in log scale 
p2<-ggplot(cen_tib, aes(x=NonB, y=forcats::fct_rev(Chr), fill=value))+
  geom_tile(color="white", alpha=0) +
  scale_fill_gradientn(name = "Fold enrichment",
                       colors=c("#075AFF", "white","#FF0000"),
                       breaks=c(0 ,max(cen_tib$value)/2, 10),
                       labels = c("0","1","10"))+
  theme(legend.position="bottom",
      legend.margin=margin(c(0,0,0,0)),
      legend.key.height = unit(0.4, 'cm'),
      legend.key.width = unit(2, 'cm'),
      legend.title.position = "top",
      legend.title = element_text(hjust=0.5),
      strip.background = element_rect(fill='white'))
p2
leg <- ggpubr::get_legend(p2, position="bottom")

pcomb<-cowplot::plot_grid(p, leg, ncol=1,
          rel_heights=c(1, 0.05))
pcomb
outfile=paste(plotdir,"FigS16_Enrichment_Alt_Centromere.png", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=9,height=12)



# SOME STATS 
# To get the total number of centromeres 
cen_tib %>% group_by(Species, Chr) %>% summarize(n=n())
# 115
# And the ones of them that have at east one non-B DNA >1 
cen_tib %>% group_by(Species, Chr) %>% filter(value>1)
# 81 groups remain
# And with a significant value 
cen_tib %>% group_by(Species, Chr) %>% filter(value>1 & Significance==TRUE)
# 61
# Chromosomes with more than 2-fold enrichment 
cen_tib %>% group_by(Species, Chr) %>% filter(value>2)
# 33 
# and significant:
cen_tib %>% group_by(Species, Chr) %>% filter(value>2 & Significance==TRUE)
# 30
# Number of centromeres with IR enrichments 
cen_tib %>% group_by(Species, Chr) %>% filter(value>1 & NonB=="IR")
#56
cen_tib %>% group_by(Species, Chr) %>% filter(value>1 & NonB=="IR"& Significance==TRUE)
#40


################################################################################
# MAKE TABLE WITH SIGNIFICANCE, INCLUDE SF AND CENPB INFO (FOR TABLE S5) 

SF_file="centromeres/alt_6sp_SF_merged.tsv"
tib_SF<-SF_file %>% read.table(header=TRUE) %>% as_tibble() %>% distinct()%>% mutate(Hap="Alt")
cenpb_file="centromeres/alt_6sp_cenpb.tsv"
tib_cenpb<-cenpb_file %>% read.table(header=TRUE) %>% as_tibble() %>% distinct()

EXPORT <- cen_tib %>%
  select(Species, Chr, NonB, Significance, PrintNum) %>% 
  pivot_wider(id_cols=c(Species,Chr), values_from=PrintNum, names_from=(NonB)) %>% 
  left_join(tib_SF) %>% left_join(tib_cenpb) %>% 
  select(Species, LongChr, SF, Hap, APR, DR, STR,IR,MR,TRI,G4,Z,CENPB)
write.table(EXPORT, row.names = FALSE, quote=FALSE, sep="\t", file = "centromeres/alt_6sp_significance_with_SF_CENPB.tsv")
