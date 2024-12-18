################################################################################
### R code for plotting enrichment in centromeres with and without CENP-B motifs.
### written by Linnéa Smeds 10-Dec-2024

################################################################################
# Setting up, loading R libraries and set working directory
rm(list=ls())
require(tidyverse)
require(patchwork)  
plotdir="plots/"


# Input files 
epfile="centromeres/pri_6sp_enrichment_merged.tsv"
eafile="centromeres/alt_6sp_enrichment_merged.tsv"
cpfile="centromeres/pri_6sp_cenpb.tsv"
cafile="centromeres/alt_6sp_cenpb.tsv"
ufpfile="centromeres/pri_6sp_upstream1Mb_enrichment_merged.tsv"
ufafile="centromeres/alt_6sp_upstream1Mb_enrichment_merged.tsv"
dfpfile="centromeres/pri_6sp_downstream1Mb_enrichment_merged.tsv"
dfafile="centromeres/alt_6sp_downstream1Mb_enrichment_merged.tsv"


# Reading in the data 
# Centromeres, fold enrichment
enrichp_tib<-epfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
      pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
      mutate(Hap="pri")
enricha_tib<-eafile %>% read.table(header=TRUE) %>% as_tibble() %>% 
      pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
      mutate(Hap="alt")
# Cenp-b annotation
cenpp_tib<-cpfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
      mutate(Hap="pri")
cenpa_tib<-cafile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  mutate(Hap="alt")

# Upstream 1kb
ufp_tib<-ufpfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(Hap="pri")
ufa_tib<-ufafile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(Hap="alt")
# Downstream 1kb
dfp_tib<-dfpfile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(Hap="pri")
dfa_tib<-dfafile %>% read.table(header=TRUE) %>% as_tibble() %>% 
  pivot_longer(c(-Chr,-Species), names_to="NonB", values_to="value") %>%
  mutate(Hap="alt")

DATA_FOLD <- bind_rows(enrichp_tib,enricha_tib) %>% 
  inner_join(bind_rows(cenpp_tib,cenpa_tib), relationship = "many-to-many") %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>% mutate(Type="B")

DATA_UP1Mb<-bind_rows(ufp_tib,ufa_tib) %>% 
  inner_join(bind_rows(cenpp_tib,cenpa_tib), relationship = "many-to-many") %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>% mutate(Type="A")

DATA_DOWN1Mb<-bind_rows(dfp_tib,dfa_tib) %>% 
  inner_join(bind_rows(cenpp_tib,cenpa_tib), relationship = "many-to-many") %>% 
  mutate(NonB = str_replace_all(NonB, c('GQ' = 'G4'))) %>% mutate(Type="C")

# And combined 
DATA <- DATA_FOLD %>% bind_rows(DATA_UP1Mb, DATA_DOWN1Mb)


nonBcol=c("#AE75A2","#7BB0DF","#008A69","#E9DC6D","#F4A637","#DB5829","#894B45")

################################################################################
# PLOT ALL TOGETHER 

DATA$NonB <- factor(DATA$NonB, levels=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"))
DATA$CENPB <- factor(DATA$CENPB, levels=c("Yes","No"))
max<-max(DATA$value)

SUB<-DATA%>%filter(Type=="A")
pA<-ggplot(SUB, aes(x=CENPB, y=value, fill=NonB, color=NonB))+
  geom_boxplot(show.legend=FALSE, alpha=0.6)+
  facet_wrap(~NonB, nrow=1)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=1)+
  labs(x="", y="Fold-enrichment") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'lightgray', colour="white"),
        strip.text = element_blank(),
        panel.grid.major = element_line(colour = 'beige'),
         panel.spacing.y=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.position="top",
        legend.justification = "right",
        legend.title = element_blank(),
        plot.title = element_text(size=18),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))+
  labs(title = bold('A') ~ '')+
  ylim(c(0,max))
pA

SUB<-DATA%>%filter(Type=="B")
pB<-ggplot(SUB, aes(x=CENPB, y=value, fill=NonB, color=NonB))+
  geom_boxplot(show.legend=TRUE, alpha=0.6)+
  facet_wrap(~NonB, nrow=1)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=1)+
  labs(x="CENP-B", y="") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'lightgray', colour="white"),
        strip.text = element_blank(),
        panel.grid.major = element_line(colour = 'beige'),
        panel.spacing.y=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.position="bottom",
        legend.justification = "center",
        legend.title = element_blank(),
        plot.title = element_text(size=18),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))+
  labs(title = bold('B') ~ '')+
  ylim(c(0,max))
pB

SUB<-DATA%>%filter(Type=="C")
pC<-ggplot(SUB, aes(x=CENPB, y=value, fill=NonB, color=NonB))+
  geom_boxplot(show.legend=FALSE, alpha=0.6)+
  facet_wrap(~NonB, nrow=1)+
  scale_fill_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  scale_color_manual(breaks=c("APR", "DR", "STR", "IR", "MR", "G4", "Z"), values=nonBcol)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=1)+
  labs(x="", y="") +
  theme(panel.background = element_rect(fill = 'white', colour="white"),
        strip.background = element_rect(fill = 'lightgray', colour="white"),
        strip.text = element_blank(),
        panel.grid.major = element_line(colour = 'beige'),
        panel.spacing.y=unit(1, "lines"),
        legend.background = element_rect(color = NA),
        legend.position="top",
        legend.justification = "right",
        legend.title = element_blank(),
        plot.title = element_text(size=18),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  guides(fill = guide_legend(nrow = 1))+
  labs(title = bold('C') ~ '')+
  ylim(c(0,max))
pC

pcomb <- pA +  pB +  pC + plot_layout(widths = c(5, 5, 5)) 
pcomb
outfile=paste(plotdir,"FigS16_ABC_CENPB.jpg", sep="")
ggsave(outfile,plot = pcomb,scale = 1,dpi = 600,limitsize = TRUE,width=13,height=6)

