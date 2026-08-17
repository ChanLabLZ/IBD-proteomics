
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(dplyr)
library(readxl)
library(ggplot2)
library(RColorBrewer)

protein_counts <- read_excel("./Allsample_ProteinIdentification.2.xlsx")
protein_counts <- as.data.frame(protein_counts)

protein_counts$Group2 <- "NotKnow"
protein_counts[protein_counts$Group == 'UNT',]$Group2 <- 'UNT'
protein_counts[protein_counts$Group == 'HC',]$Group2 <- 'HC'
protein_counts[protein_counts$Group == 'NR' & 
                   protein_counts$Description == 'anti-TNF',]$Group2 <- 'a1NR'
protein_counts[protein_counts$Group == 'R' & 
                protein_counts$Description == 'anti-TNF',]$Group2 <- 'a1R'
protein_counts[protein_counts$Group == 'NR' & 
                protein_counts$Description == 'anti-p40',]$Group2 <- 'a2NR'
protein_counts[protein_counts$Group == 'R' & 
                protein_counts$Description == 'anti-p40',]$Group2 <- 'a2R'

protein_counts <- protein_counts[order(protein_counts$Group2, protein_counts$`Identified proteins`),]
protein_counts$Rank <- c(1:82,1:66,1:53,1:48,1:101,1:153)

protein_counts <- protein_counts[order(protein_counts$Group, protein_counts$`Identified proteins`),]
protein_counts$Rank2 <- c(1:101,1:135,1:114,1:153)

color_manual = c(brewer.pal(9,"Set1"))
color_manual <- c("#377EB8","#FF7F00","#984EA3", "#A65628", "#E41A1C", "#999999", "#F781BF","#4DAF4A")

p1<-ggplot(protein_counts, aes(Rank, `Identified proteins`, col = Group2)) +
    geom_point() + geom_line() +
    theme_wf +
    scale_color_manual(values = color_manual)
    
p2<-ggplot(protein_counts, aes(Rank2, `Identified proteins`, col = Group)) +
    geom_point() + geom_line() +
    theme_wf +
    scale_color_manual(values = color_manual)

p1+p2

pdf("./XHOM_PM_sampleInfor.3.pdf",9.5,8)
print(p1)
print(p2)
dev.off()

