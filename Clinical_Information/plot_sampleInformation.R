
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(dplyr)
library(readxl)
library(ggpubr)

model_excel <- read_excel("../XHOM_PM_sampleInfor.xlsx")
model_excel <- as.data.frame(model_excel)

## counts ------------
table(model_excel$Group, model_excel$Description)

model_excel <- model_excel[,c(1:35)]


## plot ---------

sampleInfor <- model_excel

sampleInfor$Group2 <- "NotKnow"
sampleInfor[sampleInfor$Group == 'UNT',]$Group2 <- 'UNT'
sampleInfor[sampleInfor$Group == 'HC',]$Group2 <- 'HC'
sampleInfor[sampleInfor$Group == 'NR' & 
                sampleInfor$Description == 'anti-TNF',]$Group2 <- 'a1NR'
sampleInfor[sampleInfor$Group == 'R' & 
                sampleInfor$Description == 'anti-TNF',]$Group2 <- 'a1R'
sampleInfor[sampleInfor$Group == 'NR' & 
                sampleInfor$Description == 'anti-p40',]$Group2 <- 'a2NR'
sampleInfor[sampleInfor$Group == 'R' & 
                sampleInfor$Description == 'anti-p40',]$Group2 <- 'a2R'


sampleInfo <- sampleInfor[!sampleInfor$Group %in% c("HC"),]
rm(sampleInfor)
sampleInfo$Group2 <- factor(sampleInfo$Group2, levels = c("UNT","a1NR","a1R","a2NR","a2R"))


# for (i in c("UNT","a1NR","a1R","a2NR","a2R")) {
#     temp <- sampleInfo[sampleInfo$Group2==i,]
#     suT <- summary(temp$Hb)
#     cil <- suT[2]
#     ciu <- suT[5]
#     m <- suT[4]
#     print(paste0(m,' (',cil,'–',ciu,')'))
# }

(p0 <- ggplot(sampleInfo, aes(Group2, `SES-CD(<3)`)) +
    geom_boxplot(na.rm = T, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 1, size = 1, color = "darkred") +
        ggpubr::stat_compare_means(
            comparisons = list(c("a2NR","a2R"),c("a1NR", "a1R")),
            na.rm = TRUE,
            method = "wilcox.test",  
            label = "p.format",      
            label.x.npc = "center",  
            label.y.npc = "top",     
            vjust = -0.5             
        ) + 
    theme_wf2 +
    theme(axis.title.x = element_blank()))


(p1 <- ggplot(sampleInfo, aes(Group2, `CRP(<10)`)) +
    geom_boxplot(na.rm = T, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 1, size = 1, color = "darkred") +
        ggpubr::stat_compare_means(
            comparisons = list(c("a2NR","a2R"),c("a1NR", "a1R")),
            na.rm = TRUE,
            method = "wilcox.test",  
            label = "p.format",      
            label.x.npc = "center",  
            label.y.npc = "top",     
            vjust = -0.5             
        ) + 
    theme_wf2 + ylim(0,130)+
    theme(axis.title.x = element_blank()))


(p2 <- ggplot(sampleInfo, aes(Group2, `FCCP(<250)`)) +
    geom_boxplot(na.rm = T, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 1, size = 1, color = "darkred") +
        ggpubr::stat_compare_means(
            comparisons = list(c("a2NR","a2R"),c("a1NR", "a1R")),
            na.rm = TRUE,
            method = "wilcox.test",  
            label = "p.format",      
            label.x.npc = "center",  
            label.y.npc = "top",     
            vjust = -0.5             
        ) + 
    theme_wf2 + ylim(0,5000)+
    theme(axis.title.x = element_blank()))


(p3 <- ggplot(sampleInfo, aes(Group2,`ESR`)) + 
    geom_boxplot(na.rm = T, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 1, size = 1, color = "darkred") +
        ggpubr::stat_compare_means(
            comparisons = list(c("a2NR","a2R"),c("a1NR", "a1R")),
            na.rm = TRUE,
            method = "wilcox.test",  
            label = "p.format",      
            label.x.npc = "center",  
            label.y.npc = "top",     
            vjust = -0.5             
        ) + 
    theme_wf2 +
    theme(axis.title.x = element_blank()))


(p4 <- ggplot(sampleInfo, aes(Group2,`Hb`)) + 
    geom_boxplot(na.rm = T, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 1, size = 1, color = "darkred") +
        ggpubr::stat_compare_means(
            comparisons = list(c("a2NR","a2R"),c("a1NR", "a1R")),
            na.rm = TRUE,
            method = "wilcox.test",  
            label = "p.format",      
            label.x.npc = "center",  
            label.y.npc = "top",     
            vjust = -0.5             
        ) + 
    theme_wf2 +
    theme(axis.title.x = element_blank()))



pdf("./XHOM_PM_sampleInfor.2.pdf",16,8)
p0+p1+p2+p3+p4
dev.off()

