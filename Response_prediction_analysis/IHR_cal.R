
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
color_manual = c(brewer.pal(9,"Set1"))


## read protein matrix ------------- 
df_final <- read_excel('../protein_Samplematrix_imputeNA1.xlsx')
df_final <- as.data.frame(df_final)
rownames(df_final) <- df_final$protein
df_final <- df_final[,-1]
df_final[1:4,1:4]

# read sample information ------------- 
sampleInfor <- read_excel("../XHOM_PM_sampleInfor.xlsx")
sampleInfor <- as.data.frame(sampleInfor)
sampleInfor$Type

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
model_excel <- sampleInfor
rownames(model_excel) <- model_excel$Type

sample.inf <- data.frame(col = colnames(df_final))
row.names(sample.inf) <- sample.inf$col

merged_df <- merge(sample.inf, model_excel, by = "row.names", all = TRUE)
selIndex <- grep(pattern="rep", merged_df$col)
merged_df2 <- merged_df[selIndex,]
merged_df <- merged_df[!is.na(merged_df$sampleID),]

temp_col <- c()
for (i in merged_df2$Row.names) {
  temp_col <- c(temp_col, unlist(strsplit(i, "_"))[1])
}

merged_df2$Row.names <- temp_col
rm(temp_col)

rownames(merged_df2) <- merged_df2$Row.names
merged_df2 <- merged_df2[,1:2]
pp <- cbind(merged_df2, model_excel[row.names(merged_df2),])
merged_df2 <- pp
rownames(merged_df2) <- merged_df2$col

pp <- rbind(merged_df, merged_df2)
merged_df <- pp
rm(pp)

sample.inf <- merged_df
rm(merged_df)
rownames(sample.inf) <- sample.inf$col


## table S2
tableS2 <- read_excel("./tableS2.xlsx")
tableS2 <- as.data.frame(tableS2)
tableS2 <- tableS2[tableS2$`follow-up` == 'Yes',]
sample.inf <- sample.inf[sample.inf$patientID2 %in% tableS2$patientID,]
sample.inf <- sample.inf[sample.inf$Group=='UNT',]
sample.inf <- sample.inf[!sample.inf$col %in% c("No150_rep", "No267_rep", "No277_rep", "No29_rep",  
                                                "No361_rep", "No439_rep", "No465_rep","No472_rep", 
                                                "No514_rep", "No528_rep"),]
rownames(sample.inf) <- sample.inf$patientID3
rownames(tableS2) <- tableS2$sampleID
tableS2$associatedNo <- sample.inf[row.names(tableS2),'Type']
tableS2 <- tableS2[!tableS2$Group_follow_up %in% c('NotAssign'),]
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-TNF","anti-p40"))
tableS2 <- tableS2[tableS2$Paired == 'Yes',]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No14','No140','No238'),]

# tableS2 <- tableS2[tableS2$Paired == 'No',]
# tableS2 <- tableS2[!tableS2$associatedNo %in% c('No210','No234',
#                                                 'No358','No376','No373','No371','No417','No424'),]


ExtractProteinName <- function(rawProtein){
  cleanName <- c()
  for (i in c(1:length(rawProtein))) {
    cleanName <- c(cleanName, unlist(strsplit(rawProtein[i],'_'))[2])
  }
  return(cleanName)
}

proName <- data.frame(name1 = rownames(df_final),
                      name2 = ExtractProteinName(rownames(df_final)))
rownames(proName) <- proName$name2

PC2_3 <- read.delim("./pre_all_PCA_loading.genes_using4.all60.txt", sep='\t') 

protein_set2 <- PC2_3[PC2_3$type2 == 'PC2','protein']
protein_set3 <- PC2_3[PC2_3$type2 == 'PC3','protein']

protein_set2 <- proName[protein_set2,'name1']
protein_set3 <- proName[protein_set3,'name1']


df <- as.data.frame(t(df_final[c(protein_set2, protein_set3),
                               tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']


df <- df %>%
  rowwise() %>%
  mutate(Set1_Mean = mean(c_across(all_of(protein_set3)), na.rm = TRUE),
         Set2_Mean = mean(c_across(all_of(protein_set2)), na.rm = TRUE),
         Set1_Set2_Ratio = Set1_Mean - Set2_Mean) %>%
  ungroup()

p1<-ggplot(df, aes(Drug, Set1_Set2_Ratio, col = Drug)) +
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_color_manual(values = color_manual) +
  ggpubr::stat_compare_means(
    comparisons = list(c("anti-TNF", "anti-p40")),
    method = "wilcox.test",
    label = "p.format",
    label.x.npc = "center",
    label.y.npc = "top",
    vjust = -0.5
  ) + 
  labs(x = '', y = 'Inflammation-to-Homeostasis Ratio', title = 'paired') +
  theme_wf2





## table S2
tableS2 <- read_excel("./tableS2.xlsx")
tableS2 <- as.data.frame(tableS2)
tableS2 <- tableS2[tableS2$`follow-up` == 'Yes',]
sample.inf <- sample.inf[sample.inf$patientID2 %in% tableS2$patientID,]
sample.inf <- sample.inf[sample.inf$Group=='UNT',]
sample.inf <- sample.inf[!sample.inf$col %in% c("No150_rep", "No267_rep", "No277_rep", "No29_rep",  
                                                "No361_rep", "No439_rep", "No465_rep","No472_rep", 
                                                "No514_rep", "No528_rep"),]
rownames(sample.inf) <- sample.inf$patientID3
rownames(tableS2) <- tableS2$sampleID
tableS2$associatedNo <- sample.inf[row.names(tableS2),'Type']
tableS2 <- tableS2[!tableS2$Group_follow_up %in% c('NotAssign'),]
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-TNF","anti-p40"))
tableS2 <- tableS2[tableS2$Paired == 'No',]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No210','No234',
                                                'No358','No376','No373','No371','No417','No424'),]


df <- as.data.frame(t(df_final[c(protein_set2, protein_set3),
                               tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']


df <- df %>%
  rowwise() %>%
  mutate(Set1_Mean = mean(c_across(all_of(protein_set3)), na.rm = TRUE),
         Set2_Mean = mean(c_across(all_of(protein_set2)), na.rm = TRUE),
         Set1_Set2_Ratio = Set1_Mean - Set2_Mean) %>%
  ungroup()

p2<-ggplot(df, aes(Drug, Set1_Set2_Ratio, col = Drug)) +
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_color_manual(values = color_manual) +
  ggpubr::stat_compare_means(
    comparisons = list(c("anti-TNF", "anti-p40")),
    method = "wilcox.test",
    label = "p.format",
    label.x.npc = "center",
    label.y.npc = "top",
    vjust = -0.5
  ) + 
  labs(x = '', y = 'Inflammation-to-Homeostasis Ratio', title = 'unpaired') +
  theme_wf2





## table S2
tableS2 <- read_excel("./tableS2.xlsx")
tableS2 <- as.data.frame(tableS2)
tableS2 <- tableS2[tableS2$`follow-up` == 'Yes',]
sample.inf <- sample.inf[sample.inf$patientID2 %in% tableS2$patientID,]
sample.inf <- sample.inf[sample.inf$Group=='UNT',]
sample.inf <- sample.inf[!sample.inf$col %in% c("No150_rep", "No267_rep", "No277_rep", "No29_rep",  
                                                "No361_rep", "No439_rep", "No465_rep","No472_rep", 
                                                "No514_rep", "No528_rep"),]
rownames(sample.inf) <- sample.inf$patientID3
rownames(tableS2) <- tableS2$sampleID
tableS2$associatedNo <- sample.inf[row.names(tableS2),'Type']
tableS2 <- tableS2[!tableS2$Group_follow_up %in% c('NotAssign'),]
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-TNF","anti-p40"))
# tableS2 <- tableS2[tableS2$Paired == 'No',]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No14','No140','No238'),]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No210','No234',
                                                'No358','No376','No373','No371','No417','No424'),]


df <- as.data.frame(t(df_final[c(protein_set2, protein_set3),
                               tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']


df <- df %>%
  rowwise() %>%
  mutate(Set1_Mean = mean(c_across(all_of(protein_set3)), na.rm = TRUE),
         Set2_Mean = mean(c_across(all_of(protein_set2)), na.rm = TRUE),
         Set1_Set2_Ratio = Set1_Mean - Set2_Mean) %>%
  ungroup()

p3<-ggplot(df, aes(Drug, Set1_Set2_Ratio, col = Drug)) +
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_color_manual(values = color_manual) +
  ggpubr::stat_compare_means(
    comparisons = list(c("anti-TNF", "anti-p40")),
    method = "wilcox.test",
    label = "p.format",
    label.x.npc = "center",
    label.y.npc = "top",
    vjust = -0.5
  ) + 
  labs(x = '', y = 'Inflammation-to-Homeostasis Ratio', title = 'all') +
  theme_wf2


pdf("./pre_all_PCA_loading.genes_using4_testInUNT_IHR_paired_unpaired.pdf",
    8.7,4)
p1+p2+p3+plot_layout(guides='collect')
dev.off()
