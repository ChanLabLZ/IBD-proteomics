
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(readxl)
library(tidyverse)
library(sjPlot)
library(emmeans)
library(ggplot2)
library(patchwork)

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
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-p40","anti-TNF"))
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

protein_set3 <- c("PSMD3",
                  "ACLY",
                  "PSMC3",
                  "ADD1",
                  "GNL1",
                  "PSMD4",
                  "CUL3",
                  "NAP1L4",
                  "PSMC4",
                  "PSMC5",
                  "DMTN",
                  "PSMC1",
                  "PSMC2",
                  "COPS7B",
                  "RILP",
                  "PSMD13",
                  "VCP",
                  "ADSL",
                  "NCAPH",
                  "PA2G4")
protein_set3 <- proName[protein_set3,'name1']

df <- as.data.frame(t(df_final[protein_set3,tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']

df <- df %>% 
  mutate(Protein_Set3_Score = rowMeans(scale(select(., all_of(protein_set3)))))

# logit(P(Response=1)) = beta_0 + beta_1 * Protein_Score + beta_2 * Drug + beta_3 * ({Protein_Score} * {Drug})

fit_interaction <- glm(Response ~ Protein_Set3_Score * Drug, 
                       data = df, 
                       family = binomial(link = "logit"))

summary(fit_interaction)

p1 <- plot_model(fit_interaction, type = "pred", terms = c("Protein_Set3_Score", "Drug"),
                 title = "Interaction Effect of Proteasome Score and Drug Type on Response",
                 axis.title = c("Proteasome Score (Set 3)", "Predicted Response Probability")) +
  theme_wf2

slopes <- emtrends(fit_interaction, ~ Drug, var = "Protein_Set3_Score")
print(slopes)
pairs(slopes)



## for unpaired
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
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-p40","anti-TNF"))
# tableS2 <- tableS2[tableS2$Paired == 'Yes',]
# tableS2 <- tableS2[!tableS2$associatedNo %in% c('No14','No140','No238'),]

tableS2 <- tableS2[tableS2$Paired == 'No',]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No210','No234',
                                                'No358','No376','No373','No371','No417','No424'),]

proName <- data.frame(name1 = rownames(df_final),
                      name2 = ExtractProteinName(rownames(df_final)))
rownames(proName) <- proName$name2

df <- as.data.frame(t(df_final[protein_set3,tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']

df <- df %>% 
  mutate(Protein_Set3_Score = rowMeans(scale(select(., all_of(protein_set3)))))


# logit(P(Response=1)) = beta_0 + beta_1 * Protein_Score + beta_2 * Drug + beta_3 * ({Protein_Score} * {Drug})

fit_interaction <- glm(Response ~ Protein_Set3_Score * Drug, 
                       data = df, 
                       family = binomial(link = "logit"))

p2 <- plot_model(fit_interaction, type = "pred", terms = c("Protein_Set3_Score", "Drug"),
           title = "Interaction Effect of Proteasome Score and Drug Type on Response",
           axis.title = c("Proteasome Score (Set 3)", "Predicted Response Probability")) +
  theme_wf2

slopes <- emtrends(fit_interaction, ~ Drug, var = "Protein_Set3_Score")
print(slopes)
pairs(slopes)




## for all samples
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
tableS2$Drug_follow_up <- factor(tableS2$Drug_follow_up, levels = c("anti-p40","anti-TNF"))

tableS2 <- tableS2[!tableS2$associatedNo %in% c('No14','No140','No238'),]
tableS2 <- tableS2[!tableS2$associatedNo %in% c('No210','No234',
                                                'No358','No376','No373','No371','No417','No424'),]

proName <- data.frame(name1 = rownames(df_final),
                      name2 = ExtractProteinName(rownames(df_final)))
rownames(proName) <- proName$name2

df <- as.data.frame(t(df_final[protein_set3,tableS2$associatedNo]))
rownames(tableS2) <- tableS2$associatedNo
tableS2$Response <- ifelse(tableS2$Group_follow_up=='R',1,0)
df$Response <- tableS2[row.names(df),'Response']
df$Drug <- tableS2[row.names(df),'Drug_follow_up']

df <- df %>% 
  mutate(Protein_Set3_Score = rowMeans(scale(select(., all_of(protein_set3)))))

fit_interaction <- glm(Response ~ Protein_Set3_Score * Drug, 
                       data = df, 
                       family = binomial(link = "logit"))

p3 <- plot_model(fit_interaction, type = "pred", terms = c("Protein_Set3_Score", "Drug"),
           title = "Interaction Effect of Proteasome Score and Drug Type on Response",
           axis.title = c("Proteasome Score (Set 3)", "Predicted Response Probability")) +
  theme_wf2

slopes <- emtrends(fit_interaction, ~ Drug, var = "Protein_Set3_Score")
print(slopes)
pairs(slopes)


pdf("./pre_all_PCA_loading.genes_using4_testInUNT_InteractionAnalysis_paired_unpaired_all.pdf", 10,7.5)
print(p1)
print(p2)
print(p3)
dev.off()

