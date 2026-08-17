
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(Vennerable)
library(readxl)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

ExtractProteinName <- function(rawProtein){
  cleanName <- c()
  for (i in c(1:length(rawProtein))) {
    cleanName <- c(cleanName, unlist(strsplit(rawProtein[i],'_'))[2])
  }
  return(cleanName)
}

# DEPs a1 vs a2 in Responders --------------
## read the DEP1 of R ---------------
DEP1R <- read.delim("./a1_preRvsPostR_deg2.txt", sep='\t')
DEP1R$protein2 <- ExtractProteinName(DEP1R$protein)

proteinNames <- DEP1R[,c('protein', 'protein2')]
rownames(proteinNames) <- proteinNames$protein2

DEP1R <- DEP1R[DEP1R$type!='NOT',]

## read the DEP2 of R ---------------
DEP2R <- read.delim("./a2_preRvsPostR_deg2.txt", sep='\t')
DEP2R$protein2 <- ExtractProteinName(DEP2R$protein)
DEP2R <- DEP2R[DEP2R$type!='NOT',]

## Venn of DEP1R and DEP2R -----------

pdf("./a1vsa2_pre_post_Venn.pdf",6,5)
plot(Venn(list("DEP1R" = unique(DEP1R$protein2),
               "DEP2R" = unique(DEP2R$protein2))),
     doWeight=T)
plot(Venn(list("DEP1R_UP" = unique(DEP1R[DEP1R$type=='UP',]$protein2),
               "DEP2R_UP" = unique(DEP2R[DEP2R$type=='UP',]$protein2))),
     doWeight=T)
plot(Venn(list("DEP1R_DOWN" = unique(DEP1R[DEP1R$type=='DOWN',]$protein2),
               "DEP2R_DOWN" = unique(DEP2R[DEP2R$type=='DOWN',]$protein2))),
     doWeight=T)
dev.off()

MyGeneV(unique(DEP1R$protein2),
        unique(DEP2R$protein2))
# "STC2"  "RBP4"  "TTR"   "ALB"   "EPHB2" "THBS3" "FBLN7" "SFRP2"


## example to show -----------
### read protein matrix ------------- 

df_final <- read_excel('../protein_Samplematrix_imputeNA1.xlsx')
df_final <- as.data.frame(df_final)
rownames(df_final) <- df_final$protein
df_final <- df_final[,-1]
df_final[1:4,1:4]

merge.final <- df_final

### read sample information ------------- 

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

pre_a1 <- read_excel("../Sample_Pre_Post.xlsx")
pre_a1 <- as.data.frame(pre_a1)
pre_a1 <- pre_a1[pre_a1$Group=='a1',]

sample.inf.temp <- sample.inf[sample.inf$sampleID %in% pre_a1$sampleID,]
sample.inf <- sample.inf.temp[sample.inf.temp$Description == '/',]
sample.inf.temp <- sample.inf.temp[sample.inf.temp$Description != '/',]

sample.inf.temp <- sample.inf.temp[!sample.inf.temp$col %in% c('No240_rep','No519_rep'),]

rownames(sample.inf) <- sample.inf$patientID2
rownames(sample.inf.temp) <- sample.inf.temp$patientID2

sample.inf$Group_WF <- sample.inf.temp[row.names(sample.inf),'Group']
sample.inf.temp$Group_WF <- sample.inf.temp$Group

sample.inf$Group_WF2 <- paste0('pre-',sample.inf$Group_WF)
sample.inf.temp$Group_WF2 <- paste0('post-',sample.inf.temp$Group_WF)

a1_sample.inf <- rbind(sample.inf, sample.inf.temp)
a1_sample.inf$Group_WF3 <- 'a1'


### select pre samples by a2 --------
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

pre_a2 <- read_excel("../Sample_Pre_Post.xlsx")
pre_a2 <- as.data.frame(pre_a2)
pre_a2 <- pre_a2[pre_a2$Group=='a2',]

sample.inf.temp <- sample.inf[sample.inf$sampleID %in% pre_a2$sampleID,]
sample.inf <- sample.inf.temp[sample.inf.temp$Description == '/',]
sample.inf.temp <- sample.inf.temp[sample.inf.temp$Description != '/',]

sample.inf <- sample.inf[!sample.inf$col %in% c('No361_rep','No472_rep'),]
sample.inf.temp <- sample.inf.temp[!sample.inf.temp$col %in% c('No322_rep'),]

rownames(sample.inf) <- sample.inf$patientID2
rownames(sample.inf.temp) <- sample.inf.temp$patientID2

sample.inf$Group_WF <- sample.inf.temp[row.names(sample.inf),'Group']
sample.inf.temp$Group_WF <- sample.inf.temp$Group

sample.inf$Group_WF2 <- paste0('pre-',sample.inf$Group_WF)
sample.inf.temp$Group_WF2 <- paste0('post-',sample.inf.temp$Group_WF)

a2_sample.inf <- rbind(sample.inf, sample.inf.temp)
a2_sample.inf$Group_WF3 <- 'a2'

mergeSampleInfor <- rbind(a1_sample.inf, a2_sample.inf)
rownames(mergeSampleInfor) <- mergeSampleInfor$col



### plot heatmap 2: Select DEP --------------------

selectProtein <- MyGeneV(unique(DEP1R$protein2),
                         unique(DEP2R$protein2))

a1_spe_up <- MyGeneE(unique(DEP1R[DEP1R$type=='UP',]$protein2),
                     unique(DEP2R[DEP2R$type=='UP',]$protein2))
a2_spe_up <- MyGeneE(unique(DEP2R[DEP2R$type=='UP',]$protein2),
                     unique(DEP1R[DEP1R$type=='UP',]$protein2))

a1_spe_down <- MyGeneE(unique(DEP1R[DEP1R$type=='DOWN',]$protein2),
                       unique(DEP2R[DEP2R$type=='DOWN',]$protein2))
a2_spe_down <- MyGeneE(unique(DEP2R[DEP2R$type=='DOWN',]$protein2),
                       unique(DEP1R[DEP1R$type=='DOWN',]$protein2))


merge.final$a1preR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-R') &
                                                            mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1postR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-R') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1preNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-NR') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1postNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-NR') &
                                                              mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)

merge.final$a2preR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-R') &
                                                            mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2postR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-R') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2preNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-NR') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2postNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-NR') &
                                                              mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)


plot_res <- merge.final[proteinNames[c(selectProtein,
                                       a1_spe_down,
                                       a1_spe_up,
                                       a2_spe_down,
                                       a2_spe_up),'protein'],
                        c('a1preR','a1postR','a2preR','a2postR',
                          'a1preNR','a1postNR','a2preNR','a2postNR')]
rownames(plot_res) <- c(selectProtein,
                        a1_spe_down,
                        a1_spe_up,
                        a2_spe_down,
                        a2_spe_up)

geneDf <- data.frame(geneName = c(selectProtein,a1_spe_down,a1_spe_up,a2_spe_down,a2_spe_up),
                     type2 = c(rep('Common',length(selectProtein)),
                               rep('a1_spe_down', length(a1_spe_down)),
                               rep('a1_spe_up', length(a1_spe_up)),
                               rep('a2_spe_down', length(a2_spe_down)),
                               rep('a2_spe_up', length(a2_spe_up))))
rownames(geneDf) <- geneDf$geneName
row_ha <- rowAnnotation(
  df = geneDf[,c('type2'), drop = FALSE],
  col = list(type2 = c('Common' = "#F8CDC2", 
                       "a1_spe_down" = "#C6E6E1",'a1_spe_up' = "#D3C3DE", 
                       "a2_spe_down" = "#D0DAC3",'a2_spe_up'="#D1B9AF")))


genes_to_show <- c("CRP","SAA1","SAA2","C9","FCGR3A",
                   "LIPC","SKP1","BMP6",
                   "PSMD3","DNAJC13","MYH10","NF2", # R up in pre
                   "APOA1","AHSG","MYOC","CHAD","VASN","FCGBP"
)

idx <- c(which(geneDf$geneName %in% genes_to_show))
labels <- geneDf$geneName[idx]
right_ha <- rowAnnotation(
  mark = anno_mark(at = idx, labels = labels, 
                   which = "row", 
                   side = "right"))

pdf("./a1vsa2_pre_post.Heatmap.1.new.pdf",5.5,7)
Heatmap(t(scale(t(plot_res[,c('a1preR','a1postR','a2preR','a2postR')]))),
        name = "z-score",
        show_row_names = F,
        show_column_names = T, cluster_columns = F,cluster_rows = F, 
        # column_split = orderSample[,c('Group_WF2')], 
        # row_split = geneDf[,c('type2')],
        left_annotation = row_ha,
        # bottom_annotation = annotation,
        right_annotation = right_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()

MyWriteTab(t(scale(t(plot_res[,c('a1preR','a1postR','a2preR','a2postR')]))),
           row.names = T, "./a1vsa2_pre_post.Heatmap.1.new.csv", sep=',')

pdf("./a1vsa2_pre_post.Heatmap.1.pdf",4,15)
Heatmap(t(scale(t(plot_res[,c('a1preR','a1postR','a2preR','a2postR')]))),
        name = "z-score",
        show_row_names = T,
        show_column_names = T, cluster_columns = F,cluster_rows = F, 
        # column_split = orderSample[,c('Group_WF2')], 
        # row_split = geneDf[,c('type2')],
        left_annotation = row_ha,
        # bottom_annotation = annotation,
        # right_annotation = right_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()


selectProtein2 <- c("CRP","FCGR3A","SAA1","SAA2","C9","TXN","PLA2G2A","AGPAT2","ACTN4","HIP1","ARHGAP25")
selectProtein3 <- c('BMP6','COL11A2','COMP','NOTUM','ASPN','LIPC','ANGPTL8','SKP1','CAST')
selectProtein4 <- c('GOLGA3','MYH10','NF2','GIT2','DYNC1LI1','KLC4','LRPPRC','PSMC5','PSMD3','RANBP1','SHTN1','VPS26A','LPCAT1','ARK2N','DNAJC13','PARP9','PTPRF')
selectProtein5 <- c("AHSG","APOA1","CST3","PLG","SPP2","TF","SERPINF2","ALDOB","ABCA1","F13B","MASP1","CHAD","GPLD1","MGP","MYOC","SFRP4","CHRDL2","FCN3","PGLYRP2","TMEM126A","ADH1B","ADH4","PON1","SORD","DPP4","AFM","PON3","AOC3","ANPEP","BTD","GPT","HHIPL1","VASN")

plot_res <- merge.final[proteinNames[c(selectProtein,
                                       selectProtein2,
                                       selectProtein3,
                                       selectProtein4,
                                       selectProtein5),'protein'],
                        c('a1preR','a1postR','a2preR','a2postR',
                          'a1preNR','a1postNR','a2preNR','a2postNR')]
rownames(plot_res) <- c(selectProtein,
                        selectProtein2,
                        selectProtein3,
                        selectProtein4,
                        selectProtein5)

geneDf <- data.frame(geneName = c(selectProtein,selectProtein2,selectProtein3,selectProtein4,selectProtein5),
                     type2 = c(rep('Common',length(selectProtein)),
                               rep('a1_spe_down', length(selectProtein2)),
                               rep('a1_spe_up', length(selectProtein3)),
                               rep('a2_spe_down', length(selectProtein4)),
                               rep('a2_spe_up', length(selectProtein5))))
rownames(geneDf) <- geneDf$geneName
row_ha <- rowAnnotation(
  df = geneDf[,c('type2'), drop = FALSE],
  col = list(type2 = c('Common' = "#F8CDC2", 
                       "a1_spe_down" = "#C6E6E1",'a1_spe_up' = "#D3C3DE", 
                       "a2_spe_down" = "#D0DAC3",'a2_spe_up'="#D1B9AF")))


pdf("./a1vsa2_pre_post.Heatmap.pdf",4,13)
Heatmap(t(scale(t(plot_res[,c('a1preR','a1postR','a2preR','a2postR')]))),
        name = "z-score",
        show_row_names = T,
        show_column_names = T, cluster_columns = F,cluster_rows = F, 
        left_annotation = row_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()


### select the example to show -------
mergeSampleInfor$Group_WF4 <- paste0(mergeSampleInfor$Group_WF3,'_',mergeSampleInfor$Group_WF2)
mergeSampleInfor$Group_WF4 <- factor(mergeSampleInfor$Group_WF4,
                                     levels = c("a1_pre-R","a1_post-R","a2_pre-R","a2_post-R",
                                                "a1_pre-NR","a1_post-NR","a2_pre-NR","a2_post-NR"))

mergeSampleInfor <- mergeSampleInfor[order(mergeSampleInfor$Group_WF4,
                                           mergeSampleInfor$patientID2),]

plot2show <- mergeSampleInfor
## common
plot2show$P02753_RBP4 <- unlist(merge.final['P02753_RBP4', plot2show$col])
plot2show$P02766_TTR <- unlist(merge.final['P02766_TTR', plot2show$col])
plot2show$P02768_ALB <- unlist(merge.final['P02768_ALB', plot2show$col])

## a1R_DOWN
plot2show$P02741_CRP <- unlist(merge.final['P02741_CRP', plot2show$col])
plot2show$P0DJI8_SAA1 <- unlist(merge.final['P0DJI8_SAA1', plot2show$col])
plot2show$P0DJI9_SAA2 <- unlist(merge.final['P0DJI9_SAA2', plot2show$col])

## a1R_UP
plot2show$P22004_BMP6 <- unlist(merge.final['P22004_BMP6', plot2show$col])
plot2show$P11150_LIPC <- unlist(merge.final['P11150_LIPC', plot2show$col])
plot2show$Q6UXH0_ANGPTL8 <- unlist(merge.final['Q6UXH0_ANGPTL8', plot2show$col])


## a2R_DOWN
plot2show$P42704_LRPPRC <- unlist(merge.final['P42704_LRPPRC', plot2show$col])
plot2show$O75165_DNAJC13 <- unlist(merge.final['O75165_DNAJC13', plot2show$col])
plot2show$O43242_PSMD3 <- unlist(merge.final['O43242_PSMD3', plot2show$col])


## a2R_UP
plot2show$P02647_APOA1 <- unlist(merge.final['P02647_APOA1', plot2show$col])
plot2show$P02765_AHSG <- unlist(merge.final['P02765_AHSG', plot2show$col])


MyWriteTab(plot2show, "./a1vsa2_plot2show.csv", sep=',')


MyRunPlot2show <- function(plot2show, colName, selSample){
  name <- unlist(strsplit(colName,'_'))[2]
  plot2show <- plot2show[plot2show$Group_WF %in% c(selSample),]
  p1<-ggplot(plot2show, aes(x=Group_WF4, y=.data[[colName]],group = patientID2)) +
    geom_line(color = "black", linewidth = 0.5) +          # 连接线
    geom_point(color = "black", size = 2) +
    ggpubr::stat_compare_means(
      comparisons = list(c(c(levels(plot2show$Group_WF4)[1:2])), 
                         c(c(levels(plot2show$Group_WF4)[3:4]))),
      method = "wilcox.test",  # 或 "t.test", "anova", "kruskal.test" 等
      label = "p.format",      # 显示 p 值格式
      label.x.npc = "center",  # 水平位置
      label.y.npc = "top",     # 垂直位置
      vjust = -0.5             # 垂直调整
    ) + 
    labs(x = "", y = expression(log[2](intensity)), title = name) +
    theme_wf +
    theme(plot.title = element_text(size = 20, colour = "black", hjust = 0.65, face = "bold"),
          plot.title.position = "plot",
          axis.text.x.bottom = element_text(angle = 30, hjust = 1))
  return(p1)
}


MyRunPlot2show2 <- function(plot2show, colName, selSample){
  nameS <- unlist(strsplit(colName,'_'))[2]
  plot2show <- plot2show[plot2show$Group_WF %in% c(selSample),]
  plot2show$Group_WF2 <- factor(plot2show$Group_WF2, levels = c('pre-R', 'post-R'))
  p1<-ggplot(plot2show, aes(x = Group_WF2, y = .data[[colName]], color = Group_WF3, group = Group_WF3)) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          ymin = quantile(x, 0.25),
          ymax = quantile(x, 0.75),
          y = median(x))},
      geom = "errorbar",
      width = 0.25,
      linewidth = 0.8) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point",size = 3,shape = 19) +
    scale_color_manual(values = c("a1" = "#8B4513", "a2" = "#4682B4"),name = "") +
    labs(x = "", y = "log2(intensity)", title = nameS) +
    theme_wf2 +
    theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1),
          legend.key.width = unit(0.8, "cm"))
  return(p1)
}

p1 <- MyRunPlot2show2(plot2show,'P02753_RBP4',selSample='R')
p2 <- MyRunPlot2show2(plot2show,'P02766_TTR',selSample='R')
p3 <- MyRunPlot2show2(plot2show,'P02768_ALB',selSample='R')
p4 <- MyRunPlot2show2(plot2show,'P02741_CRP',selSample='R')
p5 <- MyRunPlot2show2(plot2show,'P0DJI8_SAA1',selSample='R')
p6 <- MyRunPlot2show2(plot2show,'P0DJI9_SAA2',selSample='R')
p7 <- MyRunPlot2show2(plot2show,'P22004_BMP6',selSample='R')
p8 <- MyRunPlot2show2(plot2show,'P11150_LIPC',selSample='R')
p9 <- MyRunPlot2show2(plot2show,'Q6UXH0_ANGPTL8',selSample='R')
p10 <- MyRunPlot2show2(plot2show,'P02647_APOA1',selSample='R')
p11 <- MyRunPlot2show2(plot2show,'P02765_AHSG',selSample='R')
p12 <- MyRunPlot2show2(plot2show,'P42704_LRPPRC',selSample='R')
p13 <- MyRunPlot2show2(plot2show,'O75165_DNAJC13',selSample='R')
p14 <- MyRunPlot2show2(plot2show,'O43242_PSMD3',selSample='R')


pdf("./a1vsa2_pre_post_plot2show.pdf",11,15)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p12+p13+p10+p11+plot_layout(guides='collect')
dev.off()


MyRunPvalue <- function(plot2show, colName, selSample){
  plot2show <- plot2show[plot2show$Group_WF %in% c(selSample),]
  print(paste0('a1:',wilcox.test(plot2show[plot2show$Group_WF2=='pre-R' &
                                             plot2show$Group_WF3 =='a1',colName],
                                 plot2show[plot2show$Group_WF2=='post-R' &
                                             plot2show$Group_WF3 =='a1',colName])$p.value))
  print(paste0('a2:',wilcox.test(plot2show[plot2show$Group_WF2=='pre-R' &
                                             plot2show$Group_WF3 =='a2',colName],
                                 plot2show[plot2show$Group_WF2=='post-R' &
                                             plot2show$Group_WF3 =='a2',colName])$p.value))
}

MyRunPvalue(plot2show,'P02753_RBP4',selSample='R')
MyRunPvalue(plot2show,'P02766_TTR',selSample='R')
MyRunPvalue(plot2show,'P02768_ALB',selSample='R')
MyRunPvalue(plot2show,'P02741_CRP',selSample='R')
MyRunPvalue(plot2show,'P0DJI8_SAA1',selSample='R')
MyRunPvalue(plot2show,'P0DJI9_SAA2',selSample='R')
MyRunPvalue(plot2show,'P22004_BMP6',selSample='R')
MyRunPvalue(plot2show,'P11150_LIPC',selSample='R')
MyRunPvalue(plot2show,'Q6UXH0_ANGPTL8',selSample='R')
MyRunPvalue(plot2show,'P02647_APOA1',selSample='R')
MyRunPvalue(plot2show,'P02765_AHSG',selSample='R')
MyRunPvalue(plot2show,'P42704_LRPPRC',selSample='R')
MyRunPvalue(plot2show,'Q96B23_ARK2N',selSample='R')
  




