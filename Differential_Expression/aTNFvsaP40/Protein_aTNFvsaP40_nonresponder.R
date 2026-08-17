
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
## read the DEP1 of NR ---------------
DEP1NR <- read.delim("./a1_preNRvsPostNR_deg2.txt", sep='\t')
DEP1NR$protein2 <- ExtractProteinName(DEP1NR$protein)

proteinNames <- DEP1NR[,c('protein', 'protein2')]
rownames(proteinNames) <- proteinNames$protein2

DEP1NR <- DEP1NR[DEP1NR$type!='NOT',]

## read the DEP2 of NR ---------------
DEP2NR <- read.delim("./a2_preNRvsPostNR_deg2.txt", sep='\t')
DEP2NR$protein2 <- ExtractProteinName(DEP2NR$protein)
DEP2NR <- DEP2NR[DEP2NR$type!='NOT',]


## Venn of DEP1NR and DEP2NR -----------
pdf("./a1vsa2_pre_post_Venn_inNR.pdf",6,5)
plot(Venn(list("DEP1NR" = unique(DEP1NR$protein2),
               "DEP2NR" = unique(DEP2NR$protein2))),
     doWeight=T)
plot(Venn(list("DEP1NR_UP" = unique(DEP1NR[DEP1NR$type=='UP',]$protein2),
               "DEP2NR_UP" = unique(DEP2NR[DEP2NR$type=='UP',]$protein2))),
     doWeight=T)
plot(Venn(list("DEP1NR_DOWN" = unique(DEP1NR[DEP1NR$type=='DOWN',]$protein2),
               "DEP2NR_DOWN" = unique(DEP2NR[DEP2NR$type=='DOWN',]$protein2))),
     doWeight=T)
dev.off()


MyGeneV(unique(DEP1NR$protein2), unique(DEP2NR$protein2))


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

selectProtein <- MyGeneV(unique(DEP1NR$protein2),
                         unique(DEP2NR$protein2))

a1_spe_up <- MyGeneE(unique(DEP1NR[DEP1NR$type=='UP',]$protein2),
                     unique(DEP2NR[DEP2NR$type=='UP',]$protein2))
a2_spe_up <- MyGeneE(unique(DEP2NR[DEP2NR$type=='UP',]$protein2),
                     unique(DEP1NR[DEP1NR$type=='UP',]$protein2))

a1_spe_down <- MyGeneE(unique(DEP1NR[DEP1NR$type=='DOWN',]$protein2),
                       unique(DEP2NR[DEP2NR$type=='DOWN',]$protein2))
a2_spe_down <- MyGeneE(unique(DEP2NR[DEP2NR$type=='DOWN',]$protein2),
                       unique(DEP1NR[DEP1NR$type=='DOWN',]$protein2))


merge.final$a1preNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-NR') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1postNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-NR') &
                                                              mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1preR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-R') &
                                                            mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)
merge.final$a1postR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-R') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a1'),'col']], 1, mean)

merge.final$a2preNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-NR') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2postNR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-NR') &
                                                              mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2preR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('pre-R') &
                                                            mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)
merge.final$a2postR <- apply(merge.final[,mergeSampleInfor[mergeSampleInfor$Group_WF2 %in% c('post-R') &
                                                             mergeSampleInfor$Group_WF3 %in% c('a2'),'col']], 1, mean)

a1_spe_up <- a1_spe_up[a1_spe_up != "HEBP1"]
a2_spe_down <- a2_spe_down[a2_spe_down != "HEBP1"]

plot_res <- merge.final[proteinNames[c(selectProtein,
                                       a1_spe_down,
                                       a1_spe_up,
                                       a2_spe_down,
                                       a2_spe_up),'protein'],
                        c('a1preNR','a1postNR','a2preNR','a2postNR',
                          'a1preR','a1postR','a2preR','a2postR')]
rownames(plot_res) <- c(selectProtein,
                        a1_spe_down,
                        a1_spe_up,
                        a2_spe_down,
                        a2_spe_up)

geneDf <- data.frame(geneName = c(selectProtein,a1_spe_down,a1_spe_up,a2_spe_down,a2_spe_up),
                     type2 = c(rep('Contrast',length(selectProtein)),
                               rep('a1_spe_down', length(a1_spe_down)),
                               rep('a1_spe_up', length(a1_spe_up)),
                               rep('a2_spe_down', length(a2_spe_down)),
                               rep('a2_spe_up', length(a2_spe_up))))
rownames(geneDf) <- geneDf$geneName
row_ha <- rowAnnotation(
  df = geneDf[,c('type2'), drop = FALSE],
  col = list(type2 = c('Contrast' = "#F8CDC2", 
                       "a1_spe_down" = "#C6E6E1",'a1_spe_up' = "#D3C3DE", 
                       "a2_spe_down" = "#D0DAC3",'a2_spe_up'="#D1B9AF")))


genes_to_show <- c("GPR180","TELO2",
                   "APOA4","LAMB1","TGM2","CRAD9","CDH11","CDH13",
                   "CRP","LTF","S100A8","S100A9","SAA2","HMGB2",
                   "MBP","DDX42","NUP133"
)

idx <- c(which(geneDf$geneName %in% genes_to_show))
labels <- geneDf$geneName[idx]
right_ha <- rowAnnotation(
  mark = anno_mark(at = idx, labels = labels, 
                   which = "row", 
                   side = "right"))


pdf("./a1vsa2_pre_post.Heatmap.1_inNR.new.pdf",5.5,7)
Heatmap(t(scale(t(plot_res[,c('a1preNR','a1postNR','a2preNR','a2postNR')]))),
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

MyWriteTab(t(scale(t(plot_res[,c('a1preNR','a1postNR','a2preNR','a2postNR')]))),
           row.names = T, "./a1vsa2_pre_post.Heatmap.1_inNR.new.csv", sep=',')

pdf("./a1vsa2_pre_post.Heatmap.1_inNR.pdf",4,15)
Heatmap(t(scale(t(plot_res[,c('a1preNR','a1postNR','a2preNR','a2postNR')]))),
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



### select the example to show -------
mergeSampleInfor$Group_WF4 <- paste0(mergeSampleInfor$Group_WF3,'_',mergeSampleInfor$Group_WF2)
mergeSampleInfor$Group_WF4 <- factor(mergeSampleInfor$Group_WF4,
                                     levels = c("a1_pre-R","a1_post-R","a2_pre-R","a2_post-R",
                                                "a1_pre-NR","a1_post-NR","a2_pre-NR","a2_post-NR"))

mergeSampleInfor <- mergeSampleInfor[order(mergeSampleInfor$Group_WF4,
                                           mergeSampleInfor$patientID2),]

plot2show <- mergeSampleInfor
## contrast
plot2show$Q9NRV9_HEBP1 <- unlist(merge.final['Q9NRV9_HEBP1', plot2show$col])

## a1NR_DOWN
plot2show$Q86V85_GPR180 <- unlist(merge.final['Q86V85_GPR180', plot2show$col])
plot2show$Q9Y4R8_TELO2 <- unlist(merge.final['Q9Y4R8_TELO2', plot2show$col])

## a1NR_UP
plot2show$P49746_THBS3 <- unlist(merge.final['P49746_THBS3', plot2show$col])
plot2show$P51884_LUM <- unlist(merge.final['P51884_LUM', plot2show$col])
plot2show$P55290_CDH13 <- unlist(merge.final['P55290_CDH13', plot2show$col])
plot2show$Q14623_IHH <- unlist(merge.final['Q14623_IHH', plot2show$col])

## a2NR_DOWN
plot2show$P0DJI9_SAA2 <- unlist(merge.final['P0DJI9_SAA2', plot2show$col])
plot2show$P02741_CRP <- unlist(merge.final['P02741_CRP', plot2show$col])
plot2show$P00738_HP <- unlist(merge.final['P00738_HP', plot2show$col])

## a2NR_UP
plot2show$P02686_MBP <- unlist(merge.final['P02686_MBP', plot2show$col])
plot2show$Q86XP3_DDX42 <- unlist(merge.final['Q86XP3_DDX42', plot2show$col])
plot2show$Q8IWZ3_ANKHD1 <- unlist(merge.final['Q8IWZ3_ANKHD1', plot2show$col])


MyWriteTab(plot2show, "./a1vsa2_plot2show_inNR.csv", sep=',')


MyRunPlot2show2 <- function(plot2show, colName, selSample){
  nameS <- unlist(strsplit(colName,'_'))[2]
  plot2show <- plot2show[plot2show$Group_WF %in% c(selSample),]
  plot2show$Group_WF2 <- factor(plot2show$Group_WF2, levels = c('pre-NR', 'post-NR'))
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

p1 <- MyRunPlot2show2(plot2show,'Q9NRV9_HEBP1',selSample='NR')
p2 <- MyRunPlot2show2(plot2show,'Q86V85_GPR180',selSample='NR')
p3 <- MyRunPlot2show2(plot2show,'Q9Y4R8_TELO2',selSample='NR')

p4 <- MyRunPlot2show2(plot2show,'P49746_THBS3',selSample='NR')
p5 <- MyRunPlot2show2(plot2show,'P51884_LUM',selSample='NR')
p6 <- MyRunPlot2show2(plot2show,'P55290_CDH13',selSample='NR')
p13 <- MyRunPlot2show2(plot2show,'Q14623_IHH',selSample='NR')

p7 <- MyRunPlot2show2(plot2show,'P0DJI9_SAA2',selSample='NR')
p8 <- MyRunPlot2show2(plot2show,'P02741_CRP',selSample='NR')
p9 <- MyRunPlot2show2(plot2show,'P00738_HP',selSample='NR')
p10 <- MyRunPlot2show2(plot2show,'P02686_MBP',selSample='NR')
p11 <- MyRunPlot2show2(plot2show,'Q86XP3_DDX42',selSample='NR')
p12 <- MyRunPlot2show2(plot2show,'Q8IWZ3_ANKHD1',selSample='NR')


pdf("./a1vsa2_pre_post_plot2show_inNR.pdf",11,15)
p1+p2+p3+p4+p5+p6+p13+p7+p8+p9+p10+p11+p12+plot_layout(guides='collect')
dev.off()


MyRunPvalue <- function(plot2show, colName, selSample){
  plot2show <- plot2show[plot2show$Group_WF %in% c(selSample),]
  print(paste0('a1:',wilcox.test(plot2show[plot2show$Group_WF2=='pre-NR' &
                                             plot2show$Group_WF3 =='a1',colName],
                                 plot2show[plot2show$Group_WF2=='post-NR' &
                                             plot2show$Group_WF3 =='a1',colName])$p.value))
  print(paste0('a2:',wilcox.test(plot2show[plot2show$Group_WF2=='pre-NR' &
                                             plot2show$Group_WF3 =='a2',colName],
                                 plot2show[plot2show$Group_WF2=='post-NR' &
                                             plot2show$Group_WF3 =='a2',colName])$p.value))
}

MyRunPvalue(plot2show,'Q9NRV9_HEBP1',selSample='NR')
MyRunPvalue(plot2show,'Q86V85_GPR180',selSample='NR')
MyRunPvalue(plot2show,'Q9Y4R8_TELO2',selSample='NR')
MyRunPvalue(plot2show,'P49746_THBS3',selSample='NR')
MyRunPvalue(plot2show,'P51884_LUM',selSample='NR')
MyRunPvalue(plot2show,'P55290_CDH13',selSample='NR')
MyRunPvalue(plot2show,'Q14623_IHH',selSample='NR')
MyRunPvalue(plot2show,'P0DJI9_SAA2',selSample='NR')
MyRunPvalue(plot2show,'P02741_CRP',selSample='NR')
MyRunPvalue(plot2show,'P00738_HP',selSample='NR')
MyRunPvalue(plot2show,'P02686_MBP',selSample='NR')
MyRunPvalue(plot2show,'Q86XP3_DDX42',selSample='NR')
MyRunPvalue(plot2show,'Q8IWZ3_ANKHD1',selSample='NR')


