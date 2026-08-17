
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(readxl)
library(FactoMineR)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


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


## select pre samples --------
## a1 is aTNF
pre_a1 <- read_excel("../Sample_Pre_Post.xlsx")
pre_a1 <- as.data.frame(pre_a1)
pre_a1 <- pre_a1[pre_a1$Group=='a1',]

sample.inf.temp <- sample.inf[sample.inf$sampleID %in% pre_a1$sampleID,]
sample.inf <- sample.inf.temp[sample.inf.temp$Description == '/',]
sample.inf.temp <- sample.inf.temp[sample.inf.temp$Description != '/',]

# sample.inf.temp <- sample.inf.temp[!sample.inf.temp$col %in% c('No240_rep','No519_rep'),]

rownames(sample.inf) <- sample.inf$patientID2
rownames(sample.inf.temp) <- sample.inf.temp$patientID2

sample.inf$Group_WF <- sample.inf.temp[row.names(sample.inf),'Group']
sample.inf.temp$Group_WF <- sample.inf.temp$Group

sample.inf$Group_WF2 <- paste0('pre-',sample.inf$Group_WF)
sample.inf.temp$Group_WF2 <- paste0('post-',sample.inf.temp$Group_WF)

merge.final <- df_final
merge.sample.inf <- rbind(sample.inf, sample.inf.temp)


## DEP2 post vs pre in NR -------------
sample.inf <- merge.sample.inf[merge.sample.inf$Group_WF2 %in% c("pre-NR","post-NR"),]
df_final <- merge.final[,sample.inf$col]
rownames(sample.inf) <- sample.inf$col

df_final2 <- df_final

rowMedian <- function(df){
  return(apply(df, MARGIN=1, median))
}

df_final2$log2FC <- rowMedian(df_final2[,row.names(sample.inf[sample.inf$Group_WF2 == 'post-NR',])]) - rowMedian(df_final2[,row.names(sample.inf[sample.inf$Group_WF2 == 'pre-NR',])])

u_test_function <- function(row) {
  group_a <- row[row.names(sample.inf[sample.inf$Group_WF2 == 'pre-NR',])]
  group_b <- row[row.names(sample.inf[sample.inf$Group_WF2 == 'post-NR',])]
  test_result <- wilcox.test(group_a, group_b)
  return(test_result$p.value)
}

p_values <- apply(df_final2, 1, u_test_function)
df_final2$p_value <- p_values
df_final2$FDR <- p.adjust(p_values, method = "fdr")
df_final2$type <- 'NOT'
df_final2[df_final2$log2FC >= log2(1.2) & df_final2$p_value < 0.05,]$type <- 'UP'
df_final2[df_final2$log2FC <= -log2(1.2) & df_final2$p_value < 0.05,]$type <- 'DOWN'
table(df_final2$type)

# DOWN  NOT   UP 
# 2 2747   26 

ExtractProteinName <- function(rawProtein){
  cleanName <- c()
  for (i in c(1:length(rawProtein))) {
    cleanName <- c(cleanName, unlist(strsplit(rawProtein[i],'_'))[2])
  }
  return(cleanName)
}

df_final2$protein <- row.names(df_final2)
df_final2$protein2 <- ExtractProteinName(df_final2$protein)

write.table(df_final2[,c('log2FC','p_value','FDR','type','protein')], 
            "./Protein/a1_pre_post/preNR_postNR/a1_preNRvsPostNR_deg2.txt", sep='\t',
            quote = F, col.names = T, row.names = F)



## DEP2 post vs pre in R -------------
sample.inf <- merge.sample.inf[merge.sample.inf$Group_WF2 %in% c("pre-R","post-R"),]
df_final <- merge.final[,sample.inf$col]
rownames(sample.inf) <- sample.inf$col

df_final2 <- df_final

df_final2$log2FC <- rowMedian(df_final2[,row.names(sample.inf[sample.inf$Group_WF2 == 'post-R',])]) - rowMedian(df_final2[,row.names(sample.inf[sample.inf$Group_WF2 == 'pre-R',])])

u_test_function <- function(row) {
  group_a <- row[row.names(sample.inf[sample.inf$Group_WF2 == 'pre-R',])]
  group_b <- row[row.names(sample.inf[sample.inf$Group_WF2 == 'post-R',])]
  test_result <- wilcox.test(group_a, group_b)
  return(test_result$p.value)
}

p_values <- apply(df_final2, 1, u_test_function)
df_final2$p_value <- p_values
df_final2$FDR <- p.adjust(p_values, method = "fdr")
df_final2$type <- 'NOT'
df_final2[df_final2$log2FC >= log2(1.2) & df_final2$p_value < 0.05,]$type <- 'UP'
df_final2[df_final2$log2FC <= -log2(1.2) & df_final2$p_value < 0.05,]$type <- 'DOWN'
table(df_final2$type)

# DOWN  NOT   UP 
# 19 2724   32 


df_final2$protein <- row.names(df_final2)
df_final2$protein2 <- ExtractProteinName(df_final2$protein)

write.table(df_final2[,c('log2FC','p_value','FDR','type','protein')], 
            "./Protein/a1_pre_post/preR_postR/a1_preRvsPostR_deg2.txt", sep='\t',
            quote = F, col.names = T, row.names = F)




## read DEP --------------------
library(Vennerable)

a1_preNRvsPostNR_deg2 <- MyReadDelim("./Protein/a1_pre_post/preNR_postNR/a1_preNRvsPostNR_deg2.txt", sep='\t')
# a1_preNRvsPostNR_deg2 <- a1_preNRvsPostNR_deg2[a1_preNRvsPostNR_deg2$type!='NOT',]
a1_preNRvsPostNR_deg2$type2 <- paste0('preNRvspostNR',a1_preNRvsPostNR_deg2$type)
a1_preNRvsPostNR_deg2$protein2 <- ExtractProteinName(a1_preNRvsPostNR_deg2$protein)
a1_preNRvsPostNR_deg2$Group <- 'NR'

a1_preRvsPostR_deg2 <- MyReadDelim("./Protein/a1_pre_post/preR_postR/a1_preRvsPostR_deg2.txt", sep='\t')
# a1_preRvsPostR_deg2 <- a1_preRvsPostR_deg2[a1_preRvsPostR_deg2$type != 'NOT',]
a1_preRvsPostR_deg2$type2 <- paste0('preRvspostR',a1_preRvsPostR_deg2$type)
a1_preRvsPostR_deg2$protein2 <- ExtractProteinName(a1_preRvsPostR_deg2$protein)
a1_preRvsPostR_deg2$Group <- 'R'


## plot volcanno --------------------

mergeAll <- rbind(a1_preNRvsPostNR_deg2,a1_preRvsPostR_deg2)
library(ggplot2)

pdf("./Protein/a1_pre_post/pre_post_DE_protein.2.pdf",6.2,5)
ggplot(mergeAll, aes(log2FC, -log10(p_value), col = type, shape = Group)) +
  geom_point(size = 2) +
  labs(x = 'log2(Post/Pre)', y = '-log10(p_value)') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', linewidth = 0.5) + #横向虚线
  geom_vline(xintercept = log2(1.2), linetype = 'dashed', linewidth = 0.5) +
  geom_vline(xintercept = -log2(1.2), linetype = 'dashed', linewidth = 0.5) +
  scale_color_manual(values = c("#327db7",'grey60',"#b71c2c")) +
  theme_wf2 + 
  theme(
    legend.text = element_text(color = 'black',size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15),
    axis.title = element_text(color = 'black',size = 15),
    axis.ticks = element_line(color = 'black'),
    axis.ticks.length = unit(.25, "cm"))

ggplot(mergeAll, aes(log2FC, -log10(p_value), col = type, shape = Group)) +
  geom_point(size = 2) +
  labs(x = 'log2(Post/Pre)', y = '-log10(p_value)') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', linewidth = 0.5) + #横向虚线
  geom_vline(xintercept = log2(1.2), linetype = 'dashed', linewidth = 0.5) +
  geom_vline(xintercept = -log2(1.2), linetype = 'dashed', linewidth = 0.5) +
  scale_color_manual(values = c("#327db7",'grey60',"#b71c2c")) +
  theme_wf2 + 
  theme(
    legend.text = element_text(color = 'black',size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15),
    axis.title = element_text(color = 'black',size = 15),
    axis.ticks = element_line(color = 'black'),
    axis.ticks.length = unit(.25, "cm")) + 
  ggrepel::geom_text_repel(
    data = subset(mergeAll, -log10(p_value)>2),
    aes(label = protein2),
    size = 4, # colour = "grey50",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
dev.off()


## plot venn ----------
pdf("./Protein/a1_pre_post/pre_post_Venn.pdf",6,5)
plot(Venn(list("NR" = a1_preNRvsPostNR_deg2[a1_preNRvsPostNR_deg2$type!='NOT','protein2'],
               "R" = a1_preRvsPostR_deg2[a1_preRvsPostR_deg2$type != 'NOT','protein2'])),
     doWeight=T)
dev.off()


## run PCA using deg2 --------
deg2 <- mergeAll[mergeAll$type != 'NOT',]

X <- as.data.frame(t(merge.final[deg2$protein,merge.sample.inf$Type]))
X[1:4,1:4]
pca.data <- X
pca.res <- PCA(pca.data, graph = FALSE)

head(pca.res$eig)

library(ggplot2)
library(readxl)
pcaeig <- as.data.frame(pca.res$eig)
pc.percent <- pcaeig$`percentage of variance`
pca.coord <- pca.res$ind$coord
pca.coord <- as.data.frame(pca.coord)

tempSample <- merge.sample.inf
rownames(tempSample) <- tempSample$Type

pca.df <- cbind(x.pos = pca.coord[,1],
                y.pos = pca.coord[,2],
                tempSample[row.names(pca.coord),])
library(RColorBrewer)
color_manual = c(brewer.pal(9,"Set1"))


pdf("./Protein/a1_pre_post/pre_post_PCA_usingDeg2.pdf", 5.5,4)
ggplot(pca.df, aes(x = x.pos, y = y.pos, col = Group_prepost)) +
  geom_point(aes(size = 5)) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  stat_ellipse(aes(col=Group_prepost),alpha=0.2,fill=NA,linewidth=0.7, level = 0.90, geom="polygon") +
  scale_color_manual(values = color_manual) +
  # scale_fill_manual(values = color_manual) +
  theme_wf2 +
  guides(colour = guide_legend(title = "",
                               keywidth = 1,
                               keyheight = 1,
                               override.aes = list(size = 6),
                               order = 1)) +
  guides(shape = guide_legend(title = "",
                              keywidth = 1,
                              keyheight = 1,
                              override.aes = list(size = 6),
                              order = 2)) +
  guides(size = guide_none())

ggplot(pca.df, aes(x = x.pos, y = y.pos, col = Group_prepost, shape = Group_WF)) +
  geom_point(aes(size = 5)) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  # stat_ellipse(aes(col=Group_WF2),alpha=0.2,fill=NA,linewidth=0.7, level = 0.90, geom="polygon") +
  scale_color_manual(values = color_manual) +
  scale_shape_manual(values = c(16, 1))+
  theme_wf2 +
  guides(colour = guide_legend(title = "",
                               keywidth = 1,
                               keyheight = 1,
                               override.aes = list(size = 6),
                               order = 1)) +
  guides(shape = guide_legend(title = "",
                              keywidth = 1,
                              keyheight = 1,
                              override.aes = list(size = 6),
                              order = 2)) +
  guides(size = guide_none())
dev.off()


MyWriteTab(pca.df, "./Protein/a1_pre_post/pre_post_PCA_usingDeg2.csv", sep=',')



## plot heatmap 1: all DEP --------------------
a1_preNRvsPostNR_deg2 <- a1_preNRvsPostNR_deg2[a1_preNRvsPostNR_deg2$type!='NOT',]
a1_preRvsPostR_deg2 <- a1_preRvsPostR_deg2[a1_preRvsPostR_deg2$type != 'NOT',]
vM <- unique(c(a1_preNRvsPostNR_deg2$protein, a1_preRvsPostR_deg2$protein))


orderSample <- merge.sample.inf
orderSample$Group_WF2 <- factor(orderSample$Group_WF2, levels = c("pre-R","post-R",
                                                                  "pre-NR","post-NR"))
orderSample <- orderSample[order(orderSample$Group_WF2,
                                 orderSample$sampleID),]

annotation <- HeatmapAnnotation(
  df = orderSample[,c('Group_WF2'), drop = FALSE],
  col = list(Group_WF2 = c('pre-R' = "#FFA74F", "post-R" = "#008000",
                           'pre-NR' = "#A65628", "post-NR" = "#377EB8"))
)

geneDf <- rbind(a1_preNRvsPostNR_deg2, a1_preRvsPostR_deg2)
geneDf[geneDf$protein %in% c('P49746_THBS3',
                             'P49747_COMP',
                             'Q9NVZ3_NECAP2'),'type2'] <- 'CommonUP'
geneDf$type2 <- factor(geneDf$type2, levels = c("preRvspostRDOWN","preRvspostRUP",
                                                "preNRvspostNRDOWN","preNRvspostNRUP",'CommonUP'))
geneDf <- geneDf[order(geneDf$type2, geneDf$log2FC),]
geneDf <- geneDf[!duplicated(geneDf$protein2),]
rownames(geneDf) <- geneDf$protein


row_ha <- rowAnnotation(
  df = geneDf[,c('type2'), drop = FALSE],
  col = list(type2 = c('preRvspostRUP' = "#FFA74F", "preRvspostRDOWN" = "#008000",
                       'preNRvspostNRUP' = "#A65628", "preNRvspostNRDOWN" = "#377EB8",'CommonUP'="#F781BF"))
)

# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"

genes_to_show <- c("P0DJI8_SAA1", "P0DJI9_SAA2","P02741_CRP","P08637_FCGR3A","P14555_PLA2G2A","Q13217_DNAJC3","P02748_C9", # R up in pre
                   "Q53RD9_FBLN7","P02766_TTR","P22004_BMP6","P11150_LIPC","P02753_RBP4","P02768_ALB","P20810_CAST", "P13987_CD59","Q96HF1_SFRP2",# R up in post
                   'Q9Y4R8_TELO2', # NR up in pre
                   'P07942_LAMB1',"P21980_TGM2","P51884_LUM",'P27169_PON1',"P32119_PRDX2","P30043_BLVRB","P06727_APOA4",# NR up in post
                   'P49746_THBS3','P49747_COMP','Q9NVZ3_NECAP2')

idx <- c(which(geneDf$protein %in% genes_to_show))
labels <- geneDf$protein2[idx]
right_ha <- rowAnnotation(
  mark = anno_mark(at = idx, labels = labels, 
                   which = "row", 
                   side = "right"))

pdf("./Protein/a1_pre_post/pre_post.Heatmap.pdf",8,7)
Heatmap(t(scale(t(merge.final[geneDf$protein,orderSample$col]))),
        name = "z-score",
        show_row_names = F,
        show_column_names = F, cluster_columns = F,cluster_rows = F, 
        column_split = orderSample[,c('Group_WF2')], 
        row_split = geneDf[,c('type2')],
        left_annotation = row_ha,
        bottom_annotation = annotation,
        right_annotation = right_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()

MyWriteTab(t(scale(t(merge.final[geneDf$protein,orderSample$col]))),row.names = T,
           "./Protein/a1_pre_post/pre_post.Heatmap.csv", sep=',')


## plot heatmap 2: Select DEP --------------------

selectProtein <- c('SAA1','SAA2', 'CRP','FCGR3A','PLA2G2A','C9','AGPAT2', # Acute inflammatory response
                   'ACTN4','HIP1', # actin filament-based process
                   'BMP6','CXCL12', # regulation of epithelial cell proliferation
                   'SFRP2','SKP1','NOTUM', # Wnt signaling pathway
                   'COL11A2','COMP',# extracellular matrix organization
                   'APOA4','CDH11','CDH13','LAMB1','CDH23','FAT4' # , # cell-cell adhesion
                   # 'HBA2','PARK7','KPRP', # inflammatory response
                   # 'IHH','PRDX2','PGLYRP2' # negative regulation of hemopoiesis
                   )

a1_preNRvsPostNR_deg2 <- MyReadDelim("./Protein/a1_pre_post/preNR_postNR/a1_preNRvsPostNR_deg2.txt", sep='\t')
a1_preNRvsPostNR_deg2$type2 <- paste0('preNRvspostNR',a1_preNRvsPostNR_deg2$type)
a1_preNRvsPostNR_deg2$protein2 <- ExtractProteinName(a1_preNRvsPostNR_deg2$protein)
a1_preNRvsPostNR_deg2$Group <- 'NR'

tt <- a1_preNRvsPostNR_deg2[a1_preNRvsPostNR_deg2$protein2 %in% selectProtein,]
rownames(tt) <- tt$protein2

merge.final$preR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('pre-R'),'col']], 1, mean)
merge.final$postR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('post-R'),'col']], 1, mean)
merge.final$preNR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('pre-NR'),'col']], 1, mean)
merge.final$postNR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('post-NR'),'col']], 1, mean)

plot_res <- merge.final[tt$protein,c('preR','postR','preNR','postNR')]
rownames(plot_res) <- tt$protein2
plot_res <- plot_res[selectProtein,]


pdf("./Protein/a1_pre_post/pre_post.Heatmap.2.pdf",5,8)
Heatmap(t(scale(t(plot_res))),
        name = "z-score",
        show_row_names = T,
        show_column_names = F, cluster_columns = F,cluster_rows = F, 
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()



## select the example to show -------
merge.sample.inf$Group_WF2 <- factor(merge.sample.inf$Group_WF2,
                                     levels = c("pre-R","post-R","pre-NR","post-NR" ))

merge.sample.inf <- merge.sample.inf[order(merge.sample.inf$Group_WF2,
                                           merge.sample.inf$patientID2),]

plot2show <- merge.sample.inf
plot2show$P02741_CRP <- unlist(merge.final['P02741_CRP', plot2show$col])
plot2show$P0DJI8_SAA1 <- unlist(merge.final['P0DJI8_SAA1', plot2show$col])
plot2show$P0DJI9_SAA2 <- unlist(merge.final['P0DJI9_SAA2', plot2show$col])
plot2show$P08637_FCGR3A <- unlist(merge.final['P08637_FCGR3A', plot2show$col])
plot2show$P22004_BMP6 <- unlist(merge.final['P22004_BMP6', plot2show$col])
plot2show$P48061_CXCL12 <- unlist(merge.final['P48061_CXCL12', plot2show$col])
plot2show$Q96HF1_SFRP2 <- unlist(merge.final['Q96HF1_SFRP2', plot2show$col])
plot2show$P07942_LAMB1 <- unlist(merge.final['P07942_LAMB1', plot2show$col])
plot2show$Q5T749_KPRP <- unlist(merge.final['Q5T749_KPRP', plot2show$col])
plot2show$P55290_CDH13 <- unlist(merge.final['P55290_CDH13', plot2show$col])
plot2show$P55287_CDH11 <- unlist(merge.final['P55287_CDH11', plot2show$col])
plot2show$Q9H257_CARD9 <- unlist(merge.final['Q9H257_CARD9', plot2show$col])
plot2show$P21980_TGM2 <- unlist(merge.final['P21980_TGM2', plot2show$col])


MyWriteTab(plot2show, "./Protein/a1_pre_post/plot2show_1.csv", sep=',')

# [1] "O00291_HIP1"    "O15120_AGPAT2"  "O43707_ACTN4"   "P02741_CRP"     "P02748_C9"      "P06727_APOA4"   "P07942_LAMB1"   "P08637_FCGR3A" 
# [9] "P0DJI8_SAA1"    "P0DJI9_SAA2"    "P13942_COL11A2" "P14555_PLA2G2A" "P22004_BMP6"    "P32119_PRDX2"   "P48061_CXCL12"  "P49747_COMP"   
# [17] "P55287_CDH11"   "P55290_CDH13"   "P63208_SKP1"    "P69905_HBA2"    "Q14623_IHH"     "Q5T749_KPRP"    "Q6P988_NOTUM"   "Q6V0I7_FAT4"   
# [25] "Q96HF1_SFRP2"   "Q96PD5_PGLYRP2" "Q99497_PARK7"   "Q9H251_CDH23"  

MyRunPlot2show <- function(plot2show, colName){
  name <- unlist(strsplit(colName,'_'))[2]
  p1<-ggplot(plot2show, aes(x=Group_WF2, y=plot2show[,colName],group = patientID2)) +
     geom_line(color = "black", linewidth = 0.5) +          # 连接线
     geom_point(color = "black", size = 2) +
     ggpubr::stat_compare_means(
       comparisons = list(c("pre-R", "post-R"), c("pre-NR", "post-NR")),
       method = "wilcox.test", 
       label = "p.format",
       label.x.npc = "center",
       label.y.npc = "top",
       vjust = -0.5
     ) + 
     labs(x = "", y = expression(log[2](intensity)), title = name) +
     theme_wf +
     theme(plot.title = element_text(size = 20, colour = "black", hjust = 0.65, face = "bold"),
           plot.title.position = "plot",
           axis.text.x.bottom = element_text(angle = 30, hjust = 1))
  return(p1)
}

p1 <- MyRunPlot2show(plot2show,'P02741_CRP')
p2 <- MyRunPlot2show(plot2show,'P0DJI8_SAA1')
p3 <- MyRunPlot2show(plot2show,'P0DJI9_SAA2')
p4 <- MyRunPlot2show(plot2show,'P08637_FCGR3A')
p5 <- MyRunPlot2show(plot2show,'P22004_BMP6')
p6 <- MyRunPlot2show(plot2show,'P48061_CXCL12')
p7 <- MyRunPlot2show(plot2show,'Q96HF1_SFRP2')
p8 <- MyRunPlot2show(plot2show,'P07942_LAMB1')
p9 <- MyRunPlot2show(plot2show,'Q5T749_KPRP')
p10 <- MyRunPlot2show(plot2show,'P55287_CDH11')
p11 <- MyRunPlot2show(plot2show,'P55290_CDH13')
p12 <- MyRunPlot2show(plot2show,'Q9H257_CARD9')
p13 <- MyRunPlot2show(plot2show,'P21980_TGM2')

pdf("./Protein/a1_pre_post/pre_post_plot2show.pdf",11.5,16)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13
dev.off()

