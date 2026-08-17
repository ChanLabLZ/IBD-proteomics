
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

## read protein matrix ------------- 
library(readxl)
library(FactoMineR)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


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
pre_a2 <- read_excel("../Sample_Pre_Post.xlsx")
pre_a2 <- as.data.frame(pre_a2)
pre_a2 <- pre_a2[pre_a2$Group=='a2',]


sample.inf.temp <- sample.inf[sample.inf$sampleID %in% pre_a2$sampleID,]
sample.inf <- sample.inf.temp[sample.inf.temp$Description == '/',]
sample.inf.temp <- sample.inf.temp[sample.inf.temp$Description != '/',]

# sample.inf <- sample.inf[!sample.inf$col %in% c('No361_rep','No472_rep'),]
# sample.inf.temp <- sample.inf.temp[!sample.inf.temp$col %in% c('No322_rep'),]

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
# 28 2742    5

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
            "./Protein/a2_pre_post/preNR_postNR/a2_preNRvsPostNR_deg2.txt", sep='\t',
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
# 21 2704   50 

df_final2$protein <- row.names(df_final2)
df_final2$protein2 <- ExtractProteinName(df_final2$protein)

write.table(df_final2[,c('log2FC','p_value','FDR','type','protein')], 
            "./Protein/a2_pre_post/preR_postR/a2_preRvsPostR_deg2.txt", sep='\t',
            quote = F, col.names = T, row.names = F)

 
## read DEP --------------------
library(Vennerable)

a2_preNRvsPostNR_deg2 <- MyReadDelim("./Protein/a2_pre_post/preNR_postNR/a2_preNRvsPostNR_deg2.txt", sep='\t')
a2_preNRvsPostNR_deg2$type2 <- paste0('preNRvspostNR',a2_preNRvsPostNR_deg2$type)
a2_preNRvsPostNR_deg2$protein2 <- ExtractProteinName(a2_preNRvsPostNR_deg2$protein)
a2_preNRvsPostNR_deg2$Group <- 'NR'

a2_preRvsPostR_deg2 <- MyReadDelim("./Protein/a2_pre_post/preR_postR/a2_preRvsPostR_deg2.txt", sep='\t')
a2_preRvsPostR_deg2$type2 <- paste0('preRvspostR',a2_preRvsPostR_deg2$type)
a2_preRvsPostR_deg2$protein2 <- ExtractProteinName(a2_preRvsPostR_deg2$protein)
a2_preRvsPostR_deg2$Group <- 'R'


## plot volcanno --------------------
mergeAll <- rbind(a2_preNRvsPostNR_deg2,a2_preRvsPostR_deg2)
library(ggplot2)

pdf("./Protein/a2_pre_post/pre_post_DE_protein.2.pdf",6.2,5)
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


## plot venn --------------
pdf("./Protein/a2_pre_post/pre_post_Venn.pdf",6,5)
plot(Venn(list("NR" = a2_preNRvsPostNR_deg2[a2_preNRvsPostNR_deg2$type!='NOT','protein2'],
               "R" = a2_preRvsPostR_deg2[a2_preRvsPostR_deg2$type != 'NOT','protein2'])),
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


pdf("./Protein/a2_pre_post/pre_post_PCA_usingDeg2.pdf", 5.5,4)
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

MyWriteTab(pca.df, "./Protein/a2_pre_post/pre_post_PCA_usingDeg2.csv", sep=',')


## plot heatmap 1: all DEP --------------------

a2_preNRvsPostNR_deg2 <- a2_preNRvsPostNR_deg2[a2_preNRvsPostNR_deg2$type!='NOT',]
a2_preRvsPostR_deg2 <- a2_preRvsPostR_deg2[a2_preRvsPostR_deg2$type != 'NOT',]
vM <- unique(c(a2_preNRvsPostNR_deg2$protein, a2_preRvsPostR_deg2$protein))

library(ComplexHeatmap)
library(circlize)

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

geneDf <- rbind(a2_preNRvsPostNR_deg2, a2_preRvsPostR_deg2)
geneDf[geneDf$protein %in% c('Q8IWZ3_ANKHD1'),'type2'] <- 'DOWNinR_UPinNR'
geneDf$type2 <- factor(geneDf$type2, levels = c("preRvspostRDOWN","preRvspostRUP",
                                                "preNRvspostNRDOWN","preNRvspostNRUP",'DOWNinR_UPinNR'))
geneDf <- geneDf[order(geneDf$type2, geneDf$log2FC),]
geneDf <- geneDf[!duplicated(geneDf$protein2),]
rownames(geneDf) <- geneDf$protein


row_ha <- rowAnnotation(
  df = geneDf[,c('type2'), drop = FALSE],
  col = list(type2 = c('preRvspostRUP' = "#FFA74F", "preRvspostRDOWN" = "#008000",
                       'preNRvspostNRUP' = "#A65628", "preNRvspostNRDOWN" = "#377EB8",'DOWNinR_UPinNR'="#F781BF"))
)


# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
genes_to_show <- c("O43242_PSMD3","O75165_DNAJC13","Q8NF37_LPCAT1","P35580_MYH10","P42704_LRPPRC", # R up in pre
                   "P02753_RBP4","P02768_ALB","Q96HF1_SFRP2","P02766_TTR","Q53RD9_FBLN7","P00325_ADH1B","P08319_ADH4", # R up in post
                   "P0DJI9_SAA2","P02741_CRP","P14555_PLA2G2A","P26583_HMGB2","P05109_S100A8","P06702_S100A9","P50995_ANXA11", # NR up in pre
                   "P02686_MBP",# NR up in post
                   'Q8IWZ3_ANKHD1' # DOWNinR_UPinNR
)

idx <- c(which(geneDf$protein %in% genes_to_show))
labels <- geneDf$protein2[idx]
right_ha <- rowAnnotation(
  mark = anno_mark(at = idx, labels = labels, 
                   which = "row", 
                   side = "right"))


pdf("./Protein/a2_pre_post/pre_post.Heatmap.pdf",8,7)
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
           "./Protein/a2_pre_post/pre_post.Heatmap.csv", sep=',')




## plot heatmap 2: Select DEP --------------------

selectProtein <- c("PSMD3","DNAJC13","DYNC1LI1","LPCAT1","NF2", # Neutrophil degranulation
                   'CRP','SAA2', # 'HMGB2','S100A8','S100A9','NAMPT', # inflammatory response
                   'AZU1','BPI','LTF', # defense response to bacterium
                   'SFRP2','SFRP4',"MYOC", # non-canonical Wnt signaling pathway
                   'AHSG','CHAD','RBP4','SPP2','THBS3', # skeletal system development
                   "APOA1","MASP1","FCN3" # positive regulation of endocytosis
)

a2_preNRvsPostNR_deg2 <- MyReadDelim("./Protein/a2_pre_post/preNR_postNR/a2_preNRvsPostNR_deg2.txt", sep='\t')
a2_preNRvsPostNR_deg2$type2 <- paste0('preNRvspostNR',a2_preNRvsPostNR_deg2$type)
a2_preNRvsPostNR_deg2$protein2 <- ExtractProteinName(a2_preNRvsPostNR_deg2$protein)
a2_preNRvsPostNR_deg2$Group <- 'NR'

tt <- a2_preNRvsPostNR_deg2[a2_preNRvsPostNR_deg2$protein2 %in% selectProtein,]
rownames(tt) <- tt$protein2

merge.final$preR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('pre-R'),'col']], 1, mean)
merge.final$postR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('post-R'),'col']], 1, mean)
merge.final$preNR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('pre-NR'),'col']], 1, mean)
merge.final$postNR <- apply(merge.final[,merge.sample.inf[merge.sample.inf$Group_WF2 %in% c('post-NR'),'col']], 1, mean)

plot_res <- merge.final[tt$protein,c('preR','postR','preNR','postNR')]
rownames(plot_res) <- tt$protein2
plot_res <- plot_res[selectProtein,]


pdf("./Protein/a2_pre_post/pre_post.Heatmap.2.pdf",5,8)
Heatmap(t(scale(t(plot_res))),
        name = "z-score",
        show_row_names = T,
        show_column_names = F, cluster_columns = F,cluster_rows = F, 
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
        col = colorRamp2(c(-1.5, 0, 1.5), c("#327db7", "white", "#b71c2c")))
dev.off()

MyWriteTab(t(scale(t(plot_res))), row.names = T,
           './Protein/a2_pre_post/pre_post.Heatmap.2.csv', sep=',')


## select the example to show -------
merge.sample.inf$Group_WF2 <- factor(merge.sample.inf$Group_WF2,
                                     levels = c("pre-R","post-R","pre-NR","post-NR" ))

merge.sample.inf <- merge.sample.inf[order(merge.sample.inf$Group_WF2,
                                           merge.sample.inf$patientID2),]

plot2show <- merge.sample.inf
plot2show$P02741_CRP <- unlist(merge.final['P02741_CRP', plot2show$col])
plot2show$P02788_LTF <- unlist(merge.final['P02788_LTF', plot2show$col])
plot2show$P0DJI9_SAA2 <- unlist(merge.final['P0DJI9_SAA2', plot2show$col])

plot2show$O43242_PSMD3 <- unlist(merge.final['O43242_PSMD3', plot2show$col])
plot2show$O15335_CHAD <- unlist(merge.final['O15335_CHAD', plot2show$col])
plot2show$P02765_AHSG <- unlist(merge.final['P02765_AHSG', plot2show$col])
plot2show$O75636_FCN3 <- unlist(merge.final['O75636_FCN3', plot2show$col])

plot2show$P49746_THBS3 <- unlist(merge.final['P49746_THBS3', plot2show$col])
plot2show$P02647_APOA1 <- unlist(merge.final['P02647_APOA1', plot2show$col])

plot2show$Q96HF1_SFRP2 <- unlist(merge.final['Q96HF1_SFRP2', plot2show$col])
plot2show$Q99972_MYOC <- unlist(merge.final['Q99972_MYOC', plot2show$col])
plot2show$Q6FHJ7_SFRP4 <- unlist(merge.final['Q6FHJ7_SFRP4', plot2show$col])

plot2show$P35240_NF2 <- unlist(merge.final['P35240_NF2', plot2show$col])
plot2show$Q14161_GIT2 <- unlist(merge.final['Q14161_GIT2', plot2show$col])
plot2show$Q96PD5_PGLYRP2 <- unlist(merge.final['Q96PD5_PGLYRP2', plot2show$col])
plot2show$P26583_HMGB2 <- unlist(merge.final['P26583_HMGB2', plot2show$col])




MyWriteTab(plot2show, "./Protein/a2_pre_post/plot2show.csv", sep=',')

# [1] "O15335_CHAD"     "O43242_PSMD3"    "O75165_DNAJC13"  "O75636_FCN3"     "P02647_APOA1"    "P02741_CRP"      "P02753_RBP4"    
# [8] "P02765_AHSG"     "P02788_LTF"      "P0DJI9_SAA2"     "P17213_BPI"      "P20160_AZU1"     "P35240_NF2"      "P48740_MASP1"   
# [15] "P49746_THBS3"    "Q13103_SPP2"     "Q6FHJ7_SFRP4"    "Q8NF37_LPCAT1"   "Q96HF1_SFRP2"    "Q99972_MYOC"     "Q9Y6G9_DYNC1LI1"


MyRunPlot2show <- function(plot2show, colName){
  name <- unlist(strsplit(colName,'_'))[2]
  p1<-ggplot(plot2show, aes(x=Group_WF2, y=plot2show[,colName],group = patientID2)) +
    geom_line(color = "black", linewidth = 0.5) +          # 连接线
    geom_point(color = "black", size = 2) +
    ggpubr::stat_compare_means(
      comparisons = list(c("pre-R", "post-R"), c("pre-NR", "post-NR")),
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

p1 <- MyRunPlot2show(plot2show,'P02741_CRP')
p2 <- MyRunPlot2show(plot2show,'P0DJI9_SAA2')
p3 <- MyRunPlot2show(plot2show,'O43242_PSMD3')
p4 <- MyRunPlot2show(plot2show,'O15335_CHAD')
p5 <- MyRunPlot2show(plot2show,'P02765_AHSG')
p6 <- MyRunPlot2show(plot2show,'O75636_FCN3')
p7 <- MyRunPlot2show(plot2show,'P49746_THBS3')
p8 <- MyRunPlot2show(plot2show,'P02647_APOA1')
p9 <- MyRunPlot2show(plot2show,'Q99972_MYOC')
p10 <- MyRunPlot2show(plot2show,'Q96HF1_SFRP2')
p11 <- MyRunPlot2show(plot2show,'Q6FHJ7_SFRP4')
p12 <- MyRunPlot2show(plot2show,'P35240_NF2')
p13 <- MyRunPlot2show(plot2show,'Q14161_GIT2')
p14 <- MyRunPlot2show(plot2show,'Q96PD5_PGLYRP2')
p15 <- MyRunPlot2show(plot2show,'P26583_HMGB2')
p16 <- MyRunPlot2show(plot2show,'P02788_LTF')

pdf("./Protein/a2_pre_post/pre_post_plot2show.pdf",11.5,16)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16
dev.off()


