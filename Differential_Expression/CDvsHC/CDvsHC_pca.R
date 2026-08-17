
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(readxl)
library(FactoMineR)
library(ggplot2)
library(RColorBrewer)

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
rm(pp)
rownames(merged_df2) <- merged_df2$col

pp <- rbind(merged_df, merged_df2)
merged_df <- pp
rm(pp)

sample.inf <- merged_df
rm(merged_df)
rownames(sample.inf) <- sample.inf$col


## run PCA in all samples --------
# log2
pca.data <- t(df_final)
pca.res <- PCA(pca.data, graph = FALSE)

head(pca.res$eig)

pcaeig <- as.data.frame(pca.res$eig)
pc.percent <- pcaeig$`percentage of variance`
pca.coord <- pca.res$ind$coord
pca.coord <- as.data.frame(pca.coord)

pca.df <- cbind(x.pos = pca.coord[,1],
                y.pos = pca.coord[,2],
                sample.inf[row.names(pca.coord),])

color_manual = c(brewer.pal(9,"Set1"))
color_manual <- color_manual[c(1:5,7)]

pca.df$Group_prepost_drug <- pca.df$Description
pca.df[pca.df$Group_prepost %in% 'pre','Group_prepost_drug'] <- pca.df[pca.df$Group_prepost %in% 'pre','followingDrug']

pca.df$Type2 <- factor(pca.df$Type2, levels = c("HC","CD"))

pdf("./Protein/HC_allIBD/PCA_1.pdf", 7,5.2)
ggplot(pca.df, aes(x = x.pos, y = y.pos, col = Type2)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  stat_ellipse(aes(col=Type2),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
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
dev.off()

MyWriteTab(pca.df,"./PCA_1.csv", sep=',')


pdf("./PCA.pdf", 7,5.2)
ggplot(pca.df, aes(x = x.pos, y = y.pos, col = Group2)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  stat_ellipse(aes(col=Group2),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
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

ggplot(pca.df[pca.df$Group %in% c('NR','R'),], 
       aes(x = x.pos,y = y.pos, col = Group)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  stat_ellipse(aes(col=Group),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
  guides(colour=guide_legend(title=NULL)) +
  scale_color_manual(values = color_manual) +
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


ggplot(pca.df[pca.df$Group %in% c('NR','R'),], 
       aes(x = x.pos,y = y.pos, col = Group2)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  stat_ellipse(aes(col=Group2),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
  guides(colour=guide_legend(title=NULL)) +
  scale_color_manual(values = color_manual) +
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

ggplot(pca.df[pca.df$Group %in% c('HC','UNT'),], 
       aes(x = x.pos,y = y.pos, col = Group)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  stat_ellipse(aes(col=Group),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
  scale_color_manual(values = color_manual) +
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

ggplot(pca.df[!pca.df$Group %in% c('HC','UNT'),], 
       aes(x = x.pos,y = y.pos, col = Description)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  stat_ellipse(aes(col=Description),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
  scale_color_manual(values = color_manual) +
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

ggplot(pca.df[!is.na(pca.df$Group_prepost),], 
       aes(x = x.pos,y = y.pos, col = Group_prepost_drug, shape = Group_prepost)) +
  geom_point(size = 5) +
  xlab(paste("PC",1,"(",round(pc.percent[1],1),"%)",sep = "")) +
  ylab(paste("PC",2,"(",round(pc.percent[2],1),"%)",sep = "")) +
  guides(colour=guide_legend(title=NULL)) +
  #stat_ellipse(aes(col=Group_prepost),alpha=0.2,fill=NA,linewidth=0.7, level = 0.95, geom="polygon") + 
  scale_color_manual(values = color_manual) +
  scale_shape_manual(values = c(1, 16)) +
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






