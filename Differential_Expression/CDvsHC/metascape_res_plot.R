
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251124_PrMetabolism_XHOM_WF/")

library(readxl)
library(ggplot2)


## read metascape result ---------------
df <- read_excel('./metascape_result.xlsx', sheet = "Enrichment")
df <- as.data.frame(df)

df_summary <- df[grep('Summary',df$GroupID),]


merge_res_plot <- df_summary
merge_res_plot$X_axis <- -1 * (merge_res_plot$LogP)
tempColor <- colorRampPalette(c("#E57539",'white'))(30)

pdf("./metascape_result.pdf",
    5,9)
ggplot(merge_res_plot, aes(reorder(Description, X_axis), X_axis, fill = X_axis))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=20.5))+
  geom_text(data = merge_res_plot[which(merge_res_plot$X_axis>0),],aes(x=Description, y=-0.1, label=Description),
            hjust=1, size=5)+
  geom_text(data = merge_res_plot[which(merge_res_plot$X_axis<0),],aes(x=Description, y=0.1, label=Description),
            hjust=0, size=5)+
  # scale_fill_manual(values = c("#008000","#FFA74F"))+
  scale_fill_gradient(low=tempColor[10], high=tempColor[1]) +
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  labs(x='', y='-log10(P)')
dev.off()



