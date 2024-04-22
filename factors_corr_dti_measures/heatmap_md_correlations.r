

#使用ComplexHeatmap包绘制热图；
library(BiocManager)
library(ggplot2)
#安装相关R包；

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)
library(circlize)


rm(list=ls())

setwd('ukb/wmh')

p1 <- read.csv("ukb/wmh/wmh_corr_structural_mri/results/MD_pvalues.csv",row.names = 1) 
p2=p1*27

r2 <- read.csv("ukb/wmh/wmh_corr_structural_mri/results/MD_rvalues.csv",row.names = 1)


names(r2) <- gsub("_", " ", names(r2))
names(p2) <- gsub("_", " ", names(p2))


r2<-t(r2)
p2<-t(p2)



p2[p2>=0 & p2 < 0.001] <- "***"
p2[p2>=0.001 & p2 < 0.01] <- "**"
p2[p2>=0.01 & p2 < 0.05] <- "*"
p2[p2>=0.05 ] <- ""

#预览替换后的矩阵；
#p2[1:9,1:9]

range(r2)

col_fun1 = colorRamp2(c(-0.14, 0, 0.1), c("#0f86a9", "white", "#FC8452"))


cellwidth = 0.6
cellheight = 0.4
cn = dim(r2)[2]
rn = dim(r2)[1]
w=cellwidth*cn
h=cellheight*rn




tiff( 
  filename = "ukb/wmh/wmh_corr_structural_mri/figs/MD_heatmap.tiff", # 文件名称
  width = 6000,           # 宽
  height = 3000,          # 高
  units = "px",          # 单位
  bg = "white",          # 背景颜色
  res = 300)              # 分辨率


#绘制热图显示显著性星号标记；
h_fig<-Heatmap(r2,name ="r", col = col_fun1,
        #格子大小设置；
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        border_gp = gpar(col = "#0f86a9",lty = 2,lwd = 1.2),
        #聚类树样式设置；
        column_dend_height = unit(1.5, "cm"),
        row_dend_width = unit(1.5, "cm"),
        column_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        row_dend_gp = gpar(col = "#0f86a9",lwd = 1.4),
        
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        
        #设置聚类gap数量和大小；
        #row_split = 0, column_split = 0,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        #行列标签文字样式设置；
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 12),
        row_names_side = c("left"),
        
        
        
        #图例样式设置；
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.3, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 5)),
        #显示星号标记设置；
        #vjust垂直微调星号的位置；
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(p2[i, j], x, y, vjust = 0.6,
                    gp = gpar(fontsize = 8,col="black"))
        })

print(h_fig)

dev.off()



