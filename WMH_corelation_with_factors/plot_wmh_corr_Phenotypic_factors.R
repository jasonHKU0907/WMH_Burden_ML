#require(data.table)
library(ggplot2)
library(ggrepel)
require(dplyr)
require(tidyr)
require(RColorBrewer)

rm(list=ls())

setwd('ukb/wmh')

dfres <- read.csv('ukb/wmh/Cross_section_correlation/results/wmh_corr_Phenotypic_results.csv')

dfres$Category <- dfres$category
dfres$bonferroni_p <- p.adjust(dfres$raw_p,method="bonferroni",n=length(dfres$raw_p))
dfres$P <- dfres$raw_p

dfres$plotx=seq(1,3*nrow(dfres)+10,3)[1:nrow(dfres)]

ss=brewer.pal(12, "Paired")
coldata=data.frame(category=unique(dfres$category),cols=ss[c(7,9,3,4,6,5,2,8,1,10,11,12)])
coldata=coldata[order(coldata$category),]


#thr1<-0.05/241
thr <- 0.05
# 为dfres中的每个行设置一个名为size的变量，初始值为4，当dfres$x小于阈值thr时，将size设置为6  
dfres$size <- 1.5
dfres$size[dfres$bonferroni_p<thr] <- 2.5


thr_r=dfres$rvalue[dfres$bonferroni_p<thr]

r_upper=min(thr_r[thr_r>0])
r_lower=max(thr_r[thr_r<0])





# 使用dplyr包中的group_by和summarize函数来计算每个category的plotx的均值，存储到X_axis中  
X_axis <- dfres %>% group_by(Category) %>% summarize(center=( max(plotx) +min(plotx) ) / 2 )
X_axis$Category <- unique(dfres$category)



x_max <- max(dfres$plotx)
y_max <-  round(max(thr_r)+0.02,2)
y_min <-  round(min(thr_r)-0.02,2)

p<-ggplot(dfres, aes(x=plotx, y=rvalue)) + 
  geom_hline(yintercept=r_upper, color=ss[6], size=0.4,linetype="longdash")+
  geom_hline(yintercept=r_lower, color=ss[6], size=0.4,linetype="longdash")+
  geom_hline(yintercept=0, color="black", size=0.55)+
  # geom_hline(yintercept=-log10(thr), color=ss[6], size=0.4,linetype="longdash")+
  geom_point(aes(colour = Category),size=dfres$size) + 
  #shape=Correlation_directions
  labs(color="Category", x="", y=expression('R value')) +
 
  ggrepel::geom_label_repel(data=. %>% mutate(label = ifelse(bonferroni_p < 0.05, as.character(Phenotypic_bl_names), "")), aes(label=label), 
                            size=4, box.padding = unit(0.2, "lines"),
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 550),
                            min.segment.length = 0)+
  
  
  
  scale_colour_manual(values = as.vector(coldata$cols))+
 
  theme_classic() + 
  theme(axis.title = element_text(face="bold",size=16),
        axis.line.y = element_line(color = "black",size=0.4),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black",size=10),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        #legend.key.size = unit(0.3, "inches")
        #panel.grid.minor=element_blank())
  )
print(p)

ggsave('Figs/wmh_corr_Phenotypic_regressed_disase_full_results.jpg', p, dpi = 500, width = 450, height = 180, units = "mm")
