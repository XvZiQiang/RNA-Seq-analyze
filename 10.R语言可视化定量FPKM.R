# 在R中可视化定量结果
install.packages("reshape2")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggh4x")
install.packages("ggrepel")
install.packages("ggsci")
install.packages("RColorBrewer")








#加载R包
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggrepel)
library(ggsci)
library(RColorBrewer)

# 检查所有包是否成功加载
sessionInfo()

# 清空缓存
rm(list = ls())

# 定义通用保存函数
save_plot <- function(plot_object, file_prefix) {
  # 保存为PDF
  ggsave(paste0(file_prefix, ".pdf"), plot_object, width = 9, height = 6)
  # 保存为PNG
  png_filename <- paste0(file_prefix, ".png")
  png(png_filename, width = 3000, height = 2000, res = 300)
  print(plot_object)
  dev.off()
}

# 读取样品表达矩阵
count_data <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/gene.count.FPKM.csv",header = T,sep = ",",
                         row.names = 1,quote = "",check.names = F)
head(count_data)

# 设定分组
sample_data <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/sample_desc.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(sample_data) <- c("group","sample")
sample_info <- data.frame(row.names = sample_data$sample,group = sample_data$group)
head(sample_info)

# 取对数
log_count_data <- log10(count_data)

# 宽转长
log_count_data_long <- melt(log_count_data, variable.name = "sample", value.name = "expression")

# 添加分组信息
boxplot_data <- merge(log_count_data_long, sample_data, by = "sample")

# 去除异常值
boxplot_data <- boxplot_data %>%
  filter(!is.na(expression) & !is.infinite(expression))

# 分组表达水平箱线图
boxplot <-
  ggplot(boxplot_data, aes(x = sample, y = expression, color = group)) +
  geom_boxplot(outlier.shape = NA, lineend = "round") +
  scale_y_continuous(
    limits = c(-2,4),
    expand = c(0,0)
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "The boxplot of gene expression level distribution", 
       x = NULL, y = expression("Expression Level (Log"[10]*"FPKM)"),
       color = NULL) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 18) +
  theme(panel.border = element_rect(fill=NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5, lineend = "round")) +
  theme(legend.key = element_rect(color = NA, fill = NA, linewidth = 0.4),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank()) +
  theme(aspect.ratio = 2/3)
boxplot
save_plot(boxplot, "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.expression.boxplot")









# 表达量分布堆积柱状图
# 把表达量矩阵 count_data（FPKM 表）复制为 stack_data
stack_data <- count_data
# 计算基因总数（行数）
stack_sum <- nrow(stack_data)

# 划分区间
level_bin <- c(0, 1, 10, 100, 1000, Inf)
level_label <- c("0~1", "1~10", "10~100", "100~1000", ">1000")

# 统计区间数量
# cut() 把表达量分到定义好的区间；
# table() 统计每个区间内的数量；
# apply(..., 2, ...) 表示对每一列（样本）执行
level_count <- apply(stack_data, 2, 
                     function(col) table(cut(col, breaks = level_bin, labels = level_label, right = FALSE)))
level_count <- as.data.frame(level_count)
level_count <- cbind(Level = factor(rownames(level_count),levels = level_label), level_count)

# 保存各表达量区间的基因数量
write.csv(level_count,"C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.expression.intervals.csv",row.names = F)

# 宽转长
level_melt <- melt(level_count, id.vars = 'Level', variable.name = 'Sample', value.name = 'Count')

# 计算占比
level_melt <- level_melt  %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count) * 100)
# 绘制堆积柱状图，每个样本一个柱子，不同表达量区间堆叠
stack <- 
  ggplot(level_melt, aes(x = Sample, y = Proportion, fill = Level)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75, color = NA) +
# 设置 Y 轴为百分比，刻度从 0–100，每隔 20
  scale_y_continuous(
    expand = c(0,0),
    name = "Percentage(%)",
    breaks = seq(0, 100, by = 20)
  ) +
# 允许标签超出绘图区
  coord_cartesian(clip = "off") +
# 使用ggsci的BMJ风格配色
  scale_fill_bmj() +
# 图标题和图例标签设置
  labs(x = NULL, fill = 'Expression Level') +
  ggtitle("The Stacked bar chart of gene expression level distribution") +
# 设置字体大小、坐标轴文字角度与颜色
# hjust = 0.5 → 水平居中（0=左对齐，1=右对齐）  
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black")) +
# 去掉背景网格线
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
# 控制坐标轴线和刻度线的样式（这里 Y 轴有刻度，X 轴无）
  theme(axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.length.y = unit(0.25,"cm"),
        axis.ticks.x = element_blank()) +
# guides() 控制轴显示
  guides(y = "axis_truncated") +
# 图例框样式  
  theme(legend.key = element_rect(color = NA, fill = NA, linewidth = 0.4),  # 图例条目边框
        legend.key.size = unit(0.5, "cm")
  ) +
# 图高宽比设为 1:2（较扁平）
  theme(aspect.ratio = 1/2)
# 绘图并保存为 PDF/PNG 文件
stack
save_plot(stack, "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.expression.stackchart")










# PCA计算
# pca_count 保存经过过滤后的矩阵
# 去掉表达量全为 0 的基因（没有信息的行）rowSums(count_data) > 0
pca_count = count_data[which(rowSums(count_data) >0),]
# 检测每个基因是否在所有样品中表达量相同,如果某行（基因）表达量完全一致，PCA 对它没有信息量，因此剔除,保存为pca_count文件
if ( nrow(pca_count[which(apply(pca_count,1,function(x) max(x)==min(x))),]) > 0 ){
  pca_count <- pca_count[-which(apply(pca_count,1,function(x) max(x)==min(x))),]
}
# scale(pca_count)：对pca_count里每个基因做 Z-score 标准化（均值 0，方差 1），eigen()：计算特征值和特征向量
eigen <- scale(pca_count) %>%
  cor() %>% 
  eigen()
eigen$values # 提取特征值
eigen$vectors # 提取特征向量
scale(pca_count) %*% eigen$vectors %>% head() #计算主成分得分


# t(pca_count)：转置矩阵，使样品为行，基因为列（prcomp 默认行是样品）。
# rank. = 4：只计算前 4 个主成分。
# retx = T：返回每个样品在主成分空间的得分（PC Scores）
pca_res <- prcomp(t(pca_count),
                  rank. = 4,
                  retx = T) # 计算PCA
summary(pca_res) # 查看PCA_res结果，输出每个主成分的标准差、方差贡献率（解释方差百分比）




pca_data <- as.data.frame(pca_res$x) # 提取PC主成分
pca_data$sample <- row.names(pca_data)# 添加一列 sample 作为样品名，方便后续绘图

# 用 merge 将 PCA 得分和分组信息(sample_info)对应起来。
# 最后把行名设置为样品名，便于索引
sample_info$sample <- row.names(sample_info)
pca_data <- merge(pca_data, sample_info, by = "sample", all = TRUE)
row.names(pca_data) <- pca_data$sample

# 保存 PCA 结果
pca_out <- pca_data
write.csv(pca_out,"C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.PCA.csv",row.names = F)

summ <- summary(pca_res) # 计算权重

# summ$importance[2,1] → PC1 方差贡献率;summ$importance[2,2] → PC2 方差贡献率
# paste0() 拼接字符串，用于坐标轴标签：
# x 轴 → PC1 (xx%);y 轴 → PC2 (xx%)
xlab <- paste0("PC #1 (", round(summ$importance[2, 1] * 100, 2),"%)")
ylab <- paste0("PC #2 (", round(summ$importance[2, 2] * 100, 2),"%)")


# PCA图
pca <-
  ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 5, stroke = 1,alpha = 1,
             shape = 21, show.legend = T) +
  coord_cartesian(clip = "off") +
  labs(x = xlab, y = ylab, color = NULL) +
  geom_text_repel(aes(label = sample, color = group), size = 4, show.legend = F) +
  scale_x_continuous(expand = c(0,0), limits = c(-5000,5000)) +
  scale_y_continuous(expand = c(0,0), limits = c(-2000,2000)) +
  scale_color_lancet() +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5, lineend = "round")) +
  theme(panel.border = element_rect(fill=NA, linewidth = 0.5)) +
  theme(legend.key = element_rect(color = NA, fill = NA, linewidth = 0.4),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank()) +
  theme(aspect.ratio = 1/1) +
  theme(plot.margin = margin(t = 20, r = 0, b = 20, l = 0))


pca
save_plot(pca, "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.PCA")







# 相关性图
# 计算相关性
log10_data <- log10(count_data) %>%
  filter(if_all(everything(), ~ is.finite(.)))

cor_data <- as.data.frame(cor(log10_data, method = 'pearson')) 

# 保存相关性计算结果
write.csv(cor_data,"C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.correlation.coefficient.csv",row.names = T)

# 绘制散点图矩阵

# 自定义颜色
colors <- colorRampPalette(brewer.pal(9, "YlOrRd")[3:7])((max(cor(log10_data)) - min(cor(log10_data))) * 10000 + 1)

# 自定义面板函数：显示相关系数
panel.cor <- function(x, y, digits = 4, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  # 计算相关系数
  r <- round(cor(x, y), digits)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  
  # 设置字体大小
  if (missing(cex.cor)) {
    cex.cor <- 0.8 / strwidth(txt)
  }
  
  # 绘制背景色
  marg <- par("usr")
  mcol <- floor((r - min(cor(log10_data))) * 10000) + 1
  rect(marg[1], marg[3], marg[2], marg[4], col = colors[mcol])
  
  # 添加文本
  text(0.5, 0.5, txt, cex = cex.cor * (1 + r) / 2, col = "black")
}

# 自定义面板函数：显示直方图
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(usr[1:2], 0, 1.5))
  
  # 计算直方图
  h <- hist(x, breaks = 30, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  
  # 绘制直方图
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", border = "black", ...)
}

# 自定义面板函数：显示散点图和线性拟合
panel.lm <- function(x, y, col = "blue", bg = NA, pch = par("pch"), cex = 1, col.smooth = "red", ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  rug(x, side = 1, col = col, ticksize = 0.02, lwd = 0.2)
  rug(y, side = 2, col = col, ticksize = 0.02, lwd = 0.2)
  abline(stats::lm(y ~ x), col = col.smooth, ...)
}

# 绘制相关性图
pdf(file = "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.correlation.pdf", height = 9, width = 9)
pairs(
  log10_data, 
  pch = ".", 
  upper.panel = panel.cor, 
  diag.panel = panel.hist, 
  lower.panel = panel.lm, 
  main = "Pearson Correlation of All Samples",
  gap = 0
)

#关闭设备
dev.off()


png(filename = "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.correlation.png", height = 2000, width = 2000, res = 300, units = "px")
pairs(
  log10_data, 
  pch = ".", 
  upper.panel = panel.cor, 
  diag.panel = panel.hist, 
  lower.panel = panel.lm, 
  main = "Pearson Correlation of All Samples",
  gap = 0
)
#关闭设备
dev.off()
#R中显示图片
pairs(
  log10_data, 
  pch = ".", 
  upper.panel = panel.cor, 
  diag.panel = panel.hist, 
  lower.panel = panel.lm, 
  main = "Pearson Correlation of All Samples",
  gap = 0
)


R.version.string
# 使用官方 Bioconductor 镜像
options(BioC_mirror = "https://bioconductor.org")
# 安装 BiocManager 并设置版本
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version="3.19", ask=FALSE)
# 安装 ComplexHeatmap
BiocManager::install(c("ComplexHeatmap","circlize"), type="binary")

install.packages("circlize")



# 聚类热图
library(ComplexHeatmap)
library(circlize)
library(grid)
# 绘制热图

cor_mat <- as.matrix(cor_data)

col_anno = HeatmapAnnotation(Condition = sample_info$group,
                             col = list(Condition = c("WT" = "#00DAE0", "cpl2-2" = "#FF9289")),
                             show_legend = F)
hm <- Heatmap(
  cor_mat,  # 输入矩阵
  column_title = "All Samples Correlation Heatmap",
  name = "Correlation",  # 图例名称
  col = colorRampPalette(c("white", "red"))(100),  # 颜色映射
  cluster_rows = T,  # 关闭行聚类
  cluster_columns = T,  # 关闭列聚类
  show_row_names = T,  # 显示行名
  show_column_names = T,  # 显示列名
  row_names_side = "right",  
  column_names_side = "bottom",  
  show_heatmap_legend = F,
  top_annotation = col_anno,
  width = ncol(cor_mat)*unit(50, "pt"), 
  height = nrow(cor_mat)*unit(50, "pt"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.3f", cor_mat[i, j]),
      x, y, gp = gpar(fontsize = 15))
  }
)

lgd_cor <- Legend(
  title = "Cor",  # 图例标题
  at = seq(min(cor_mat), max(cor_mat), length.out = 5), 
  labels = sprintf("%.3f", seq(min(cor_mat), max(cor_mat), length.out = 5)),
  col_fun = colorRamp2(
    breaks = seq(min(cor_mat), max(cor_mat), length.out = 100),
    colors = colorRampPalette(c("white", "red"))(100)
  )
)

lgd_gp <- Legend(
  title = "Condition",  # 图例标题
  at = c("WT", "cpl2-2"),  # 分组名称
  legend_gp = gpar(fill = c("#00DAE0","#FF9289"))
)

hlgd <- packLegend(list = list(lgd_gp,lgd_cor),
                   ncol= 1, direction = "vertical")
heatmap <- 
  draw(hm, annotation_legend_list = hlgd)

# 显示热图
heatmap



# 保存为 PNG
png("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.correlation.heatmap.png", height = 2000, width = 2000, res = 300, units = "px")
ComplexHeatmap::draw(heatmap)
dev.off()

# 保存为 PDF
pdf("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/All.samples.correlation.heatmap.pdf", width = 8, height = 8)
ComplexHeatmap::draw(heatmap)
dev.off()

# 显示在RStudio或外部窗口#R中无法显示热图ComplexHeatmap原因
dev.new(width = 10, height = 10)
draw(hm, annotation_legend_list = hlgd)
