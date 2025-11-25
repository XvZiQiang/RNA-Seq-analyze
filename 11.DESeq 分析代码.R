###################### 差异表达分析


# 查看当前有效库路径（新库在前，原库在后）
.libPaths()

# 设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库


# 检查2个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找

# 安装 Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
BiocManager::install("edgeR", update = FALSE, ask = FALSE)

# 安装 CRAN 包
install.packages("ggplot2", lib = "C:/R packages")
install.packages("ggrepel", lib = "C:/R packages")
install.packages("dplyr", lib = "C:/R packages")
install.packages("openxlsx", lib = "C:/R packages")

# 加载R包
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(openxlsx)
})

# 清空缓存
rm(list = ls())



# 设置一个保存表格函数
save_excel_with_style <- function(data_frame, file_name) {
  # 创建一个Arial字体的样式
  arial_style <- createStyle(fontName = "Arial", fontSize = 10,
                             halign = "center")
  # 创建一个Arial加粗字体的样式
  bold_style <- createStyle(fontName = "Arial", fontSize = 10, 
                            fontColour = "#FFFFFF", 
                            fgFill = "#4682b4",
                            halign = "center",
                            valign = "center",
                            textDecoration = "bold"
  )
  # 创建工作簿和工作表
  wb <- createWorkbook(creator = "")
  addWorksheet(wb, "Sheet 1")
  # 写入数据
  writeData(wb, "Sheet 1", data_frame)
  # 计算并设置每一列的宽度
  for (col in seq_along(data_frame)) {
    max_width <- max(nchar(as.character(data_frame[[col]])), na.rm = TRUE)
    setColWidths(wb, "Sheet 1", col, max(max_width, nchar(colnames(data_frame)[col])) + 2)
  }
  # 应用样式到所有的单元格
  addStyle(wb, "Sheet 1", style = arial_style, 
           rows = 1:(nrow(data_frame) + 1), cols = 1:ncol(data_frame), gridExpand = T)
  # 应用加粗样式到标题行
  addStyle(wb, "Sheet 1", style = bold_style, 
           rows = 1, cols = 1:ncol(data_frame), gridExpand = T)
  # 设置标题行行高为20
  setRowHeights(wb, sheet = 1, rows = 1, height = 20)
  # 设置冻结窗格
  freezePane(wb, "Sheet 1", firstRow = T) 
  # 保存工作簿
  saveWorkbook(wb, file_name, overwrite = T)
}






# 设置一个保存图像函数
save_plot <- function(plot_object, file_prefix) {
  # 保存为PDF
  ggsave(paste0(file_prefix, ".pdf"), plot_object, width = 9, height = 6)
  # 保存为PNG
  png_filename <- paste0(file_prefix, ".png")
  png(png_filename, width = 3000, height = 2000, res = 300)
  print(plot_object)
  dev.off()
}













############## 数据读取和预处理

# 读取样品表达矩阵
count_mat <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/gene.readcount.csv",header = T,sep = ",",
                        row.names = 1,quote = "",check.names = F)

head(count_mat)

norm_mat <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/counts to FPKM/gene.count.FPKM.csv",header = T,sep = ",",
                       row.names = 1,quote = "",check.names = F)

head(norm_mat)










# 读取基因信息文件
gene_info <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/reference/TAIR10.62_genes.csv",header = T,sep = ",",
                        quote = "",check.names = F)

rownames(gene_info) <- gene_info$featureID
#只要这三列gene_ID symbol biotype
gene_info_subset <- cbind(gene_ID = rownames(gene_info), gene_info[, c("symbol", "biotype")])

head(gene_info_subset)











# 设定分组
sample_data <- read.table("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/sample_desc.txt", header = F, sep = "\t", stringsAsFactors = F)
sample_list <- data.frame(row.names = sample_data$V2,group = sample_data$V1)
#取V1组为group，V2组row.names
head(sample_list)





# 设定分组信息：这里需要指定对照组信息
#设置对照组名称为 WT；设置处理组名称为 cpl2-2
ctr <- "WT"
exp <- "cpl2-2"
#计算对照组样本数量 ctr_n；计算处理组样本数量 exp_n
ctr_n <- sum(sample_list$group == ctr)
exp_n <- sum(sample_list$group == exp)




# 设置比较分组
#创建一个因子向量 group 用于差异分析；先重复对照组名称 ctr ctr_n 次，再重复处理组名称 exp exp_n 次
group <- factor(c(rep(ctr, ctr_n), rep(exp, exp_n))) #指定重复数 前为对照，后为处理
#选择 sample_list 中只包含对照组或处理组的行；转换为数据框，并将字符列转换为因子
#得到最终用于差异分析的 sample_info
sample_info <- as.data.frame(subset(sample_list, group == ctr | group == exp),
                             stringsAsFactors = T)






# 设置输出文件夹名称为固定名称 "DESeq"
folder_name <- DESeq
# 设置完整路径
folder_path <- file.path("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/DESeq")
# 创建目录（如果不存在）
dir.create(folder_path, recursive = T, showWarnings = F)
# 设置当前工作目录....../rna-seq-project/DESeq
setwd(folder_path)









# 表达矩阵提取
count_data <- count_mat[rownames(sample_info)]
# 表达矩阵取整
count_data_int <- round(count_data)

# 原始表达矩阵和标准化矩阵提取
#read_data：同 count_data_int，整数化后的原始表达量
read_data <- round(count_mat[rownames(sample_info)])
#norm_data：标准化后的表达矩阵（例如FPKM已经经过某种归一化处理）
#方便后续可视化或聚类分析
norm_data <- norm_mat[rownames(sample_info)]




# 移除所有样品中表达量为0的基因
count_data_int <- count_data_int[rowSums(count_data_int) > 0,]

# 移除低表达量的基因（这里以至少在3个样品中表达量大于1为例）
count_data_filt <- count_data_int[rowSums(count_data_int >= 1) >= ctr_n + 1,]








############## DESeq2差异分析

# 构建DESeq2数据
#使用 DESeqDataSetFromMatrix 将筛选后的计数矩阵 count_data_filt 构建成 DESeq2 的对象
#colData = sample_info 提供每个样品的分组信息
#design = ~ group 指定分析设计：以 group 为因素比较差异
dds <- DESeqDataSetFromMatrix(countData = count_data_filt, colData = sample_info, design= ~ group )
#relevel 将对照组 ctr 设置为基准水平，这样差异分析结果是处理组 vs 对照组
dds$group <- relevel(dds$group, ref = ctr)




# 样本数小于30，rlog法归一
rld <- rlog(dds,blind = F)
# 样本数大于30，vst法归一
#vld <- vst(dds,blind = F)

# 差异分析
dds <- DESeq(dds)

# 获取标准化数据
normal_dds <- as.data.frame(counts(dds, normalized=T))
#计算每组样品的平均表达量
baseMean <- data.frame(sapply(c(ctr, exp), 
                              function(cond) rowMeans(normal_dds[, colData(dds)$group == cond])))



# 由于DESeq2识别不了'-'字符，这里再做一下替换
#DESeq2 默认用 . 替换了分组名中的 -，这里把 . 换回 -。
#给列名加上 baseMean(...) 前缀，方便标识。
colnames(baseMean) <- gsub("\\.", "-", colnames(baseMean))
names(baseMean) <- paste("baseMean(", names(baseMean), ")", sep="")



# 获取差异结果
#contrast = c("group", exp, ctr) 指定比较方向：处理组 vs 对照组。
#na.omit() 去掉缺失值。
deseq2_res <- results(dds,contrast = c("group",exp,ctr))
deseq2_res <- na.omit(deseq2_res)

# 设定阈值
logfc = 1
p = 0.05

# 查看和筛选差异基因
#deseq2_res创建一个新列 diffType。根据 log2FoldChange 和 padj 标记：
#"Up"：上调基因;"Down"：下调基因;"NS"：非显著
deseq2_res$diffType <- "NS"
deseq2_res[which(deseq2_res$log2FoldChange >= logfc & deseq2_res$padj < p),"diffType"] <- "Up"
deseq2_res[which(deseq2_res$log2FoldChange <= -logfc & deseq2_res$padj < p),"diffType"] <- "Down"
deseq2_res[which(abs(deseq2_res$log2FoldChange) < logfc & deseq2_res$padj < p),"diffType"] <- "NS"
deseq2_res[,c("padj","diffType")]
#保存和提取差异基因
deseq2_res_data <- as.data.frame(deseq2_res)
#deseq2_diff：筛选出显著差异的基因。
deseq2_diff <- deseq2_res[which(deseq2_res$diffType != "NS"),]

# 获得筛选结果
deseq2_diff_data <- as.data.frame(deseq2_diff)
#nrow(deseq2_diff) 查看差异基因数量。
nrow(deseq2_diff)

# 提取需要的结果
deseq2_res_slt <- deseq2_res_data[,c("log2FoldChange","pvalue","padj","diffType")]

# 重命名FC列和FDR列
#将 log2FoldChange 改名为 Log2FC(处理/对照)，更直观。
#将 padj 改名为 p.adj，符合常用命名习惯。
names(deseq2_res_slt)[names(deseq2_res_slt) == "log2FoldChange"] <- 
  paste("Log2FC(", exp, "/", ctr, ")", sep="")
names(deseq2_res_slt)[names(deseq2_res_slt) == "padj"] <- "p.adj"

# DESeq2结果添加basemean值/基因平均表达量
#将前面计算的每组基因平均表达量 (baseMean) 与差异分析结果deseq2_res_slt按行合并。
deseq2_res_out <- cbind(baseMean[rownames(baseMean) %in% rownames(deseq2_res_slt), ], 
                        deseq2_res_slt)



#保存完整结果为res_out
res_out <- deseq2_res_out
#保存完整结果为sig_out
#sig_out保存 diffType 为 "Up" 或 "Down" 的基因，即显著差异基因。
sig_out <- res_out[which(res_out$diffType != "NS"),]










############## 保存差异分析结果

# 原始表达矩阵或标准化矩阵合并,给列名加后缀，方便区分来源，导出文件时不会混淆。
#norm_data：标准化表达矩阵FPKM
norm_file <- `colnames<-`(norm_data, paste0(colnames(norm_data), ":FPKM"))
#read_data：原始 read count
read_file <- `colnames<-`(read_data, paste0(colnames(read_data), ":read_count"))
#按行名（基因ID）合并标准化矩阵和原始矩阵，保留所有基因。
expr_file <- merge(norm_file, read_file, by = "row.names", all = TRUE)




# 第一列作为列名并删除
expr_file <- `rownames<-`(expr_file[, -1], expr_file[, 1])

#根据行名匹配,合并差异分析结果、基因信息和表达矩阵
#gene_info_subset：基因注释信息（如名称、描述）
#res_out：DESeq2差异分析结果
#expr_file：标准化 + 原始表达矩阵
#cbind按行名匹配，合并成一个完整的数据框，便于导出。

# 取出 res_out 的基因名
genes <- rownames(res_out)
# 按 res_out 的基因顺序提取基因注释和表达矩阵
gene_info_matched <- gene_info_subset[genes, , drop = FALSE]
expr_matched <- expr_file[genes, , drop = FALSE]
# 合并
res_file <- cbind(gene_info_matched, res_out, expr_matched)






#res_file 是DESeq2差异分析后完整合并后的表格，包含基因信息、差异分析结果和表达矩阵。
#不管是否结果差异显著
# res_file保存为xlsx
res_xlsx <- "DESeq.res.xlsx"
save_excel_with_style(res_file, res_xlsx)

# 保存为csv
res_csv <- "DESeq.res.csv"
write.csv(res_file,res_csv,row.names = F)




# 保存显著差异基因结果
sig_file <- res_file[which(res_file$diffType != "NS"),]

#sig_file 是仅包含显著差异基因的表格（diffType != "NS"）。
#结果是差异显著的几百个基因
# 保存为xlsx
sig_xlsx <- "DESeq.SigDE.res.xlsx"
save_excel_with_style(sig_file, sig_xlsx)

# 保存为csv
sig_csv <- "DESeq.SigDE.res.csv"
write.csv(sig_file, sig_csv,row.names = F)








############## 可视化

# 绘制火山图
install.packages("ggrastr")
library(ggrastr)

de_data <- read.csv("DESeq.res.csv",header = T,sep = ",",
                    check.names = F)

# 提取指定基因
mark_gene <- subset(de_data, symbol %in% c("FT","FLC","MYB47","GBSS1","TRX5"))

volplot <- 
  ggplot(de_data, mapping = aes(x = de_data[,6], # x轴：log2 Fold Change（第6列）
                                y = -log10(de_data[, 8]), # y轴：-log10(p.adj)（第8列）
                                color = factor(diffType)))+# 点的颜色按差异类型分类（Up/Down/NS）
  # 绘制散点，使用 ggrastr 加速大数据点绘制
  geom_point_rast(size = 3.5)+
  # x轴范围固定为[-15,15]
  scale_x_continuous(limits = c(-15, 15)) +
  # y轴留一点空白边距
  scale_y_continuous(expand = c(0.05,0.05), limits = c(0, 500)) +# 设置y轴上限为500
  # 显著性阈值线（p.adj=0.05）
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey") +
  # FC阈值线（log2FC = ±1
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey") +
  # 上调基因数量/下调基因数量/非显著
  scale_color_manual(labels=c(paste0("Up (",length(de_data$diffType[de_data$diffType == "Up"]),')'),
                              paste0("Down (",length(de_data$diffType[de_data$diffType == "Down"]),')' ),
                              'NS'),
                     # 自定义颜色
                     values=c("#e13222", "#75a8d7", "#d4cfd0"),
                     # 显示顺序
                     breaks=c("Up","Down","NS")) +
  # x轴标签+ y轴标签+
  labs(x = expression(Log[2]*"FC (WT vs. cpl2-2)"), 
       y = expression("-Log"[10]*"("*italic(P)["adj"]*")"),
       #color="Significance",# 去掉图例标题
       color = NULL
  ) +
  # 基本主题字体大小
  theme_minimal(base_size = 18) +
  # 保持宽高比1:1
  theme(aspect.ratio = 1) +
  # 面板边框
  theme(panel.border = element_rect(fill=NA, linewidth = 0.75),
        panel.grid.major = element_blank(),# 去掉主网格线
        panel.grid.minor = element_blank()) +# 去掉次网格线
  # 坐标轴文字和坐标轴刻度线，颜色
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  # 图例在左上角
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),# 图例背景透明
        legend.key.spacing = unit(0.01,"cm")) +# 图例间距
  # 图例点大小
  guides(color=guide_legend(override.aes = list(size=5,alpha=1))) +
  # 图边距
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
  # 标注特定基因
  geom_label_repel(
    # 需要标注的基因数据
    data = mark_gene,
    aes(x = mark_gene[,6], 
        y = -log10(mark_gene[,8]),
        label = symbol),
    # 不在图例显示
    show.legend = F
  )


volplot
#保存
save_plot(volplot, "All DESeq analyse gene.vocalno")












# 绘制散点图

scatterplot <- 
  ggplot(de_data, mapping = aes(x = log10(de_data[, 4]), # X轴：WT样品平均表达（取第4列），取log10
                                y = log10(de_data[, 5]), # Y轴：cpl2-2样品平均表达（取第5列），取log10
                                color = factor(diffType)))+# 点颜色：根据差异类型（Up/Down/NS）区分
  # 绘制散点，使用 raster 加速大数据点绘制
  geom_point_rast(size = 3.5) +
  # 对角线 y=x，参考线，黑色虚线
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, lty = 2, color = "black", lineend = "round") +
  # X轴范围和刻度
  scale_x_continuous(limits = c(-1, 6), breaks = seq(-1, 6, 1), expand = c(0,0)) +
  # y轴范围和刻度
  scale_y_continuous(limits = c(-1, 6), breaks = seq(-1, 6, 1), expand = c(0,0)) +
  # 超出坐标轴的标签不被剪裁
  coord_cartesian(clip = "off") +
  # 上调基因数量/下调基因数量/非显著基因
  scale_color_manual(labels=c(paste0("Up (",length(de_data$diffType[de_data$diffType == "Up"]),')'),
                              paste0("Down (",length(de_data$diffType[de_data$diffType == "Down"]),')' ),
                              'NS'),
                     values=c("#e13222", "#75a8d7", "#d4cfd0"),
                     breaks=c("Up","Down","NS")) +
  #X+Y轴标签
  labs(x = "Mean Expression of WT", 
       y = "Mean Expression of cp2-2",
       #color="Significance",
       color = NULL
  ) +
  # 使用简洁主题，基础字体18
  theme_minimal(base_size = 18) +
  
  theme(panel.border = element_rect(fill=NA, linewidth = 0.5),# 面板边框
        panel.grid.major = element_blank(),# 去掉主网格
        panel.grid.minor = element_blank()) +# 去掉次网格
  # 坐标轴文字颜色# 坐标轴刻度
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),# 图例背景透明
        legend.key.spacing = unit(0.01,"cm")) +# 图例间距
  # XY轴比例1:1
  theme(aspect.ratio = 1/1) +
  # 图例点大小
  guides(color=guide_legend(override.aes = list(size=5,alpha=1))) +
  # 图边距
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
  # 特定标记基因
  geom_label_repel(
    data = mark_gene,
    aes(x = log10(mark_gene[, 4]), 
        y = log10(mark_gene[, 5]),
        label = symbol),
    show.legend = F
  )

scatterplot
#保存
save_plot(scatterplot, "All DESeq analyse gene.scatterplot")











# 基因表达热图

library(ComplexHeatmap)
library(circlize)

# 从火山图的de_data提取上调或下调的基因数据
sig_data <- de_data[which(de_data$diffType != "NS"),]

# 去掉 symbol 为 NA 或 '-' 的行
sig_data <- sig_data[!is.na(sig_data$symbol) & sig_data$symbol != "-", ]

# 将基因名称 symbol 设置为行名，方便后续热图绘制和索引
rownames(sig_data) <- sig_data$symbol

# 提取表达量列，这里假设第10到第15列是样品的 FPKM 表达量
sig_mat <- sig_data[,c(10:15)]

# 去除FPKM字符
names(sig_mat) <- gsub(":FPKM", "", names(sig_mat))

# 对sig_mat数据标准化
# 1. log2(FPKM+1) 转换表达量
# 2. scale() 将每行（每个基因）归一化为均值为0、标准差为1
for (i in 1:nrow(sig_mat)) sig_mat[i, ] <- scale(log(unlist(sig_mat[i, ] + 1), 2))


# 保存sig_mat表格
write.csv(sig_mat, "DESeq analyse Significant differences gene.heatmap.csv", row.names = F)









# 绘制热图
# 创建列注释（顶部显示分组信息）
col_anno = HeatmapAnnotation(Condition = sample_info$group,# 每个样品的分组
                             col = list(Condition = c("WT" = "#00DAE0", "cpl2-2" = "#FF9289")),# 分组颜色
                             show_legend = F)# 不显示默认图例
# 创建行注释（右侧标注特定基因）
row_anno = rowAnnotation(link = anno_mark(at = which(rownames(sig_mat) %in% mark_gene$symbol), # 需要标注的行索引+标注的基因名
                                          labels = mark_gene$symbol, labels_gp = gpar(fontsize = 10)))# 标签字体大小
# 构建热图对象
hm <- Heatmap(
  as.matrix(sig_mat),  # 输入矩阵
  column_title = "Heatmap of Differentially Expressed Genes",# 图标题
  col = colorRampPalette(c("blue","white","red"))(100), # 颜色映射
  cluster_rows = T,  # 开启行聚类
  cluster_columns = T,  # 开启列聚类
  show_row_names = F,  # 显示行名
  show_column_names = T,  # 显示列名
  row_names_side = "right",  
  column_names_side = "bottom",  
  show_heatmap_legend = F,# 不显示默认热图图例
  top_annotation = col_anno,# 添加顶部列注释
  right_annotation = row_anno,# 添加右侧行注释
)
# 构建热图图例（颜色图例）
lgd_exp <- Legend(
  title = "Normalized Express",  # 图例标题
  at = seq(-2, 2, length.out = 5), # 图例刻度
  col_fun = colorRamp2(
    breaks = seq(min(sig_mat), max(sig_mat), length.out = 100),# 对应数据范围
    colors = colorRampPalette(c("blue","white","red"))(100)# 对应颜色
  )
)
# 构建分组图例（列注释图例）
lgd_gp <- Legend(
  title = "Condition",  # 图例标题
  at = c("WT", "cpl2-2"),  # 分组名称
  legend_gp = gpar(fill = c("#00DAE0","#FF9289"))# 对应颜色
)
# 合并两个图例，纵向排列
hlgd <- packLegend(list = list(lgd_gp,lgd_exp),
                   ncol= 1, direction = "vertical")
# 绘制热图，并添加自定义图例
heatmap <- 
  draw(hm, annotation_legend_list = hlgd)

# 显示热图
heatmap

# 保存为 PNG
png("DESeq analyse Significant differences gene.heatmap.png", height = 2000, width = 2000, res = 300, units = "px")
ComplexHeatmap::draw(heatmap)
dev.off()

# 保存为 PDF
pdf("DESeq analyse Significant differences gene.heatmap.pdf", width = 6, height = 6)
ComplexHeatmap::draw(heatmap)
dev.off()

