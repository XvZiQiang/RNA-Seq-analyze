# 查看当前有效库路径（新库在前，原库在后）
.libPaths()

# 设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库


# 检查3个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找

# 安装 Bioconductor 包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
BiocManager::install("edgeR", update = FALSE, ask = FALSE)





# 加载R包
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(ggh4x)
library(tidyr)
library(dplyr)
library(aplot)
library(forcats)
library(viridis)
library(gridExtra)
library(openxlsx)

# 首先设定在根目录下
setwd("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/enrich")

# 清空缓存
rm(list = ls())

# 设置一个保存表格函数
save_excel_with_style <- function(data_frame, file_name) {
  # 创建一个Arial加粗字体的样式
  bold_style <- createStyle(fontName = "Arial", fontSize = 10, 
                            fontColour = "#FFFFFF", 
                            fgFill = "#4682b4",
                            halign = "center",
                            valign = "center",
                            textDecoration = "bold"
  )
  # 创建一个居中的数值样式
  center_style <- createStyle(halign = "center",fontName = "Arial", fontSize = 10)
  # 创建一个左对齐的文本样式
  left_style <- createStyle(halign = "left",fontName = "Arial", fontSize = 10)
  # 创建工作簿和工作表
  wb <- createWorkbook(creator = "")
  addWorksheet(wb, "Sheet 1")
  # 写入数据
  writeData(wb, "Sheet 1", data_frame)
  # 最大宽度
  max_width <- 25
  # 应用样式到单元格
  for (row in 1:nrow(data_frame)) {
    for (col in 1:ncol(data_frame)) {
      # 检查单元格的类型
      if (is.numeric(data_frame[row, col])) {
        # 如果是数值类型，应用居中样式
        addStyle(wb, "Sheet 1", style = center_style, rows = row + 1, cols = col)  # row + 1 因为第一行是标题
      } else {
        # 否则，应用左对齐样式
        addStyle(wb, "Sheet 1", style = left_style, rows = row + 1, cols = col)  # row + 1 因为第一行是标题
      }
    }
  }
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





#设置一个保存图片函数
save_plot <- function(plot_object, file_prefix) {
  # 保存为PDF
  ggsave(paste0(file_prefix, ".pdf"), plot_object, width = 10, height = 10)
  # 保存为PNG
  png_filename <- paste0(file_prefix, ".png")
  png(png_filename, width = 3000, height = 3000, res = 300)
  print(plot_object)
  dev.off()
}






# 设置要查找文件的目录
dir_path <- "C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/DESeq"  # 替换为实际目录的路径


# 找出目录中带有比较组合前缀的DESeq.res.csv文件
file_pattern <- paste0("DESeq.res.csv")
file <- list.files(path = dir_path, pattern = file_pattern, full.names = T, recursive = T)
# 读取差异分析结果表
de_file <- read.csv(file, check.names = F)




# 选取重要的几列
# 这里假设第1列是Geneid，第2列是Symbol，第6-9列是LogFC、pvalue、FDR、diffType
de_data <- de_file[, c(1:2,6:9)]
colnames(de_data) <- c("Geneid", "Symbol", "LogFC", "pvalue", "FDR", "diffType")


# 去掉 de_data 中 Geneid 前面的 "gene:"
de_data$Geneid <- gsub("^gene:", "", de_data$Geneid)

#查看de_data
head(de_data)






#设置输出目录
enrich_path <- file.path("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/enrich")
#如果没有，新建一个文件夹
if (!dir.exists(enrich_path)) {
  dir.create(enrich_path, recursive = TRUE)
}

# 设置当前工作目录
setwd("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/enrich")

# 读取基因列表de_data
# 仅保留上调和下调的基因（去掉 diffType == "NS" 的）
# 从差异分析结果中提取非 NS 的基因 ID
gene_select <- de_data$Geneid[de_data$diffType != "NS"]

#查看gene_select
head(gene_select)







########################################################## KEGG富集分析
# enrichKEGG 需要 organism = 'ath'（Arabidopsis thaliana）
kegg_rich <- enrichKEGG(gene = gene_select,
                        organism = 'ath',
                        pvalueCutoff = 1, # 暂时不过滤，后面手动筛选
                        pAdjustMethod = 'BH',# 多重校正方法 Benjamini-Hochberg
                        qvalueCutoff = 1,
                        maxGSSize = 500# 每个 KEGG 通路最多允许包含的基因数
)

# 将结果转换为数据框
kegg_rich_result <- as.data.frame(kegg_rich)

# 去除 Description 列中 "-" 后面的物种信息
# 清理 KEGG 通路描述（去掉后缀物种名）
# 比如 "Photosynthesis - Arabidopsis thaliana" 变成 "Photosynthesis"
kegg_rich_result$Description <- str_replace(kegg_rich_result$Description, " - .*", "")

# 准备上下调基因列表
up_genes <- de_data$Geneid[de_data$diffType == "Up"]   # 上调基因
down_genes <- de_data$Geneid[de_data$diffType == "Down"] # 下调基因

# 创建新列以存储上调和下调基因的列表及其数量（均为NA）
kegg_rich_result$up_genes <- NA
kegg_rich_result$down_genes <- NA
kegg_rich_result$UpNum <- NA
kegg_rich_result$DownNum <- NA









# 遍历每个 KEGG 通路，统计上/下调基因数量   并放入上一步生成的4列里
for (i in seq_len(nrow(kegg_rich_result))) {
  # 获取当前路径的基因ID
  genes_in_pathway <- unlist(strsplit(kegg_rich_result$geneID[i], "/"))
  
  # 获取上调和下调基因
  upregulated_in_pathway <- intersect(genes_in_pathway, up_genes)
  downregulated_in_pathway <- intersect(genes_in_pathway, down_genes)
  
  # 保存结果
  kegg_rich_result$up_genes[i] <- paste(upregulated_in_pathway, collapse = ", ")
  kegg_rich_result$down_genes[i] <- paste(downregulated_in_pathway, collapse = ", ")
  kegg_rich_result$UpNum[i] <- length(upregulated_in_pathway)
  kegg_rich_result$DownNum[i] <- length(downregulated_in_pathway)
}

# 将数据框中所有不含值的单元格用"-"代替
kegg_rich_result[kegg_rich_result == ""] <- "-"










##########非必须步骤：给富集结果加上Symbol信息

# 创建基因 ID 到 Symbol 的映射
# 例如：AT1G01010 -> LHY
# 这里假设 de_data 里有一列 “Symbol”，表示对应的基因简称
gene_symbol_map <- setNames(de_data$Symbol, de_data$Geneid)

# 格式化 geneID 列
# kegg_rich_result的geneID 列加 Symbol
kegg_rich_result$geneID <- sapply(strsplit(kegg_rich_result$geneID, "/"), function(genes) {
  paste(sapply(genes, function(gene) {
    if (gene %in% names(gene_symbol_map)) {
      paste0(gene, "(", gene_symbol_map[[gene]], ")")
    } else {
      gene
    }
  }), collapse = "/")
})

# kegg_rich_result的 up_genes 和 down_genes 列
# 如果为“-”，表示该通路中没有上调或下调基因，不做修改
# 否则在每个基因后面加上Symbol
kegg_rich_result$up_genes <- sapply(kegg_rich_result$up_genes, function(genes) {
  if (genes == "-") {
    "-"
  } else {
    paste(sapply(strsplit(genes, ", ")[[1]], function(gene) {
      if (gene %in% names(gene_symbol_map)) {
        paste0(gene, "(", gene_symbol_map[[gene]], ")")
      } else {
        gene
      }
    }), collapse = ", ")
  }
})

kegg_rich_result$down_genes <- sapply(kegg_rich_result$down_genes, function(genes) {
  if (genes == "-") {
    "-"
  } else {
    paste(sapply(strsplit(genes, ", ")[[1]], function(gene) {
      if (gene %in% names(gene_symbol_map)) {
        paste0(gene, "(", gene_symbol_map[[gene]], ")")
      } else {
        gene
      }
    }), collapse = ", ")
  }
})




#kegg_rich_result变成kegg_rich_out
kegg_rich_out <- kegg_rich_result

# 根据列位置重命名（确保列号对应你自己的表结构）
# 一般列意义如下：
#   9: Zscore / 通路富集得分
#   13: GeneList / 所有富集基因
#   15: UpList / 上调基因列表
#   16: DownList / 下调基因列表
# 指定列的表头重命名
colnames(kegg_rich_out)[c(9, 13, 15, 16)] <- c("Zscore", "GeneList", "UpList", "DownList")

#按照 pvalue 升序排列，显著性高的排在前面
kegg_rich_out <- arrange(kegg_rich_out, pvalue) #按照pvalue排序

# 保存结果为 Excel 文件
kegg_rich_xlsx <- paste("KEGG.enrich.res.xlsx", sep="")
# 注意：save_excel_with_style() 是自定义函数
save_excel_with_style(kegg_rich_out, kegg_rich_xlsx)














# 可视化
# 取前20个富集通路（如果总数少于20，则取实际行数）
kegg_plot_data <- kegg_rich_out[1:min(20, nrow(kegg_rich_out)), ] #取前20行
# 创建图中显示的标签，将通路名称和KEGG ID组合显示
kegg_plot_data$label <- sprintf("%s\n(%s)", kegg_plot_data$Description, kegg_plot_data$ID)









# 绘制气泡图
# 第一种气泡图：以富集因子和Pathway名称为X轴和Y轴
plot1 <- 
  ggplot(kegg_plot_data, aes(FoldEnrichment, label)) + 
  geom_point(aes(size = Count, color = pvalue),# 气泡大小对应Count，颜色对应pvalue
             alpha = 1,# 气泡透明度
             shape = 19,# 圆形
             stroke = 0.5# 边框宽度
  ) + 
  # 设置气泡大小范围
  scale_size_continuous(range = c(5, 15)) +
  # 设置颜色梯度：pvalue小（显著）颜色为红色，大为蓝色
  scale_color_gradient(low = "red",high = "blue",
                       name = expression(italic(P)["value"])) + #颜色变化
  # 图标题
  ggtitle("KEGG Pathway Enrichment")+
  # 坐标轴和图例标签
  labs(color = expression(pvalue),
       size = "Count",
       x = "FoldEnrichment",
       y = NULL) +
  # 调整图例顺序
  guides(fill = NULL,
         size = guide_legend(order = 1))+ #修改图例顺序
  # 基本主题
  theme_linedraw(base_size = 18) +
  # 标题大小相对缩放
  theme (plot.title=element_text (size = rel(1))) +
  # 坐标轴文字颜色
  theme(axis.text=element_text(color = "black"))+
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))+
  
  
  # 去掉y轴刻度线
  theme(axis.ticks.y = element_blank()) +
  # 网格线、边框和图例背景设置
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black",linetype = 2, linewidth = 0.25),
        panel.border = element_rect(fill = NA,linewidth = 0.75),
        legend.background = element_blank()
  )
# 展示气泡图
plot1
# 保存图形
plot1_file <- paste( "KEGG.enrich.bubbleplot1", sep="")
save_plot(plot1, plot1_file)












# 双向柱状图可视化

# KEGG 通路左侧柱状图

# 获取所有通路分类类型
category_type <- unique(kegg_plot_data$category)

# 为每个类别生成颜色，使用 viridis 调色板
category_color <- viridis(n = length(category_type), option = "D")

#添加显著性标记（pvalue: *** <0.001, ** <0.01, * <0.05, 其他空）
kegg_plot_data <- kegg_plot_data %>%
  mutate(significance = case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))
# 重新排序标签（label），确保柱状图按分类顺序显示
kegg_plot_data <- kegg_plot_data %>%
  mutate(label = fct_reorder(label, category)) # 确保label按顺序显示




# 绘制左侧柱状图
plotl <- ggplot(kegg_plot_data, aes(x = log10(pvalue), y = label, fill = category)) +
  geom_bar(stat="identity",color = "white",linewidth = 0.4,
           alpha = 0.8,width = 0.75) +
  # 添加显著性标记
  geom_text(aes(label = significance),
            nudge_x = -0.1,    # 水平方向调整，负值表示左移
            hjust = 1,         # 左侧对齐
            nudge_y = -0.15,     # 调整垂直方向距离
            size = 5, 
            color = "black") +
  # x 轴范围设置为负的 log10(pvalue)
  labs(x = expression("-Log"[10]*"("*italic(P)["value"]*")"), y= NULL, fill = "Pathway Category") +
  scale_x_continuous(limits = c(min(log10(kegg_plot_data$pvalue)) * 1.1, 0), labels = abs) + 
  scale_fill_manual(values = setNames(category_color, category_type)) + 
  ggtitle("KEGG Pathway Enrichment",
          subtitle = ("Order by significance on the left and by count on the right"))+
  theme_minimal(base_size = 14) +
  theme (plot.title=element_text (size = rel(1)),
         plot.subtitle = element_text(size = rel(0.75))) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm")) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 3,
                                                 color = "white")))

plotl







# 右侧柱状图
# KEGG 通路右侧柱状图（Up/Down基因数量）
# 转换数据为长格式，列 Direction 包含 UpNum/DownNum，值为 RegulateNum
kegg_plot_data_long <- kegg_plot_data %>%
  pivot_longer(cols = c(UpNum, DownNum), 
               names_to = "Direction", 
               values_to = "RegulateNum")
# 绘制右侧柱状图
plotr <- ggplot(kegg_plot_data_long, aes(x = RegulateNum, y = label, fill = Direction)) +
  geom_bar(stat="identity",color = "white",linewidth = 0.4,
           alpha = 0.5,width = 0.75) +
  scale_x_continuous(limits = c(0, max(kegg_plot_data_long$RegulateNum) * 1.1)) +
  scale_fill_manual(values = c("UpNum" = "#d73032", "DownNum" = "#4575b4"),
                    labels = c("Down", "Up")) +
  labs(x = "Count", y= NULL,fill = "Feature Type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm")) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 3,
                                                 color = "white")))
plotr





#将左右图合并
plot3 <- plotl %>% insert_right(plotr, width = 0.5)
# 显示图形
plot3
# 保存图形
plot3_file <- paste("KEGG.enrich.barchat", sep="")
save_plot(plot3, plot3_file)









########################################################## GO富集分析
BiocManager::install("org.At.tair.db", update = FALSE, ask = FALSE)

library(org.At.tair.db)


#类似KEGG，选定物种数据库
go_rich <- enrichGO(gene = gene_select,
                    OrgDb = org.At.tair.db,
                    keyType = 'TAIR',
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 1,
                    maxGSSize = 100,
                    ont = "all"
)

go_rich_result <- as.data.frame(go_rich)

# 添加上下调结果
up_genes <- de_data$Geneid[de_data$diffType == "Up"]   # 上调基因
down_genes <- de_data$Geneid[de_data$diffType == "Down"] # 下调基因

# 创建新列以存储上调和下调基因的列表及其数量
go_rich_result$up_genes <- NA
go_rich_result$down_genes <- NA
go_rich_result$UpNum <- NA
go_rich_result$DownNum <- NA

# 遍历每个GO功能
for (i in seq_len(nrow(go_rich_result))) {
  # 获取当前路径的基因ID
  genes_in_term <- unlist(strsplit(go_rich_result$geneID[i], "/"))
  
  # 获取上调和下调基因
  upregulated_in_term<- intersect(genes_in_term, up_genes)
  downregulated_in_term <- intersect(genes_in_term, down_genes)
  
  # 保存结果
  go_rich_result$up_genes[i] <- paste(upregulated_in_term, collapse = ", ")
  go_rich_result$down_genes[i] <- paste(downregulated_in_term, collapse = ", ")
  go_rich_result$UpNum[i] <- length(upregulated_in_term)
  go_rich_result$DownNum[i] <- length(downregulated_in_term)
}

# 将数据框中所有不含值的单元格用"-"代替
go_rich_result[go_rich_result == ""] <- "-"







##########非必须步骤：给富集结果加上Symbol信息

# 创建基因 ID 到 Symbol 的映射
gene_symbol_map <- setNames(de_data$Symbol, de_data$Geneid)

# 格式化 geneID 列
go_rich_result$geneID <- sapply(strsplit(go_rich_result$geneID, "/"), function(genes) {
  paste(sapply(genes, function(gene) {
    if (gene %in% names(gene_symbol_map)) {
      paste0(gene, "(", gene_symbol_map[[gene]], ")")
    } else {
      gene
    }
  }), collapse = "/")
})

# 格式化 up_genes 和 down_genes 列
go_rich_result$up_genes <- sapply(go_rich_result$up_genes, function(genes) {
  if (genes == "-") {
    "-"
  } else {
    paste(sapply(strsplit(genes, ", ")[[1]], function(gene) {
      if (gene %in% names(gene_symbol_map)) {
        paste0(gene, "(", gene_symbol_map[[gene]], ")")
      } else {
        gene
      }
    }), collapse = ", ")
  }
})

go_rich_result$down_genes <- sapply(go_rich_result$down_genes, function(genes) {
  if (genes == "-") {
    "-"
  } else {
    paste(sapply(strsplit(genes, ", ")[[1]], function(gene) {
      if (gene %in% names(gene_symbol_map)) {
        paste0(gene, "(", gene_symbol_map[[gene]], ")")
      } else {
        gene
      }
    }), collapse = ", ")
  }
})

go_rich_out <- go_rich_result

# 指定列重命名
colnames(go_rich_out)[c(8, 12, 14, 15)] <- c("Zscore", "GeneList", "UpList", "DownList")

# 排序
go_rich_out <- arrange(go_rich_out,pvalue) #按照pvalue排序


#保存表格
go_rich_xlsx <- paste("GO.enrich.res.xlsx", sep="")
save_excel_with_style(go_rich_out, go_rich_xlsx)









# 可视化

go_plot_data <- go_rich_out[1:min(20, nrow(go_rich_out)), ] #取前20行
go_plot_data$label <- sprintf("%s\n(%s)", go_plot_data$Description, go_plot_data$ID)

# 第一种气泡图：以富集因子和Pathway名称为X轴和Y轴
plot4 <- 
  ggplot(go_plot_data,aes(FoldEnrichment,label)) + #以富集因子和Pathway名称为X轴和Y轴
  geom_point(aes(size = Count, color = pvalue),
             alpha = 1,
             shape = 19,
             stroke=0.5
  ) + #点图大小和颜色数据
  scale_size_continuous(range = c(5, 15)) +
  scale_color_gradient(low = "red",high = "blue",
                       name = expression(italic(P)["value"])) + #颜色变化
  ggtitle("GO Term Enrichment")+
  labs(color = expression(pvalue),
       size = "Count",
       x = "FoldEnrichment",
       y = NULL) +
  guides(fill = NULL,
         size = guide_legend(order=1))+ #修改图例顺序
  theme_minimal(base_size = 18) +
  theme (plot.title=element_text (size = rel(1)),
         plot.subtitle = element_text(size = rel(0.75))) +
  theme(axis.text=element_text(color = "black")) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 30))+
  theme(text = element_text(family = "sans"),
        axis.ticks.y = element_blank()) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", linetype = 2, linewidth = 0.25),
        panel.border = element_rect(fill = NA,linewidth = 0.75),
        legend.background = element_blank())
plot4
#保存
plot4_file <- paste("GO.enrich.bubbleplot1", sep="")
save_plot(plot4, plot4_file)






# 双向柱状图

# 左侧柱状图

# 获取所有 ONTOLOGY 种类
ontology_type <- unique(go_plot_data$ONTOLOGY)

# 动态生成足够多的颜色，使用 viridis 
ontology_color <- viridis(n = length(ontology_type), option = "D")

go_plot_data <- go_plot_data %>%
  mutate(significance = case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))

go_plot_data <- go_plot_data %>%
  mutate(label = fct_reorder(label, ONTOLOGY)) # 确保label按顺序显示

plotl_go <- ggplot(go_plot_data, aes(x = log10(pvalue), y = label, fill = ONTOLOGY)) +
  geom_bar(stat="identity",color = "white",linewidth = 0.4,
           alpha = 0.8,width = 0.75) +
  geom_text(aes(label = significance),
            nudge_x = -0.1,    # 水平方向调整，负值表示左移
            hjust = 1,         # 左侧对齐
            nudge_y = -0.15,     # 调整垂直方向距离
            size = 5, 
            color = "black") +
  labs(x = expression("-Log"[10]*"("*italic(P)["value"]*")"), y= NULL, fill = "Ontology Category") +
  scale_x_continuous(limits = c(min(log10(go_plot_data$pvalue)) * 1.1, 0), labels = abs) + 
  scale_fill_manual(values = ontology_color) +
  ggtitle("GO Term Enrichment",
          subtitle = ("Order by significance on the left and by count on the right"))+
  theme_minimal(base_size = 14) +
  theme (plot.title=element_text (size = rel(1)),
         plot.subtitle = element_text(size = rel(0.75))) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm")) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 5,
                                                 color = "white")))

plotl_go

# 右侧柱状图
go_plot_data_long <- go_plot_data %>%
  pivot_longer(cols = c(UpNum, DownNum), 
               names_to = "Direction", 
               values_to = "RegulateNum")

plotr_go <- ggplot(go_plot_data_long, aes(x = RegulateNum, y = label, fill = Direction)) +
  geom_bar(stat="identity",color = "white",linewidth = 0.4,
           alpha = 0.5,width = 0.75) +
  scale_x_continuous(limits = c(0, max(go_plot_data_long$RegulateNum) * 1.1)) +
  scale_fill_manual(values = c("UpNum" = "#d73032", "DownNum" = "#4575b4"), 
                    labels = c("Down", "Up")) +
  labs(x = "Count", y= NULL,fill = "Feature Type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm")) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(fill = guide_legend(override.aes = list(shape = 21,size = 5,
                                                 color = "white")))
plotr_go



#合并5图+6图
plot6 <- plotl_go %>% insert_right(plotr_go, width = 0.75)
plot6
#保存图片
plot6_file <- paste("GO.enrich.barchat", sep="")
save_plot(plot6, plot6_file)














########################################################## GSEA KEGG富集分析
# GSEA 使用的是连续值：通常是 LogFC、t 值等
# 下面用差异分析结果中的 LogFC 构建排序基因列表
# 提取 LogFC
gene_lfc <- de_data$LogFC

# 设置基因名
names(gene_lfc) <- de_data$Geneid

# --- 关键清洗步骤（解决你的报错） ---
# 去除含 NA 基因名、空字符基因名、含 NA 数值
valid_idx <- !is.na(names(gene_lfc)) & 
  names(gene_lfc) != "" &
  !is.na(gene_lfc)

gene_lfc <- gene_lfc[valid_idx]

# 查看是否仍有 NA
sum(is.na(names(gene_lfc)))     # 必须是 0
sum(is.na(gene_lfc))            # 必须是 0

# 去除重复基因名（GSEA 不允许重复）
gene_lfc <- gene_lfc[!duplicated(names(gene_lfc))]

# 最后排序
gene_lfc <- sort(gene_lfc, decreasing = TRUE)

# 检查前 10 个
head(gene_lfc)




########## KEGG GSEA富集

gsea_kegg <- gseKEGG(
  geneList = gene_lfc,      # 排序后的基因列表
  organism = 'ath',          # 拟南芥对应 KEGG 物种代码
  pvalueCutoff = 1,          # 不设过滤，让结果全部输出
  pAdjustMethod = 'BH',      # FDR 校正方法
  verbose = F
)

head(gsea_kegg)




# 去除结果 Description 中 "-" 后面的物种信息
gsea_kegg@result$Description <- str_replace(gsea_kegg@result$Description, " - .*", "")
#转换为数据表
gsea_kegg_result <- as.data.frame(gsea_kegg)




# 添加 core enrichment 的基因数量（每个通路贡献的核心基因数）
gsea_kegg_result$core_gene_num <- sapply(strsplit(gsea_kegg_result$core_enrichment, "/"), length) 


# 根据 NES 判断通路激活/抑制
gsea_kegg_result[which(gsea_kegg_result$NES > 0),'stats'] <- 'Activated'
gsea_kegg_result[which(gsea_kegg_result$NES < 0),'stats'] <- 'Suppressed'


# 按 FDR 排序（p.adjust）
gsea_kegg_result <- arrange(gsea_kegg_result,p.adjust) #按照p.adjust即FDR排序

# 输出对象
gsea_kegg_out <- gsea_kegg_result


# 输出为 Excel
gsea_kegg_xlsx <- paste("KEGG.GSEA.res.xlsx", sep="")
save_excel_with_style(gsea_kegg_out, gsea_kegg_xlsx)





########## GO GSEA富集

gsea_go <- gseGO(
  geneList       = gene_lfc,           # 排序基因
  OrgDb          = org.At.tair.db,     # 拟南芥注释数据库
  keyType        = 'TAIR',             # 基因类型
  pvalueCutoff   = 0.05,               # GSEA 默认 FDR < 0.05
  pAdjustMethod  = 'BH',
  maxGSSize      = 500,                # GO 类别较多，可适当放大
  ont            = "all",              # BP + CC + MF
  verbose        = FALSE
)

# 查看 GO GSEA 的前几行
head(gsea_go)


# 转为 data.frame表格
gsea_go_result <- as.data.frame(gsea_go)



# 添加 core enrichment 基因数量
gsea_go_result$core_gene_num <- sapply(strsplit(gsea_go_result$core_enrichment, "/"), length) 

# 根据 NES 判断激活/抑制
gsea_go_result[which(gsea_go_result$NES > 0),'stats'] <- 'Activated'
gsea_go_result[which(gsea_go_result$NES < 0),'stats'] <- 'Suppressed'

# 按照 FDR 排序
gsea_go_result <- arrange(gsea_go_result,p.adjust) #按照p.adjust即FDR排序

# 输出对象
gsea_go_out <- gsea_go_result

# 写入 Excel
gsea_go_xlsx <- paste("GO.GSEA.res.xlsx", sep="")
save_excel_with_style(gsea_go_out, gsea_go_xlsx)






########## GSEA可视化

# GSEA柱状图可视化

# 筛选显著性前 10 的KEGG通路，并定义透明度
gsea_kegg_bardata <-
  gsea_kegg_out %>%
  filter(stats %in% c("Suppressed", "Activated")) %>%
  mutate(
    transparency = case_when(
      p.adjust < 0.001 ~ 1,          # 高度显著
      p.adjust < 0.01 ~ 0.8,        # 显著
      p.adjust < 0.05 ~ 0.6,        # 中等显著
      TRUE ~ 0.4                    # 不显著
    )
  ) %>%
  arrange(stats, desc(abs(NES))) %>% # 按 NES 绝对值排序
  group_by(stats) %>%
  slice_head(n = 10) %>% # 每组取前10
  ungroup()

# 将透明度转换为离散型变量
gsea_kegg_bardata <- gsea_kegg_bardata %>%
  mutate(transparency_discrete = as.factor(transparency))

# 绘制 KEGG GSEA 杆状图
gsea_kegg_barplot <-
  ggplot(gsea_kegg_bardata, aes(x = abs(NES), y = fct_reorder(Description, stats))) +
  geom_col(aes(fill = stats, alpha = transparency_discrete),
           position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.45) +
  geom_text(aes(x = 0.05, y = Description, label = Description),
            size = 5, hjust = 0) +
  scale_x_continuous(limits = c(0,2.5), 
                     expand = c(0,0)) +
  scale_fill_manual(
    values = c("Suppressed" = "#4ca763", "Activated" = "#497fb2"), # 颜色设置
    name = "Status"
  ) +
  scale_alpha_manual(
    values = c("0.4" = 0.4, "0.6" = 0.6, "0.8" = 0.8, "1" = 1),
    labels = c("0.4" = "NS", "0.6" = "< 0.05", "0.8" = "< 0.01", "1" = "< 0.001"),
    name = expression(italic(P)["adj"])
  ) +
  labs(
    title = "GSEA Enrichment of KEGG Pathways",
    x = "Absolute value of Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.key.spacing.y = unit(0.01,"cm"),
        legend.position = "right",
        legend.direction = "vertical") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(
    fill = guide_legend(override.aes = list(color = NA)),
    alpha = guide_legend(override.aes = list(color = NA))) +
  theme(aspect.ratio = 3/2)
#显示图片
gsea_kegg_barplot

#保存图片
save_plot(gsea_kegg_barplot, "KEGG.GSEA.barplot")









# 筛选显著性前 10 的GO功能，并定义透明度
gsea_go_bardata <-
  gsea_go_out %>%
  filter(stats %in% c("Suppressed", "Activated")) %>%
  mutate(
    transparency = case_when(
      p.adjust < 0.001 ~ 1,          # 高度显著
      p.adjust < 0.01 ~ 0.8,        # 显著
      p.adjust < 0.05 ~ 0.6,        # 中等显著
      TRUE ~ 0.4                    # 不显著
    )
  ) %>%
  arrange(stats, desc(abs(NES))) %>% # 按 NES 绝对值排序
  group_by(stats) %>%
  slice_head(n = 10) %>% # 每组取前10
  ungroup()

# 将透明度转换为离散型变量
gsea_go_bardata <- gsea_go_bardata %>%
  mutate(transparency_discrete = as.factor(transparency))

#GO的GSEA可视化绘图
gsea_go_barplot <-
  ggplot(gsea_go_bardata, aes(x = abs(NES), y = fct_reorder(Description, stats))) +
  geom_col(aes(fill = stats, alpha = transparency_discrete),
           position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.45) +
  geom_text(aes(x = 0.05, y = Description, label = Description),
            size = 5, hjust = 0) +
  scale_x_continuous(limits = c(0,2.5), 
                     expand = c(0,0)) +
  scale_fill_manual(
    values = c("Suppressed" = "#4ca763", "Activated" = "#497fb2"), # 颜色设置
    name = "Status"
  ) +
  scale_alpha_manual(
    values = c("0.4" = 0.4, "0.6" = 0.6, "0.8" = 0.8, "1" = 1),
    labels = c("0.4" = "NS", "0.6" = "< 0.05", "0.8" = "< 0.01", "1" = "< 0.001"),
    name = expression(italic(P)["adj"])
  ) +
  labs(
    title = "GSEA Enrichment of GO Terms",
    x = "Absolute value of Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.length.x = unit(0.3,"cm")) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.key.spacing.y = unit(0.01,"cm"),
        legend.position = "right",
        legend.direction = "vertical") +
  guides(x = "axis_truncated", y = "axis_truncated") +
  guides(
    fill = guide_legend(override.aes = list(color = NA)),
    alpha = guide_legend(override.aes = list(color = NA))) +
  theme(aspect.ratio = 3/2)

#显示图片
gsea_go_barplot
#保存图片
save_plot(gsea_go_barplot, "GO.GSEA.barplot")







# 定义 gsInfo 函数：用于提取某个基因集（pathway）的 GSEA 运行轨迹和相关信息
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  # 如果传入的是数字索引，则根据 result 表选取对应的基因集 ID
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  # 取出该基因集包含的基因
  geneSet <- object@geneSets[[geneSetID]]
  # GSEA 权重参数 exponent（默认 p=1）
  exponent <- object@params[["exponent"]]
  # 计算 GSEA 轨迹（running score 等数据），并转为数据框
  df <- gseaScores(geneList, geneSet, exponent, fortify=T)
  df$ymin <- 0# 初始化用于标记 peak 的最小值
  df$ymax <- 0# 初始化 peak 的最大值
  
  # position == 1 的位置是“命中基因”（hit）
  pos <- df$position == 1
  # h 为 peak 线条的高度，用于可视化
  h <- diff(range(df$runningScore))/20
  
  
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList# 保存排序后的 geneList 值
  
  # 如果没有基因对应 symbol，则使用基因 ID
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  # 添加 GSEA 结果中的详细信息
  df$Description <- object@result[geneSetID, "Description"]
  df$ID <- object@result[geneSetID, "ID"]
  df$NES <- object@result[geneSetID, "NES"]
  df$pvalue <- object@result[geneSetID, "pvalue"]
  df$p.adjust <- object@result[geneSetID, "p.adjust"]
  df$setSize <- object@result[geneSetID, "setSize"]
  # 计算核心富集基因数目
  df$core_gene_num <- sapply(strsplit(object@result[geneSetID, "core_enrichment"], "/"), length) 
  
  return(df)
}
# 导入 gseaScores 函数（DOSE 包内部函数）
gseaScores <- getFromNamespace("gseaScores", "DOSE")

# 定义提取基因集信息的函数
# 定义函数用于提取多个或单个 gene set 的轨迹数据
getGeneSetData <- function(x, geneSetID) {
  if (length(geneSetID) == 1) {
    # 如果 geneSetID 是单个值，直接调用 gsInfo
    # 单个 pathway
    gsdata <- gsInfo(x, geneSetID)
  } else {
    # 如果 geneSetID 是多个值，使用 lapply 循环调用 gsInfo，并将结果合并
    # 多个 pathway：循环 gsInfo 再合并
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  return(gsdata)
}





# KEGG GSEA 显著性前5通路可视化

# 提取 KEGG 前 5 个显著 pathway 的 GSEA 数据
gsea_kegg_data <- getGeneSetData(gsea_kegg, geneSetID = c(1:5))

# 生成 label（Description + ID）
gsea_kegg_data$label <- sprintf("%s (%s)", gsea_kegg_data$Description, gsea_kegg_data$ID)




# 可视化

#上半部分：Running Score 曲线
gsea_kegg_1 <-
  ggplot(gsea_kegg_data) +
  geom_line(aes(x = x, y = runningScore, color = label), linewidth = 0.75, lineend = "round", show.legend = F) +
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.5, lineend = "round") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  scale_colour_viridis_d() +
  labs(x = NULL,
       y = "RunningScore") +
  theme_linedraw(base_size = 18) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1)) +
  theme(axis.text.y = element_text(color = "black", angle = 90, hjust = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.border = element_rect(fill = NA,linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  annotate("text", x = -Inf, y = -Inf, label = "WT", size = 9, hjust = -0.2, vjust = -0.5) +  
  annotate("text", x = Inf, y = -Inf, label = "cp2-2", size = 8, hjust = 1.2, vjust = -0.5)
gsea_kegg_1

#下半部分：hit gene 的位置
gsea_kegg_2 <-
  ggplot(gsea_kegg_data) +
  geom_linerange(aes(x = x, ymin = ymin / 2, ymax = ymax / 2, color = label),show.legend = F) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 21000), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_viridis_d() +
  facet_wrap(~ label, ncol = 1) +
  labs(x = "Position in the Rank List",
       y = NULL) +
  theme(strip.text = element_text(hjust = 0, size = rel(1.5)),
        panel.spacing.y = unit(0, "cm")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = rel(1.5)),
        axis.title.x = element_text(color = "black", size = rel(1.5)),
        axis.line.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.length.x = unit(0.2,"cm")
  ) +
  guides(x = "axis_truncated")
gsea_kegg_2


# 拼接上下图
gsea_kegg_plot <- gsea_kegg_2 %>% insert_top(gsea_kegg_1, height = 0.95)
#显示图片
gsea_kegg_plot

# 保存 KEGG 图
save_plot(gsea_kegg_plot, "KEGG.GSEAplot")












# GO GSEA 显著性前5通路可视化

gsea_go_data <- getGeneSetData(gsea_go, geneSetID = c(1:5))
gsea_go_data$label <- sprintf("%s (%s)", gsea_go_data$Description, gsea_go_data$ID)

# 可视化

#上半部分：Running Score 曲线
gsea_go_1 <-
  ggplot(gsea_go_data) +
  geom_line(aes(x = x, y = runningScore, color = label), linewidth = 0.75, lineend = "round", show.legend = F) +
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.5, lineend = "round") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  scale_colour_viridis_d() +
  labs(x = NULL,
       y = "RunningScore") +
  theme_linedraw(base_size = 18) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1)) +
  theme(axis.text.y = element_text(color = "black", angle = 90, hjust = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.border = element_rect(fill = NA,linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  annotate("text", x = -Inf, y = -Inf, label = "WT", size = 9, hjust = -0.2, vjust = -0.5) +  
  annotate("text", x = Inf, y = -Inf, label = "cp2-2", size = 8, hjust = 1.2, vjust = -0.5)
gsea_go_1

#下半部分：hit gene 的位置
gsea_go_2 <-
  ggplot(gsea_go_data) +
  geom_linerange(aes(x = x, ymin = ymin / 2, ymax = ymax / 2, color = label),show.legend = F) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 21000), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_viridis_d() +
  facet_wrap(~ label, ncol = 1) +
  labs(x = "Position in the Rank List",
       y = NULL) +
  theme(strip.text = element_text(hjust = 0, size = rel(1.5)),
        panel.spacing.y = unit(0, "cm")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = rel(1.5)),
        axis.title.x = element_text(color = "black", size = rel(1.5)),
        axis.line.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5, lineend = "round"),
        axis.ticks.length.x = unit(0.2,"cm")
  ) +
  guides(x = "axis_truncated")
gsea_go_2


# 拼接上下图
gsea_go_plot <- gsea_go_2 %>% insert_top(gsea_go_1, height = 0.95)
#显示图片
gsea_go_plot
#保存图片
save_plot(gsea_go_plot, "GO.GSEAplot")








library(ggh4x)

# 生成单个通路或者功能的GSEA折线图

#图片函数
generate_gsea_plot <- function(object, geneSetIDs, type = "KEGG") {
  
  # 内部函数：生成图表并保存
  save_plot_gsea <- function(plot_object, file_name, type) {
    # 确保文件夹存在，根据type创建不同的文件夹
    folder_name <- ifelse(type == "KEGG", "GSEA_KEGG", "GSEA_GO")
    if (!dir.exists(folder_name)) {
      dir.create(folder_name)
    }
    
    # 设置文件路径
    file_path <- file.path("C:/Users/XZQ/Desktop/NCBI/Plant communications_RNA_seq/rna-seq-project/enrich",
                           folder_name,
                           paste0(file_name))
    
    # 保存为PDF
    ggsave(paste0(file_path, ".pdf"), plot_object, width = 9, height = 6)
    # 保存为PNG
    png_filename <- paste0(file_path, ".png")
    png(png_filename, width = 2250, height = 1500, res = 300)
    print(plot_object)
    dev.off()
  }
  
  # 获取指定通路的数据，并生成图表
  plot_for_geneSetID <- function(object, geneSetID, type) {
    # 获取相关信息
    gsea_data <- gsInfo(object, geneSetID)
    
    # 生成数据框，用于表格显示
    pd <- data.frame(
      Value = c(round(unique(gsea_data$NES), 4),
                format(unique(gsea_data$pvalue), digits = 4), 
                format(unique(gsea_data$p.adjust), digits = 4)),
      row.names = c("NES", "pvalue", "p.adjust"))
    
    # 生成标题
    gsea_title <- paste0(unique(gsea_data$Description), " (", unique(gsea_data$ID), ")")
    
    # 创建图表
    gsea_plot <- 
      ggplot(gsea_data) +
      geom_linerange(aes(x = x, ymin = ymin / 2, ymax = ymax / 2), color = "black", show.legend = F) +
      geom_line(aes(x = x, y = runningScore), color = "green", linewidth = 1, show.legend = F) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      ggtitle(gsea_title) +
      labs(x = "Position in the Rank List", y = "Enrichment Score") +
      theme_minimal(base_size = 18) +
      theme(plot.title = element_text(size = rel(1))) +
      theme(axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black")) +
      theme(axis.ticks.length = unit(2, "mm"),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            ggh4x.axis.ticks.length.minor = rel(0.5)) +
      theme(panel.border = element_rect(fill = NA, linewidth = 0.75)) +
      guides(x = "axis_minor", y = "axis_minor") +
      theme(aspect.ratio = 2/3)
    
    # 添加表格注释
    gsea_plot_with_table <- 
      gsea_plot + 
      annotation_custom(tableGrob(pd, cols = NULL, 
                                  theme = ttheme_default(rowhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontface = "bold"),
                                                                        bg_params = list(fill = "grey80")))),
                        xmin = quantile(gsea_plot$data$x, 0.05),
                        xmax = quantile(gsea_plot$data$x, 0.5),
                        ymin = quantile(gsea_plot$data$runningScore, 0.05),
                        ymax = quantile(gsea_plot$data$runningScore, 0.5))
    
    # 确定文件名
    if (type == "GO") {
      # GO类型的文件名为 GO_xxxxxxxx
      file_name <- gsub(":", "_", paste0(unique(gsea_data$ID)))
    } else {
      # 其他类型直接使用ID
      file_name <- unique(gsea_data$ID)
    }
    
    # 保存图表，并将其放入指定的文件夹中
    save_plot_gsea(gsea_plot_with_table, file_name, type)
  }
  
  # 如果传入的是数字（编号范围或单个编号）
  if (is.numeric(geneSetIDs)) {
    for (geneSetID in geneSetIDs) {
      plot_for_geneSetID(object, geneSetID, type)
    }
  }
  
  # 如果传入的是通路ID（字符串）
  if (is.character(geneSetIDs)) {
    for (geneSetID in geneSetIDs) {
      geneSetIndex <- which(object@result$ID == geneSetID)
      plot_for_geneSetID(object, geneSetIndex, type)
    }
  }
}

#创造并保存图片
generate_gsea_plot(gsea_kegg, 1:5, type = "KEGG") 
#创造并保存图片
generate_gsea_plot(gsea_go, 1:5, type = "GO") 
