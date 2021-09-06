########################################################################
# Functions for analysis of COVID-19 adult severity manuscript
########################################################################

load_packages <- function(){
  library(ggplot2)
  library(viridis)
  library(ggpubr)
  library(stringr)
  library(DESeq2)
  library(biomaRt)
  library(grid)
  library(gridExtra)
  library(PCAtools)
  library(tidyverse)
  library(forcat)
  library(pheatmap)
}

paths_fun <- function(deg_list){
  up_genes <- deg_list[deg_list$log2FoldChange > 0 & deg_list$padj < 0.05,]
  down_genes <- deg_list[deg_list$log2FoldChange < 0 & deg_list$padj < 0.05,]
  # identify pathways from genes up regulated
  paths_up <- gost(up_genes$gene,
                   organism = 'hsapiens',
                   correction_method = 'fdr',
                   ordered_query = T,
                   source = "GO:BP")
  # identify pathways from genes down regulated
  paths_down <- gost(down_genes$gene,
                     organism = 'hsapiens',
                     correction_method = 'fdr',
                     ordered_query = T,
                     source = "GO:BP")
  # filter so only significant 
  paths_up <- paths_up$result[paths_up$result$p_value < 0.05,]
  paths_down <- paths_down$result[paths_down$result$p_value < 0.05,]
  return(list(paths_up = paths_up,
              paths_down = paths_down))
}

read_revigo <- function(file_name){
  f <- read.csv(file = file_name, header = T)
  f <- f[f$Eliminated==' False',]
  f <- f[order(f$Value),]
  # add on the intersection size from the original results file 
  #original_pathways <- original_pathways[match(f$term_ID, original_pathways$term_id),]
  #f$intersection <- original_pathways$intersection_size/original_pathways$term_size
  return(f)
}

s_cell <- function(df, cell_name, cell_name_col){
  p_df <- data.frame(Severity = df$Severity.Sample, cell = df[,cell_name_col])
  m <- max(p_df$cell)+0.05
  p <- ggplot(p_df, aes(x = as.factor(Severity),
                        y = cell,
                        color = Severity))+
    geom_boxplot()+
    geom_jitter(width = 0.2)+
    theme_bw()+
    labs(x = "", y = paste("Proportion of ", cell_name, sep = ""))+
    stat_compare_means(label = 'p.format',
                       comparisons = comparisons,
                       p.adjust.method = "bonferroni")+
    scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(color = 'black'), 
          legend.position = "none")
  return(p)
}

pca_plot_fun <- function(pc1, pc2, color, num1, num2, col_vec, importance, col_name){
  df <- data.frame(dim_1 = pc1,
                   dim_2 = pc2,
                   col = color)
  ggplot(df, aes(x = dim_1, y = dim_2, color = col))+
    geom_point()+
    theme_bw()+
    stat_ellipse()+
    scale_color_manual(values = col_vec)+
    labs(x =paste("PC", num1, ": ", round(importance[2,num1]*100, 2), "%", sep = ""),
         y = paste("PC", num2, ": ", round(importance[2,num2]*100, 2), "%", sep = ""),
         color = col_name)+
    theme(legend.position="bottom", 
          legend.title = element_blank())
}

treatment_plot <- function(df, cell_name, cell_name_col, treatment, treatment_col, col_vec){
  p_df <- data.frame(severity = df$Severity.Sample,
                     t = df[,treatment_col], 
                     cell = df[,cell_name_col])
  m <- max(p_df$cell)+0.05
  p <- ggplot(p_df, aes(x = as.factor(severity),
                        y = cell,
                        color = t))+
    geom_boxplot()+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2))+
    theme_bw()+
    labs(x = "", 
         y = paste("Proportion of ", cell_name, sep = ""),
         color = treatment)+
    stat_compare_means(label = 'p.format', label.y = m)+
    scale_color_manual(values = col_vec)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(color = 'black'))
  return(p)
}


cross_plot <- function(df1, df2, lims, title){
  df <- data.frame()
  df <- data.frame(ID = rownames(df1), 
                   Gene = df1$gene, 
                   Cell_P = df1$padj,
                   Cell_LFC = df1$log2FoldChange, 
                   Treatment_P = df2$padj[match(rownames(df1), rownames(df2))], 
                   Treatment_LFC = df2$log2FoldChange[match(rownames(df1), rownames(df2))])
  
  df$Significance <- "NS"
  df$Significance[df$Cell_P < 0.05] <- 'Cell corrected'
  df$Significance[df$Treatment_P < 0.05] <- 'Treatment corrected'
  df$Significance[df$Cell_P < 0.05 & df$Treatment_P < 0.05] <- 'Both'
  df$Significance <- as.factor(df$Significance)
  df <- df[order(df$Significance),]
  df$alpha <- 0.1
  df$alpha[df$Significance %in% c('Both', 'Cell corrected', 'Treatment corrected')]<- 1
  
  # plot function 
  p <- ggplot(df, aes(x = Cell_LFC, 
                      y = Treatment_LFC, 
                      color = fct_relevel(Significance, c('Both', 'Cell corrected', 'Treatment corrected', 'NS')),
                      alpha = alpha))+
    theme_bw()+
    geom_point()+
    scale_color_manual(values = c("red", "orange", '#2ca25f', 'darkgrey'))+
    labs(x = "Cell corrected LFC", y = "Treatment corrected LFC", color = "Significance")+
    theme(legend.text.align = 0, 
          legend.position = 'bottom', 
          plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(limits = lims)+
    scale_y_continuous(limits = lims)+
    geom_label_repel(
      data = subset(df, Significance == 'Both' ),
      aes(label = Gene),
      size = 3,
      max.overlaps = 20,
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines"), 
      fill = 'white',
      color = 'black')+
    guides(alpha = FALSE)+
    ggtitle(title)
  print(p)
  return(df)
}

deseq.results <- function(res, gene.names, outfile){
  res <- res[order(res$padj),]
  res.df <- data.frame(res)
  rownames(res.df) <- rownames(res)
  
  gene.names <- gene.names[match(rownames(res.df), gene.names$ensembl_gene_id),]
  res.df$gene <- gene.names$external_gene_name
  res.df$Significant <- 0
  res.df$Significant[res.df$padj < 0.05] <- 1
  
  
  # create volcano plot column
  res.df$Color <- 'black'
  res.df$Color[abs(res.df$log2FoldChange) > 1] <- 'limegreen'
  res.df$Color[res.df$padj < 0.05] <- 'gold'
  res.df$Color[res.df$padj < 0.05 & 
                 abs(res.df$log2FoldChange) > 1] <- 'red'
  
  res.df$Color <- factor(res.df$Color, levels = c("black", 'limegreen', "gold", 'red'))
  
  res.df <- res.df[!(is.na(res.df$log2FoldChange)),]
  res.df <- res.df[!(is.na(res.df$padj)),]
  print(range(res.df$log2FoldChange))
  write.csv(file = paste(outfile, '.csv', sep = ""), res.df)
  return(res.df)
}
