##################################################################################################
# Analysis for SARS-CoV-2 severity paper
##################################################################################################

source(file = 'functions.R')
load_packages()

########################################################################
# Read in data
########################################################################

# read in metadata and counts
load("counts.RData", verbose = T)
counts <- P200163_gene_counts
meta.interest <- read.csv(file = 'metadata.csv', header = T)
counts <- counts[,match(meta.interest$ID, colnames(counts))]

########################################################################
# Plot levels of cells estimated from CIBERSORT
########################################################################

cell_types <- meta.interest[,c(8, 12, 14, 16:19)]

# immune cell levels v severity
comparisons <- list(c("Severe", "Moderate"),
                    c("Severe", "Mild"),
                    c("Moderate", "Mild"))

s_cd8 <- s_cell(cell_types,  "CD8 T cells", "T.cells.CD8")
s_mono <- s_cell(cell_types,  "Monocytes", "Monocytes")
s_cd4 <- s_cell(cell_types,  "CD4 T cells", "CD4.T.cells")
s_neut <- s_cell(cell_types,  "Neutrophils", "Neutrophils")
s_b <- s_cell(cell_types,  "B cells", "B.cells")
s_nk <- s_cell(cell_types,  "NK cells", "NK.cells")

grid.arrange(s_cd4, s_cd8, s_b, s_neut, s_mono, s_nk, ncol = 3)

########################################################################
#  DESeq Analysis
########################################################################

# create DESeq object
des.obj <- DESeqDataSetFromMatrix(countData = counts, 
                                  colData = meta.interest, 
                                  design = ~1)  
colnames(des.obj) <- colnames(counts)

# remove lowly expressed
keep <- rowSums(counts(des.obj) >= 20 ) >= 3
des.obj <- des.obj[keep,]
des.obj <- estimateSizeFactors(des.obj)

# create gene.data object containing gene names and ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene.names <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), 
                    filters = 'ensembl_gene_id', 
                    values = rownames(des.obj) , 
                    mart = ensembl)
gene.names <- gene.names[match(rownames(des.obj), gene.names$ensembl_gene_id),]

# extract normalised counts
norm <- counts(des.obj, normalized = TRUE)
colData(des.obj)$Group_Severity <- colData(des.obj)$Severity.Sample

########################################################################
# PCA
########################################################################

# principal component analysis 
pc <- prcomp(t(norm), scale. = T)
p <- summary(pc)
pc <- data.frame(pc$x,
                 group_severity = colData(des.obj)$Group_Severity,
                 sex = colData(des.obj)$Sex,
                 age = colData(des.obj)$Age)

pc1_2_group_sev <- pca_plot_fun(pc1 = pc$PC1, pc2 = pc$PC2,
                                color = pc$group_severity,
                                num1 = 1, num2 = 2,
                                col_vec = c("#fdae61", "#f46d43", '#d73027'),
                                importance = p$importance,
                                col_name = "Disease/severity group")

pc3_4_group_sev <- pca_plot_fun(pc1 = pc$PC3, pc2 = pc$PC4,
                                color = pc$group_severity,
                                num1 = 3, num2 = 4,
                                col_vec = c("#fdae61", "#f46d43", '#d73027'),
                                importance = p$importance,
                                col_name = "Disease/severity group")

pc5_6_group_sev <- pca_plot_fun(pc1 = pc$PC5, pc2 = pc$PC6,
                                color = pc$group_severity,
                                num1 = 5, num2 = 6,
                                col_vec = c("#fdae61", "#f46d43", '#d73027'),
                                importance = p$importance,
                                col_name = "Disease/severity group")

pc1_2_sex <- pca_plot_fun(pc1 = pc$PC1, pc2 = pc$PC2,
                          color = pc$sex,
                          num1 = 1, num2 = 2,
                          col_vec = c("#a6611a", "#018571"),
                          importance = p$importance,
                          col_name = "Sex")
#
pc3_4_sex <- pca_plot_fun(pc1 = pc$PC3, pc2 = pc$PC4,
                          color = pc$sex,
                          num1 = 3, num2 = 4,
                          col_vec = c("#a6611a", "#018571"),
                          importance = p$importance,
                          col_name = "Sex")
#
pc5_6_sex <- pca_plot_fun(pc1 = pc$PC5, pc2 = pc$PC6,
                          color = pc$sex,
                          num1 = 5, num2 = 6,
                          col_vec = c("#a6611a", "#018571"),
                          importance = p$importance,
                          col_name = "Sex")

pc1_2_age <- ggplot(pc, aes(x = PC1, y = PC2, color = age))+
  geom_point()+
  theme_bw()+
  labs(x =paste("PC1: ", round(p$importance[2,1]*100, 2), "%", sep = ""),
       y = paste("PC2: ", round(p$importance[2,2]*100, 2), "%", sep = ""),
       color = "Age")+
  scale_color_viridis_c(direction= -1)+
  theme(legend.position = 'bottom')


pc3_4_age <- ggplot(pc, aes(x = PC3, y = PC4, color = age))+
  geom_point()+
  theme_bw()+
  labs(x =paste("PC3: ", round(p$importance[2,3]*100, 2), "%", sep = ""),
       y = paste("PC4: ", round(p$importance[2,4]*100, 2), "%", sep = ""),
       color = "Age")+
  scale_color_viridis_c(direction= -1)+
  theme(legend.position = 'bottom')


pc5_6_age <- ggplot(pc, aes(x = PC5, y = PC6, color = age))+
  geom_point()+
  theme_bw()+
  labs(x =paste("PC5: ", round(p$importance[2,5]*100, 2), "%", sep = ""),
       y = paste("PC6: ", round(p$importance[2,6]*100, 2), "%", sep = ""),
       color = "Age")+
  scale_color_viridis_c(direction= -1)+
  theme(legend.position = 'bottom')

# plot all together
grid.arrange(pc1_2_group_sev,
             pc3_4_group_sev,
             pc5_6_group_sev,
             pc1_2_sex,
             pc3_4_sex,
             pc5_6_sex,
             pc1_2_age,
             pc3_4_age,
             pc5_6_age,
             ncol = 3)

########################################################################
# 1) Evaluate impacts of immunomodulatory treatment on moderate COVID
########################################################################

# take out only COVID 
des_covid_treatment_comp <- des.obj
# remove samples receiving steroids and immunmoodulatory 
des_covid_treatment_comp <- des_covid_treatment_comp[,!(colnames(des_covid_treatment_comp) %in% des_covid_treatment_comp$ID[des_covid_treatment_comp$steroids=="Y" & des_covid_treatment_comp$monoclonal=="Y"])]

# only keep moderate samples 
des_covid_treatment_comp <- des_covid_treatment_comp[,colnames(des_covid_treatment_comp) %in% des_covid_treatment_comp$ID[des_covid_treatment_comp$Severity.Sample=="Moderate"]]
des_covid_treatment_comp$Sex <- as.factor(des_covid_treatment_comp$Sex)
des_covid_treatment_comp$steroids <- as.factor(des_covid_treatment_comp$steroids)

# update design
design(des_covid_treatment_comp) <- formula(~0 + steroids + Age + Sex)

# run DESeq
des_covid_treatment_comp <- DESeq(des_covid_treatment_comp)
des_covid_treatment_comp <- des_covid_treatment_comp[which(mcols(des_covid_treatment_comp)$betaConv),]

# extract results comparing severe to moderate
res_mod_steroids <- results(des_covid_treatment_comp,
                            independentFiltering=TRUE,
                            alpha=0.05,
                            pAdjustMethod="BH",
                            parallel=FALSE,
                            contrast=c('steroids', "Y", "N"))

res_mod_steroids <- deseq.results(res_mod_steroids, gene.names, "Moderate_Steroids")
table(res_mod_steroids$Significant)
dim(res_mod_steroids[res_mod_steroids$Significant==1 & res_mod_steroids$log2FoldChange > 0,])
dim(res_mod_steroids[res_mod_steroids$Significant==1 & res_mod_steroids$log2FoldChange < 0,])

# results are to be exported and pathway analysis ran using IPA

########################################################################
# 2) Severity of COVID
# Pairwise comparisons correcting for immune cells or treatment
########################################################################

# First make DESeq two models
# this one will correct for immunomodulatory treatment
des_covid_t <- des.obj
des_covid_t$Severity.Sample <- as.factor(des_covid_t$Severity.Sample)

# this one will correct for cell proportions
des_covid_cells <- des.obj
des_covid_cells$Severity.Sample <- as.factor(des_covid_cells$Severity.Sample)

# update design for treatment
des_covid_t$steroids <- as.factor(des_covid_t$steroids)
des_covid_t$monoclonal <- as.factor(des_covid_t$monoclonal)
des_covid_t$interferon <- as.factor(des_covid_t$interferon)
des_covid_t$Sex <- as.factor(des_covid_t$Sex)
des_covid_cells$Sex <- as.factor(des_covid_cells$Sex)

design(des_covid_t) <- formula(~0 + Age + 
                                 Sex + 
                                 Severity.Sample + 
                                 steroids + 
                                 monoclonal +
                                 interferon)
# run DESeq
des_covid_t <- DESeq(des_covid_t)
des_covid_t <- des_covid_t[which(mcols(des_covid_t)$betaConv),]

# update design for cells
design(des_covid_cells) <- formula(~0 + Age + 
                                     Sex + 
                                     Severity.Sample + 
                                     Monocytes +
                                     Neutrophils +
                                     B.cells +
                                     CD4.T.cells +
                                     NK.cells +
                                     T.cells.CD8)
# run DESeq
des_covid_cells <- DESeq(des_covid_cells)
des_covid_cells <- des_covid_cells[which(mcols(des_covid_cells)$betaConv),]

########################################################################
# First comparison is moderate COVID-19 vs mild COVID-19 

# extract results comparing mild to moderate correcting for treatment
res_mild_moderate_treatment <-  results(des_covid_t, 
                                        independentFiltering=TRUE,
                                        alpha=0.05, 
                                        pAdjustMethod="BH", 
                                        parallel=FALSE, 
                                        contrast=c("Severity.Sample","Moderate", "Mild"))

res_mild_moderate_treatment <- deseq.results(res_mild_moderate_treatment, 
                                             gene.names, 
                                             "Moderate_Mild_Treatment")

table(res_mild_moderate_treatment$Significant)
dim(res_mild_moderate_treatment[res_mild_moderate_treatment$Significant==1 & res_mild_moderate_treatment$log2FoldChange > 0,])
dim(res_mild_moderate_treatment[res_mild_moderate_treatment$Significant==1 & res_mild_moderate_treatment$log2FoldChange < 0,])

# extract results comparing mild to moderate correcting for cells
res_mild_moderate_cells <-  results(des_covid_cells, 
                                    independentFiltering=TRUE,
                                    alpha=0.05, 
                                    pAdjustMethod="BH", 
                                    parallel=FALSE, 
                                    contrast=c("Severity.Sample", "Moderate", "Mild"))

res_mild_moderate_cells <- deseq.results(res_mild_moderate_cells, 
                                         gene.names, 
                                         "Moderate_Mild_Cells")

table(res_mild_moderate_cells$Significant)
dim(res_mild_moderate_cells[res_mild_moderate_cells$Significant==1 & res_mild_moderate_cells$log2FoldChange > 0,])
dim(res_mild_moderate_cells[res_mild_moderate_cells$Significant==1 & res_mild_moderate_cells$log2FoldChange < 0,])

# cross plot comparing mild vs moderate with and without cell correction (baseline treatment correction)
cross_mild_mod <- cross_plot(res_mild_moderate_cells, res_mild_moderate_treatment, c(-6,6), "Moderate COVID-19 vs Mild COVID-19")
table(cross_mild_mod$Significance)

# export results and run IPA pathway analysis

########################################################################
# Next comparison is severe COVID-19 vs mild COVID-19 

# results from treatment correction
res_mild_severe_treatment <-  results(des_covid_t, 
                                      independentFiltering=TRUE,
                                      alpha=0.05, 
                                      pAdjustMethod="BH", 
                                      parallel=FALSE, 
                                      contrast=c("Severity.Sample", "Severe", "Mild"))

res_mild_severe_treatment <- deseq.results(res_mild_severe_treatment, 
                                           gene.names, 
                                           "Severe_Mild_Treatment")
table(res_mild_severe_treatment$Significant)
dim(res_mild_severe_treatment[res_mild_severe_treatment$Significant==1 & res_mild_severe_treatment$log2FoldChange > 0,])
dim(res_mild_severe_treatment[res_mild_severe_treatment$Significant==1 & res_mild_severe_treatment$log2FoldChange < 0,])

# results from cell correction 
res_mild_severe_cells <-  results(des_covid_cells, 
                                  independentFiltering=TRUE,
                                  alpha=0.05, 
                                  pAdjustMethod="BH", 
                                  parallel=FALSE, 
                                  contrast=c("Severity.Sample", "Severe", "Mild"))

res_mild_severe_cells <- deseq.results(res_mild_severe_cells, 
                                       gene.names, 
                                       "Severe_Mild_Cells")
table(res_mild_severe_cells$Significant)
dim(res_mild_severe_cells[res_mild_severe_cells$Significant==1 & res_mild_severe_cells$log2FoldChange > 0,])
dim(res_mild_severe_cells[res_mild_severe_cells$Significant==1 & res_mild_severe_cells$log2FoldChange < 0,])

# cross plot comparing results 
cross_mild_mod <- cross_plot(res_mild_severe_cells, res_mild_severe_treatment, c(-10,10), "Severe COVID-19 vs Mild COVID-19")
table(cross_mild_mod$Significance)

# export results and run IPA 

########################################################################
# Final pairwise comparison is severe COVID-19 vs moderate COVID-19 

# get results correcting for treatment
res_mod_sev_treatment <-  results(des_covid_t, 
                                  independentFiltering=TRUE,
                                  alpha=0.05, 
                                  pAdjustMethod="BH", 
                                  parallel=FALSE, 
                                  contrast=c("Severity.Sample", "Severe","Moderate"))

res_mod_sev_treatment <- deseq.results(res_mod_sev_treatment, 
                                       gene.names, 
                                       "Severe_Moderate_Treatment")

table(res_mod_sev_treatment$Significant)
dim(res_mod_sev_treatment[res_mod_sev_treatment$Significant==1 & res_mod_sev_treatment$log2FoldChange > 0,])
dim(res_mod_sev_treatment[res_mod_sev_treatment$Significant==1 & res_mod_sev_treatment$log2FoldChange < 0,])

# get results correcting for cells
res_mod_sev_cells <-  results(des_covid_cells, 
                              independentFiltering=TRUE,
                              alpha=0.05, 
                              pAdjustMethod="BH", 
                              parallel=FALSE, 
                              contrast=c("Severity.Sample", "Severe","Moderate"))

res_mod_sev_cells <- deseq.results(res_mod_sev_cells, 
                                   gene.names, 
                                   "Severe_Moderate_Cells")

table(res_mod_sev_cells$Significant)

# cross plot comparing mild vs severe with and without cell correction (baseline treatment correction)
cross_mod_sev <- cross_plot(res_mod_sev_cells, res_mod_sev_treatment, c(-6,6), "Moderate COVID-19 vs Severe COVID-19")
p <- ggplot(cross_mod_sev, aes(x = Cell_LFC, 
                               y = Treatment_LFC, 
                               color = fct_relevel(Significance, c('Both', 'Treatment corrected', 'NS')),
                               alpha = alpha))+
  theme_bw()+
  geom_point()+
  scale_color_manual(values = c('red', '#2ca25f', 'darkgrey'))+
  labs(x = "Cell corrected LFC", y = "Treatment corrected LFC", color = "Significance")+
  theme(legend.text.align = 0, 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(-8.8, 8.8))+
  scale_y_continuous(limits =  c(-8.8, 8.8))+
  geom_label_repel(
    data = subset(cross_mod_sev, abs(Treatment_LFC) > 3 & Treatment_P < 0.05 | 
                    Significance == "Both"),
    aes(label = Gene),
    size = 3,
    max.overlaps = 20,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"), 
    fill = 'white',
    color = 'black')+
  guides(alpha = FALSE)+
  ggtitle("Severe COVID-19 vs Moderate COVID-19")

print(p)
table(cross_mod_sev$Significance)

# run IPA for pathway analysis 

########################################################################
# 4) Severity of COVID, using an additive model (supplementary results)
########################################################################

#  are any genes SDE with age when accounting for severity?
des.covid <- des.obj
des.covid$Sex <- as.factor(des.covid$Sex)

# update design 
design(des.covid) <- formula(~ 0 + Age + Sex + severity.code:Age + severity.code)

# run deseq
des.covid <- DESeq(des.covid, betaPrior=FALSE)
des.covid <- des.covid[which(mcols(des.covid)$betaConv),]

# extract results 
# which genes are SDE with age whilst correcting for additive severity? 
res.covid.age <- results(des.covid, 
                         independentFiltering=TRUE,
                         alpha=0.05, 
                         pAdjustMethod="BH", 
                         parallel=FALSE, 
                         name='Age')

res.covid.age <- deseq.results(res.covid.age, gene.names, "COVID_Age")
table(res.covid.age$Significant)
dim(res.covid.age[res.covid.age$Significant==1 & res.covid.age$log2FoldChange > 0,])
dim(res.covid.age[res.covid.age$Significant==1 & res.covid.age$log2FoldChange < 0,])

######## create DESeq models with additive severity correcting for cells or treatment

# create DESeq model with additive severity and cells 
des.covid.cells <- des.obj
des.covid.cells$Sex <- as.factor(des.covid.cells$Sex)

# update design 
design(des.covid.cells) <- formula(~ 0 + Age + Sex + 
                                     Monocytes +
                                     Neutrophils +
                                     B.cells +
                                     CD4.T.cells +
                                     NK.cells +
                                     T.cells.CD8 + 
                                     severity.code)

# run deseq
des.covid.cells <- DESeq(des.covid.cells, betaPrior=FALSE)
des.covid.cells <- des.covid.cells[which(mcols(des.covid.cells)$betaConv),]

# genes SDE with additive severity whilst correcting for cell proportions
res.covid.sev.cells <- results(des.covid.cells, 
                               independentFiltering=TRUE,
                               alpha=0.05, 
                               pAdjustMethod="BH", 
                               parallel=FALSE, 
                               name='severity.code')

res.covid.sev.cells <- deseq.results(res.covid.sev.cells, gene.names, "COVID_Severity_Cell_Corrected")
table(res.covid.sev.cells$Significant)
dim(res.covid.sev.cells[res.covid.sev.cells$Significant==1 & res.covid.sev.cells$log2FoldChange > 0,])
dim(res.covid.sev.cells[res.covid.sev.cells$Significant==1 & res.covid.sev.cells$log2FoldChange < 0,])

# create DESeq model with additive severity and treatment 
des.covid.treatment <- des.obj
des.covid.treatment$Sex <- as.factor(des.covid.treatment$Sex)
des.covid.treatment$steroids <- as.factor(des.covid.treatment$steroids)
des.covid.treatment$monoclonal <- as.factor(des.covid.treatment$monoclonal)
des.covid.treatment$interferon <- as.factor(des.covid.treatment$interferon)

design(des.covid.treatment) <- formula(~ 0 + Age + 
                                         Sex + 
                                         severity.code + 
                                         steroids +
                                         monoclonal + 
                                         interferon)
# run deseq
des.covid.treatment <- DESeq(des.covid.treatment, betaPrior=FALSE)
des.covid.treatment <- des.covid.treatment[which(mcols(des.covid.treatment)$betaConv),]

#  results from correcting for treatment 
res.covid.treatment <- results(des.covid.treatment, 
                               independentFiltering=TRUE,
                               alpha=0.05, 
                               pAdjustMethod="BH", 
                               parallel=FALSE, 
                               name='severity.code')

res.covid.treatment <- deseq.results(res.covid.treatment, gene.names, "COVID_Severity_Treatment_Corrected")
table(res.covid.treatment$Significant)
dim(res.covid.treatment[res.covid.treatment$Significant==1 & res.covid.treatment$log2FoldChange > 0,])
dim(res.covid.treatment[res.covid.treatment$Significant==1 & res.covid.treatment$log2FoldChange < 0,])


####################################################################################################
# heatmap just with additive genes 
norm_additive_sde <- counts(des.covid.treatment, normalized = T)
# filter so top genes 
include <- rownames(res.covid.treatment)[res.covid.treatment$Significant==1 & abs(res.covid.treatment$log2FoldChange) > 2 & res.covid.treatment$padj < 0.0001]
norm_additive_sde <- norm_additive_sde[rownames(norm_additive_sde) %in% include,]
rownames(norm_additive_sde) <- res.covid.treatment$gene[match(rownames(norm_additive_sde), rownames(res.covid.treatment))]

treatment_df <- data.frame(steroids = des.covid.treatment$steroids, 
                           monoclonal = des.covid.treatment$monoclonal, 
                           interferon = des.covid.treatment$interferon, 
                           ID = des.covid.treatment$ID)

treatment_df$summary <- NA
treatment_df$summary[treatment_df$steroids=="Y"] <- 'Steroids'
treatment_df$summary[treatment_df$monoclonal=="Y"] <- 'Monoclonal'
treatment_df$summary[treatment_df$interferon=="Y"] <- 'Interferon'
treatment_df$summary[treatment_df$interferon=="Y" &
                       treatment_df$steroids=="Y" &
                       treatment_df$monoclonal=="N"] <- 'Interferon + Steroids'
treatment_df$summary[treatment_df$interferon=="Y" &
                       treatment_df$monoclonal=="Y" &
                       treatment_df$steroids=="N"] <- 'Interferon + Monoclonal'
treatment_df$summary[treatment_df$monoclonal=="Y" &
                       treatment_df$steroids=="Y" &
                       treatment_df$interferon=="N"] <- 'Steroids + Monoclonal'

annotation_additive <- data.frame(Treatment = treatment_df$summary[match(colnames(norm_additive_sde), treatment_df$ID)], 
                                  Age = des.covid.treatment$Age[match(colnames(norm_additive_sde), des.covid.treatment$ID)],
                                  Sex = des.covid.treatment$Sex[match(colnames(norm_additive_sde), des.covid.treatment$ID)], 
                                  Severity =  des.covid.treatment$Severity.Sample[match(colnames(norm_additive_sde), des.covid.treatment$ID)], 
                                  ID = des.covid.treatment$ID[match(colnames(norm_additive_sde), des.covid.treatment$ID)])

annotation_additive <- annotation_additive[order(annotation_additive$Severity),]
norm_additive_sde <- norm_additive_sde[,match(annotation_additive$ID, colnames(norm_additive_sde))]

rownames(annotation_additive) <- annotation_additive$ID

ann_colors <- list(Treatment = c("Steroids" = "#004529",
                                 "Steroids + Monoclonal" = "#41ab5d",
                                 "Interferon + Steroids" = "#addd8e", 
                                 "Monoclonal" = "#f7fcb9"), 
                   Age = c(low="white", middle = 'grey', high="black"),
                   Sex = c('Male' = "#dadaeb", 
                           'Female' = '#54278f'),
                   Severity = c('Mild' = '#f1b6da',
                                'Moderate' = '#de77ae',
                                'Severe' = '#8e0152'))

norm_additive_sde <- data.frame(t(norm_additive_sde))

newnames_1 <- lapply(
  colnames(norm_additive_sde),
  function(x) bquote(italic(.(x))))


t_additive <- pheatmap(t(log2(norm_additive_sde+1)), 
                       annotation_col = annotation_additive[,1:4],
                       #viridis(100, direction = -1),
                       legend =TRUE, 
                       main = "", 
                       annotation_colors = ann_colors, 
                       cluster_cols = F, 
                       cluster_rows = T,
                       show_colnames = F, 
                       annotation_legend_side = "bottom",
                       labels_row = as.expression(newnames_1))


grid::grid.newpage()
grid::grid.draw(t_additive$gtable)


##### cross plot comparing SDE genes with cell and with treatment correction 
additive_cross <- cross_plot(res.covid.sev.cells, res.covid.treatment, c(-5, 5), "Additive severity model")

table(additive_cross$Significance)

########################################################################
# 5) Genes SDE in additive and pairwise models 
########################################################################

norm_treatment <- counts(des_covid_t, normalized = T)

# make venn diagram looking at genes SDE with treatment correction

library(venn)
venn(x = list(moderate_mild_treatment = rownames(res_mild_moderate_treatment)[res_mild_moderate_treatment$Significant==1], 
              severe_mild_treatment = rownames(res_mild_severe_treatment)[res_mild_severe_treatment$Significant==1], 
              severe_moderate_treatment = rownames(res_mod_sev_treatment)[res_mod_sev_treatment$Significant==1], 
              additive_treatment = rownames(res.covid.treatment)[res.covid.treatment$Significant==1]), 
     zcolor = "style")

# which genes are SDE in all treatment corrections 
t_overlapping_all <- intersect(rownames(res.covid.treatment)[res.covid.treatment$Significant==1], 
                               intersect(rownames(res_mild_moderate_treatment)[res_mild_moderate_treatment$Significant==1], 
                                         intersect(rownames(res_mild_severe_treatment)[res_mild_severe_treatment$Significant==1], 
                                                   rownames(res_mod_sev_treatment)[res_mod_sev_treatment$Significant==1])))

##### add on LFC to check concordance
t_overlapping_all_df <- data.frame(ID = t_overlapping_all, 
                                   additive = res.covid.treatment$log2FoldChange[match(t_overlapping_all, rownames(res.covid.treatment))], 
                                   additive_p = res.covid.treatment$padj[match(t_overlapping_all, rownames(res.covid.treatment))],
                                   mild_mod = res_mild_moderate_treatment$log2FoldChange[match(t_overlapping_all, rownames(res_mild_moderate_treatment))], 
                                   mild_sev = res_mild_severe_treatment$log2FoldChange[match(t_overlapping_all, rownames(res_mild_severe_treatment))],
                                   mod_sev = res_mod_sev_treatment$log2FoldChange[match(t_overlapping_all, rownames(res_mod_sev_treatment))])



t_overlapping_all_df$disco <- apply(t_overlapping_all_df, 1, function(x){
  if(x[2] > 0 & x[4] > 0 & x[5] > 0 & x[6] > 0){
    disco = 'concordant'
  }
  else if(x[2] < 0 & x[4] < 0 & x[5] < 0 & x[6] < 0){
    disco = 'concordant'
  }
  else{
    disco = 'discordant'
  }
})

dim(t_overlapping_all_df[t_overlapping_all_df$additive > 0,])
dim(t_overlapping_all_df[t_overlapping_all_df$additive < 0,])

######## heatmap for the overlapping genes 

# filter so top genes 
t_overlapping_all_df <- t_overlapping_all_df[abs(t_overlapping_all_df$additive) > 2 & t_overlapping_all_df$additive_p < 0.0001,]


exp_treatment <- norm_treatment[t_overlapping_all_df$ID,]
exp_treatment <- exp_treatment[,match(colnames(des_covid_t), colnames(exp_treatment))]
exp_treatment <- data.frame(Severity = des_covid_t$Severity.Sample[match(colnames(exp_treatment), des_covid_t$ID)], 
                            t(exp_treatment))

exp_treatment <- exp_treatment[!(exp_treatment$Severity==""),]

colnames(exp_treatment)[2:ncol(exp_treatment)] <- gene.names$external_gene_name[match(colnames(exp_treatment)[2:ncol(exp_treatment)], gene.names$ensembl_gene_id)]

exp_treatment <- exp_treatment[order(exp_treatment$Severity),]

treatment_df <- data.frame(steroids = des_covid_t$steroids, 
                           monoclonal = des_covid_t$monoclonal, 
                           interferon = des_covid_t$interferon, 
                           ID = des_covid_t$ID)

treatment_df$summary <- NA
treatment_df$summary[treatment_df$steroids=="Y"] <- 'Steroids'
treatment_df$summary[treatment_df$monoclonal=="Y"] <- 'Monoclonal'
treatment_df$summary[treatment_df$interferon=="Y"] <- 'Interferon'
treatment_df$summary[treatment_df$interferon=="Y" &
                       treatment_df$steroids=="Y" &
                       treatment_df$monoclonal=="N"] <- 'Interferon + Steroids'
treatment_df$summary[treatment_df$interferon=="Y" &
                       treatment_df$monoclonal=="Y" &
                       treatment_df$steroids=="N"] <- 'Interferon + Monoclonal'
treatment_df$summary[treatment_df$monoclonal=="Y" &
                       treatment_df$steroids=="Y" &
                       treatment_df$interferon=="N"] <- 'Steroids + Monoclonal'

annotation <- data.frame(Treatment = treatment_df$summary[match(rownames(exp_treatment), treatment_df$ID)], 
                         Age = des.obj$Age[match(rownames(exp_treatment), des.obj$ID)],
                         Sex = des.obj$Sex[match(rownames(exp_treatment), des.obj$ID)], 
                         #  Final_Severity = des.obj$Severity.Final[match(rownames(exp_treatment), des.obj$ID)],
                         Severity = exp_treatment[,1])

rownames(annotation) <- rownames(exp_treatment)

ann_colors <- list(Treatment = c("Steroids" = "#004529",
                                 "Steroids + Monoclonal" = "#41ab5d",
                                 "Interferon + Steroids" = "#addd8e", 
                                 "Monoclonal" = "#f7fcb9"), 
                   Age = c(low="white", middle = 'grey', high="black"),
                   Sex = c('Male' = "#dadaeb", 
                           'Female' = '#54278f'),
                   Severity = c('Mild' = '#f1b6da',
                                'Moderate' = '#de77ae',
                                'Severe' = '#8e0152'))
newnames_1 <- lapply(
  colnames(exp_treatment)[2:ncol(exp_treatment)],
  function(x) bquote(italic(.(x))))

t_p <- pheatmap(t(log2(exp_treatment[,2:ncol(exp_treatment)]+1)), 
                annotation_col = annotation,
                legend =TRUE, 
                main = "", 
                annotation_colors = ann_colors, 
                cluster_cols = F, 
                cluster_rows = T,
                show_colnames = F, 
                annotation_legend_side = "bottom",
                labels_row = as.expression(newnames_1))

grid::grid.newpage()
grid::grid.draw(t_p$gtable)

# Boxplots of the 10 top genes highlighted in the heatmap 
CD177 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                 y = log2(CD177+1),
                 color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "CD177")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")

OLFM4 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                   y = log2(OLFM4+1),
                                   color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "OLFM4")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")

ABCA13 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                   y = log2(ABCA13+1),
                                   color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "ABCA13")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


MMP8 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                    y = log2(MMP8+1),
                                    color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "MMP8")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


METTL7B <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                  y = log2(METTL7B+1),
                                  color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "METTL7B")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


PCOLCE2 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                     y = log2(PCOLCE2+1),
                                     color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "PCOLCE2")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


COL17A1 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                     y = log2(COL17A1+1),
                                     color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "COL17A1")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


SERPINB10 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                     y = log2(SERPINB10+1),
                                     color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "SERPINB10")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


CD177P1 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                       y = log2(CD177P1+1),
                                       color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "CD177P1")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")


PCSK9 <- ggplot(exp_treatment, aes(x = as.factor(Severity),
                                     y = log2(PCSK9+1),
                                     color = Severity))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  theme_bw()+
  labs(x = "", y = "PCSK9")+
  scale_color_manual(values = c("#de77ae", "#c51b7d", '#8e0152'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'), 
        legend.position = "none")

grid.arrange(CD177, OLFM4, ABCA13, MMP8, METTL7B, PCOLCE2, COL17A1, SERPINB10, CD177P1, PCSK9, ncol = 2)
