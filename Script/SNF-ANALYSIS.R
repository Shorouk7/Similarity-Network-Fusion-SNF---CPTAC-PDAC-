# *** Similarity Network Fusion (SNF) *** > PAN-CAN-PROJ
#============================================
# 
#                  BY *** SHOROUK ALDEYARBI ***
#
# Description:
# *** Integration of 8 omics layers
# ***  Proteomics
# *** Phosphoproteomics (Gene / Site / Multi-site)
# *** RNA-seq
# *** DNA Methylation
# *** Glycoproteomics (Gene / Site)
# *** Use of MOFA latent factor (Factor1) as tumor program signature
# *** Network-based clustering using SNF
# *** Identification of 3 tumor molecular subtypes
#============================================


#*** REQUIRED PACKAGERS==================================================
# install.packages("devtools")
# devtools::install_github("branchlab/metasnf")
# install.packages("SNFtool")# install.packages("devtools")
# devtools::install_github("branchlab/metasnf")
# install.packages("SNFtool")
#*** LOADIG LIBEREARIES====================================================

library(MOFA2)
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(randomForest)
library(e1071)
library(caret)
library(pROC)
library(reshape2)

#*** SNF ***
library(metasnf)
library(SNFtool)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
#*** LOADING FILES=========================================================

# *** MOFAOBJECT+ MODEL *** 
MOFAobject <- readRDS("Data/MOFA-OBJECT.rds")
factor_mat <- read.csv("Data/MOFA-MODEL-FACTORS.csv")
model <- load_model("Data/MOFA-MODEL.hdf5")
factors_df <- get_factors(MOFAobject, factors = 1, groups = "all", as.data.frame = TRUE)

# ***  META-DATA+ FACTORES ***
# merge factor values with stage info
meta <- samples_metadata(MOFAobject)
factors_df <- merge(factors_df, meta[, c("sample", "state","tumor_size", "tumor_metastasis", "lymph_Metastasis", "tumor_stage_pathological")], by = "sample")

head(factors_df)
# *** TOP FEATURE LOADING *** 

top_files <- list(
  Proteomic = "TOP-FEATURES/Top300_Positive_Factor1_Proteomic.csv",
  PhosphoGene = "TOP-FEATURES/Top300_Positive_Factor1_PhosphoGene.csv",
  PhosphoSite = "TOP-FEATURES/Top300_Positive_Factor1_PhosphoSite.csv",
  PhosphoMultiSite = "TOP-FEATURES/Top300_Positive_Factor1_PhosphoMultiSite.csv",
  RNAseq = "TOP-FEATURES/Top300_Positive_Factor1_RNAseq.csv",
  Methylation = "TOP-FEATURES/Top300_Positive_Factor1_Methylation.csv",
  GlycoGene = "TOP-FEATURES/Top300_Positive_Factor1_GlycoGene.csv",
  GlycoSite = "TOP-FEATURES/Top300_Positive_Factor1_GlycoSite.csv"
)

feature_lists <- lapply(top_files, function(f){
  df <- read_csv(f)
  unique(df$feature)})

###

data_list <- get_data(MOFAobject)

# *** Subset matrices using your top 300 features

omics_factor1 <- list()

for(v in names(feature_lists)){
  
  # extract view and group matrix
  mat <- data_list[[v]][["group1"]]
  
  # ensure it is matrix
  mat <- as.matrix(mat)
  
  # intersect with top features
  features <- intersect(rownames(mat), feature_lists[[v]])
  
  # subset matrix
  omics_factor1[[v]] <- mat[features, , drop = FALSE]
  
}

# STRUCTURE? 
str(data_list, max.level = 2)
# ** transpose matrices

omics_factor1 <- lapply(omics_factor1, function(x){
  
  t(x)
  
})


# ** scale data

omics_factor1 <- lapply(omics_factor1, scale)

# ** compute distances

dist_list <- lapply(omics_factor1, dist)

# ** Similarity Network: 

affinity_list <- lapply(dist_list, function(d){
  
  affinityMatrix(
    as.matrix(d),
    K = 20,
    sigma = 0.5
  )
  
})

# ** Fuse Network 

fused_network <- SNF(
  affinity_list,
  K = 20,
  t = 20
)

# ** Cluster Numbering: 
estimateNumberOfClustersGivenGraph(
  fused_network,
  NUMC = 2:6
)

# ** Clustring: 
sapply(omics_factor1, nrow)
sapply(omics_factor1, ncol)

clusters <- spectralClustering(
  fused_network,
  K = 3
)


# ** cluster with metadtata:

meta <- samples_metadata(MOFAobject)

meta$SNF_cluster <- clusters

# ** visualization:

ggplot(meta, aes(x = factor(SNF_cluster), fill = state)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    x = "SNF Cluster",
    y = "Proportion",
    title = "Cluster composition by Tumor/Normal state"
  )

# Does the Factor1 multi-omics program change with tumor progression?


meta_clean <- meta[!is.na(meta$tumor_stage_pathological), ]
ggplot(meta_clean,
       aes(x = factor(SNF_cluster),
           fill = tumor_stage_pathological)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    x = "SNF Cluster",
    y = "Proportion",
    title = "Tumor Stage Distribution across SNF clusters"
  )



# ** correlate Factor1 score with stage
factors_df <- get_factors(MOFAobject, factors = 1, as.data.frame = TRUE)

factors_df <- merge(
  factors_df,
  meta[, c("sample","tumor_stage_pathological")],
  by="sample"
)

factors_df <- factors_df[!is.na(factors_df$tumor_stage_pathological), ]


ggplot(factors_df,
       aes(x = tumor_stage_pathological,
           y = value,
           fill = tumor_stage_pathological)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x="Tumor Stage",
    y="Factor1 score",
    title="Association between Factor1 and Tumor Stage"
  )

# * Compute layer activity per sample

layer_scores <- data.frame(sample = rownames(omics_factor1[[1]]))

for(layer in names(omics_factor1)){
  
  mat <- omics_factor1[[layer]]
  
  layer_scores[[layer]] <- rowMeans(mat)
  
}


layer_scores$SNF_cluster <- clusters

# layer activity per cluster
cluster_layer_activity <- layer_scores %>%
  group_by(SNF_cluster) %>%
  summarise(across(-sample, mean))

# *** plottiing

cluster_layer_long <- cluster_layer_activity %>%
  pivot_longer(
    cols = -SNF_cluster,
    names_to = "Layer",
    values_to = "Activity"
  )

ggplot(cluster_layer_long,
       aes(x = Layer,
           y = Activity,
           fill = factor(SNF_cluster))) +
  geom_bar(stat="identity",
           position="dodge") +
  theme_minimal() +
  labs(
    title="Omics layer activity across SNF clusters",
    x="Omics layer",
    y="Mean activity",
    fill="Cluster"
  ) +
  theme(axis.text.x = element_text(angle=45,hjust=1))



## *** VISUILIZE CLUSTERS VIA 'Heatmap' : 

#library(pheatmap)
cluster_matrix <- as.matrix(cluster_layer_activity[,-1])
rownames(cluster_matrix) <- paste("Cluster",cluster_layer_activity$SNF_cluster)

pheatmap(cluster_matrix,
         scale="row",
         main="Omics layer activity across SNF clusters")


## Statistics: 

chisq.test(
  table(meta_clean$SNF_cluster,
        meta_clean$tumor_stage_pathological)
)


# ** calculate mean tumor program score per sample
# Assuming omics_factor1 is a list: each element is samples × features matrix with tumor features only

# ** Combine all tumor-positive features across layers for each sample (column bind features)
combined_features <- do.call(cbind, omics_factor1)

# ** Calculate average tumor score per sample
sample_scores <- rowMeans(combined_features, na.rm = TRUE)

# ** Add cluster info and metadata
meta$tumor_program_score <- sample_scores
meta$SNF_cluster <- factor(meta$SNF_cluster)  # make sure factor

# ** Visualize tumor score by cluster
library(ggplot2)
ggplot(meta, aes(x = SNF_cluster, y = tumor_program_score, fill = SNF_cluster)) +
  geom_boxplot() +
  labs(title = "Tumor Program Score across Clusters",
       x = "SNF Cluster",
       y = "Mean Tumor Positive Feature Expression") +
  theme_minimal()


#Which genes / proteins / phosphosites are enriched in each cluster?


#1* feature matrix: 

combined_features <- do.call(cbind, omics_factor1)

combined_df <- as.data.frame(combined_features)
combined_df$SNF_cluster <- clusters  # make sure clusters is a vector of cluster labels, ordered to match rows
#combined_df$cluster <- clusters

#2* features associated with each cluster:

length(clusters)
nrow(combined_df)
combined_df$SNF_cluster <- clusters
head(combined_df$SNF_cluster)
# check column exists
tail(colnames(combined_df), 10)
head(colnames(combined_df), 10)
"SNF_cluster" %in% colnames(combined_df)
# *** 


cluster_feature_means <- combined_df %>%
  group_by(SNF_cluster) %>%
  summarise(across(all_of(setdiff(colnames(combined_df), "SNF_cluster")), mean))

# ** Score cluster “tumor program strength”


# Calculate mean tumor-positive program per cluster
cluster_scores <- cluster_feature_means %>%
  rowwise() %>%
  mutate(tumor_program_score = mean(c_across(-SNF_cluster))) %>%
  arrange(desc(tumor_program_score))

cluster_scores

# ** categorization 

cluster_scores <- cluster_scores %>%
  mutate(tumor_program_strength = case_when(
    tumor_program_score > quantile(tumor_program_score, 0.66) ~ "Strong",
    tumor_program_score > quantile(tumor_program_score, 0.33) ~ "Intermediate",
    TRUE ~ "Weak"
  ))

# ** Identify which features drive each cluster

top_features_per_cluster <- cluster_feature_means %>%
  pivot_longer(-SNF_cluster, names_to = "feature", values_to = "mean_expression") %>%
  group_by(SNF_cluster) %>%
  slice_max(mean_expression, n = 20)  # top 20 features
top_features_per_cluster
# ** top features plotting=========================================== 

top_features_per_cluster %>%
  mutate(feature = reorder(feature, mean_expression)) %>%
  ggplot(aes(feature, mean_expression, fill=factor(SNF_cluster))) +
  geom_col() +
  coord_flip() +
  facet_wrap(~SNF_cluster, scales="free") +
  theme_minimal() +
  labs(
    title="Top molecular features driving each SNF cluster",
    x="Feature",
    y="Mean activity",
    fill="Cluster"
  )

# heatmap CLUSTERS VS STAGES GROUPING 

cluster_matrix <- cluster_feature_means %>%
  column_to_rownames("SNF_cluster") %>%
  as.matrix()

pheatmap(
  cluster_matrix,
  scale="row",
  show_colnames=FALSE,
  main="Multi-omics tumor program across SNF clusters"
)

ggplot(meta, aes(SNF_cluster, fill=tumor_stage_pathological)) +
  geom_bar(position="fill") +
  theme_minimal()


# ** Pathway / biological interpretation======================================

# ** PATHWAY ENRCHMENT PER CLUSTER 

# library(clusterProfiler)
# library(org.Hs.eg.db)
# gene_list <- gsub("_.*$", "", top_features_per_cluster$feature)
# gene_list <- unlist(strsplit(gene_list, ";"))
# gene_list <- unique(gene_list)

cluster_genes <- top_features_per_cluster %>%
  mutate(gene = gsub("_.*","",feature)) %>%
  separate_rows(gene, sep=";")

# ** gene list in the top cluster ** 
cluster3_genes <- cluster_genes %>%
  filter(SNF_cluster == 3) %>%
  pull(gene) %>%
  unique()
cluster3_genes
length(cluster3_genes)


# ** Enrichment
ego_c3 <- enrichGO(
  gene = cluster3_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

as.data.frame(ego_c3)
dotplot(ego_c3,
        showCategory = 15) +
  ggtitle("GO Biological Process enrichment of SNF tumor program")

barplot(ego_c3,
        showCategory = 15,
        title = "Top enriched GO Biological Processes")


# KEGG PATHWAYS
## Entrez IDs convertion

gene_df <- bitr(
  cluster3_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_df)
# ** KEGG pathway enrichment
kegg_c3 <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

as.data.frame(kegg_c3)

# ** PLOTTING: 

dotplot(kegg_c3, showCategory = 15) +
  ggtitle("KEGG pathway enrichment of Cluster 3 tumor program")

kegg_c3 <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 1
)

as.data.frame(kegg_c3)















