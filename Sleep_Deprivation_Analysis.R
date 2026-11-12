############################################################
# PROJECT: Analysis of Gene Expression Changes During
#          Prolonged Sleep Deprivation
#
# DATASET: GSE98564
# ORGANISM: Homo sapiens
# PLATFORM: Affymetrix Human Gene 1.0 ST
#
# OBJECTIVE:
# To investigate gene expression changes in human blood
# during baseline, sleep deprivation, and recovery phases.
#
# ANALYSIS:
# 1. Data loading and preparation
# 2. Exploratory data analysis
# 3. Differential gene expression analysis
# 4. Gene annotation
# 5. Principal Component Analysis (PCA)
# 6. Hierarchical clustering and heatmap
# 7. Gene Set Enrichment Analysis (GSEA)
#
############################################################



############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Load required packages
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)




############################################################
# 2. LOAD GENE EXPRESSION DATASET
############################################################

# Download the GEO dataset
gse <- getGEO("GSE98564", GSEMatrix = TRUE)

# Extract the ExpressionSet
eset <- gse[[1]]

# Extract sample information
sample_info <- pData(eset)

# Extract gene expression values
expr <- exprs(eset)

# Check dataset information
cat("Dataset: GSE98564\n")
cat("Platform:", annotation(eset), "\n")
cat("Number of transcript clusters:", nrow(expr), "\n")
cat("Number of samples:", ncol(expr), "\n")

# Check expression matrix
dim(expr)
summary(as.vector(expr))



############################################################
# 3. SAMPLE METADATA PREPARATION
############################################################

# Extract important sample information
meta <- data.frame(
  GSM = sample_info$geo_accession,
  
  subject_id = sub(
    "^[^:]+:\\s*", "",
    sample_info$characteristics_ch1
  ),
  
  fatigue_category = sub(
    "^[^:]+:\\s*", "",
    sample_info$characteristics_ch1.1
  ),
  
  subject_group = sub(
    "^[^:]+:\\s*", "",
    sample_info$characteristics_ch1.2
  ),
  
  study_phase = sub(
    "^[^:]+:\\s*", "",
    sample_info$characteristics_ch1.3
  ),
  
  timepoint = sub(
    "^[^:]+:\\s*", "",
    sample_info$characteristics_ch1.4
  )
)

# Convert variables into factors
meta$subject_group <- factor(
  meta$subject_group,
  levels = c("Control", "Sleep Deprived")
)

meta$study_phase <- factor(
  meta$study_phase,
  levels = c("baseline", "experimental", "recovery")
)

meta$subject_id <- factor(meta$subject_id)

# Check metadata
head(meta)

# Count participants in each group
table(meta$subject_group)

# Count samples in each study phase
table(meta$subject_group, meta$study_phase)

# Check that expression samples match metadata
identical(colnames(expr), meta$GSM)




############################################################
# 4. EXPLORATORY DATA ANALYSIS (EDA)
############################################################

# Check dimensions of expression matrix
dim(expr)

# Check missing values
sum(is.na(expr))

# Check minimum and maximum expression values
range(expr, na.rm = TRUE)

# Get summary statistics
summary(as.vector(expr))

# Check number of unique participants
length(unique(meta$subject_id))

# Check sample distribution across groups and phases
table(meta$subject_group, meta$study_phase)

# Check number of samples per participant
table(meta$subject_id)



############################################################
# 5. DIFFERENTIAL GENE EXPRESSION ANALYSIS
############################################################

# Create the statistical design matrix
design <- model.matrix(
  ~ 0 + subject_group:study_phase,
  data = meta
)

# Assign meaningful column names
colnames(design) <- c(
  "Control_Baseline",
  "Sleep_Baseline",
  "Control_Experimental",
  "Sleep_Experimental",
  "Control_Recovery",
  "Sleep_Recovery"
)

# Check the design matrix
head(design)

# Check design matrix dimensions
dim(design)



############################################################
# 6. ACCOUNT FOR REPEATED MEASUREMENTS
############################################################

# Estimate within-participant correlation
corfit <- duplicateCorrelation(
  expr,
  design,
  block = meta$subject_id
)

# Display estimated consensus correlation
corfit$consensus


############################################################
# 7. FIT THE LINEAR MODEL
############################################################

# Fit the linear model while accounting for
# repeated measurements from the same participants

fit <- lmFit(
  expr,
  design,
  block = meta$subject_id,
  correlation = corfit$consensus
)

# Check the fitted model
fit



############################################################
# 8. DEFINE BIOLOGICAL COMPARISONS
############################################################

contrast_matrix <- makeContrasts(
  
  # Sleep-deprived group: experimental vs baseline
  SD_Experimental_vs_Baseline =
    Sleep_Experimental - Sleep_Baseline,
  
  # Sleep-deprived group: recovery vs experimental
  SD_Recovery_vs_Experimental =
    Sleep_Recovery - Sleep_Experimental,
  
  # Sleep-deprived group: recovery vs baseline
  SD_Recovery_vs_Baseline =
    Sleep_Recovery - Sleep_Baseline,
  
  # Control group: experimental vs baseline
  Control_Experimental_vs_Baseline =
    Control_Experimental - Control_Baseline,
  
  # Difference in changes between sleep-deprived and control groups
  Experimental_Difference =
    (Sleep_Experimental - Sleep_Baseline) -
    (Control_Experimental - Control_Baseline),
  
  # Difference in recovery changes between the two groups
  Recovery_Difference =
    (Sleep_Recovery - Sleep_Baseline) -
    (Control_Recovery - Control_Baseline),
  
  levels = design
)

# Display the contrast matrix
contrast_matrix



############################################################
# 9. APPLY CONTRASTS AND EMPIRICAL BAYES
############################################################

# Apply the contrast matrix to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Check the results
fit2

# Display available comparisons
colnames(fit2$coefficients)



############################################################
# 10. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
############################################################

# Extract results for all six comparisons
all_results <- lapply(
  colnames(contrast_matrix),
  function(comparison) {
    topTable(
      fit2,
      coef = comparison,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
  }
)

# Assign comparison names to the results
names(all_results) <- colnames(contrast_matrix)

# Save each comparison as a CSV file
for (comparison in names(all_results)) {
  
  write.csv(
    all_results[[comparison]],
    file = paste0(
      "data/DE_",
      comparison,
      ".csv"
    ),
    row.names = TRUE
  )
}

# Count statistically significant transcript clusters
significant_counts <- sapply(
  all_results,
  function(results) {
    sum(results$adj.P.Val < 0.05, na.rm = TRUE)
  }
)

# Create a summary table
de_summary <- data.frame(
  Comparison = names(significant_counts),
  Significant_Transcript_Clusters =
    as.integer(significant_counts)
)

# Display the summary
print(de_summary)

exists("fit2")


# Extract results for all six comparisons
all_results <- lapply(
  colnames(contrast_matrix),
  function(comparison) {
    topTable(
      fit2,
      coef = comparison,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
  }
)

# Assign comparison names
names(all_results) <- colnames(contrast_matrix)


names(all_results)


# Save each comparison as a CSV file
for (comparison in names(all_results)) {
  
  write.csv(
    all_results[[comparison]],
    file = paste0(
      "data/DE_",
      comparison,
      ".csv"
    ),
    row.names = TRUE
  )
}


# Count statistically significant transcript clusters
significant_counts <- sapply(
  all_results,
  function(results) {
    sum(results$adj.P.Val < 0.05, na.rm = TRUE)
  }
)

# Create summary table
de_summary <- data.frame(
  Comparison = names(significant_counts),
  Significant_Transcript_Clusters =
    as.integer(significant_counts)
)

# Display summary
print(de_summary)



############################################################
# 11. SAVE DIFFERENTIAL EXPRESSION SUMMARY
############################################################

write.csv(
  de_summary,
  "data/Differential_Expression_Summary.csv",
  row.names = FALSE
)

# Display saved summary
print(de_summary)



############################################################
# 12. GENE ANNOTATION
############################################################

# Extract significant results from recovery vs baseline
significant_recovery <- all_results[[
  "SD_Recovery_vs_Baseline"
]]

significant_recovery <- significant_recovery[
  significant_recovery$adj.P.Val < 0.05,
]

# Extract transcript cluster IDs
probe_ids <- rownames(significant_recovery)

# Map transcript cluster IDs to gene information
gene_annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

# Display gene annotation
print(gene_annotation)

# Save annotation results
write.csv(
  gene_annotation,
  "data/Significant_Gene_Annotation.csv",
  row.names = FALSE
)


############################################################
# 13. TAF4 EXPRESSION ANALYSIS
############################################################

# Extract TAF4 expression values
taf4_expression <- as.numeric(expr["8067361", ])

# Create a data frame for plotting
taf4_data <- data.frame(
  Expression = taf4_expression,
  Group = meta$subject_group,
  Phase = meta$study_phase
)

# Calculate average expression by group and phase
taf4_summary <- aggregate(
  Expression ~ Group + Phase,
  data = taf4_data,
  FUN = mean
)

# Display summary
print(taf4_summary)

# Create TAF4 expression boxplot
taf4_plot <- ggplot(
  taf4_data,
  aes(x = Phase, y = Expression, fill = Group)
) +
  geom_boxplot() +
  labs(
    title = "TAF4 Expression Across Study Phases",
    x = "Study Phase",
    y = "Expression Level",
    fill = "Group"
  ) +
  theme_minimal()

# Display plot
print(taf4_plot)

# Save plot
ggsave(
  "figures/TAF4_Expression_Across_Study_Phases.jpeg",
  plot = taf4_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Save summary
write.csv(
  taf4_summary,
  "data/TAF4_Expression_Summary.csv",
  row.names = FALSE
)



############################################################
# 14. PRINCIPAL COMPONENT ANALYSIS (PCA)
############################################################

# Calculate variance for each transcript cluster
gene_variance <- apply(expr, 1, var)

# Select the 500 most variable transcript clusters
top_variable <- names(
  sort(gene_variance, decreasing = TRUE)[1:500]
)

# Create expression matrix using selected clusters
pca_data <- t(expr[top_variable, ])

# Perform PCA
pca_result <- prcomp(
  pca_data,
  center = TRUE,
  scale. = TRUE
)

# Calculate percentage of variance explained
variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

# Display variance explained by first two components
pca_summary <- data.frame(
  Component = c("PC1", "PC2"),
  Variance_Explained = variance_explained[1:2],
  Cumulative_Variance = c(
    variance_explained[1],
    sum(variance_explained[1:2])
  )
)

print(pca_summary)

# Create PCA plotting data
pca_plot_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = meta$subject_group,
  Phase = meta$study_phase
)

# Create PCA plot
pca_plot <- ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Phase
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Gene Expression Profiles",
    x = paste0(
      "PC1 (", round(variance_explained[1], 2), "%)"
    ),
    y = paste0(
      "PC2 (", round(variance_explained[2], 2), "%)"
    )
  ) +
  theme_minimal()

# Display PCA plot
print(pca_plot)

# Save PCA plot
ggsave(
  "figures/PCA_of_Gene_Expression_Profiles.jpeg",
  plot = pca_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Save PCA summary
write.csv(
  pca_summary,
  "data/PCA_Variance_Summary.csv",
  row.names = FALSE
)


############################################################
# 15. HIERARCHICAL CLUSTERING AND HEATMAP
############################################################

# Install pheatmap if it is not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# Calculate variance for each transcript cluster
gene_variance <- apply(expr, 1, var)

# Select the 30 most variable transcript clusters
top_30 <- names(
  sort(gene_variance, decreasing = TRUE)[1:30]
)

# Extract expression values
heatmap_data <- expr[top_30, ]

# Prepare sample annotations
sample_annotation <- data.frame(
  Group = meta$subject_group,
  Phase = meta$study_phase
)

rownames(sample_annotation) <- colnames(heatmap_data)

# Generate heatmap
heatmap_plot <- pheatmap(
  heatmap_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = sample_annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Heatmap of Top 30 Variable Transcript Clusters"
)

# Save heatmap
png(
  "figures/Heatmap_Top_30_Variable_Clusters.png",
  width = 1600,
  height = 1200,
  res = 150
)

pheatmap(
  heatmap_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = sample_annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Heatmap of Top 30 Variable Transcript Clusters"
)

dev.off()


# ==========================================
# 21. PREPARE RANKED GENE LIST FOR GSEA
# ==========================================

# Extract results from the primary comparison
primary_results <- all_results[["Experimental_Difference"]]

# Add transcript cluster IDs
primary_results$PROBEID <- rownames(primary_results)

# Map transcript clusters to gene symbols and Entrez IDs
primary_annotations <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = primary_results$PROBEID,
  keytype = "PROBEID",
  columns = c("SYMBOL", "ENTREZID")
)

# Combine differential expression results and annotations
primary_annotated <- merge(
  primary_results,
  primary_annotations,
  by = "PROBEID"
)

# Remove entries without Entrez IDs
primary_annotated <- primary_annotated[
  !is.na(primary_annotated$ENTREZID),
]

# Keep the transcript with the largest absolute moderated t-statistic
primary_annotated <- primary_annotated[
  order(abs(primary_annotated$t), decreasing = TRUE),
]

primary_unique <- primary_annotated[
  !duplicated(primary_annotated$ENTREZID),
]

# Create the ranked gene list
gene_list <- primary_unique$t
names(gene_list) <- primary_unique$ENTREZID

# Sort from highest to lowest
gene_list <- sort(gene_list, decreasing = TRUE)

# Check the result
length(gene_list)
head(gene_list)

# ==========================================
# 22. KEGG GENE SET ENRICHMENT ANALYSIS
# ==========================================

# Set seed for reproducibility
set.seed(123)

# Add a very small random value to resolve tied ranking scores
gene_list_checked <- gene_list +
  runif(length(gene_list), -1e-10, 1e-10)

# Sort the adjusted gene list
gene_list_checked <- sort(
  gene_list_checked,
  decreasing = TRUE
)

# Perform KEGG GSEA
gsea_checked <- gseKEGG(
  geneList = gene_list_checked,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0
)

# Convert results into a data frame
gsea_checked_table <- as.data.frame(gsea_checked)

# Save all pathway results
write.csv(
  gsea_checked_table,
  "data/KEGG_GSEA_All_Results.csv",
  row.names = FALSE
)

# Extract statistically significant pathways
verified_significant <- gsea_checked_table[
  !is.na(gsea_checked_table$p.adjust) &
    gsea_checked_table$p.adjust < 0.05,
]

# Save significant pathways
write.csv(
  verified_significant,
  "data/KEGG_GSEA_Significant_Pathways.csv",
  row.names = FALSE
)

# Display summary
cat("Total pathways:", nrow(gsea_checked_table), "\n")
cat("Significant pathways:", nrow(verified_significant), "\n")

# Display the top 10 pathways
head(
  gsea_checked_table[
    order(gsea_checked_table$p.adjust),
    c("Description", "NES", "pvalue", "p.adjust")
  ],
  10
)


# ==========================================
# 23. KEGG GSEA DOT PLOT
# ==========================================

# Select the top 15 significant pathways
top_pathways <- verified_significant[
  order(verified_significant$p.adjust),
]

top_pathways <- head(top_pathways, 15)

# Create the dot plot
gsea_dotplot <- ggplot(
  top_pathways,
  aes(
    x = NES,
    y = reorder(Description, NES),
    size = setSize,
    color = p.adjust
  )
) +
  geom_point() +
  scale_color_gradient(
    low = "red",
    high = "blue"
  ) +
  labs(
    title = "KEGG Gene Set Enrichment Analysis",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "Gene Set Size",
    color = "Adjusted P-value"
  ) +
  theme_minimal()

# Display the plot
print(gsea_dotplot)

# Save the plot
ggsave(
  "figures/KEGG_GSEA_Dotplot.png",
  plot = gsea_dotplot,
  width = 11,
  height = 8,
  dpi = 300
)

# ==========================================
# 24. RIBOSOME PATHWAY ENRICHMENT PLOT
# ==========================================

# Generate the enrichment plot
ribosome_plot <- gseaplot2(
  gsea_checked,
  geneSetID = "hsa03010",
  title = "Ribosome Pathway Enrichment"
)

# Display the plot
print(ribosome_plot)

# Save the plot
ggsave(
  "figures/GSEA_Ribosome_Enrichment.png",
  plot = ribosome_plot,
  width = 10,
  height = 7,
  dpi = 300
)

list.files("figures")

browseURL(
  normalizePath("figures/KEGG_GSEA_Dotplot.png")
)

browseURL(
  normalizePath("figures/GSEA_Ribosome_Enrichment.png")
)


# ==========================================
# 25. EXTRACT RIBOSOME LEADING-EDGE GENES
# ==========================================

# Extract Ribosome pathway results
ribosome_row <- as.data.frame(gsea_checked@result)

ribosome_row <- ribosome_row[
  ribosome_row$ID == "hsa03010",
]

# Extract Entrez gene IDs
ribosome_entrez <- unlist(
  strsplit(ribosome_row$core_enrichment[1], "/")
)

# Count leading-edge genes
cat(
  "Number of Ribosome leading-edge genes:",
  length(ribosome_entrez),
  "\n"
)

# Convert Entrez IDs into gene symbols and names
ribosome_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ribosome_entrez,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

# Save annotations
write.csv(
  ribosome_annotation,
  "data/Ribosome_Leading_Edge_Genes.csv",
  row.names = FALSE
)

# Display the first 10 genes
head(ribosome_annotation, 10)


# ==========================================
# 26. RIBOSOME GENE-EXPRESSION STATISTICS
# ==========================================

# Extract statistics for Ribosome leading-edge genes
ribosome_gene_stats <- primary_unique[
  primary_unique$ENTREZID %in% ribosome_entrez,
  c("ENTREZID", "SYMBOL", "logFC", "t",
    "P.Value", "adj.P.Val")
]

# Remove duplicate Entrez IDs
ribosome_gene_stats <- ribosome_gene_stats[
  !duplicated(ribosome_gene_stats$ENTREZID),
]

# Sort by moderated t-statistic
ribosome_gene_stats <- ribosome_gene_stats[
  order(ribosome_gene_stats$t, decreasing = TRUE),
]

# Save results
write.csv(
  ribosome_gene_stats,
  "data/Ribosome_Gene_Expression_Statistics.csv",
  row.names = FALSE
)

# Display the first 10 genes
head(ribosome_gene_stats, 10)

# Count positive and negative statistics
cat(
  "Positive t-statistics:",
  sum(ribosome_gene_stats$t > 0),
  "\n"
)

cat(
  "Negative t-statistics:",
  sum(ribosome_gene_stats$t < 0),
  "\n"
)

# ==========================================
# 27. OXIDATIVE PHOSPHORYLATION ANALYSIS
# ==========================================

# Extract Oxidative Phosphorylation pathway
oxidative_row <- as.data.frame(gsea_checked@result)

oxidative_row <- oxidative_row[
  oxidative_row$ID == "hsa00190",
]

# Extract leading-edge Entrez IDs
oxidative_entrez <- unlist(
  strsplit(oxidative_row$core_enrichment[1], "/")
)

# Count leading-edge genes
cat(
  "Number of Oxidative Phosphorylation leading-edge genes:",
  length(oxidative_entrez),
  "\n"
)

# Annotate genes
oxidative_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = oxidative_entrez,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

# Save gene annotations
write.csv(
  oxidative_annotation,
  "data/Oxidative_Phosphorylation_Gene_Annotations.csv",
  row.names = FALSE
)

# Display the first 10 genes
head(oxidative_annotation, 10)


# ==========================================
# 28. OXIDATIVE PHOSPHORYLATION STATISTICS
# ==========================================

# Extract gene-expression statistics
oxidative_gene_stats <- primary_unique[
  primary_unique$ENTREZID %in% oxidative_entrez,
  c("ENTREZID", "SYMBOL", "logFC", "t",
    "P.Value", "adj.P.Val")
]

# Remove duplicate Entrez IDs
oxidative_gene_stats <- oxidative_gene_stats[
  !duplicated(oxidative_gene_stats$ENTREZID),
]

# Sort by moderated t-statistic
oxidative_gene_stats <- oxidative_gene_stats[
  order(oxidative_gene_stats$t, decreasing = TRUE),
]

# Save results
write.csv(
  oxidative_gene_stats,
  "data/Oxidative_Phosphorylation_Gene_Statistics.csv",
  row.names = FALSE
)

# Display the first 10 genes
head(oxidative_gene_stats, 10)

# Count positive and negative statistics
cat(
  "Positive t-statistics:",
  sum(oxidative_gene_stats$t > 0),
  "\n"
)

cat(
  "Negative t-statistics:",
  sum(oxidative_gene_stats$t < 0),
  "\n"
)


# ==========================================
# 29. FINAL RESULTS SUMMARY
# ==========================================

results_summary <- data.frame(
  Analysis = c(
    "Dataset",
    "Differential Expression",
    "TAF4",
    "KEGG GSEA",
    "Ribosome",
    "Oxidative Phosphorylation"
  ),
  
  Findings = c(
    "199 blood samples from 17 participants",
    "0 significant transcript clusters in primary comparison",
    "1 significant transcript cluster in SD recovery vs baseline",
    paste(
      nrow(verified_significant),
      "significant pathways out of",
      nrow(gsea_checked_table)
    ),
    paste(
      "NES =",
      round(gsea_checked_table$NES[
        gsea_checked_table$ID == "hsa03010"
      ], 3),
      "; adjusted p =",
      signif(gsea_checked_table$p.adjust[
        gsea_checked_table$ID == "hsa03010"
      ], 3)
    ),
    paste(
      "NES =",
      round(gsea_checked_table$NES[
        gsea_checked_table$ID == "hsa00190"
      ], 3),
      "; adjusted p =",
      signif(gsea_checked_table$p.adjust[
        gsea_checked_table$ID == "hsa00190"
      ], 3)
    )
  )
)

# Display summary
print(results_summary)

# Save summary
write.csv(
  results_summary,
  "data/Final_Results_Summary.csv",
  row.names = FALSE
)


# ==========================================
# 30. CHECK PROJECT OUTPUT FILES
# ==========================================

# List all generated data files
cat("DATA FILES:\n")
list.files("data")

# List all generated figures
cat("\nFIGURE FILES:\n")
list.files("figures")

# Check the final results summary
file.exists("data/Final_Results_Summary.csv")

# Check the main GSEA outputs
file.exists("data/KEGG_GSEA_All_Results.csv")
file.exists("data/KEGG_GSEA_Significant_Pathways.csv")

# Check the important gene statistics
file.exists("data/Ribosome_Gene_Expression_Statistics.csv")
file.exists("data/Oxidative_Phosphorylation_Gene_Statistics.csv")