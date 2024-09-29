library(Biobase)
library(hypeR)
library(readr)
library(dplyr)
library(ggplot2)
library(preprocessCore)
library(reshape2)
library(DESeq2)
library(pheatmap)
library(randomForest)
library(caret)
library(glmnet)
library(e1071)
library(pROC)
library(MLmetrics)
library(fgsea)
library(limma)
library(msigdbr)

###### Data Preprocessing



# load expression data
expdata <- readRDS("/projectnb/bs831/projects/data/HNSC_htseq_raw_counts_AEvsG1vsG3.RDS")
table(expdata$grade)
# gene expression data
dim(exprs(expdata))
# participant data, phenotype data
dim(pData(expdata))
# feature data;gene data
dim(fData(expdata))
#current distribution of the data
hist(exprs(expdata))  

# Assuming `expdata` is your ExpressionSet object
exprs_data <- exprs(expdata)
sample_info <- pData(expdata)
feature_info <- fData(expdata)

# Filter out 'ae' samples and keep only 'g1' and 'g3'
samples_to_keep <- sample_info$grade %in% c("g1", "g3")
filtered_exprs_data <- exprs_data[, samples_to_keep]
filtered_sample_info <- sample_info[samples_to_keep, ]

# Update feature data to match the filtered expression data
filtered_feature_info <- feature_info[rowSums(filtered_exprs_data) > 0, ]

# Create a new ExpressionSet object with the filtered data
expdata_filtered <- ExpressionSet(
  assayData = filtered_exprs_data[rowSums(filtered_exprs_data) > 0, ],
  phenoData = AnnotatedDataFrame(filtered_sample_info),
  featureData = AnnotatedDataFrame(filtered_feature_info)
)

table(expdata_filtered$grade)
# gene expression data
dim(exprs(expdata_filtered))
# participant data, phenotype data
dim(pData(expdata_filtered))
# feature data;gene data
dim(fData(expdata_filtered))
#current distribution of the data
hist(exprs(expdata_filtered)) 


# Function to create violin plots by group
violin_plot_exprs <- function(exprs_data, pheno_data, title) {
  melted_data <- melt(exprs_data)
  colnames(melted_data) <- c("Gene", "Sample", "Expression")
  melted_data <- merge(melted_data, pheno_data, by.x = "Sample", by.y = "row.names")
  
  ggplot(melted_data, aes(x = grade, y = Expression, fill = grade)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
    theme_minimal() +
    labs(title = title, x = "Grade", y = "Expression") +
    scale_fill_manual(values = c("g1" = "#1f77b4", "g3" = "#ff7f0e")) +
    theme(legend.position = "none")
}

# Boxplot of median expression per group
boxplot_median_exprs <- function(exprs_data, pheno_data, title) {
  melted_data <- melt(exprs_data)
  colnames(melted_data) <- c("Gene", "Sample", "Expression")
  melted_data <- merge(melted_data, pheno_data, by.x = "Sample", by.y = "row.names")
  
  median_exprs <- aggregate(Expression ~ grade + Gene, data = melted_data, median)
  
  ggplot(median_exprs, aes(x = grade, y = Expression, fill = grade)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(title = title, x = "Grade", y = "Expression") +
    scale_fill_manual(values = c("g1" = "#1f77b4", "g3" = "#ff7f0e")) +
    theme(legend.position = "none")
}

# Violin plot before normalization
print(violin_plot_exprs(exprs(expdata_filtered), pData(expdata_filtered), "Expression Data Before Normalization"))

# Boxplot of median expression before normalization
# print(boxplot_median_exprs(exprs(expdata_filtered), pData(expdata_filtered), "Median Expression Data Before Normalization"))

# Check for missing values
missing_count <- sum(is.na(exprs(expdata_filtered)))
cat("Number of missing values in the dataset:", missing_count, "\n")

# Log2 Transform the Data
log2_exprs_data <- log2(exprs(expdata_filtered) + 1)

# Perform Quantile Normalization
quantile_normalized_exprs <- normalize.quantiles(as.matrix(log2_exprs_data))
rownames(quantile_normalized_exprs) <- rownames(log2_exprs_data)
colnames(quantile_normalized_exprs) <- colnames(log2_exprs_data)

variance_filter <- function(exprs_data, threshold = 0.1) {
  variances <- apply(exprs_data, 1, var)
  high_variance_genes <- variances >= quantile(variances, probs = threshold)
  exprs_data[high_variance_genes, ]
}

# Filter out low variance genes
filtered_quantile_normalized_exprs <- variance_filter(quantile_normalized_exprs, threshold = 0.25)

# Violin plot after normalization
print(violin_plot_exprs(filtered_quantile_normalized_exprs, pData(expdata_filtered), "Expression Data After Quantile Normalization"))

# Boxplot of median expression after normalization
# print(boxplot_median_exprs(filtered_quantile_normalized_exprs, pData(expdata_filtered), "Median Expression Data After Quantile Normalization"))

hist(filtered_quantile_normalized_exprs)

# Filter the featureData to include only those features that are present in the new expression data
updated_feature_data <- featureData(expdata_filtered)[rownames(filtered_quantile_normalized_exprs), ]


# Check if the rownames match now
all(rownames(filtered_quantile_normalized_exprs) == rownames(updated_feature_data))
new_expdata <- ExpressionSet(
  assayData = filtered_quantile_normalized_exprs,
  phenoData = phenoData(expdata_filtered),
  featureData = updated_feature_data
)


# Validate the ExpressionSet
validObject(new_expdata)


# Check dimensions to ensure they match
dim(exprs(new_expdata))  
dim(pData(new_expdata))
dim(fData(new_expdata))

# Optionally check the first few rows of each component
head(exprs(new_expdata))
head(pData(new_expdata))
head(fData(new_expdata))




# Verify that the expression data has been updated
hist(exprs(new_expdata))











#### Differential Expression Analysis

complete_cases_index <- complete.cases(pData(expdata_filtered)[, c("patient.number_pack_years_smoked",
                                                              "patient.gender",
                                                              "patient.age_at_initial_pathologic_diagnosis",
                                                              "grade")])



# Subset new_expdata to remove samples with NA in any of the specified columns
new_expdata_filtered <- expdata_filtered[, complete_cases_index]

# Check the dimensions of the new dataset to confirm samples were removed
dim(new_expdata_filtered)
# Verify no NAs in the critical columns
sum(is.na(pData(new_expdata_filtered)$patient.number_pack_years_smoked))
sum(is.na(pData(new_expdata_filtered)$patient.gender))
sum(is.na(pData(new_expdata_filtered)$patient.age_at_initial_pathologic_diagnosis))
sum(is.na(pData(new_expdata_filtered)$grade))

# Perform Differential Expression Analysis using DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(exprs(new_expdata_filtered)),
                              colData = pData(new_expdata_filtered),
                              design = ~ patient.number_pack_years_smoked + patient.age_at_initial_pathologic_diagnosis + grade)


# Relevel the grade factor to have 'g1' as the reference level
dds$grade <- relevel(dds$grade, ref = "g1")

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
res

# View summary of results
summary(res)


# List significant genes
significant_genes <- res[which(res$padj < 0.05), ]
significant_genes <- significant_genes[order(significant_genes$padj), ]
length(significant_genes)


# Print top significant genes
print(head(significant_genes))

# Plot MA plot
plotMA(res, main = "MA Plot of Differential Expression", ylim = c(-5, 5))

# Create a list of subsets for each level of the 'grade' variable
grade_levels <- levels(dds$grade)
grade_subsets <- lapply(grade_levels, function(x) {
  subset(dds, grade == x)
})

# Rerun the DESeq2 analysis for each subset and store the results in a list
grade_results <- lapply(grade_subsets, function(x) {
  x <- DESeq(x)
  results(x)
})

# Create an MA-plot for each subset using the results in the 'grade_results' list
lapply(seq_along(grade_results), function(i) {
  res <- grade_results[[i]]
grade_level <- grade_levels[i]
plotMA(res, main = paste0("MA-plot for the 'grade == ", grade_level, "' subset"))
})


volcano_data$significant <- ifelse(volcano_data$pvalue < 0.05, "Significant", "Not Significant")
volcano_data <- data.frame(log2FoldChange = res$log2FoldChange, 
                           negLog10pval = res$pvalue,
                           pvalue = res$pvalue)
volcano_data$significant <- ifelse(volcano_data$pvalue < 0.05, "Significant", "Not Significant")
# Convert pvalue to -log10(pvalue)
volcano_data$negLog10pval <- -log10(volcano_data$pvalue)

# Remove the original pvalue column
volcano_data <- volcano_data[, -which(names(volcano_data) == "pvalue")]

volcano_data$significant <- ifelse(volcano_data$pvalue < 0.05, "Significant", "Not Significant")

ggplot(volcano_data, aes(x = log2FoldChange, y = negLog10pval, color = significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("Significant" = "olivedrab4", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot of Differential Expression", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme(legend.position = "right")





### Hierarchical Clustering of only top 50 genes

# Order genes by significance
resOrdered <- res[order(res$padj), ]

# Select the top 50 most significant genes
topGenes <- head(rownames(resOrdered), 50)

# Extract normalized counts for these top genes
normalized_counts <- counts(dds, normalized = TRUE)
# Subset to the top genes
sig_gene_data <- normalized_counts[topGenes, ]

# Extract sample information and filter for 'g1' and 'g3'
sample_info <- colData(dds)
filtered_sample_info <- sample_info[sample_info$grade %in% c('g1', 'g3'), ]
filtered_sample_ids <- rownames(filtered_sample_info)

# Filter columns of sig_gene_data based on these sample IDs
sig_gene_data <- sig_gene_data[, filtered_sample_ids]

# Annotations should match the filtered sample IDs
anno_col <- data.frame(Grade = filtered_sample_info$grade)
rownames(anno_col) <- filtered_sample_ids

# Create a color mapping for the annotations
ann_colors <- list(Grade = c(g1 = "skyblue", g3 = "salmon"))

# Color palette from blue to red
my_color <- colorRampPalette(c("navy", "white", "firebrick4"))(100)

# Prepare the data for the heatmap
data_for_heatmap <- as.matrix(sig_gene_data)

# Create the heatmap with annotations
pheatmap(data_for_heatmap,
         color = my_color,
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",
         fontsize_row = 10,
         fontsize = 12,
         main = "Top 50 Significant Genes Heatmap")


#### Clustered by Grade 

# Ensure the samples are ordered by grade
ordered_sample_ids <- filtered_sample_ids[order(filtered_sample_info$grade)]

# Reorder the data matrix according to the grade
data_for_heatmap <- sig_gene_data[, ordered_sample_ids]

# Annotations should match the sample IDs
anno_col <- data.frame(Grade = filtered_sample_info$grade[order(filtered_sample_info$grade)])
rownames(anno_col) <- ordered_sample_ids

# Create a color mapping for the annotations
ann_colors <- list(Grade = c(g1 = "skyblue", g3 = "salmon"))

# Color palette from blue to red
my_color <- colorRampPalette(c("navy", "white", "firebrick4"))(100)

# Create the heatmap with annotations and ordered samples
pheatmap(data_for_heatmap,
         color = my_color,
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "row",
         fontsize_row = 10,
         fontsize = 12,
         main = "Top 50 Significant Genes Heatmap")


#### Pearson Correlation Version

de_genes <- rownames(subset(res, padj < 0.05))

# Select the top 50 significant genes
topGenes <- head(de_genes, 50)

# Subset count data to top 50 significant genes
counts_data_top50 <- exprs(new_expdata_filtered)[topGenes, ]

vst_data <- varianceStabilizingTransformation(dds, blind=FALSE)

# Subset to top 50 significant genes
vst_top50 <- assay(vst_data)[topGenes, ]

# Calculate Pearson correlation coefficient
correlation_matrix <- cor(t(vst_top50), method = "pearson")

# Perform hierarchical clustering based on Pearson correlation coefficient
hclust_res <- hclust(as.dist(1 - correlation_matrix), method = "complete")

# Convert hclust_res to dendrogram
dend <- as.dendrogram(hclust_res)

# Plot dendrogram
plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")

# Reorder counts data based on dendrogram
heatmap_data <- vst_top50[order.dendrogram(dend), ]

pheatmap(heatmap_data,
         color = my_color,
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         fontsize = 12,
         main = "Top 50 Significant Genes Heatmap")





##### Geneset Enrichment Analysis



# Load all "HALLMARK" gene sets for "Homo sapiens"
msigdbr_species <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_genesets <- split(msigdbr_species$gene_symbol, msigdbr_species$gs_name)
# Filter to keep only 'HALLMARK' gene sets
hallmark_genesets <- hallmark_genesets[grep("^HALLMARK_", names(hallmark_genesets))]
# Using log fold change for ranking
group_labels <- pData(expdata_filtered)$grade
group_labels <- factor(group_labels, levels = c("g1", "g3"))

# Make design matrix
design_matrix <- model.matrix(~ 0 + group_labels)
colnames(design_matrix) <- levels(group_labels)  # Naming as per your groups

fit <- lmFit(filtered_quantile_normalized_exprs, design = design_matrix)
# Fit to model
contrast_matrix <- makeContrasts(g3vsG1 = g3 - g1, levels = design_matrix)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
summary(fit2)
ranked_genes <- sort(logFC, decreasing = TRUE)

# Extract HGNC symbols
hgnc_symbols <- fData(expdata_filtered)$hgnc_symbol
# Create a named vector where names are Ensembl IDs and values are HGNC symbols
ensembl_to_hgnc <- setNames(hgnc_symbols, rownames(fData(expdata_filtered)))
hgnc_ranked_genes <- ranked_genes[names(ranked_genes) %in% names(ensembl_to_hgnc)]
names(hgnc_ranked_genes) <- ensembl_to_hgnc[names(hgnc_ranked_genes)]

# Ensure no NAs in the names
hgnc_ranked_genes <- hgnc_ranked_genes[!is.na(names(hgnc_ranked_genes))]

# Run FGSEA
fgsea_results <- fgsea(pathways = hallmark_genesets, stats = hgnc_ranked_genes, minSize = 15, maxSize = 500)

# Filter for significant pathways
significant_pathways <- fgsea_results[fgsea_results$pval < 0.05, ]

# Plot the results
ggplot(significant_pathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0), show.legend = TRUE) +
  scale_fill_manual(name = "NES Direction",
                    values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Positive", "Negative")) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", title = "Top Enriched Pathways in G1 vs G3") +
  theme_minimal() +
  theme(legend.position = "right")







##### Classifiers



# Extract expression data and grade labels
exprs_data <- exprs(new_expdata)
sample_info <- pData(new_expdata)

# Convert grade labels to a factor
sample_info$grade <- as.factor(sample_info$grade)
sample_info$grade <- factor(sample_info$grade, levels = levels(sample_info$grade)[levels(sample_info$grade) != "AE"])

# Train-test split using caret
set.seed(42)
trainIndex <- createDataPartition(sample_info$grade, p = 0.7, list = FALSE)
train_data <- exprs_data[, trainIndex]
test_data <- exprs_data[, -trainIndex]
train_labels <- sample_info$grade[trainIndex]
test_labels <- sample_info$grade[-trainIndex]

# Remove constant features
non_constant_features <- apply(train_data, 1, function(x) length(unique(x)) > 1)
train_data <- train_data[non_constant_features, ]
test_data <- test_data[non_constant_features, ]

# Train RandomForest model
rf_model <- randomForest(x = t(train_data), y = train_labels, importance = TRUE)
rf_predictions <- predict(rf_model, t(test_data))

# Train Elastic Net model
set.seed(42)
cv_enet <- cv.glmnet(t(train_data), train_labels, family = "binomial", alpha = 0.5)
enet_model <- glmnet(t(train_data), train_labels, family = "binomial", alpha = 0.5, lambda = cv_enet$lambda.min)
enet_predictions <- predict(enet_model, t(test_data), type = "class")

# Train SVM model (linear)
set.seed(42)
svm_linear_model <- svm(t(train_data), train_labels, kernel = "linear", probability = TRUE, scale = FALSE)
svm_predictions <- predict(svm_linear_model, t(test_data))

# Evaluate model performance
rf_confusion_mtx <- confusionMatrix(rf_predictions, test_labels)
enet_confusion_mtx <- confusionMatrix(as.factor(enet_predictions), test_labels)
svm_confusion_mtx <- confusionMatrix(svm_predictions, test_labels)

# Print confusion matrices
print("Random Forest Confusion Matrix:")
print(rf_confusion_mtx)

print("Elastic Net Confusion Matrix:")
print(enet_confusion_mtx)

print("SVM Confusion Matrix:")
print(svm_confusion_mtx)

# Calculate AUC for each model
rf_probs <- predict(rf_model, t(test_data), type = "prob")
enet_probs <- predict(enet_model, t(test_data), type = "response")
svm_probs <- attr(predict(svm_linear_model, t(test_data), probability = TRUE), "probabilities")

# ROC for each model
roc_rf <- roc(test_labels, rf_probs[, "g3"], levels = c("g1", "g3"), direction = ">")
roc_enet <- roc(test_labels, as.numeric(enet_probs), levels = c("g1", "g3"), direction = ">")
roc_svm <- roc(test_labels, svm_probs[, "g3"], levels = c("g1", "g3"), direction = ">")

# Plot ROC
plot(roc_rf, col = "steelblue", main = "ROC Curves for Classification Models", lwd = 2)
lines(roc_enet, col = "darkred", lwd = 2)
lines(roc_svm, col = "darkgreen", lwd = 2)
legend("topleft", legend = c("Random Forest", "Elastic Net", "SVM"), col = c("steelblue", "darkred", "darkgreen"), lwd = 2)

# Calculate F1 scores
f1_rf <- F1_Score(test_labels, rf_predictions, positive = "g3")
f1_enet <- F1_Score(test_labels, as.factor(enet_predictions), positive = "g3")
f1_svm <- F1_Score(test_labels, svm_predictions, positive = "g3")

# Summary of performance metrics
metric_summary <- data.frame(
  Model = c("Random Forest", "Elastic Net", "SVM"),
  Accuracy = c(rf_confusion_mtx$overall["Accuracy"], enet_confusion_mtx$overall["Accuracy"], svm_confusion_mtx$overall["Accuracy"]),
  Precision = c(rf_confusion_mtx$byClass["Pos Pred Value"], enet_confusion_mtx$byClass["Pos Pred Value"], svm_confusion_mtx$byClass["Pos Pred Value"]),
  Recall = c(rf_confusion_mtx$byClass["Sensitivity"], enet_confusion_mtx$byClass["Sensitivity"], svm_confusion_mtx$byClass["Sensitivity"]),
  F1_Score = c(f1_rf, f1_enet, f1_svm),
  AUC = c(auc(roc_rf), auc(roc_enet), auc(roc_svm))
)

print(metric_summary)

# Feature importance for Random Forest
rf_importance_df <- data.frame(Feature = rownames(importance(rf_model)), Importance = importance(rf_model)[, 1])
rf_importance_df <- rf_importance_df[order(rf_importance_df$Importance, decreasing = TRUE), ]
rf_top_features <- head(rf_importance_df, 20)

ggplot(rf_top_features, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Important Features - Random Forest", x = "Feature", y = "Importance") +
  theme_minimal()

# Feature importance for Elastic Net
enet_coeffs <- as.data.frame(as.matrix(coef(enet_model)))
colnames(enet_coeffs) <- "Coefficient"
enet_coeffs$Feature <- rownames(enet_coeffs)
enet_coeffs <- enet_coeffs[order(abs(enet_coeffs$Coefficient), decreasing = TRUE), ]
enet_top_features <- head(enet_coeffs, 20)

ggplot(enet_top_features, aes(x = reorder(Feature, abs(Coefficient)), y = abs(Coefficient))) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  labs(title = "Top 20 Important Features - Elastic Net", x = "Feature", y = "Coefficient (Absolute Value)") +
  theme_minimal()

# Feature importance for SVM
svm_coeffs <- abs(t(svm_linear_model$coefs) %*% svm_linear_model$SV)
svm_importance_df <- data.frame(Feature = rownames(train_data), Importance = svm_coeffs[1, ])
svm_importance_df <- svm_importance_df[order(svm_importance_df$Importance, decreasing = TRUE), ]
svm_top_features <- head(svm_importance_df, 20)

ggplot(svm_top_features, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top 20 Important Features - SVM", x = "Feature", y = "Importance") +
  theme_minimal()



# PCA Plot with RandomForest 
pca_rf <- prcomp(t(test_data))
pca_rf_data <- data.frame(Sample = colnames(test_data),
                          PC1 = pca_rf$x[, 1],
                          PC2 = pca_rf$x[, 2],
                          Grade = test_labels,
                          Predicted = rf_predictions)

ggplot(pca_rf_data, aes(x = PC1, y = PC2, color = Predicted, shape = Grade)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Test Data with Random Forest Predicted Labels", x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("g1" = "#1f77b4", "g3" = "#ff7f0e"))

# PCA Plot with Elastic Net 
pca_enet <- prcomp(t(test_data))
pca_enet_data <- data.frame(Sample = colnames(test_data),
                            PC1 = pca_enet$x[, 1],
                            PC2 = pca_enet$x[, 2],
                            Grade = test_labels,
                            Predicted = enet_predictions)

ggplot(pca_enet_data, aes(x = PC1, y = PC2, color = s0, shape = Grade)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Test Data with Elastic Net Predicted Labels", x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("g1" = "#1f77b4", "g3" = "#ff7f0e"))

# PCA Plot with SVM 
pca_svm <- prcomp(t(test_data))
pca_svm_data <- data.frame(Sample = colnames(test_data),
                           PC1 = pca_svm$x[, 1],
                           PC2 = pca_svm$x[, 2],
                           Grade = test_labels,
                           Predicted = svm_predictions)

ggplot(pca_svm_data, aes(x = PC1, y = PC2, color = Predicted, shape = Grade)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Test Data with SVM Predicted Labels", x = "PC1", y = "PC2") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("g1" = "#1f77b4", "g3" = "#ff7f0e"))




