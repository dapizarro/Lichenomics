#!/usr/bin/env Rscript

# Lichenomics metabolomics statistics
# Input:
#   docs/metabolomics/example_feature_table.csv
#   docs/metabolomics/example_metadata.tsv
#
# Output:
#   results/metabolomics/

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

dir.create("results/metabolomics", recursive = TRUE, showWarnings = FALSE)

feature_table <- read.csv("docs/metabolomics/example_feature_table.csv",
                          check.names = FALSE)

metadata <- read.delim("docs/metabolomics/example_metadata.tsv")

feature_matrix <- feature_table %>%
  column_to_rownames("feature_id") %>%
  t() %>%
  as.data.frame()

feature_matrix$sample_id <- rownames(feature_matrix)

data <- metadata %>%
  inner_join(feature_matrix, by = "sample_id")

numeric_features <- data %>%
  select(where(is.numeric))

# PCA
pca <- prcomp(numeric_features, scale. = TRUE)

pca_scores <- as.data.frame(pca$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata, by = "sample_id")

write.csv(pca_scores,
          "results/metabolomics/pca_scores.csv",
          row.names = FALSE)

pdf("results/metabolomics/PCA_metabolomics.pdf")
plot(pca_scores$PC1, pca_scores$PC2,
     xlab = "PC1",
     ylab = "PC2",
     main = "Metabolomics PCA")
text(pca_scores$PC1, pca_scores$PC2,
     labels = pca_scores$sample_id,
     pos = 3)
dev.off()

# Bray-Curtis distance
bray <- vegdist(numeric_features, method = "bray")

# PERMANOVA
if ("chemotype" %in% colnames(metadata)) {
  permanova <- adonis2(bray ~ chemotype, data = metadata)
  write.csv(as.data.frame(permanova),
            "results/metabolomics/PERMANOVA_chemotype.csv")
}

# Feature-wise Wilcoxon test
if ("chemotype" %in% colnames(data)) {
  features <- colnames(numeric_features)

  wilcox_results <- map_dfr(features, function(f) {
    test <- wilcox.test(data[[f]] ~ data$chemotype)

    tibble(
      feature_id = f,
      p_value = test$p.value
    )
  }) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR)

  write.csv(wilcox_results,
            "results/metabolomics/differential_features_wilcoxon.csv",
            row.names = FALSE)
}

message("Metabolomics statistics completed.")
