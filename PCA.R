# PCA on fits_all (no Treatment column needed)
library(tidyverse)
library(cluster)
library(ggplot2)

#  scale data
pca_dat <- fits_all %>%
  select(Test, rep, Tcrit, T50) %>%
  filter(!is.na(Tcrit), !is.na(T50))

pca_res <- prcomp(scale(pca_dat[, c("Tcrit", "T50")]), center = TRUE, scale. = FALSE)
var_expl <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2], 1)
cat("PC1 =", var_expl[1], "%, PC2 =", var_expl[2], "%\n")

# PCA scores with metadata
scores <- as_tibble(pca_res$x) %>%
  mutate(rep = pca_dat$rep, Test = pca_dat$Test)

# k-means with K = 2
set.seed(42)
km <- kmeans(scores[, c("PC1", "PC2")], centers = 2, nstart = 50)
scores$cluster <- factor(km$cluster)

# Silhouette
sil <- silhouette(km$cluster, dist(scores[, c("PC1", "PC2")]))
cat("Mean silhouette width =", round(mean(sil[, "sil_width"]), 3), "\n")

# Cluster summaries
cluster_summary <- pca_dat %>%
  mutate(cluster = scores$cluster) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_Tcrit = mean(Tcrit), sd_Tcrit = sd(Tcrit),
    mean_T50 = mean(T50), sd_T50 = sd(T50),
    .groups = "drop"
  )
print(cluster_summary)

scores <- scores %>% mutate(cluster = factor(cluster))

p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_text_repel(aes(label = Test),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.size = 0.25,
                  show.legend = FALSE) +
  stat_ellipse(aes(fill = cluster), type = "norm", geom = "polygon", alpha = 0.12, show.legend = FALSE) +
  stat_summary(aes(group = cluster), fun = mean, geom = "point", shape = 21, size = 4, color = "black", stroke = 0.6, show.legend = FALSE) +
  geom_label(
    data = cluster_summary %>%
      mutate(cluster = factor(cluster)) %>%
      left_join(
        scores %>% group_by(cluster) %>% summarise(cx = mean(PC1), cy = mean(PC2)),
        by = "cluster"
      ) %>%
      mutate(
        label = paste0("mean Tcrit = ", round(mean_Tcrit, 1), "°C\nmean T50 = ", round(mean_T50, 1), "°C")
      ),
    aes(x = cx, y = cy, label = label, fill = NULL),
    color = "black",
    size = 3,
    label.size = 0.25,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("1" = "#e66101", "2" = "#5ab4ac")) +   # adjust colors if desired
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)"),
    title = "PCA of thermal thresholds (Tcrit, T50)",
    color = "Cluster"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray92"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  )

print(p_pca)
ggsave("Figure6_PCA_clusters.tiff", p_pca, width = 9, height = 6, dpi = 300)

 
 
 # create a simple grouping: T4 tests -> "Drought", all others -> "Control"
 scores <- scores %>%
   mutate(Drought = if_else(grepl("^T4", Test), "Drought", "Control"))
 
 # percent variance values from your PCA results
 pc1_pct <- pc1_pct  # e.g. 80.41  (numeric)

 p_box <- ggplot(scores, aes(x = Drought, y = PC1, fill = Drought)) +
   geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black", alpha = 0.9) +
   geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.6, aes(color = Drought)) +
   stat_summary(fun = median, geom = "point", shape = 95, size = 14, color = "black") + 
   scale_fill_manual(values = c("Control" = "#f48b7b", "Drought" = "#33bfc3")) +
   scale_color_manual(values = c("Control" = "#f48b7b", "Drought" = "#33bfc3")) +
   labs(
     x = "",
     y = paste0("PC1 (", pc1_pct, "%)"),
     title = "Effect of drought stress on PC1"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     legend.position = "right",
     panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     axis.title.x = element_blank()
   )
 
 print(p_box)
ggsave("PC1_drought_boxplot.tiff", p_box, width = 6, height = 5, dpi = 300)
 