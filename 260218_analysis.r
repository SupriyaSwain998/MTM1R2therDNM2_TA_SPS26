#volcano plot

##repeat analysis
log2FC_cutoff <- 0.5
neglog10_cutoff <- 1.3

#volcano plot

library(tidyverse)
library(ggrepel)

plot_volcano_proteomics <- function(df, title, log2fc_col, neglog10_col, gene_col) {
  
  df <- df %>%
    mutate(
      Significance = case_when(
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] >  log2FC_cutoff  ~ "Up",
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] < -log2FC_cutoff  ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  top_up <- df %>%
    filter(Significance == "Up") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  top_down <- df %>%
    filter(Significance == "Down") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  top_proteins <- bind_rows(top_up, top_down)
  
  ggplot(df, aes(x = .data[[log2fc_col]], y = .data[[neglog10_col]], color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_classic(base_size = 14) +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10(p-value)",
         color = "Regulation") +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = neglog10_cutoff,
               linetype = "dashed") +
    geom_text_repel(
      data = top_proteins,
      aes(label = .data[[gene_col]]),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 30
    )
}
getwd ()
setwd ("../Desktop/IGBMC/Mtmr2_Dnm2/")
disease_prot <- read.csv("260218_disease_control.csv", check.names = FALSE)
treated_prot <- read.csv("260218_treated_disease.csv", check.names = FALSE)
treated__WT_prot <- read.csv("260218_treated_control.csv", check.names = FALSE)
head (treated__WT_prot)
head (treated__WT_prot)
colnames (treated_prot$`-log10(p-value)`)
volcano_disease <- plot_volcano_proteomics(
  df            = disease_prot,
  title         = expression(italic("Mtmr2-/-") ~ "vs WT"),
  log2fc_col    = "Difference",
  neglog10_col  = "-log10(p-value)",
  gene_col      = "Gene Name"
)

volcano_disease
colnames (treated_prot)
volcano_treated <- plot_volcano_proteomics(
  df            = treated_prot,
  title         = expression(italic("Mtmr2-/-TgDnm2") ~ "vs WT"),
  log2fc_col    = "Difference",       
  neglog10_col  = "-log10(p-value)", 
  gene_col      = "Gene Name"
)

volcano_treated
colnames (treated__WT_prot$`-log10(p-value)`)
volcano_WT_prot <- plot_volcano_proteomics(
  df            = treated__WT_prot,
  title         = expression(italic("Mtmr2-/-TgDnm2") ~ "vs" ~ italic("Mtmr2-/-")),
  log2fc_col    = "Difference",       
  neglog10_col  = "-log10(p-value)", 
  gene_col      = "Gene Name"
)
volcano_WT_prot
ggsave("260429_output/260218_volcanoplot_disease.png", plot = volcano_disease, dpi = 300)
ggsave("260429_output/260218_volcanoplot_treated_vs_disease.png", plot = volcano_treated, dpi = 300)
ggsave("260429_output/260218_volcanoplot_treated_vs_control.png", plot = volcano_WT_prot, dpi = 300)

#specifically highlight "DNM2"

plot_volcano_proteomics_1 <- function(df, title, log2fc_col, neglog10_col, gene_col) {
  
  df <- df %>%
    mutate(
      Significance = case_when(
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] >  log2FC_cutoff  ~ "Up",
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] < -log2FC_cutoff  ~ "Down",
        TRUE ~ "NS"
      ),
      Highlight = ifelse(tolower(.data[[gene_col]]) == "dnm2", "dnm2", "Other")
    )
  
  # Top significant proteins
  top_up <- df %>%
    filter(Significance == "Up") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  top_down <- df %>%
    filter(Significance == "Down") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  # Extract dnm2 row (if present)
  dnm2_row <- df %>%
    filter(tolower(.data[[gene_col]]) == "dnm2")
  
  top_proteins <- bind_rows(top_up, top_down, dnm2_row)
  
  ggplot(df, aes(x = .data[[log2fc_col]], y = .data[[neglog10_col]])) +
    
    # All points
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    
    # Highlight dnm2 with bigger black circle
    geom_point(
      data = dnm2_row,
      color = "black",
      size = 4
    ) +
    
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    
    theme_classic(base_size = 14) +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10(p-value)",
         color = "Regulation") +
    
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = neglog10_cutoff,
               linetype = "dashed") +
    
    geom_text_repel(
      data = top_proteins,
      aes(label = .data[[gene_col]]),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 30
    )
}
volcano_disease_1 <- plot_volcano_proteomics_1(
  df            = disease_prot,
  title         = "Proteomics: Diseased vs Control",
  log2fc_col    = "Difference",
  neglog10_col  = "-log10(p-value)",
  gene_col      = "Gene Name"
)

volcano_disease_1

volcano_treated_1 <- plot_volcano_proteomics_1(
  df            = treated_prot,
  title         = "Proteomics: Treated vs Diseased",
  log2fc_col    = "Difference",       
  neglog10_col  = "-log10(p-value)", 
  gene_col      = "Gene Name"
)

volcano_treated_1

volcano_WT_prot_1 <- plot_volcano_proteomics_1(
  df            = treated__WT_prot,
  title         = "Proteomics: Treated vs Control",
  log2fc_col    = "Difference",       
  neglog10_col  = "-log10(p-value)", 
  gene_col      = "Gene Name"
)
volcano_WT_prot_1
ggsave("260429_output/260218_volcanoplot_disease_highlighted_DNM2.png", plot = volcano_disease_1, dpi = 300)
ggsave("260429_output/260218_volcanoplot_treated_vs_disease_highlighted_DNM2.png", plot = volcano_treated_1, dpi = 300)
ggsave("260429_output/260218_volcanoplot_treated_vs_control_highlighted_DNM2.png", plot = volcano_WT_prot_1, dpi = 300)

#highlighting both mtmr2 and dnm2:

# Volcano plot function highlighting DNM2 and MTMR2
plot_volcano_proteomics_highlight2 <- function(df, title, log2fc_col, neglog10_col, gene_col) {
  
  # Add significance column
  df <- df %>%
    mutate(
      Significance = case_when(
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] >  log2FC_cutoff  ~ "Up",
        .data[[neglog10_col]] >= neglog10_cutoff & .data[[log2fc_col]] < -log2FC_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Top significant proteins
  top_up <- df %>%
    filter(Significance == "Up") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  top_down <- df %>%
    filter(Significance == "Down") %>%
    arrange(desc(.data[[neglog10_col]])) %>%
    head(10)
  
  # Extract DNM2 and MTMR2 rows
  dnm2_row <- df %>%
    filter(tolower(.data[[gene_col]]) == "dnm2")
  
  mtmr2_row <- df %>%
    filter(tolower(.data[[gene_col]]) == "mtmr2")
  
  # Combine for labeling
  top_proteins <- bind_rows(top_up, top_down, dnm2_row, mtmr2_row) %>%
    distinct()
  
  # Plot
  ggplot(df, aes(x = .data[[log2fc_col]], y = .data[[neglog10_col]])) +
    
    # All points colored by significance
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    
    # Highlight DNM2 (black)
    geom_point(data = dnm2_row, color = "black", size = 4) +
    
    # Highlight MTMR2 (dark green)
    geom_point(data = mtmr2_row, color = "darkgreen", size = 4) +
    
    # Manual color for Up, Down, NS
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    
    # Labels
    geom_text_repel(
      data = top_proteins,
      aes(label = .data[[gene_col]]),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 30
    ) +
    
    # Threshold lines
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
    geom_hline(yintercept = neglog10_cutoff, linetype = "dashed") +
    
    # Theme and labels
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "Regulation"
    )
}

volcano_disease_2 <- plot_volcano_proteomics_highlight2(
  df = disease_prot,
  title = expression(italic("Mtmr2-/-") ~ "vs WT"),
  log2fc_col = "Difference",
  neglog10_col = "-log10(p-value)",
  gene_col = "Gene Name"
)
volcano_disease_2

volcano_treated_2 <- plot_volcano_proteomics_highlight2(
  df = treated_prot,
  title = expression(italic("Mtmr2-/-TgDnm2") ~ "vs" ~ italic("Mtmr2-/-")),
  log2fc_col = "Difference",
  neglog10_col = "-log10(p-value)",
  gene_col = "Gene Name"
)

volcano_WT_prot_2 <- plot_volcano_proteomics_highlight2(
  df = treated__WT_prot,
  title = expression(italic("Mtmr2-/-TgDnm2") ~ "vs WT"),
  log2fc_col = "Difference",
  neglog10_col = "-log10(p-value)",
  gene_col = "Gene Name"
)
# Save plot
ggsave("260429_output/260223_volcanoplot_disease_highlighted_DNM2_MTM2.png",
       plot = volcano_disease_2, dpi = 300)
ggsave("260429_output/260223_volcanoplot_treated_vs_disease_highlighted_DNM2_MTM2.png",
       plot = volcano_treated_2, dpi = 300)
ggsave("260429_output/260223_volcanoplot_treated_vs_control_highlighted_DNM2_MTM2.png",
       plot = volcano_WT_prot_2, dpi = 300)

#GO enrichment

disease_sig <- disease_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff)

disease_genes <- unique(disease_sig$`Gene Name`)
length(disease_genes)
treated_sig <- treated_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff)

treated_genes <- unique(treated_sig$`Gene Name`)
length(treated_genes)
treated_WT_sig <- treated__WT_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff)

treated_WT_genes <- unique(treated_WT_sig$`Gene Name`)
length(treated_WT_genes)
library(clusterProfiler)
library(org.Mm.eg.db)

disease_entrez <- mapIds(org.Mm.eg.db, keys = disease_genes,
                         column = "ENTREZID", keytype = "SYMBOL",
                         multiVals = "first")

treated_entrez <- mapIds(org.Mm.eg.db, keys = treated_genes,
                         column = "ENTREZID", keytype = "SYMBOL",
                         multiVals = "first")

treated_WT_entrez <- mapIds(org.Mm.eg.db, keys = treated_WT_genes,
                            column = "ENTREZID", keytype = "SYMBOL",
                            multiVals = "first")
#remove NAs
disease_entrez   <- na.omit(disease_entrez)
treated_entrez   <- na.omit(treated_entrez)
treated_WT_entrez <- na.omit(treated_WT_entrez)
run_GO_all_ont <- function(entrez_ids, label) {
  
  ego_BP <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  ego_MF <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  ego_CC <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
  
  list(BP = ego_BP, MF = ego_MF, CC = ego_CC, label = label)
}


go_disease <- run_GO_all_ont(
  disease_entrez,
  expression("Disease Signature " * italic("Mtmr2-/-") * " vs WT")
)
go_treated   <- run_GO_all_ont(treated_entrez, "Treated vs Diseased")
go_treatedWT <- run_GO_all_ont(treated_WT_entrez, "Treated vs WT")

save_go_combined <- function(go_list, prefix) {
  
  bp <- as.data.frame(go_list$BP)
  mf <- as.data.frame(go_list$MF)
  cc <- as.data.frame(go_list$CC)
  
  bp$Ontology <- "BP"
  mf$Ontology <- "MF"
  cc$Ontology <- "CC"
  
  combined <- rbind(bp, mf, cc)
  
  write.csv(combined,
            paste0(prefix, "_GO_All_Ontologies.csv"),
            row.names = FALSE)
}

save_go_combined(go_disease, "260429_output/Diseased_vs_Control")
save_go_combined(go_treated, "260429_output/Treated_vs_Diseased")
save_go_combined(go_treatedWT, "260429_output/Treated_vs_WT")

nrow (go_disease)

library(dplyr)

combine_GO_results <- function(go_list_object) {
  
  bp_df <- as.data.frame(go_list_object$BP) %>% 
    head(20) %>% 
    mutate(Ontology = "BP")
  
  mf_df <- as.data.frame(go_list_object$MF) %>% 
    head(20) %>% 
    mutate(Ontology = "MF")
  
  cc_df <- as.data.frame(go_list_object$CC) %>% 
    head(20) %>% 
    mutate(Ontology = "CC")
  
  bind_rows(bp_df, mf_df, cc_df)
}
go_disease_df   <- combine_GO_results(go_disease)
go_treated_df   <- combine_GO_results(go_treated)
go_treatedWT_df <- combine_GO_results(go_treatedWT)
library(ggplot2)

plot_GO_combined <- function(go_df, title_text) {
  ggplot(go_df, aes(x = reorder(Description, Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~Ontology, , ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("BP" = "steelblue", 
                                 "MF" = "darkorange", 
                                 "CC" = "forestgreen")) +
    labs(title = title_text,
         x = "GO Term",
         y = "Gene Count") +
    theme_classic(base_size = 12)
}
p_disease_combined <- plot_GO_combined(go_disease_df, 
                                       expression("Disease signature " * italic("(Mtmr2-/-") * " vs WT)"))

p_treated_combined <- plot_GO_combined(go_treated_df, 
                                       "GO Enrichment (BP/MF/CC): Treated vs Diseased")

p_treatedWT_combined <- plot_GO_combined(go_treatedWT_df, 
                                         "GO Enrichment (BP/MF/CC): Treated vs WT")

p_disease_combined
p_treated_combined
p_treatedWT_combined
ggsave("260429_output/260218_GO_combined_disease.png", p_disease_combined, width = 12, height = 10, dpi = 300)
ggsave("260429_output/260218_GO_combined_treated.png", p_treated_combined, width = 12, height = 10, dpi = 300)
ggsave("260429_output/260218_GO_combined_treatedWT.png", p_treatedWT_combined, width = 12, height = 10, dpi = 300)
write.csv(as.data.frame(go_disease_df), "260429_output/260218_GO_disease.csv", row.names = FALSE)
write.csv(as.data.frame(go_treated_df), "260429_output/260218_GO_treated.csv", row.names = FALSE)
write.csv(as.data.frame(go_treatedWT_df), "260429_output/260218_GO_treated_WT.csv", row.names = FALSE)

#rescue matrix


colnames(disease_prot) <- make.unique(colnames(disease_prot))
colnames(treated_prot) <- make.unique(colnames(treated_prot))

diseased_sig <- disease_prot[
  abs(disease_prot$Difference) >= log2FC_cutoff &
    disease_prot$`-log10(p-value)` >= neglog10_cutoff,]
cat("Number of significant disease genes:", nrow(diseased_sig), "\n")
nrow (diseased_sig)
common_genes <- intersect(diseased_sig$`Gene Name`, treated_prot$`Gene Name`)
diseased_matched <- diseased_sig[diseased_sig$`Gene Name` %in% common_genes, ]
treated_matched  <- treated_prot[treated_prot$`Gene Name` %in% common_genes, ]
nrow (treated_matched)

treated_matched <- treated_matched[match(diseased_matched$`Gene Name`,
                                         treated_matched$`Gene Name`), ] #same order
nrow (diseased_matched)
log2FC_disease <- diseased_matched$Difference
log2FC_treated <- treated_matched$Difference   # leading space in treated file!
log2FC_treated_vs_WT <- log2FC_disease + log2FC_treated
nrow (log2FC_treated_vs_WT)
rescue_metric <- 100 * (log2FC_disease - log2FC_treated_vs_WT) / log2FC_disease
rescue_table_prot$rescue_metric <- rescue_metric
rescue_table_prot <- data.frame(
  Gene = diseased_matched$`Gene Name`,
  log2FC_disease = log2FC_disease,
  log2FC_treated = log2FC_treated,
  neglog10_disease = diseased_matched$`-log10(p-value)`,
  neglog10_treated = treated_matched$`-log10(p-value)`,
  rescue_metric = rescue_metric
)
head (rescue_table_prot)

rescue_table_prot$rescue_status <- cut(
  rescue_table_prot$rescue_metric,
  breaks = c(-Inf, 0, 30, 80, 120, Inf),
  labels = c("Worsened", "Not rescued", "Partially rescued", "Rescued", "over-rescued")
)
table(rescue_table_prot$rescue_status)
print (rescue_table_prot)
write.csv(rescue_table_prot,
          file = "260429_output/260224_rescue_analysis.csv",
          row.names = FALSE)

library(ggplot2)

ggplot(rescue_table_prot, aes(x = rescue_metric, fill = rescue_status)) +
  geom_histogram(bins = 40, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Rescue Metric",
       x = "Rescue Metric",
       y = "Protein Count")
rescued_prots <- subset(rescue_table_prot, rescue_status %in% c("Partially rescued","Rescued", "over-rescued"))

nrow(rescued_prots)
print(rescued_prots)
background_genes <- unique(rescue_table_prot$Gene)

cat("treatment worked", median(rescue_table_prot$rescue_metric), "\n")

# ============================
# Plot distribution
# ============================

png(
  "260429_output/260218_rescue_metric_distribution.png",
  width = 1200,
  height = 900,
  res = 150
)

hist(
  rescue_table_prot$rescue_metric,
  breaks = 100,
  main = "Distribution of rescue metric",
  xlab = "Rescue metric"
)

abline(v = 0, col = "red", lwd = 2)

dev.off()

#percentage rescue

rescue_summary_prot <- rescue_table_prot %>%
  filter(rescue_status %in% c("Rescued", "Partially rescued", "Not rescued", "over-rescued")) %>%
  count(rescue_status) %>%
  mutate(
    percent = n / sum(n) * 100,
    Treatment = "Dnm2 overexpression"
  )

print(rescue_summary_prot)

plot_rescue <- ggplot(rescue_summary_prot,
                      aes(x = Treatment, y = percent, fill = rescue_status)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(
    values = c(
      "Rescued" = "darkgreen",
      "Partially rescued" = "yellow",
      "Not rescued" = "darkgrey", "over-rescued" = "darkblue"
    )
  ) +
  labs(
    x = NULL,
    y = "% of Dysregulated Proteins in Mtmr2-/-",
    fill = "Rescue Status"
  ) +
  ylim(0, 100) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(plot_rescue)
ggsave("260429_output/260219_Final_stackplot_rescue.png", plot = plot_rescue, dpi = 300)
log2FC_cutoff <- 0.5
neglog10_cutoff <- 1.3

# Disease vs Control
disease_sig_genes <- disease_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff) %>%
  pull(`Gene Name`) %>%
  unique()

# Treated vs Diseased
treated_sig_genes <- treated_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff) %>%
  pull(`Gene Name`) %>%
  unique()

# Treated vs WT
treated_WT_sig_genes <- treated__WT_prot %>%
  filter(abs(`Difference`) >= log2FC_cutoff,
         `-log10(p-value)` >= neglog10_cutoff) %>%
  pull(`Gene Name`) %>%
  unique()
library(VennDiagram)
library(grid)   # needed to draw
venn_plot <- venn.diagram(
  x = list(
    Disease_vs_Control = disease_sig_genes,
    Treated_vs_Disease = treated_sig_genes,
    Treated_vs_WT      = treated_WT_sig_genes
  ),
  filename = NULL,  # important: keeps plot in R
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.4,
  cat.pos = c(-20, 20, 180),
  cat.dist = 0.05,
  main = "Overlap of Significantly Dysregulated Proteins"
)

grid.draw(venn_plot)
png("260429_output/260219_venn_protein_overlap.png", width = 2000, height = 2000, res = 300)
grid.draw(venn_plot)
dev.off()


#GO of the rescue
library(clusterProfiler)
library(org.Mm.eg.db)  # mouse annotation
library(dplyr)
head (rescue_table_prot)
length (background_genes)
rescued_prots <- rescue_table_prot %>%
  dplyr::filter(rescue_status %in% 
                  c("Partially rescued","Rescued","over-rescued")) %>%
  dplyr::select(Gene)

print (rescued_prots)
rescued_prot_entrez <- mapIds(
  org.Mm.eg.db,
  keys = rescued_prots$Gene,   # <-- this is the fix
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
background_prot_entrez = mapIds(
  org.Mm.eg.db,
  keys = background_genes,   # <-- this is the fix
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

length (rescued_prot_entrez)
ego_rescued_prot <- enrichGO(
  gene = rescued_prot_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",       # BP = Biological Process, you can also use "MF" or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head (ego_rescued_prot)
p1 = barplot(ego_rescued_prot, showCategory = 50, title = "Enriched GO-BP terms of proteins rescued by the treatment", font.size = 8)
print (p1)
write.csv (ego_rescued_prot, file = "260429_output/260219_final_rescued_0.05_GO_BP.csv")

ggsave (filename = "260429_output/260219_final_rescued_0.05_GO_BP.png", plot = p1, width = 15, height = 10, dpi = 300)
print(
  rescue_table_prot %>%
    dplyr::filter(rescue_status == "over-rescued")
)


#All components
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

# Step 1: Prepare Entrez IDs for rescued proteins
rescued_prots <- rescue_table_prot %>%
  dplyr::filter(rescue_status %in% c("Partially rescued", "Rescued", "over-rescued")) %>%
  dplyr::select(Gene)

rescued_prot_entrez <- mapIds(
  org.Mm.eg.db,
  keys = rescued_prots$Gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
) %>% na.omit()  # remove NAs

# Step 2: Function to run GO for all 3 ontologies
run_GO_all_ont_rescued <- function(entrez_ids, label) {
  
  ego_BP <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_MF <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_CC <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  list(BP = ego_BP, MF = ego_MF, CC = ego_CC, label = label)
}

# Step 3: Run GO enrichment for rescued proteins
go_rescued <- run_GO_all_ont_rescued(rescued_prot_entrez, "Rescued Proteins")
save_go_combined <- function(go_list, prefix) {
  
  bp <- as.data.frame(go_list$BP)
  mf <- as.data.frame(go_list$MF)
  cc <- as.data.frame(go_list$CC)
  
  bp$Ontology <- "BP"
  mf$Ontology <- "MF"
  cc$Ontology <- "CC"
  
  combined <- rbind(bp, mf, cc)
  
  write.csv(combined,
            paste0(prefix, "_GO_All_Ontologies.csv"),
            row.names = FALSE)
}

save_go_combined(go_rescued, "260429_output/260224_rescue_proteins")



# Step 4: Combine GO results into one dataframe for plotting
combine_GO_results <- function(go_list_object) {
  bp_df <- as.data.frame(go_list_object$BP) %>% head(20) %>% mutate(Ontology = "BP")
  mf_df <- as.data.frame(go_list_object$MF) %>% head(20) %>% mutate(Ontology = "MF")
  cc_df <- as.data.frame(go_list_object$CC) %>% head(20) %>% mutate(Ontology = "CC")
  bind_rows(bp_df, mf_df, cc_df)
}

go_rescued_df <- combine_GO_results(go_rescued)

# Step 5: Plot combined GO (BP/MF/CC)
plot_GO_combined <- function(go_df, title_text) {
  ggplot(go_df, aes(x = reorder(Description, Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~Ontology, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("BP" = "steelblue", "MF" = "darkorange", "CC" = "forestgreen")) +
    labs(title = title_text, x = "GO Term", y = "Gene Count") +
    theme_classic(base_size = 12)
}

p_rescued_combined <- plot_GO_combined(
  go_rescued_df, 
  expression("Rescued proteins (" * italic("Mtmr2-/-TgDnm2") * " vs " * italic("Mtmr2-/-") * ")")
)

# Step 6: Visualize and save
p_rescued_combined
ggsave("260429_output/260223_GO_combined_rescued.png",
       p_rescued_combined, width = 12, height = 10, dpi = 300)

# Step 7: Save the GO table
write.csv(as.data.frame(go_rescued_df), 
          "260429_output/260223_GO_combined_rescued.csv", 
          row.names = FALSE)


library(ReactomePA)
length (rescued_prot_entrez)
reactome_rescue <- enrichPathway(gene          = rescued_prot_entrez,
                                 organism      = "mouse",
                                 pvalueCutoff  = 0.05,
                                 readable      = TRUE)

library(ggplot2)

# Assign dotplot to a variable
p_reactome_rescue <- dotplot(reactome_rescue, showCategory = 15, title = "Rescued Reactome Pathways")

# Print the plot
print(p_reactome_rescue)
ggsave (filename = "260429_output/260219_reactome_rescued_0.05_GO_BP.png", plot = p_reactome_rescue, width = 15, height = 10, dpi = 300)
print (reactome_rescue_mtrix)
reactome_rescue_mtrix = data.frame(reactome_rescue)
write.csv (reactome_rescue_mtrix, file = "260429_output/260219_reactome_rescued_0.05_GO_BP.csv")


#GSEA




library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)

# Function to convert gene symbols to Entrez
convert_ids <- function(genes) {
  df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  unique(df$ENTREZID)
}
head (disease_sig_genes)
# Disease ranked list
disease_ranked <- disease_prot %>%
  dplyr::select(`Gene Name`, `Difference`) %>%
  distinct() %>%
  mutate(ENTREZID = mapIds(org.Mm.eg.db,
                           keys = `Gene Name`,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")) %>%
  na.omit()

geneList_disease <- disease_ranked$`Difference`
names(geneList_disease) <- disease_ranked$ENTREZID
geneList_disease <- sort(geneList_disease, decreasing = TRUE)

# Treated ranked list
treated_ranked <- treated_prot %>%
  dplyr::select(`Gene Name`, `Difference`) %>%
  distinct() %>%
  mutate(ENTREZID = mapIds(org.Mm.eg.db,
                           keys = `Gene Name`,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")) %>%
  na.omit()

geneList_treated <- treated_ranked$`Difference`
names(geneList_treated) <- treated_ranked$ENTREZID
geneList_treated <- sort(geneList_treated, decreasing = TRUE)
gsea_disease <- gsePathway(geneList     = geneList_disease,
                           organism     = "mouse",
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)

gsea_treated <- gsePathway(geneList     = geneList_treated,
                           organism     = "mouse",
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)

treated_WT_ranked <- treated__WT_prot %>%
  dplyr::select(`Gene Name`, `Difference`) %>%
  distinct() %>%
  mutate(ENTREZID = mapIds(org.Mm.eg.db,
                           keys = `Gene Name`,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")) %>%
  na.omit()

geneList_treated_WT <- treated_WT_ranked$`Difference`
names(geneList_treated_WT) <- treated_WT_ranked$ENTREZID
geneList_treated_WT<- sort(geneList_treated_WT, decreasing = TRUE)
gsea_treated_WT <- gsePathway(geneList     = geneList_treated_WT,
                           organism     = "mouse",
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)


disease_df <- as.data.frame(gsea_disease)[, c("ID", "Description", "NES", "p.adjust")]
treated_df <- as.data.frame(gsea_treated)[, c("ID", "NES", "p.adjust")]
colnames(treated_df)[2:3] <- c("NES_treated", "p.adjust_treated")

# Merge by pathway ID
merged_gsea <- merge(disease_df, treated_df, by = "ID")
merged_gsea$NES_change <- merged_gsea$NES_treated - merged_gsea$NES

# Pathways rescued (reversal of dysregulation)
rescued_pathways <- merged_gsea %>%
  filter((NES > 0 & NES_treated < 0) | (NES < 0 & NES_treated > 0))
head(rescued_pathways)
library(ggplot2)
library(grid)

# Top 20 rescued pathways
top20_rescued <- rescued_pathways[1:20, ]

ggplot(top20_rescued, aes(y = reorder(Description, NES))) +
  geom_segment(aes(x = NES, xend = NES_treated, yend = Description),
               arrow = arrow(length = unit(0.2,"cm")), color = "grey") +
  geom_point(aes(x = NES), color = "red", size = 3) +
  geom_point(aes(x = NES_treated), color = "blue", size = 3) +
  labs(title = "Reactome Pathways Reversed by Treatment",
       x = "NES (Red = Disease, Blue = Treated)",
       y = "Pathway") +
  theme_classic()
# Convert to data frames
gsea_disease_df <- as.data.frame(gsea_disease)
gsea_treated_df <- as.data.frame(gsea_treated)
gsea_treated_WT_df = as.data.frame(gsea_treated_WT)

# Save as CSV
write.csv(gsea_disease_df, "260429_output/260219_GSEA_disease.csv", row.names = FALSE)
write.csv(gsea_treated_df, "260429_output/260219_GSEA_treated.csv", row.names = FALSE)
library(enrichplot)  # required for gseaplot2
library(ggplot2)

# Top 5 enriched pathways in disease
top_disease <- head(gsea_disease_df$ID, 5)
for (pid in top_disease) {
  p <- gseaplot2(gsea_disease, geneSetID = pid, title = gsea_disease_df$Description[gsea_disease_df$ID == pid])
  print(p)
  ggsave(filename = paste0("260429_output/260219_GSEA_disease_", pid, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

# Top 5 enriched pathways in treated
top_treated <- head(gsea_treated_df$ID, 5)
for (pid in top_treated) {
  p <- gseaplot2(gsea_treated, geneSetID = pid, title = gsea_treated_df$Description[gsea_treated_df$ID == pid])
  print(p)
  ggsave(filename = paste0("260429_output/260219_GSEA_treated_", pid, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

library(ggplot2)
library(reshape2)

# Prepare data for heatmap
heatmap_df <- merged_gsea[, c("Description", "NES", "NES_treated")]
colnames(heatmap_df) <- c("Pathway", "Disease", "Treated")

# Optional: select top pathways by absolute change in NES
heatmap_df$NES_change <- abs(heatmap_df$Treated - heatmap_df$Disease)
top_heatmap <- heatmap_df[order(-heatmap_df$NES_change), ][1:20, ]

# Melt for ggplot
melted <- melt(top_heatmap[, c("Pathway","Disease","Treated")], id.vars = "Pathway")
melted$Pathway <- factor(melted$Pathway, levels = rev(top_heatmap$Pathway))  # for ordering

# Heatmap
p5 = ggplot(melted, aes(x = variable, y = Pathway, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  labs(title = "Comparative GSEA NES: Disease vs Treated",
       x = "",
       y = "Pathway",
       fill = "NES") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave (filename = "260429_output/260219_GSEA_pathway.png", plot = p5, width = 15, height = 10, dpi = 300)


#biological process gsea
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# GSEA for GO BP
gseaGO_disease <- gseGO(geneList     = geneList_disease,
                        OrgDb        = org.Mm.eg.db,
                        ont          = "BP",
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)

gseaGO_treated <- gseGO(geneList     = geneList_treated,
                        OrgDb        = org.Mm.eg.db,
                        ont          = "BP",
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)

# Convert to data frames
gseaGO_disease_df <- as.data.frame(gseaGO_disease)
gseaGO_treated_df <- as.data.frame(gseaGO_treated)

# Save as CSV
write.csv(gseaGO_disease_df, "260429_output/260219_GSEA_GO_disease.csv", row.names = FALSE)
write.csv(gseaGO_treated_df, "260429_output/260219_GSEA_GO_treated.csv", row.names = FALSE)

# Merge by pathway ID
merged_gseaGO <- merge(gseaGO_disease_df, gseaGO_treated_df, by = "ID")
merged_gseaGO$NES_change <- merged_gsea$NES_treated - merged_gsea$NES

# Rename columns after merge
colnames(merged_gseaGO)[which(colnames(merged_gseaGO) == "Description.x")] <- "Description"
colnames(merged_gseaGO)[which(colnames(merged_gseaGO) == "NES.x")] <- "NES_disease"
colnames(merged_gseaGO)[which(colnames(merged_gseaGO) == "NES.y")] <- "NES_treated"

# Recalculate NES change
merged_gseaGO$NES_change <- merged_gseaGO$NES_treated - merged_gseaGO$NES_disease

# Identify rescued GO BP terms
rescued_GO <- merged_gseaGO %>%
  filter((NES_disease > 0 & NES_treated < 0) | (NES_disease < 0 & NES_treated > 0))

# Top 20 rescued GO terms
top20_rescued_GO <- rescued_GO[1:20, ]

# Comparative NES plot
p6 = ggplot(top20_rescued_GO, aes(y = reorder(Description, NES_disease))) +
  geom_segment(aes(x = NES_disease, xend = NES_treated, yend = Description),
               arrow = arrow(length = unit(0.2, "cm")), color = "grey") +
  geom_point(aes(x = NES_disease), color = "red", size = 3) +
  geom_point(aes(x = NES_treated), color = "blue", size = 3) +
  labs(title = "GO BP Terms Reversed by Treatment",
       x = "NES (Red = Disease, Blue = Treated)",
       y = "GO Term") +
  theme_classic()
ggsave (filename = "260429_output/260219_GSEA_GO.png", plot = p6, width = 15, height = 10, dpi = 300)


# 1. Load libraries
# -------------------------
library(tidyverse)

# -------------------------
# 2. Sample annotation
# -------------------------
coldata <- tibble(
  Sample = c("XIC S03", "XIC S03bis","XIC S04", "XIC S05",
             "XIC S06","XIC S06bis","XIC S07","XIC S07bis",
             "XIC S08", "XIC S09", "XIC S10",
             "XIC S11","XIC S11bis","XIC S12","XIC S12bis",
             "XIC S13","XIC S14","XIC S15"),
  Condition = c("WT","WT","WT","WT",
                "KO","KO","KO","KO","KO","KO","KO",
                "KO_Tg","KO_Tg","KO_Tg","KO_Tg","KO_Tg","KO_Tg","KO_Tg")
)

# -------------------------
# 3. Load proteomics matrix
# -------------------------
prot <- read.csv("260429_output/your_proteomics_matrix_1.csv",
                 check.names = FALSE)

prot <- prot[, colnames(prot) != ""]
colnames(prot)[1] <- "Protein"

# -------------------------
# 4. Pivot to long format
# -------------------------
long_df <- prot %>%
  pivot_longer(
    cols = -Protein,
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  left_join(coldata, by = "Sample") %>%
  mutate(Intensity = as.numeric(Intensity))

# -------------------------
# 5. Proteins of interest
# -------------------------
proteins_of_interest <- c("ARHGAP18", "CAPZA1", "CAPZA2", "CDC42EP2", "CFL1",
                          "MYADM", "PRKCD", "SPIRE2", "TPM1", "ATP2A1",
                          "CAV1", "EMP2", "ROCK1", "SRI", "CORO2B",
                          "CTTNBP2NL", "MYH11", "MYL9", "MYLK", "PALLD",
                          "PDLIM5", "PPP1R12B", "TPM3")

plot_df <- long_df %>%
  filter(Protein %in% proteins_of_interest)

# -------------------------
# 6. ANOVA + Tukey
# -------------------------
anova_results <- list()
tukey_results <- list()
library(ggpubr)


# -------------------------
# 7. Plotting with Tukey stars (p < 0.05 only)
# -------------------------

# -------------------------
# 6â€“7. ANOVA + Stars Plotting
# -------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)  # for stat_compare_means

for(prot_name in proteins_of_interest){
  
  df_prot <- plot_df %>% 
    filter(Protein == prot_name) %>%
    drop_na(Intensity)
  
  df_prot$Condition <- factor(df_prot$Condition,
                              levels = c("WT","KO","KO_Tg"))
  
  if(length(unique(df_prot$Condition)) < 2) next
  
  # ANOVA for subtitle
  model <- aov(Intensity ~ Condition, data = df_prot)
  anova_table <- summary(model)[[1]]
  anova_p <- anova_table$`Pr(>F)`[1]
  
  # -------------------------
  # Define pairwise comparisons
  # -------------------------
  comparisons <- list(c("WT", "KO"),
                      c("WT", "KO_Tg"),
                      c("KO", "KO_Tg"))
  
  # -------------------------
  # Plot with stars for each comparison
  # -------------------------
  p <- ggplot(df_prot,
              aes(x = Condition, y = Intensity, fill = Condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    theme_classic(base_size = 14) +
    labs(title = paste0("Protein: ", prot_name),
         #subtitle = paste0("ANOVA p = ", signif(anova_p, 3)),
         y = "log2 Intensity",
         x = NULL) +
    scale_fill_manual(values = c("WT"="white",
                                 "KO"="red",
                                 "KO_Tg"="blue")) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold")) +
    
    stat_compare_means(comparisons = comparisons,
                       method = "t.test",  # or "wilcox.test" if non-parametric
                       label = "p.signif",
                       tip.length = 0.02, # length of line
                       size = 6)          # text size
  
  # Save
  ggsave(filename = paste0("260429_output/", prot_name,
                           "_boxplot.png"),
         plot = p,
         width = 6,
         height = 5)
}

#Pearson correlation between two wts


# ----------------------------
# 1. Load data
# ----------------------------
file_path <- "260429_output/260218_disease_control.csv"

data <- read.csv(file_path, header = TRUE, check.names = FALSE)
head (data)
# ----------------------------
# 2. Extract WT columns
# ----------------------------
WT1 <- data$`XIC S03`
WT2= data$`XIC S03bis`
WT3 <- data$`XIC S04`
WT4 = data$`XIC S05`

# Remove rows with NA in WT
valid_rows <- complete.cases(WT1, WT2, WT3, WT4)
WT1 <- WT1[valid_rows]
WT2 <- WT2[valid_rows]
WT3 = WT3[valid_rows]
WT4 = WT4[valid_rows]

# Combine all WT samples into matrix
WT_matrix <- cbind(WT1, WT2, WT3, WT4)

# Compute correlation matrix
correlation_matrix <- cor(WT_matrix, method = "pearson")

print(correlation_matrix)
# Mean pairwise correlation (excluding diagonal)
mean_cor <- mean(correlation_matrix[upper.tri(correlation_matrix)])
print(paste("Mean pairwise WT correlation:", mean_cor))
pairs(WT_matrix,
      main = "Scatterplot Matrix of WT Replicates",
      pch = 16)

# ----------------------------
# Calculate CV for 4 WT samples
# ----------------------------

WT_linear <- 2^WT_matrix

means_linear <- rowMeans(WT_linear)
sds_linear <- apply(WT_linear, 1, sd)

CV <- (sds_linear / means_linear) * 100

median_CV <- median(CV, na.rm = TRUE)

print(paste("Corrected Median WT CV (4 replicates):", median_CV))

# Now plot
hist(CV,
     breaks = 50,
     main = "Distribution of WT CV% (4 replicates)",
     xlab = "CV (%)")

#hierarchail clustering

library(readxl)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
# Read file (Excel)
prot <- read.csv("260429_output/260218_disease_control.csv", header = TRUE, check.names = FALSE)
head (prot_mat)

# Extract intensity columns only
intensity_cols <- grep("^XIC", colnames(prot), value = TRUE)
print (intensity_cols)

prot_mat <- prot[, intensity_cols]

# Use Gene Name as rownames (better for interpretation)
rownames(prot_mat) <- prot$`Gene Name`


# Convert to matrix
prot_mat <- as.matrix(prot_mat)
head (prot_mat)
coldata <- data.frame(
  sample = colnames(prot_mat),
  condition = c("WT","WT", "WT", "WT",
                "KO","KO","KO","KO","KO","KO", "KO",
                "KO_Tg","KO_Tg","KO_Tg","KO_Tg","KO_Tg","KO_Tg", "KO_Tg")
)

rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(prot_mat), , drop = FALSE]
# Distance between samples
dist_mat <- dist(t(prot_mat), method = "euclidean")

# Annotation
annotation_col <- data.frame(condition = coldata$condition)
rownames(annotation_col) <- rownames(coldata)

condition_colors <- c(
  "WT" = "black",
  "KO" = "red",
  "KO_Tg" = "blue"
)

library (pheatmap)
pheatmap(
  as.matrix(dist_mat),
  annotation_col = annotation_col,
  annotation_colors = list(condition = condition_colors),
  show_colnames = TRUE,
  labels_col = rownames(coldata),
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  clustering_method = "complete",
  main = "Hierarchical clustering â€“ Proteomics (All proteins)"
)

dev.off()
# Calculate variance
protein_vars <- apply(prot_mat, 1, var, na.rm = TRUE)

pca_res <- prcomp(t(prot_mat), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- rownames(pca_df)
pca_df$condition <- coldata[rownames(pca_df), "condition"]

percentVar <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

p <- ggplot(pca_df,
            aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = condition_colors) +
  xlab(paste0("PC1 (", round(percentVar[1],1), "%)")) +
  ylab(paste0("PC2 (", round(percentVar[2],1), "%)")) +
  theme_classic() +
  ggtitle("PCA â€“ Proteomics")
print (p)

ggsave("260429_output/Proteomics_PCA.png",
       plot = p, width = 8, height = 5, dpi = 300)