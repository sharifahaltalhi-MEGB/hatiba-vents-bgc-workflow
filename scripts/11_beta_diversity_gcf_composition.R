# =========================================================
# GCF beta diversity analysis by habitat type
#
# Inputs:
#   - mix_clustering_c0.7_v2.tsv
#   - bgc_level_antismash_v8_with_genome_metadata.csv
#
# Outputs:
#   - Table_PERMANOVA_GCF_Jaccard_TypeOnly.csv
#   - GCF_betadisper_distances_to_centroid_TypeOnly.csv
#   - GCF_PCoA_coordinates_TypeOnly.csv
#   - Fig_GCF_PCoA_byType_TypeOnly.pdf
#   - Fig_GCF_Betadisper_byType_TypeOnly.pdf
#   - GCF_sample_presence_absence_full_TypeOnly.csv
#   - GCF_sample_jaccard_distance_TypeOnly.tsv
#
# Notes:
#   - Run from a directory containing the input files
# =========================================================


rm(list = ls())

library(tidyverse)
library(vegan)
library(ggrepel)

clean_bgc <- function(x) {
    x %>%
        as.character() %>%
        str_trim() %>%
        str_remove("\\.gbk$") %>%
        str_remove("^\\./") %>%
        str_replace_all("\\s+", "")
}

# -----------------------------
# 1) Read BiG-SCAPE clustering (BGC -> GCF)
#    Expected columns include: RedSea_BGC_ID and Family (e.g., FAM_00001)
# -----------------------------
clust_raw <- read_tsv("mix_clustering_c0.7_v2.tsv", show_col_types = FALSE)

stopifnot("RedSea_BGC_ID" %in% names(clust_raw))
stopifnot("Family" %in% names(clust_raw))

clust <- clust_raw %>%
    transmute(
        BGC = clean_bgc(RedSea_BGC_ID),
        GCF = as.character(Family) %>% str_trim()
    ) %>%
    filter(!is.na(GCF), GCF != "", !is.na(BGC), BGC != "") %>%
    distinct()

cat("Input rows in clustering file:", nrow(clust_raw), "\n")
cat("Unique BGCs in clustering:", n_distinct(clust$BGC), "\n")
cat("Unique GCFs in clustering:", n_distinct(clust$GCF), "\n")

# -----------------------------
# 2) Map BGC to genome_folder + sample metadata (ONLY Type and Site)
# -----------------------------
bgc_level <- read_csv("bgc_level_antismash_v8_with_genome_metadata.csv", show_col_types = FALSE)

stopifnot("RedSea_BGC_ID" %in% names(bgc_level))
stopifnot("user_genome" %in% names(bgc_level))
stopifnot("Sample.ID" %in% names(bgc_level))
stopifnot("Type" %in% names(bgc_level))
stopifnot("Site" %in% names(bgc_level))

bgc_map <- bgc_level %>%
    transmute(
        BGC = clean_bgc(RedSea_BGC_ID),
        genome_folder = user_genome,
        SampleID = `Sample.ID`,
        Site = Site,
        Type = Type
    ) %>%
    filter(!is.na(BGC), BGC != "") %>%
    distinct()

cat("BGCs in bgc_level table:", n_distinct(bgc_map$BGC), "\n")
cat("Genomes in bgc_level table:", n_distinct(bgc_map$genome_folder), "\n")
cat("Samples in bgc_level table:", n_distinct(bgc_map$SampleID), "\n")

# -----------------------------
# 3) Join clustering with map (keeps only your BGCs)
# -----------------------------
bgc_gcf <- clust %>%
    inner_join(bgc_map, by = "BGC")

cat("BGCs after join:", n_distinct(bgc_gcf$BGC), "\n")
cat("Rows after join:", nrow(bgc_gcf), "\n")
cat("Unique genomes represented:", n_distinct(bgc_gcf$genome_folder), "\n")
cat("Unique samples represented:", n_distinct(bgc_gcf$SampleID), "\n")
cat("Unique GCFs represented:", n_distinct(bgc_gcf$GCF), "\n")

# -----------------------------
# 4) Genome-level GCF presence/absence matrix
# -----------------------------
gcf_pa_genome <- bgc_gcf %>%
    distinct(genome_folder, GCF) %>%
    mutate(presence = 1) %>%
    pivot_wider(
        names_from = GCF,
        values_from = presence,
        values_fill = 0
    )

cat("Genome x GCF matrix dimensions:", paste(dim(gcf_pa_genome), collapse = " x "), "\n")

# -----------------------------
# 5) Genome-level metadata
# -----------------------------
meta_genome <- bgc_map %>%
    distinct(genome_folder, SampleID, Site, Type)

gcf_pa_genome_meta <- gcf_pa_genome %>%
    left_join(meta_genome, by = "genome_folder")

cat("Genomes with missing SampleID:", sum(is.na(gcf_pa_genome_meta$SampleID)), "\n")

# -----------------------------
# 6) Collapse to sample level: GCF present if any genome in sample has it
# -----------------------------
gcf_pa_sample <- gcf_pa_genome_meta %>%
    pivot_longer(
        cols = starts_with("FAM_"),
        names_to = "GCF",
        values_to = "presence"
    ) %>%
    filter(presence == 1, !is.na(SampleID)) %>%
    distinct(SampleID, GCF) %>%
    mutate(presence = 1) %>%
    pivot_wider(
        names_from = GCF,
        values_from = presence,
        values_fill = 0
    )

meta_sample <- meta_genome %>%
    distinct(SampleID, Site, Type)

dat <- gcf_pa_sample %>%
    left_join(meta_sample, by = "SampleID")

cat("Sample x (GCF + metadata) dimensions:", paste(dim(dat), collapse = " x "), "\n")
cat("Samples per Type:\n"); print(table(dat$Type, useNA = "ifany"))
cat("Samples per Site:\n"); print(table(dat$Site, useNA = "ifany"))

# -----------------------------
# 7) Jaccard distance on presence/absence
# -----------------------------
gcf_mat <- dat %>% select(starts_with("FAM_")) %>% as.data.frame()
rownames(gcf_mat) <- dat$SampleID

dist_j <- vegdist(gcf_mat, method = "jaccard", binary = TRUE)

# -----------------------------
# 8) PERMANOVA (Type only)
# -----------------------------
set.seed(1)
adon_type <- adonis2(dist_j ~ Type, data = dat, permutations = 999)
print(adon_type)

perm_table <- as.data.frame(adon_type) %>%
    rownames_to_column("term") %>%
    filter(!term %in% c("Residual", "Total")) %>%
    transmute(
        model = "Type",
        term = term,
        Df = Df,
        SumOfSqs = SumOfSqs,
        R2 = R2,
        F = F,
        p = `Pr(>F)`
    )

write_csv(perm_table, "Table_PERMANOVA_GCF_Jaccard_TypeOnly.csv")
cat("Wrote: Table_PERMANOVA_GCF_Jaccard_TypeOnly.csv\n")

r2_type <- perm_table$R2[1]
p_type  <- perm_table$p[1]

# -----------------------------
# 9) betadisper for Type
# -----------------------------
bd_type <- betadisper(dist_j, dat$Type)
bd_perm <- permutest(bd_type, permutations = 999)

print(anova(bd_type))
print(bd_perm)

p_disp <- bd_perm$tab$`Pr(>F)`[1]

bd_tbl <- tibble(
    SampleID = names(bd_type$distances),
    Type = dat$Type[match(names(bd_type$distances), dat$SampleID)],
    dist_to_centroid = bd_type$distances
)

write_csv(bd_tbl, "GCF_betadisper_distances_to_centroid_TypeOnly.csv")
cat("Wrote: GCF_betadisper_distances_to_centroid_TypeOnly.csv\n")

cat("Median dist_to_centroid by Type:\n")
print(tapply(bd_tbl$dist_to_centroid, bd_tbl$Type, median))

# -----------------------------
# 10) PCoA + export coordinates
# -----------------------------
pcoa <- cmdscale(dist_j, k = 2, eig = TRUE)
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)[1:2]
cat("Variance explained PCoA1/2:", paste(var_exp, collapse = ", "), "\n")

ord <- as.data.frame(pcoa$points)
colnames(ord) <- c("PCoA1", "PCoA2")
ord$SampleID <- rownames(ord)

ord_plot <- ord %>%
    left_join(dat %>% select(SampleID, Type, Site), by = "SampleID")

write_csv(ord_plot, "GCF_PCoA_coordinates_TypeOnly.csv")
cat("Wrote: GCF_PCoA_coordinates_TypeOnly.csv\n")

sub_pcoa <- paste0(
    "PERMANOVA (Jaccard, Type): R2 = ",
    round(r2_type, 3),
    ", p = ",
    signif(p_type, 3)
)

sub_disp <- paste0(
    "betadisper permutation test: p = ",
    signif(p_disp, 3)
)

# -----------------------------
# 11) Plot PCoA (Type)
# -----------------------------
p_pcoa <- ggplot(ord_plot, aes(PCoA1, PCoA2, color = Type)) +
    geom_point(size = 3, alpha = 0.95) +
    ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = Inf) +
    labs(
        x = paste0("PCoA1 (", var_exp[1], "%)"),
        y = paste0("PCoA2 (", var_exp[2], "%)"),
        title = "GCF composition across habitat types",
        subtitle = sub_pcoa,
        color = "Habitat type"
    ) +
    coord_equal() +
    theme_classic(base_size = 13)

ggsave("Fig_GCF_PCoA_byType_TypeOnly.pdf", p_pcoa, width = 7, height = 5, device = cairo_pdf)
cat("Saved: Fig_GCF_PCoA_byType_TypeOnly.pdf\n")

# -----------------------------
# 12) Plot betadisper distances (Type)
# -----------------------------
p_disp_plot <- ggplot(bd_tbl, aes(Type, dist_to_centroid, fill = Type)) +
    geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.85) +
    labs(
        x = NULL,
        y = "Distance to centroid (Jaccard)",
        title = "Dispersion of GCF composition by habitat type",
        subtitle = sub_disp
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")

ggsave("Fig_GCF_Betadisper_byType_TypeOnly.pdf", p_disp_plot, width = 6, height = 4.5, device = cairo_pdf)
cat("Saved: Fig_GCF_Betadisper_byType_TypeOnly.pdf\n")

# -----------------------------
# 13) Save sample-level presence/absence and distance matrix
# -----------------------------
write_csv(dat, "GCF_sample_presence_absence_full_TypeOnly.csv")
cat("Wrote: GCF_sample_presence_absence_full_TypeOnly.csv\n")

write.table(
    as.matrix(dist_j),
    "GCF_sample_jaccard_distance_TypeOnly.tsv",
    sep = "\t",
    quote = FALSE
)
cat("Wrote: GCF_sample_jaccard_distance_TypeOnly.tsv\n")
