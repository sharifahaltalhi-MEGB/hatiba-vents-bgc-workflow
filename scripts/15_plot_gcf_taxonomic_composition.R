# ============================================
# Taxonomic composition of consistent GCF sets
# ============================================
#
# Purpose:
# Visualize phylum- and class-level taxonomic composition of consistent
# shared and habitat-specific gene cluster families (GCFs).
#
# Approach:
# - Read the merged table of BGCs assigned to consistent GCF sets
# - Collapse to distinct GCF-category-classification combinations
# - Extract phylum and class directly from GTDB classification strings
# - Summarize relative fractions across shared, mats-only, and precipitate-only GCFs
# - Plot top phyla and top classes across categories
#
# Outputs:
# - Fig_top10_phyla_across_GCF_categories.pdf
# - phylum_color_key_top10.csv
# - Fig_top10_classes_no_other_GCF_categories.pdf
# - class_color_key_top10_no_other.csv
# ============================================

library(tidyverse)
library(stringr)
library(forcats)

infile <- "BGCs_from_consistent_GCF_sets_with_all_annotations_v2.csv"

top_n_phylum <- 10
top_n_class  <- 10

impact10_class <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
    "#9C755F", "#EDC948", "#76B7B2", "#FF9DA7", "#2D3A4A"
)

phylum_cols <- c(
    Acidobacteriota  = "#3A0CA3",
    Chloroflexota    = "#009E73",
    Desulfobacterota = "#8B5CF6",
    Gemmatimonadota  = "#F4A261",
    Nitrospinota     = "#56B4E9",
    Nitrospirota     = "#0072B2",
    Planctomycetota  = "#CC79A7",
    Pseudomonadota   = "#6C757D",
    Thermoproteota   = "#D55E00",
    Zhuqueibacterota = "#8C1C13"
)

df <- read_csv(infile, show_col_types = FALSE)

gcf_df <- df %>%
    select(GCF, category, classification) %>%
    distinct() %>%
    mutate(
        phylum = str_match(classification, "p__([^;]+)")[, 2],
        class  = str_match(classification, "c__([^;]+)")[, 2]
    )

# ----------------------------
# Phylum composition plot
# ----------------------------
top_phyla <- gcf_df %>%
    filter(!is.na(phylum)) %>%
    count(phylum, name = "n_gcf") %>%
    arrange(desc(n_gcf)) %>%
    slice_head(n = top_n_phylum) %>%
    pull(phylum)

tbl_phylum <- gcf_df %>%
    mutate(phylum_plot = ifelse(phylum %in% top_phyla, phylum, "Other")) %>%
    count(category, phylum_plot, name = "n_gcf") %>%
    group_by(category) %>%
    mutate(frac = n_gcf / sum(n_gcf)) %>%
    ungroup()

phylum_order <- tbl_phylum %>%
    group_by(phylum_plot) %>%
    summarise(total = sum(n_gcf), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(phylum_plot)

tbl_phylum$phylum_plot <- factor(tbl_phylum$phylum_plot, levels = phylum_order) %>%
    fct_relevel("Other", after = Inf)

phy_levels <- levels(tbl_phylum$phylum_plot)

pal_phylum <- setNames(rep(NA_character_, length(phy_levels)), phy_levels)
pal_phylum[names(phylum_cols)] <- phylum_cols[names(phylum_cols)]
pal_phylum["Other"] <- "#F2F2F2"

p_phylum <- ggplot(tbl_phylum, aes(x = category, y = frac, fill = phylum_plot)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal_phylum) +
    labs(
        x = NULL,
        y = "Fraction of GCFs",
        fill = "Phylum",
        title = "Dominant phyla contributing to shared and habitat-specific GCFs"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

ggsave(
    "Fig_top10_phyla_across_GCF_categories.pdf",
    p_phylum,
    width = 9,
    height = 5,
    device = cairo_pdf
)

write_csv(
    tibble(phylum = names(pal_phylum), color = pal_phylum),
    "phylum_color_key_top10.csv"
)

# ----------------------------
# Class composition plot
# ----------------------------
top_classes <- gcf_df %>%
    filter(!is.na(class)) %>%
    count(class, name = "n_gcf") %>%
    arrange(desc(n_gcf)) %>%
    slice_head(n = top_n_class) %>%
    pull(class)

tbl_class_top10 <- gcf_df %>%
    filter(class %in% top_classes) %>%
    count(category, class, name = "n_gcf") %>%
    group_by(category) %>%
    mutate(frac = n_gcf / sum(n_gcf)) %>%
    ungroup() %>%
    rename(class_plot = class)

class_order <- tbl_class_top10 %>%
    group_by(class_plot) %>%
    summarise(total = sum(n_gcf), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(class_plot)

tbl_class_top10$class_plot <- factor(tbl_class_top10$class_plot, levels = class_order)

pal_class <- setNames(impact10_class[seq_along(class_order)], class_order)

p_class <- ggplot(
    tbl_class_top10,
    aes(x = category, y = frac, fill = class_plot)
) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal_class) +
    labs(
        x = NULL,
        y = "Fraction within top 10 classes",
        fill = "Class",
        title = "Top 10 classes contributing to shared and habitat-specific GCFs"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

ggsave(
    "Fig_top10_classes_no_other_GCF_categories.pdf",
    p_class,
    width = 9,
    height = 5,
    device = cairo_pdf
)

write_csv(
    tibble(class = names(pal_class), color = pal_class),
    "class_color_key_top10_no_other.csv"
)
