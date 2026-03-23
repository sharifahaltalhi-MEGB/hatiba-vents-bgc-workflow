# ============================================
# Lineage enrichment of major biosynthetic categories
# ============================================
#
# Purpose:
# Identify gene cluster families (GCFs) that are significantly enriched
# within dominant phyla and summarize enrichment across major
# biosynthetic categories.
#
# Approach:
# - Read the merged BGC-level table with GCF assignments
# - Assign taxonomy from GTDB classification strings
# - Define the dominant lineage for each GCF
# - Test enrichment of each GCF within its dominant lineage using Fisher's exact test
# - Adjust p values using Benjamini-Hochberg correction
# - Retain significant GCFs from major categories only (NRPS, PKS, RiPP, terpene)
#
# Outputs:
# - GCF_lineage_enrichment_all.csv
# - GCF_lineage_enrichment_significant_major_categories_only.csv
# - GCF_lineage_enrichment_summary_major_categories_only.csv
# - GCF_enrichment_by_major_category_and_phylum.png
# - GCF_enrichment_by_major_category_and_phylum.pdf
# ============================================

library(tidyverse)
library(stringr)
library(purrr)
library(forcats)

infile <- "analysis_ready_novel_GCF.csv"

extract_rank <- function(x, rank = "p__") {
    str_match(x, paste0("(", rank, "[^;]+)"))[, 2] %>%
        str_remove(rank) %>%
        replace_na("Unassigned")
}

fisher_one <- function(gcf_id, lin, gcf_genome, genome_lineage, all_genomes) {
    g_in  <- gcf_genome %>% filter(GCF == gcf_id) %>% pull(user_genome) %>% unique()
    g_lin <- genome_lineage %>% filter(lineage == lin) %>% pull(user_genome) %>% unique()
    
    a <- length(intersect(g_in, g_lin))
    b <- length(setdiff(g_in, g_lin))
    c <- length(setdiff(g_lin, g_in))
    d <- length(setdiff(all_genomes, union(g_in, g_lin)))
    
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    
    tibble(
        a = a,
        total_lineage_genomes = a + c,
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value
    )
}

tax_level <- "phylum"
min_total_gcf <- 5
min_a <- 4
fdr_cutoff <- 0.05

type_colors <- c(
    NRPS = "#e41a1c",
    PKS = "#377eb8",
    RiPP = "#984ea3",
    terpene = "#a6cee3"
)

bgc <- read_csv(infile, show_col_types = FALSE)

bgc2 <- bgc %>%
    mutate(
        classification = as.character(classification),
        GCF = str_trim(as.character(GCF)),
        user_genome = str_trim(as.character(user_genome)),
        phylum = extract_rank(classification, "p__"),
        class  = extract_rank(classification, "c__"),
        order  = extract_rank(classification, "o__"),
        genus  = extract_rank(classification, "g__"),
        Category = str_trim(as.character(Category))
    ) %>%
    filter(!is.na(GCF), GCF != "", !is.na(user_genome), user_genome != "")

gcf_genome <- bgc2 %>%
    distinct(GCF, user_genome, lineage = .data[[tax_level]])

genome_lineage <- gcf_genome %>%
    distinct(user_genome, lineage)

all_genomes <- genome_lineage$user_genome

gcf_by_tax <- gcf_genome %>%
    count(GCF, lineage, name = "n_genomes") %>%
    group_by(GCF) %>%
    mutate(total_gcf_genomes = sum(n_genomes)) %>%
    ungroup()

dominant_per_gcf <- gcf_by_tax %>%
    group_by(GCF) %>%
    slice_max(order_by = n_genomes, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
        GCF,
        dominant_lineage = lineage,
        dominant_n = n_genomes,
        total_gcf_genomes = total_gcf_genomes,
        dominant_frac = dominant_n / total_gcf_genomes
    )

enrich <- dominant_per_gcf %>%
    mutate(
        res = map2(
            GCF, dominant_lineage,
            ~ fisher_one(.x, .y, gcf_genome, genome_lineage, all_genomes)
        )
    ) %>%
    unnest(res) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj, desc(total_gcf_genomes))

gcf_category <- bgc2 %>%
    mutate(
        Category = if_else(is.na(Category) | Category == "", "other", Category)
    ) %>%
    count(GCF, Category, name = "n") %>%
    group_by(GCF) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(GCF, Category) %>%
    mutate(
        category_simple = case_when(
            Category == "NRPS" ~ "NRPS",
            Category == "PKS" ~ "PKS",
            Category == "RiPP" ~ "RiPP",
            Category == "terpene" ~ "terpene",
            TRUE ~ NA_character_
        )
    )

enrich2 <- enrich %>%
    left_join(gcf_category, by = "GCF")

sig_gcf_clean <- enrich2 %>%
    filter(
        total_gcf_genomes >= min_total_gcf,
        a >= min_a,
        p_adj < fdr_cutoff,
        !is.na(category_simple)
    ) %>%
    mutate(
        category_simple = factor(category_simple, levels = c("NRPS", "PKS", "RiPP", "terpene"))
    )

gcf_summary <- sig_gcf_clean %>%
    count(category_simple, dominant_lineage, name = "n_GCFs") %>%
    ungroup()

phylum_order <- gcf_summary %>%
    group_by(dominant_lineage) %>%
    summarise(total = sum(n_GCFs), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(dominant_lineage)

gcf_summary <- gcf_summary %>%
    mutate(
        dominant_lineage = factor(dominant_lineage, levels = phylum_order),
        category_simple = factor(category_simple, levels = c("NRPS", "PKS", "RiPP", "terpene"))
    )

write_csv(enrich, "GCF_lineage_enrichment_all.csv")
write_csv(sig_gcf_clean, "GCF_lineage_enrichment_significant_major_categories_only.csv")
write_csv(gcf_summary, "GCF_lineage_enrichment_summary_major_categories_only.csv")

p_bar <- ggplot(
    gcf_summary,
    aes(x = dominant_lineage, y = n_GCFs, fill = category_simple)
) +
    geom_col(width = 0.95) +
    scale_fill_manual(values = type_colors, drop = FALSE) +
    labs(
        x = "Dominant phylum",
        y = "Number of lineage enriched GCFs",
        fill = "Category"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )

ggsave(
    "GCF_enrichment_by_major_category_and_phylum.png",
    p_bar,
    width = 10,
    height = 7,
    dpi = 300
)

ggsave(
    "GCF_enrichment_by_major_category_and_phylum.pdf",
    p_bar,
    width = 10,
    height = 7
)
