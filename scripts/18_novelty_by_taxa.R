# ============================================
# Novelty by taxonomic group
# ============================================
#
# Purpose:
# Quantify the fraction of novel gene cluster families (GCFs)
# across paired phylum-class taxonomic groups.
#
# Approach:
# - Read the merged BGC-level table with BiG-SLiCE distances
# - Collapse novelty to the GCF level using the 75th percentile
#   of Min_Distance across member BGCs
# - Assign each GCF to a dominant paired phylum-class label
# - Restrict analysis to GCFs with at least 2 BGCs from at least
#   2 distinct genomes
# - Calculate novelty fractions and exact binomial 95% confidence intervals
#
# Outputs:
# - novelty_by_taxa/GCF_level_q75_novelty_and_taxonomy.csv
# - novelty_by_taxa/GCF_level_q75_novelty_and_taxonomy_strict.csv
# - novelty_by_taxa/Novelty_by_taxa_all_with_CI.csv
# - novelty_by_taxa/Novelty_by_taxa_min5_with_CI.csv
# ============================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(stringr)
    library(purrr)
})

infile <- "analysis_ready_novel_GCF.csv"
out_dir <- "novelty_by_taxa"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dist_cutoff <- 900
strict_n_bgcs <- 2
strict_n_genome <- 2
min_gcfs_per_taxon <- 5

get_mode_pair <- function(phylum_vec, class_vec) {
    phylum_vec <- as.character(phylum_vec)
    class_vec  <- as.character(class_vec)
    
    ok <- !is.na(phylum_vec) & phylum_vec != "" & !is.na(class_vec) & class_vec != ""
    phylum_vec <- phylum_vec[ok]
    class_vec  <- class_vec[ok]
    
    if (length(phylum_vec) == 0) {
        return(list(top_phylum = NA_character_, top_class = NA_character_))
    }
    
    pair_vec <- paste(phylum_vec, class_vec, sep = "||")
    tab <- sort(table(pair_vec), decreasing = TRUE)
    top_pair <- names(tab)[1]
    split_pair <- str_split_fixed(top_pair, "\\|\\|", 2)
    
    list(
        top_phylum = split_pair[, 1],
        top_class  = split_pair[, 2]
    )
}

analysis_df <- read_csv(infile, show_col_types = FALSE)

req_in <- c("Min_Distance", "GCF", "user_genome", "classification")
missing_in <- setdiff(req_in, names(analysis_df))
if (length(missing_in) > 0) {
    stop("Missing required columns: ", paste(missing_in, collapse = ", "))
}

analysis_df <- analysis_df %>%
    mutate(
        classification = as.character(classification),
        phylum = str_match(classification, "p__([^;]+)")[, 2],
        class  = str_match(classification, "c__([^;]+)")[, 2],
        phylum = if_else(is.na(phylum) | phylum == "", "unclassified", phylum),
        class  = if_else(is.na(class)  | class  == "", "unclassified", class),
        phylum = str_squish(phylum),
        class  = str_squish(class),
        Min_Distance = suppressWarnings(as.numeric(Min_Distance)),
        GCF = as.character(GCF),
        user_genome = as.character(user_genome)
    )

gcf_level <- analysis_df %>%
    filter(!is.na(GCF), GCF != "") %>%
    group_by(GCF) %>%
    summarise(
        n_bgcs_total = n(),
        n_genomes = n_distinct(user_genome),
        q75_dist = as.numeric(quantile(Min_Distance, 0.75, na.rm = TRUE)),
        novel_gcf = q75_dist > dist_cutoff,
        tmp = list(get_mode_pair(phylum, class)),
        top_phylum = tmp[[1]]$top_phylum,
        top_class = tmp[[1]]$top_class,
        .groups = "drop"
    ) %>%
    select(-tmp)

write_csv(
    gcf_level,
    file.path(out_dir, "GCF_level_q75_novelty_and_taxonomy.csv")
)

gcf_for_taxa <- gcf_level %>%
    filter(n_bgcs_total >= strict_n_bgcs, n_genomes >= strict_n_genome)

write_csv(
    gcf_for_taxa,
    file.path(out_dir, "GCF_level_q75_novelty_and_taxonomy_strict.csv")
)

novelty_by_taxa <- gcf_for_taxa %>%
    filter(
        !is.na(top_phylum), top_phylum != "",
        !is.na(top_class), top_class != ""
    ) %>%
    group_by(top_phylum, top_class) %>%
    summarise(
        n_gcfs = n(),
        n_gcfs_novel = sum(novel_gcf, na.rm = TRUE),
        frac_gcfs_novel = mean(novel_gcf, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        ci_low = map2_dbl(n_gcfs_novel, n_gcfs, ~ binom.test(.x, .y)$conf.int[1]),
        ci_high = map2_dbl(n_gcfs_novel, n_gcfs, ~ binom.test(.x, .y)$conf.int[2])
    ) %>%
    arrange(desc(frac_gcfs_novel), desc(n_gcfs_novel), desc(n_gcfs))

write_csv(
    novelty_by_taxa,
    file.path(out_dir, "Novelty_by_taxa_all_with_CI.csv")
)

novelty_by_taxa_min5 <- novelty_by_taxa %>%
    filter(n_gcfs >= min_gcfs_per_taxon)

write_csv(
    novelty_by_taxa_min5,
    file.path(out_dir, "Novelty_by_taxa_min5_with_CI.csv")
)
