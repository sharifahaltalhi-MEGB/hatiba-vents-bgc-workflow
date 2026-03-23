# ============================================
# Archaea-only novelty by biosynthetic category
# ============================================
#
# Archaea only novelty by Big-SCAPE/combined Category.
# Novel GCF definition: q75(Min_Distance within GCF) > 900
# Strict GCF filter: n_bgcs_total >= 2 AND n_genomes >= 2
#
# Key choices:
# 1) Accept whatever is labeled d__Archaea in the classification table
# 2) Assign GCF Category using high vote (highest BGC count within the GCF)
# 3) Assign GCF taxonomy using modal paired phylum||class
# 4) Keep mixedness metrics so weak high-vote assignments remain visible
#
# Outputs folder:
#   archaea_only_novelty_by_CATEGORY_PAIRED_TAXONOMY_V8_HIGHVOTE_full_datasets
# ============================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(stringr)
    library(purrr)
})

infile  <- "bgc_level_antismash_v8_with_bigscapeV2_and_novelty_WITH_FAMILY.csv"
out_dir <- "archaea_only_novelty_by_CATEGORY_PAIRED_TAXONOMY_V8_HIGHVOTE_full_datasets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dist_cutoff <- 900
strict_n_bgcs   <- 2
strict_n_genome <- 2

exact_ci_low  <- function(x, n) binom.test(x, n)$conf.int[1]
exact_ci_high <- function(x, n) binom.test(x, n)$conf.int[2]

get_rank <- function(x, rank_prefix) {
    x <- as.character(x)
    pat <- paste0(rank_prefix, "__[^;]+")
    out <- str_extract(x, pat)
    out <- str_remove(out, paste0("^", rank_prefix, "__"))
    out
}

clean_genome_id <- function(x) {
    x <- as.character(x)
    x <- str_squish(x)
    x <- str_remove(x, "\\.fa(sta)?$|\\.fna$|\\.gz$|\\.bz2$|\\.xz$")
    x
}

split_pair <- function(x) {
    m <- str_split_fixed(as.character(x), fixed("||"), 2)
    tibble(phylum = m[, 1], class = m[, 2])
}

normalize_category <- function(x) {
    x <- str_squish(as.character(x))
    x[x == ""] <- NA_character_
    map_chr(x, function(z) {
        if (is.na(z)) return(NA_character_)
        parts <- unlist(str_split(z, "\\."))
        parts <- parts[!is.na(parts) & parts != ""]
        if (length(parts) == 0) return(NA_character_)
        parts <- sort(unique(parts))
        paste(parts, collapse = ".")
    })
}

high_vote_label <- function(labels) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels) & labels != ""]
    if (length(labels) == 0) return(NA_character_)
    tab <- table(labels)
    mx <- max(tab)
    winners <- sort(names(tab)[tab == mx])
    winners[1]
}

high_vote_frac <- function(labels) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels) & labels != ""]
    if (length(labels) == 0) return(NA_real_)
    tab <- table(labels)
    as.numeric(max(tab) / sum(tab))
}

df <- read_csv(infile, show_col_types = FALSE)

if (!("Category" %in% names(df))) {
    stop("Missing column: Category (your file must have a Category column).")
}

req <- c("GCF", "Min_Distance", "user_genome", "classification", "Category")
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

df0 <- df %>%
    mutate(
        user_genome_raw = user_genome,
        user_genome = clean_genome_id(user_genome),
        classification = str_squish(as.character(classification)),
        domain = str_squish(get_rank(classification, "d")),
        phylum_raw = str_squish(get_rank(classification, "p")),
        class_raw  = str_squish(get_rank(classification, "c")),
        Category_raw = str_squish(as.character(Category)),
        Category = normalize_category(Category)
    )

df0 %>%
    distinct(Category_raw, Category) %>%
    arrange(Category, Category_raw) %>%
    write_csv(file.path(out_dir, "QC_Category_normalization_mapping.csv"))

df0 %>%
    filter(!is.na(Category), Category != "") %>%
    count(Category, sort = TRUE, name = "n_bgcs_all_domains") %>%
    write_csv(file.path(out_dir, "QC_Category_counts_all_domains_BGC_level.csv"))

archaea_bgc <- df0 %>%
    filter(!is.na(domain), domain == "Archaea") %>%
    filter(!is.na(user_genome), user_genome != "") %>%
    filter(!is.na(classification), classification != "")

archaea_bgc %>%
    filter(!is.na(Category), Category != "") %>%
    count(Category, sort = TRUE, name = "n_bgcs_archaea") %>%
    write_csv(file.path(out_dir, "QC_Category_counts_Archaea_BGC_level.csv"))

genome_tax_conflicts <- archaea_bgc %>%
    distinct(user_genome, classification) %>%
    group_by(user_genome) %>%
    summarise(
        n_classifications = n(),
        classifications = paste(classification, collapse = " | "),
        .groups = "drop"
    ) %>%
    filter(n_classifications > 1) %>%
    arrange(desc(n_classifications))

write_csv(genome_tax_conflicts, file.path(out_dir, "QC_genome_taxonomy_conflicts.csv"))

genome_tax <- archaea_bgc %>%
    group_by(user_genome) %>%
    summarise(
        classification = high_vote_label(classification),
        .groups = "drop"
    ) %>%
    mutate(
        domain = str_squish(get_rank(classification, "d")),
        phylum = str_squish(get_rank(classification, "p")),
        class  = str_squish(get_rank(classification, "c"))
    )

write_csv(genome_tax, file.path(out_dir, "Genome_taxonomy_lookup.csv"))

archaea_bgc_fixed <- archaea_bgc %>%
    select(-classification) %>%
    left_join(
        genome_tax %>% select(user_genome, classification, domain, phylum, class),
        by = "user_genome"
    )

lost_tax <- archaea_bgc_fixed %>%
    filter(is.na(classification) | classification == "") %>%
    distinct(user_genome, user_genome_raw) %>%
    arrange(user_genome)

write_csv(lost_tax, file.path(out_dir, "QC_missing_taxonomy_after_fix.csv"))

archaea_bgc_fixed %>%
    filter(!is.na(phylum), phylum != "") %>%
    count(phylum, classification, sort = TRUE, name = "n_bgcs") %>%
    group_by(phylum) %>%
    slice_head(n = 25) %>%
    ungroup() %>%
    write_csv(file.path(out_dir, "EVIDENCE_phylum_classification_examples_top25.csv"))

archaea_bgc_fixed %>%
    filter(!is.na(class), class != "") %>%
    count(class, classification, sort = TRUE, name = "n_bgcs") %>%
    group_by(class) %>%
    slice_head(n = 25) %>%
    ungroup() %>%
    write_csv(file.path(out_dir, "EVIDENCE_class_classification_examples_top25.csv"))

gcf_tax_counts <- archaea_bgc_fixed %>%
    filter(!is.na(GCF), GCF != "") %>%
    filter(!is.na(phylum), phylum != "", !is.na(class), class != "") %>%
    mutate(phylum_class = paste(phylum, class, sep = "||")) %>%
    count(GCF, phylum_class, name = "n_bgcs") %>%
    arrange(GCF, desc(n_bgcs))

write_csv(gcf_tax_counts, file.path(out_dir, "QC_GCF_taxonomy_composition_counts.csv"))

gcf_tax_summary <- gcf_tax_counts %>%
    group_by(GCF) %>%
    mutate(
        gcf_total_bgcs = sum(n_bgcs),
        frac = n_bgcs / gcf_total_bgcs
    ) %>%
    summarise(
        gcf_total_bgcs = max(gcf_total_bgcs),
        n_phylum_class = n_distinct(phylum_class),
        phylum_class_mode = phylum_class[which.max(n_bgcs)][1],
        tax_top_frac = max(frac),
        shannon = -sum(frac * log(frac)),
        .groups = "drop"
    ) %>%
    arrange(desc(n_phylum_class), desc(gcf_total_bgcs), desc(tax_top_frac))

write_csv(gcf_tax_summary, file.path(out_dir, "QC_GCF_taxonomy_composition_summary.csv"))

archaea_bgc_fixed %>%
    semi_join(gcf_tax_summary %>% filter(n_phylum_class > 1) %>% select(GCF), by = "GCF") %>%
    filter(!is.na(phylum), phylum != "", !is.na(class), class != "") %>%
    mutate(phylum_class = paste(phylum, class, sep = "||")) %>%
    select(
        GCF, phylum, class, phylum_class,
        user_genome, user_genome_raw,
        Category, Min_Distance,
        classification
    ) %>%
    arrange(GCF, phylum_class, user_genome) %>%
    write_csv(file.path(out_dir, "EVIDENCE_rows_for_mixed_GCFs_taxonomy.csv"))

gcf_cat_counts <- archaea_bgc_fixed %>%
    filter(!is.na(GCF), GCF != "") %>%
    filter(!is.na(Category), Category != "") %>%
    count(GCF, Category, name = "n_bgcs") %>%
    group_by(GCF) %>%
    mutate(
        gcf_total_bgcs = sum(n_bgcs),
        frac = n_bgcs / gcf_total_bgcs
    ) %>%
    ungroup() %>%
    arrange(GCF, desc(n_bgcs))

write_csv(gcf_cat_counts, file.path(out_dir, "QC_GCF_Category_composition_counts.csv"))

gcf_cat_summary <- gcf_cat_counts %>%
    group_by(GCF) %>%
    summarise(
        n_categories = n_distinct(Category),
        category_mode = Category[which.max(n_bgcs)][1],
        category_top_frac = max(frac),
        .groups = "drop"
    ) %>%
    arrange(desc(n_categories), desc(category_top_frac))

write_csv(gcf_cat_summary, file.path(out_dir, "QC_GCF_Category_composition_summary.csv"))

archaea_bgc_fixed %>%
    semi_join(gcf_cat_summary %>% filter(n_categories > 1) %>% select(GCF), by = "GCF") %>%
    select(
        GCF, Category,
        user_genome, user_genome_raw,
        phylum, class, Min_Distance,
        classification
    ) %>%
    arrange(GCF, Category, user_genome) %>%
    write_csv(file.path(out_dir, "EVIDENCE_rows_for_mixed_GCFs_category.csv"))

gcf_archaea_full <- archaea_bgc_fixed %>%
    filter(!is.na(GCF), GCF != "", !is.na(Min_Distance)) %>%
    mutate(
        phylum = str_squish(phylum),
        class  = str_squish(class),
        phylum_class = if_else(
            !is.na(phylum) & phylum != "" & !is.na(class) & class != "",
            paste(phylum, class, sep = "||"),
            NA_character_
        )
    ) %>%
    group_by(GCF) %>%
    summarise(
        n_bgcs_total = n(),
        n_genomes    = n_distinct(user_genome),
        q75_dist     = as.numeric(quantile(Min_Distance, 0.75, na.rm = TRUE)),
        novel_gcf    = q75_dist > dist_cutoff,
        n_phylum_class = n_distinct(phylum_class),
        phylum_class_mode = high_vote_label(phylum_class),
        tax_top_frac = high_vote_frac(phylum_class),
        n_categories = n_distinct(Category),
        category_mode = high_vote_label(Category),
        category_top_frac = high_vote_frac(Category),
        .groups = "drop"
    ) %>%
    mutate(
        phylum = split_pair(phylum_class_mode)$phylum,
        class  = split_pair(phylum_class_mode)$class,
        pass_strict = (n_bgcs_total >= strict_n_bgcs) & (n_genomes >= strict_n_genome),
        fail_reason = case_when(
            pass_strict ~ "PASS",
            n_bgcs_total < strict_n_bgcs & n_genomes < strict_n_genome ~ "FAIL: n_bgcs<2 AND n_genomes<2",
            n_bgcs_total < strict_n_bgcs ~ "FAIL: n_bgcs<2",
            n_genomes < strict_n_genome ~ "FAIL: n_genomes<2",
            TRUE ~ "FAIL: other"
        )
    ) %>%
    select(
        GCF, n_bgcs_total, n_genomes, q75_dist, novel_gcf,
        phylum, class, n_phylum_class, tax_top_frac,
        category_mode, n_categories, category_top_frac,
        pass_strict, fail_reason
    ) %>%
    arrange(desc(pass_strict), desc(novel_gcf), desc(q75_dist), desc(n_bgcs_total), desc(n_genomes))

write_csv(gcf_archaea_full, file.path(out_dir, "Archaea_GCF_q75_novelty_FULL.csv"))

gcf_archaea_strict <- gcf_archaea_full %>%
    filter(pass_strict)

write_csv(gcf_archaea_strict, file.path(out_dir, "Archaea_GCF_q75_novelty_STRICT.csv"))

novelty_by_category <- gcf_archaea_strict %>%
    filter(!is.na(category_mode), category_mode != "") %>%
    group_by(category_mode) %>%
    summarise(
        n_gcfs = n(),
        n_novel = sum(novel_gcf),
        frac_novel = mean(novel_gcf),
        ci_low  = map2_dbl(n_novel, n_gcfs, exact_ci_low),
        ci_high = map2_dbl(n_novel, n_gcfs, exact_ci_high),
        .groups = "drop"
    ) %>%
    arrange(desc(frac_novel), desc(n_novel), desc(n_gcfs))

write_csv(novelty_by_category, file.path(out_dir, "Archaea_novelty_by_Category_STRICT.csv"))

gcf_archaea_strict %>%
    count(category_mode, sort = TRUE, name = "n_gcfs_strict") %>%
    write_csv(file.path(out_dir, "QC_Category_counts_STRICT_GCF_level.csv"))
