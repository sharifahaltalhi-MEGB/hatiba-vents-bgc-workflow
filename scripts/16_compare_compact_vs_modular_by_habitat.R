

suppressPackageStartupMessages({
    library(tidyverse)
    library(stringr)
})

infile <- "BGCs_from_consistent_GCF_sets_with_all_annotations_v2.csv"
df <- readr::read_csv(infile, show_col_types = FALSE)

stopifnot(all(c("BGC","GCF") %in% names(df)))

# Category (Big-SCAPE class) MUST exist
if (!("Category" %in% names(df))) stop("Column 'Category' not found. This must contain BiG-SCAPE class strings.")

# Rename GCF set category if present (Shared, Mats_only, Precip_only)
df2 <- df %>%
    rename(category_gcfset = any_of(c("category", "Category_gcfset"))) %>%
    mutate(
        class_raw = Category %>% str_trim(),
        class_raw_lower = str_to_lower(class_raw)
    )

# Convert Big-SCAPE class string to core types by token presence
to_core_class <- function(x_lower) {
    toks <- unlist(str_split(x_lower, "\\."))
    has_nrps <- "nrps" %in% toks
    has_pks  <- "pks"  %in% toks
    has_ripp <- "ripp" %in% toks
    has_terp <- "terpene" %in% toks
    
    # Multi-label presence flags
    # We will output ONE primary label for statistics, with a deterministic priority
    if (has_nrps) return("NRPS")
    if (has_pks)  return("PKS")
    if (has_ripp) return("RiPP")
    if (has_terp) return("terpene")
    return("other")
}

df2 <- df2 %>%
    mutate(class_core_bgc = map_chr(class_raw_lower, to_core_class))

# One class per GCF by majority vote on BGCs within the GCF
# Tie-break priority can be adjusted
priority <- c("NRPS","PKS","RiPP","terpene","other")

resolve_gcf_class <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_character_)
    tab <- table(v)
    max_n <- max(tab)
    winners <- names(tab)[tab == max_n]
    if (length(winners) == 1) return(winners)
    winners[order(match(winners, priority))][1]
}

gcf_level <- df2 %>%
    group_by(GCF) %>%
    summarise(
        n_bgcs = n(),
        category_gcfset = if ("category_gcfset" %in% names(df2)) {
            vals <- unique(na.omit(category_gcfset))
            if (length(vals) == 0) NA_character_
            else if (length(vals) == 1) vals
            else paste(vals, collapse = ";")
        } else NA_character_,
        class_gcf = resolve_gcf_class(class_core_bgc),
        has_conflict_category = if ("category_gcfset" %in% names(df2)) {
            length(unique(na.omit(category_gcfset))) > 1
        } else FALSE,
        has_conflict_class_bgc = length(unique(na.omit(class_core_bgc))) > 1,
        class_votes = {
            tab <- sort(table(class_core_bgc), decreasing = TRUE)
            paste0(names(tab), ":", as.integer(tab), collapse = ", ")
        },
        .groups = "drop"
    )

readr::write_csv(gcf_level, "gcf_level_table.csv")

gcf_class_summary <- gcf_level %>%
    filter(!is.na(class_gcf)) %>%
    count(class_gcf, sort = TRUE) %>%
    mutate(frac = n / sum(n))

readr::write_csv(gcf_class_summary, "gcf_class_summary.csv")



cat("\nGCF level rows:", nrow(gcf_level), "\n")
cat("GCFs with conflicting category assignment:", sum(gcf_level$has_conflict_category, na.rm = TRUE), "\n")
cat("GCFs where BGC-level core classes disagree within the GCF:", sum(gcf_level$has_conflict_class_bgc, na.rm = TRUE), "\n\n")
print(gcf_class_summary)

# Export the mixed-class GCFs for inspection
gcf_level %>%
    filter(has_conflict_class_bgc) %>%
    arrange(desc(n_bgcs)) %>%
    write_csv("GCFs_with_mixed_core_class_calls.csv")





# STEP 1: assign strategy per GCF
gcf_level2 <- gcf_level %>%
    mutate(
        strategy = case_when(
            class_gcf %in% c("RiPP", "terpene") ~ "compact",
            class_gcf %in% c("NRPS", "PKS") ~ "modular",
            TRUE ~ "other"
        )
    )

# Optional: remove "other" if you want a clean comparison
gcf_level_filtered <- gcf_level2 %>%
    filter(strategy %in% c("compact", "modular"))

# STEP 2: clean category labels
gcf_level_filtered <- gcf_level_filtered %>%
    mutate(
        category_clean = str_to_lower(str_trim(category_gcfset)),
        category_clean = case_when(
            category_clean %in% c("mats_only","mat_only") ~ "Mats_only",
            category_clean %in% c("precip_only","precipitate_only","precipitates_only") ~ "Precip_only",
            TRUE ~ NA_character_
        )
    ) %>%
    filter(category_clean %in% c("Mats_only","Precip_only"))

# STEP 3: build contingency table
table_strategy <- gcf_level_filtered %>%
    count(category_clean, strategy) %>%
    pivot_wider(names_from = strategy, values_from = n, values_fill = 0)

print(table_strategy)

# STEP 4: fisher test
mat_compact  <- table_strategy %>% filter(category_clean == "Mats_only") %>% pull(compact)
mat_modular  <- table_strategy %>% filter(category_clean == "Mats_only") %>% pull(modular)
prec_compact <- table_strategy %>% filter(category_clean == "Precip_only") %>% pull(compact)
prec_modular <- table_strategy %>% filter(category_clean == "Precip_only") %>% pull(modular)

fisher_res <- fisher.test(matrix(
    c(mat_compact, mat_modular,
      prec_compact, prec_modular),
    nrow = 2, byrow = TRUE
))

# STEP 5: output
result <- tibble(
    mat_compact = mat_compact,
    mat_modular = mat_modular,
    prec_compact = prec_compact,
    prec_modular = prec_modular,
    odds_ratio_mat_vs_precip = as.numeric(fisher_res$estimate),
    p_value = fisher_res$p.value
)

print(result)

readr::write_csv(result, "strategy_enrichment_gcf_level.csv")

