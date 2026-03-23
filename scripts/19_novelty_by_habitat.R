# ============================================
# Novelty by habitat category
# ============================================
#
# Uses stable input file (already stable GCFs)
# Novel GCF definition: q75_dist (75th percentile of Min_Distance within GCF) > 900
# Strict inclusion for comparisons: n_bgcs_total >= 2 AND n_genomes >= 2
# Outputs novelty fraction with exact binomial 95% CI for habitat category
# Adds stats:
# - Overall test: Chi-square when expected counts are adequate, else Fisher
# - Pairwise 2x2 Fisher tests with BH correction and odds ratios
# ============================================

library(tidyverse)
library(stringr)
library(purrr)

infile  <- "analysis_ready_novel_GCF.csv"
out_dir <- "stable_gcf_general_novelty_outputs_primary_only/habitat_category_only"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dist_cutoff <- 900
strict_n_bgcs   <- 2
strict_n_genome <- 2
min_groupsize <- 5
fisher_B <- 100000

top1 <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
}

exact_ci_low  <- function(x, n) binom.test(x, n)$conf.int[1]
exact_ci_high <- function(x, n) binom.test(x, n)$conf.int[2]

summarise_group <- function(df, group_var, out_name_prefix) {
    stopifnot(group_var %in% names(df))
    
    tab_all <- df %>%
        filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
        group_by(.data[[group_var]]) %>%
        summarise(
            n_gcfs = n(),
            n_gcfs_novel_primary = sum(novel_gcf),
            frac_gcfs_novel_primary = mean(novel_gcf),
            .groups = "drop"
        ) %>%
        mutate(
            ci_low  = map2_dbl(n_gcfs_novel_primary, n_gcfs, exact_ci_low),
            ci_high = map2_dbl(n_gcfs_novel_primary, n_gcfs, exact_ci_high)
        ) %>%
        rename(group = .data[[group_var]]) %>%
        arrange(desc(frac_gcfs_novel_primary), desc(n_gcfs_novel_primary), desc(n_gcfs))
    
    tab_min <- tab_all %>% filter(n_gcfs >= min_groupsize)
    
    write_csv(tab_all, file.path(out_dir, paste0(out_name_prefix, "_all_with_CI.csv")))
    write_csv(tab_min, file.path(out_dir, paste0(out_name_prefix, "_min", min_groupsize, "_with_CI.csv")))
    
    list(all = tab_all, min = tab_min)
}

safe_fisher <- function(tab, B = 100000) {
    tryCatch(
        fisher.test(tab),
        error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = B)
    )
}

run_group_stats <- function(df, group_var, label, out_dir, fisher_B = 100000) {
    stopifnot(group_var %in% names(df))
    
    dat <- df %>%
        filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
        mutate(group = as.character(.data[[group_var]])) %>%
        mutate(novel = as.integer(novel_gcf))
    
    tab <- table(dat$group, dat$novel)
    
    if (!("0" %in% colnames(tab))) tab <- cbind(tab, `0` = 0)
    if (!("1" %in% colnames(tab))) tab <- cbind(tab, `1` = 0)
    tab <- tab[, c("0","1"), drop = FALSE]
    
    ct_ok <- TRUE
    ct <- suppressWarnings(chisq.test(tab))
    if (any(is.na(ct$expected))) ct_ok <- FALSE
    if (ct_ok && any(ct$expected < 5)) ct_ok <- FALSE
    
    overall_test <- if (ct_ok) {
        tibble(
            test = "Chi-square",
            p_value = ct$p.value,
            method_note = "Used because expected counts >= 5"
        )
    } else {
        ft <- safe_fisher(tab, B = fisher_B)
        tibble(
            test = ifelse(!is.null(ft$method) && str_detect(ft$method, "simulation"),
                          "Fisher exact (simulated)",
                          "Fisher exact"),
            p_value = ft$p.value,
            method_note = "Used because some expected counts < 5 (or Fisher required)"
        )
    }
    
    group_summary <- dat %>%
        group_by(group) %>%
        summarise(
            n_gcfs = n(),
            n_novel = sum(novel),
            frac_novel = mean(novel),
            .groups = "drop"
        ) %>%
        arrange(desc(frac_novel), desc(n_novel), desc(n_gcfs))
    
    groups <- sort(unique(dat$group))
    if (length(groups) >= 2) {
        pairwise_list <- list()
        k <- 1
        
        for (i in 1:(length(groups)-1)) {
            for (j in (i+1):length(groups)) {
                g1 <- groups[i]
                g2 <- groups[j]
                
                sub <- dat %>% filter(group %in% c(g1, g2))
                tab2 <- table(sub$group, sub$novel)
                
                if (!("0" %in% colnames(tab2))) tab2 <- cbind(tab2, `0` = 0)
                if (!("1" %in% colnames(tab2))) tab2 <- cbind(tab2, `1` = 0)
                tab2 <- tab2[, c("0","1"), drop = FALSE]
                
                ft2 <- safe_fisher(tab2, B = fisher_B)
                or_val <- if (is.null(ft2$estimate)) NA_real_ else as.numeric(ft2$estimate)
                
                pairwise_list[[k]] <- tibble(
                    group1 = g1,
                    group2 = g2,
                    odds_ratio = or_val,
                    ci_low = ft2$conf.int[1],
                    ci_high = ft2$conf.int[2],
                    p_value = ft2$p.value
                )
                k <- k + 1
            }
        }
        
        pairwise_df <- bind_rows(pairwise_list) %>%
            mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
            arrange(p_adj_BH, p_value)
    } else {
        pairwise_df <- tibble()
    }
    
    write_csv(group_summary, file.path(out_dir, paste0("STATS_", label, "_group_summary.csv")))
    write_csv(overall_test,  file.path(out_dir, paste0("STATS_", label, "_overall_test.csv")))
    if (nrow(pairwise_df) > 0) {
        write_csv(pairwise_df, file.path(out_dir, paste0("STATS_", label, "_pairwise_fisher_OR.csv")))
    }
    
    invisible(list(tab = tab, overall = overall_test, summary = group_summary, pairwise = pairwise_df))
}

df <- read_csv(infile, show_col_types = FALSE)

req <- c("Min_Distance","GCF","user_genome","category")
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

df <- df %>%
    mutate(
        category = str_squish(as.character(category))
    )

gcf_level <- df %>%
    filter(!is.na(GCF), GCF != "", !is.na(Min_Distance)) %>%
    group_by(GCF) %>%
    summarise(
        n_bgcs_total = n(),
        n_genomes    = n_distinct(user_genome),
        q75_dist  = as.numeric(quantile(Min_Distance, 0.75, na.rm = TRUE)),
        novel_gcf = q75_dist > dist_cutoff,
        habitat_category = top1(category),
        .groups = "drop"
    )

write_csv(gcf_level, file.path(out_dir, "GCF_level_q75_novelty_and_habitat.csv"))

gcf_use <- gcf_level %>%
    filter(n_bgcs_total >= strict_n_bgcs, n_genomes >= strict_n_genome)

write_csv(gcf_use, file.path(out_dir, "GCF_level_q75_novelty_and_habitat_strict.csv"))

res_category <- summarise_group(gcf_use, "habitat_category", "Novelty_by_habitatCategory_q75_gt_900")

stats_category <- run_group_stats(gcf_use, "habitat_category", "habitat_category", out_dir, fisher_B = fisher_B)
