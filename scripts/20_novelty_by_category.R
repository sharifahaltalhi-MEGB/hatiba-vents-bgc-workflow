# ============================================
# Novelty by biosynthetic category
# ============================================
#
# Novelty by Category (BigScape Category)
# Input: analysis_ready_novel_GCF.csv
# Novel GCF definition: q75_dist (75th percentile of Min_Distance within GCF) > 900
# Inclusion criteria (strict): n_bgcs_total >= 2 AND n_genomes >= 2
# Outputs:
# - outputs/GCF_level_q75_novelty_and_assignments.csv
# - outputs/GCF_level_q75_novelty_and_assignments_strict.csv
# - outputs/Novelty_by_Category_q75_gt_900_all_with_CI.csv
# - outputs/Novelty_by_Category_q75_gt_900_min5_with_CI.csv
# - outputs/STATS_category_overall_glm_LR.csv
# - outputs/STATS_category_pairwise_RISKDIFF_BH_all.csv
# - outputs/STATS_category_pairwise_RISKDIFF_BH_min5.csv
# - plots/Novelty_fraction_by_Category_q75_gt_900.png
# ============================================

library(tidyverse)
library(stringr)
library(purrr)
library(forcats)

infile <- "analysis_ready_novel_GCF.csv"

out_dir  <- "outputs"
plot_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

dist_cutoff <- 900
do_strict <- TRUE
strict_n_bgcs   <- 2
strict_n_genome <- 2
min_groupsize <- 5

top1 <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
}

exact_ci_low  <- function(x, n) binom.test(x, n)$conf.int[1]
exact_ci_high <- function(x, n) binom.test(x, n)$conf.int[2]

summarise_group_with_ci <- function(df, group_var, out_prefix) {
    stopifnot(group_var %in% names(df))
    
    tab_all <- df %>%
        filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
        group_by(group = .data[[group_var]]) %>%
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
        arrange(desc(frac_gcfs_novel_primary), desc(n_gcfs_novel_primary), desc(n_gcfs))
    
    tab_min <- tab_all %>% filter(n_gcfs >= min_groupsize)
    
    write_csv(tab_all, file.path(out_dir, paste0(out_prefix, "_all_with_CI.csv")))
    write_csv(tab_min, file.path(out_dir, paste0(out_prefix, "_min", min_groupsize, "_with_CI.csv")))
    
    list(all = tab_all, min = tab_min)
}

fit_overall_glm_lr <- function(gcf_df, group_var) {
    dat <- gcf_df %>%
        filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
        mutate(
            group = factor(.data[[group_var]]),
            novel = as.integer(novel_gcf)
        )
    
    m <- glm(novel ~ group, data = dat, family = binomial())
    lr <- anova(m, test = "Chisq")
    
    as_tibble(lr, rownames = "term") %>%
        transmute(
            term,
            df = Df,
            deviance = Deviance,
            p_value = `Pr(>Chi)`
        )
}

pairwise_riskdiff <- function(gcf_df, group_var) {
    dat <- gcf_df %>%
        filter(!is.na(.data[[group_var]]), .data[[group_var]] != "") %>%
        transmute(
            group = as.character(.data[[group_var]]),
            novel = as.integer(novel_gcf)
        )
    
    groups <- sort(unique(dat$group))
    if (length(groups) < 2) return(tibble())
    
    out <- vector("list", length = choose(length(groups), 2))
    k <- 1
    
    for (i in 1:(length(groups) - 1)) {
        for (j in (i + 1):length(groups)) {
            g1 <- groups[i]
            g2 <- groups[j]
            
            s1 <- dat %>% filter(group == g1)
            s2 <- dat %>% filter(group == g2)
            
            x1 <- sum(s1$novel); n1 <- nrow(s1)
            x2 <- sum(s2$novel); n2 <- nrow(s2)
            
            p1 <- x1 / n1
            p2 <- x2 / n2
            rd <- p1 - p2
            
            pt <- suppressWarnings(prop.test(x = c(x1, x2), n = c(n1, n2), correct = FALSE))
            ci <- unname(pt$conf.int)
            
            out[[k]] <- tibble(
                group1 = g1,
                group2 = g2,
                novel1 = x1,
                n1 = n1,
                frac1 = p1,
                novel2 = x2,
                n2 = n2,
                frac2 = p2,
                risk_diff = rd,
                ci_low = ci[1],
                ci_high = ci[2],
                p_value = pt$p.value
            )
            k <- k + 1
        }
    }
    
    bind_rows(out[seq_len(k - 1)]) %>%
        mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
        arrange(p_adj_BH, p_value)
}

df <- read_csv(infile, show_col_types = FALSE)

category_col <- if ("Category" %in% names(df)) {
    "Category"
} else if ("category" %in% names(df)) {
    "category"
} else {
    NA_character_
}
if (is.na(category_col)) stop("Could not find Category column")

req <- c("Min_Distance", "GCF", "user_genome", category_col)
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

df <- df %>%
    mutate(category = str_squish(as.character(.data[[category_col]])))

gcf_level <- df %>%
    filter(!is.na(GCF), GCF != "", !is.na(Min_Distance)) %>%
    group_by(GCF) %>%
    summarise(
        n_bgcs_total = n(),
        n_genomes = n_distinct(user_genome),
        q75_dist = as.numeric(quantile(Min_Distance, 0.75, na.rm = TRUE)),
        novel_gcf = q75_dist > dist_cutoff,
        category_gcf = top1(category),
        .groups = "drop"
    )

write_csv(gcf_level, file.path(out_dir, "GCF_level_q75_novelty_and_assignments.csv"))

gcf_use <- gcf_level
if (isTRUE(do_strict)) {
    gcf_use <- gcf_use %>%
        filter(n_bgcs_total >= strict_n_bgcs, n_genomes >= strict_n_genome)
}
write_csv(gcf_use, file.path(out_dir, "GCF_level_q75_novelty_and_assignments_strict.csv"))

res_category <- summarise_group_with_ci(
    gcf_use,
    "category_gcf",
    "Novelty_by_Category_q75_gt_900"
)

glm_lr <- fit_overall_glm_lr(gcf_use, "category_gcf")
write_csv(glm_lr, file.path(out_dir, "STATS_category_overall_glm_LR.csv"))

pairwise_all <- pairwise_riskdiff(gcf_use, "category_gcf")
write_csv(pairwise_all, file.path(out_dir, "STATS_category_pairwise_RISKDIFF_BH_all.csv"))

valid_groups <- res_category$all %>%
    filter(n_gcfs >= min_groupsize) %>%
    pull(group)

gcf_use_min <- gcf_use %>%
    filter(category_gcf %in% valid_groups)

pairwise_min <- pairwise_riskdiff(gcf_use_min, "category_gcf")
write_csv(pairwise_min, file.path(out_dir, "STATS_category_pairwise_RISKDIFF_BH_min5.csv"))

plot_df <- res_category$all %>%
    mutate(group = fct_reorder(group, frac_gcfs_novel_primary))

p <- ggplot(plot_df, aes(x = group, y = frac_gcfs_novel_primary)) +
    geom_col() +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
    coord_flip() +
    labs(
        x = "Category",
        y = "Fraction novel GCFs (q75 Min_Distance > 900)"
    )

ggsave(
    filename = file.path(plot_dir, "Novelty_fraction_by_Category_q75_gt_900.png"),
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
)
