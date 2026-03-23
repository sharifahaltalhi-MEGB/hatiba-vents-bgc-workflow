# ============================================
# Balanced subsampling of GCF presence to quantify
# shared and habitat-specific gene cluster families
# ============================================
#
# Purpose:
# Estimate the number of shared, mat-specific, and precipitate-specific
# gene cluster families (GCFs) while correcting for unequal sample sizes
# between habitats.
#
# Approach:
# - Subsample microbial mat and precipitate samples to equal size
# - Repeat across 1000 iterations
# - For each iteration:
#     * Identify GCFs present in each habitat
#     * Count shared, mat-only, and precipitate-only GCFs
# - Summarize results as mean ± standard deviation across iterations
#
# Outputs:
# - GCF_shared_unique_subsample_summary.csv
#     Mean and standard deviation of shared and habitat-specific GCF counts
# - GCF_shared_unique_subsample_iterations.csv
#     Per-iteration counts for all subsampling runs
#
# Notes:
# - This step supports the Results section reporting mean ± SD values
# - Stable GCF sets (frequency ≥ 0.8) are defined in a separate script
# ============================================

library(tidyverse)

set.seed(123)

df <- read_csv("GCF_sample_presence_absence_full_TypeOnly.csv", show_col_types = FALSE)
gcf_cols <- names(df)[str_detect(names(df), "^FAM_")]

n_mats <- sum(df$Type == "microbial mat")
n_prec <- sum(df$Type == "precipitate")
n_min  <- min(n_mats, n_prec)

n_iter <- 1000

res_subsample <- map_dfr(seq_len(n_iter), function(i) {
    
    mats_sub <- df %>%
        filter(Type == "microbial mat") %>%
        slice_sample(n = n_min, replace = FALSE)
    
    prec_sub <- df %>%
        filter(Type == "precipitate") %>%
        slice_sample(n = n_min, replace = FALSE)
    
    gcf_mats <- mats_sub %>%
        select(all_of(gcf_cols)) %>%
        summarise(across(everything(), ~ any(. == 1))) %>%
        pivot_longer(everything(), names_to = "GCF", values_to = "present") %>%
        filter(present) %>%
        pull(GCF)
    
    gcf_prec <- prec_sub %>%
        select(all_of(gcf_cols)) %>%
        summarise(across(everything(), ~ any(. == 1))) %>%
        pivot_longer(everything(), names_to = "GCF", values_to = "present") %>%
        filter(present) %>%
        pull(GCF)
    
    tibble(
        iter = i,
        shared = length(intersect(gcf_mats, gcf_prec)),
        mats_only = length(setdiff(gcf_mats, gcf_prec)),
        precip_only = length(setdiff(gcf_prec, gcf_mats))
    )
})

summary_subsample <- res_subsample %>%
    summarise(
        shared_mean = mean(shared),
        shared_sd = sd(shared),
        mats_only_mean = mean(mats_only),
        mats_only_sd = sd(mats_only),
        precip_only_mean = mean(precip_only),
        precip_only_sd = sd(precip_only)
    )

write_csv(summary_subsample, "GCF_shared_unique_subsample_summary.csv")
write_csv(res_subsample, "GCF_shared_unique_subsample_iterations.csv")
