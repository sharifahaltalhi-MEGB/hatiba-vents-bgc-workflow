# ============================================
# Identification of stable shared and habitat-specific GCFs
# ============================================
#
# Purpose:
# Identify gene cluster families (GCFs) that are consistently shared
# or habitat-specific across balanced subsampling iterations.
#
# Approach:
# - Perform balanced subsampling of microbial mat and precipitate samples
#   (equal sample size per iteration)
# - Repeat across 1000 iterations
# - For each GCF, calculate frequency of occurrence as:
#     * shared (present in both habitats)
#     * mats only
#     * precipitates only
# - Define stable GCF sets using a frequency threshold (≥ 0.8)
#
# Outputs:
# - GCF_stability_frequencies_v2.csv
#     Frequency of each GCF across categories
# - GCFs_consistently_shared_v2.csv
# - GCFs_consistently_mats_only_v2.csv
# - GCFs_consistently_precip_only_v2.csv
#
# Notes:
# - These stable GCF sets are used for all downstream analyses
#   (taxonomy, novelty, and enrichment)
# ============================================

library(tidyverse)

set.seed(123)

df <- read_csv("GCF_sample_presence_absence_full_TypeOnly.csv", show_col_types = FALSE)

gcf_cols <- names(df)[str_detect(names(df), "^FAM_")]

n_mats <- sum(df$Type == "microbial mat")
n_prec <- sum(df$Type == "precipitate")
n_min  <- min(n_mats, n_prec)

n_iter <- 1000

gcf_iter <- map_dfr(seq_len(n_iter), function(i) {
    
    mats_sub <- df %>%
        filter(Type == "microbial mat") %>%
        slice_sample(n = n_min)
    
    prec_sub <- df %>%
        filter(Type == "precipitate") %>%
        slice_sample(n = n_min)
    
    mats_pres <- mats_sub %>%
        summarise(across(all_of(gcf_cols), ~ any(. == 1))) %>%
        pivot_longer(everything(), names_to = "GCF", values_to = "in_mats")
    
    prec_pres <- prec_sub %>%
        summarise(across(all_of(gcf_cols), ~ any(. == 1))) %>%
        pivot_longer(everything(), names_to = "GCF", values_to = "in_precip")
    
    mats_pres %>%
        left_join(prec_pres, by = "GCF") %>%
        mutate(iter = i)
})

gcf_stability <- gcf_iter %>%
    group_by(GCF) %>%
    summarise(
        freq_shared = mean(in_mats & in_precip),
        freq_mats_only = mean(in_mats & !in_precip),
        freq_precip_only = mean(!in_mats & in_precip),
        .groups = "drop"
    )

threshold <- 0.8

shared_stable <- gcf_stability %>% filter(freq_shared >= threshold)
mats_only_stable <- gcf_stability %>% filter(freq_mats_only >= threshold)
precip_only_stable <- gcf_stability %>% filter(freq_precip_only >= threshold)

write_csv(gcf_stability, "GCF_stability_frequencies_v2.csv")
write_csv(shared_stable, "GCFs_consistently_shared_v2.csv")
write_csv(mats_only_stable, "GCFs_consistently_mats_only_v2.csv")
write_csv(precip_only_stable, "GCFs_consistently_precip_only_v2.csv")
