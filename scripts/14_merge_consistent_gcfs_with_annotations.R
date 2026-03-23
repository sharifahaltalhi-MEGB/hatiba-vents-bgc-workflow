# ============================================
# Merge stable GCF sets with BGC annotations
# ============================================
#
# Purpose:
# Link consistently shared and habitat-specific gene cluster families (GCFs)
# to their underlying BGC annotations and metadata.
#
# Approach:
# - Read stable GCF sets defined from balanced subsampling
# - Read BiG-SCAPE clustering output to map BGCs to GCFs
# - Read annotated BGC table
# - Keep only BGCs assigned to stable GCF sets
# - Merge habitat category, stability frequency, and annotation metadata
#
# Outputs:
# - BGCs_from_consistent_GCF_sets_with_all_annotations_v2.csv
#
# Notes:
# - This merged table is the foundation for downstream analyses of
#   taxonomy, novelty, and habitat-specific biosynthetic patterns
# ============================================

library(tidyverse)
library(stringr)

freq_cutoff_shared <- 0.80
freq_cutoff_unique <- 0.80

bgc_file <- "bgc_level_antismash_v8_with_bigscapeV2_annotations.csv"
clust_file <- "mix_clustering_c0.7_v2.tsv"

gcf_shared_file <- "GCFs_consistently_shared_v2.csv"
gcf_mats_file   <- "GCFs_consistently_mats_only_v2.csv"
gcf_prec_file   <- "GCFs_consistently_precip_only_v2.csv"

needed <- c(bgc_file, clust_file, gcf_shared_file, gcf_mats_file, gcf_prec_file)
missing <- needed[!file.exists(needed)]

if (length(missing) > 0) {
    stop("Missing files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

clean_bgc <- function(x) {
    x %>%
        as.character() %>%
        str_trim() %>%
        str_remove("\\.gbk$") %>%
        str_remove("^\\./") %>%
        str_replace_all("\\s+", "")
}

gcf_sets_raw <- bind_rows(
    read_csv(gcf_shared_file, show_col_types = FALSE) %>%
        filter(freq_shared >= freq_cutoff_shared) %>%
        transmute(GCF, category = "Shared", freq = freq_shared),
    
    read_csv(gcf_mats_file, show_col_types = FALSE) %>%
        filter(freq_mats_only >= freq_cutoff_unique) %>%
        transmute(GCF, category = "Mats_only", freq = freq_mats_only),
    
    read_csv(gcf_prec_file, show_col_types = FALSE) %>%
        filter(freq_precip_only >= freq_cutoff_unique) %>%
        transmute(GCF, category = "Precip_only", freq = freq_precip_only)
)

overlaps <- gcf_sets_raw %>%
    count(GCF) %>%
    filter(n > 1)

if (nrow(overlaps) > 0) {
    stop(
        "Some GCFs appear in multiple sets. Fix upstream files.\n",
        paste(overlaps$GCF, collapse = ", "),
        call. = FALSE
    )
}

gcf_sets <- gcf_sets_raw %>%
    distinct(GCF, .keep_all = TRUE)

clust_map <- read_tsv(clust_file, show_col_types = FALSE) %>%
    rename(
        BGC = RedSea_BGC_ID,
        GCF = Family
    ) %>%
    mutate(
        BGC = clean_bgc(BGC),
        GCF = str_trim(GCF)
    ) %>%
    filter(!is.na(GCF), GCF != "") %>%
    select(BGC, GCF) %>%
    distinct()

clust_map_keep <- clust_map %>%
    semi_join(gcf_sets, by = "GCF")

bgc_df <- read_csv(bgc_file, show_col_types = FALSE) %>%
    mutate(BGC = clean_bgc(RedSea_BGC_ID))

out <- clust_map_keep %>%
    left_join(gcf_sets, by = "GCF") %>%
    left_join(bgc_df, by = "BGC") %>%
    relocate(category, freq, GCF, .after = RedSea_BGC_ID)

write_csv(out, "BGCs_from_consistent_GCF_sets_with_all_annotations_v2.csv")
