
library(tidyverse)

# Load files
novelty <- read_tsv("novelty_summary_bgc_ids_v8.tsv", show_col_types = FALSE)

bgc_full <- read_csv("BGCs_from_consistent_GCF_sets_with_all_annotations_v2.csv", show_col_types = FALSE)

# Check column consistency
# Ensure no hidden spaces
novelty <- novelty %>%
    mutate(RedSea_BGC_ID = trimws(RedSea_BGC_ID))

bgc_full <- bgc_full %>%
    mutate(RedSea_BGC_ID = trimws(RedSea_BGC_ID))

# Merge novelty info into annotated BGC table
analysis_ready_novel_GCF <- bgc_full %>%
    left_join(novelty, by = "RedSea_BGC_ID")

# Optional sanity checks
cat("Total BGCs:", nrow(analysis_ready_novel_GCF), "\n")
cat("BGCs with novelty annotation:", sum(!is.na(analysis_ready_novel_GCF$Is_Novel)), "\n")

# Save output
write_csv(analysis_ready_novel_GCF, "analysis_ready_novel_GCF.csv")
