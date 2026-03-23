library(tidyverse)

df <- read_csv(
    "bgc_level_antismash_v8_with_genome_metadata.csv",
    show_col_types = FALSE
)

bigscape_df <- read_tsv(
    "record_annotations.tsv",
    show_col_types = FALSE
)

novelty_df <- read_tsv(
    "novelty_summary_bgc_ids_v8.tsv",
    show_col_types = FALSE
)

merged_df <- df %>%
    left_join(
        bigscape_df %>%
            transmute(
                RedSea_BGC_ID,
                Class,
                Category,
                Description
            ) %>%
            distinct(RedSea_BGC_ID, .keep_all = TRUE),
        by = "RedSea_BGC_ID"
    ) %>%
    left_join(
        novelty_df %>%
            distinct(RedSea_BGC_ID, .keep_all = TRUE),
        by = "RedSea_BGC_ID"
    )

write_csv(
    merged_df,
    "bgc_level_antismash_v8_with_bigscapeV2_and_novelty.csv"
)
