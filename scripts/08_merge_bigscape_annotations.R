library(tidyverse)

df <- read_csv(
    "bgc_level_antismash_v8_with_genome_metadata.csv",
    show_col_types = FALSE
)

bigscape_df <- read_tsv(
    "record_annotations.tsv",
    show_col_types = FALSE
)

merged_df <- df %>%
    mutate(RedSea_BGC_ID = str_remove(RedSea_BGC_ID, "\\.gbk$")) %>%
    left_join(
        bigscape_df %>%
            transmute(
                RedSea_BGC_ID,
                Description,
                Class,
                Category
            ) %>%
            distinct(RedSea_BGC_ID, .keep_all = TRUE),
        by = "RedSea_BGC_ID"
    )

write_csv(
    merged_df,
    "bgc_level_antismash_v8_with_bigscapeV2_annotations.csv"
)
