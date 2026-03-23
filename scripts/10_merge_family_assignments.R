suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
})

mix_file <- "mix_clustering_c0.7_v2.tsv"
bgc_file <- "bgc_level_antismash_v8_with_bigscapeV2_and_novelty.csv"

mix <- read_tsv(mix_file, show_col_types = FALSE)
bgc <- read_csv(bgc_file, show_col_types = FALSE)

stopifnot("RedSea_BGC_ID" %in% names(mix), "Family" %in% names(mix))
stopifnot("RedSea_BGC_ID" %in% names(bgc))

mix2 <- mix %>%
    mutate(
        RedSea_BGC_ID = str_trim(RedSea_BGC_ID),
        Family = str_trim(Family)
    )

bgc2 <- bgc %>%
    mutate(RedSea_BGC_ID = str_trim(RedSea_BGC_ID))

dup_check <- mix2 %>%
    count(RedSea_BGC_ID) %>%
    filter(n > 1)

if (nrow(dup_check) > 0) {
    conflicts <- mix2 %>%
        semi_join(dup_check, by = "RedSea_BGC_ID") %>%
        distinct(RedSea_BGC_ID, Family) %>%
        count(RedSea_BGC_ID) %>%
        filter(n > 1)
    
    if (nrow(conflicts) > 0) {
        print(
            mix2 %>%
                semi_join(conflicts, by = "RedSea_BGC_ID") %>%
                distinct(RedSea_BGC_ID, Family) %>%
                arrange(RedSea_BGC_ID)
        )
        stop("Same RedSea_BGC_ID maps to multiple Family values in mix_clustering file. Fix upstream before merging.")
    }
    
    mix2 <- mix2 %>%
        group_by(RedSea_BGC_ID) %>%
        summarise(
            Family = first(Family),
            .groups = "drop"
        )
} else {
    mix2 <- mix2 %>%
        select(RedSea_BGC_ID, Family)
}

merged <- bgc2 %>%
    left_join(mix2, by = "RedSea_BGC_ID")

out_file <- "bgc_level_antismash_v8_with_bigscapeV2_and_novelty_WITH_FAMILY.csv"
write_csv(merged, out_file)
