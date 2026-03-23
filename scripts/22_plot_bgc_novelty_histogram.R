suppressPackageStartupMessages({
    library(tidyverse)
    library(stringr)
    library(scales)
})

infile <- "bgc_level_antismash_v8_with_bigscapeV2_and_novelty_WITH_FAMILY.csv"
df <- read_csv(infile, show_col_types = FALSE)

req <- c("Category", "Min_Distance")
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

split_tokens <- function(x) {
    x <- tolower(x)
    x <- str_replace_all(x, "\\s+", "")
    str_split(x, "\\.", simplify = FALSE)
}

pick_first_present <- function(tokens, prefer) {
    tokens <- unique(tokens)
    hit <- prefer[prefer %in% tokens]
    if (length(hit) > 0) hit[1] else "Other"
}

normalize_primary_from_category <- function(cat) {
    toks <- split_tokens(cat)
    prefer <- c("nrps", "pks", "ripp", "terpene", "other")
    map_chr(toks, ~ pick_first_present(.x, prefer)) %>%
        recode(
            nrps = "NRPS",
            pks = "PKS",
            ripp = "RiPP",
            terpene = "terpene",
            other = "other",
            Other = "other"
        )
}

primary_levels <- c("NRPS", "PKS", "RiPP", "terpene", "other")

type_colors <- c(
    NRPS = "#e41a1c",
    PKS = "#377eb8",
    RiPP = "#984ea3",
    terpene = "#a6cee3",
    other = "#999999"
)

df_plot <- df %>%
    mutate(
        Category = as.character(Category),
        Min_Distance = as.numeric(Min_Distance),
        BGC_primary = normalize_primary_from_category(Category),
        BGC_primary = factor(BGC_primary, levels = primary_levels),
        Novel = Min_Distance > 900
    ) %>%
    filter(!is.na(Min_Distance))

p <- ggplot(df_plot, aes(x = Min_Distance, fill = BGC_primary)) +
    geom_histogram(binwidth = 25, boundary = 0, alpha = 0.9) +
    geom_vline(xintercept = 900, linewidth = 0.6) +
    scale_fill_manual(values = type_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
        x = "Min_Distance",
        y = "Number of BGCs",
        fill = "BGC class",
        title = "BGC novelty distribution",
        subtitle = "Novel threshold: Min_Distance > 900"
    ) +
    theme_classic(base_size = 12)

ggsave("HIST_novelty_MinDistance_by_primary_class.png", p, width = 8, height = 5, dpi = 300)
ggsave("HIST_novelty_MinDistance_by_primary_class.pdf", p, width = 8, height = 5)
