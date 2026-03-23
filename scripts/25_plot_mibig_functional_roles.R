
library(tidyverse)

df <- read.delim("MAG_best_hits_with_products_le0.45_header.tsv")

df <- df %>%
    mutate(category = case_when(
        
        product %in% c("ectoine") ~ "Osmoprotection",
        
        product %in% c("isorenieratene","carotenoid","zeaxanthin",
                       "APE Ec","APE Vf","indigoidine","dkxanthene 530") ~
            "Pigments and antioxidants",
        
        product %in% c("hopene","branched-chain fatty acids",
                       "eicosapentaenoic acid-like compound",
                       "heterocyst glycolipids") ~
            "Membrane adaptation",
        
        product %in% c("capsular polysaccharide") ~
            "Surface and biofilm",
        
        product %in% c("chloramphenicol","armeniaspirol A",
                       "fabclavine-polyamine","griseusin A",
                       "concanamycin A","cylindrocyclophane D",
                       "pyxidicycline A","fusaricidin B",
                       "Nocuolin A","chlorosphaerolactylate D") ~
            "Interaction metabolites",
        
        product %in% c("PreQ0 Base","murayaquinone","biotin") ~
            "Cofactors and metabolism",
        
        TRUE ~ "Other specialized compounds"
    ))

sum_df <- df %>%
    count(category, sort = TRUE) %>%
    mutate(category = factor(category, levels = rev(category)))

# muted consistent palette
cat_cols <- c(
    "Interaction metabolites" = "#8C6D62",
    "Pigments and antioxidants" = "#7A8F6B",
    "Membrane adaptation" = "#6F7F99",
    "Osmoprotection" = "#7FA6A3",
    "Surface and biofilm" = "#A08AA3",
    "Cofactors and metabolism" = "#9A8F6A",
    "Other specialized compounds" = "#9A9A9A"
)

p <- ggplot(sum_df,
            aes(x = n,
                y = category,
                fill = category)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = n), hjust = -0.2, size = 4) +
    scale_fill_manual(values = cat_cols) +
    labs(
        x = "Number of best MiBIG hits",
        y = "",
        title = "Functional roles inferred from closest MiBIG neighbors",
        subtitle = "MAG BGCs with BiG-SCAPE distance ≤ 0.45"
    ) +
    expand_limits(x = max(sum_df$n) * 1.12) +
    theme_bw(base_size = 13) +
    theme(
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11)
    )

ggsave("MiBIG_functional_roles_barplot_le0.45.pdf",
       p,
       width = 7,
       height = 5)

ggsave("MiBIG_functional_roles_barplot_le0.45.png",
       p,
       width = 7,
       height = 5,
       dpi = 300)
