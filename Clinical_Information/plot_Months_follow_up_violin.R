# Plot the distribution of follow-up duration from Table S2.

library(readxl)
library(dplyr)
library(ggplot2)

# ---- File settings ---------------------------------------------------------
input_file <- "/Users/sy_wei/Desktop/IBD_proteme/00_Data/Table S2.xlsx"
input_sheet <- "XHOM_PM_sampleInfor_following.x"
output_dir <- "/Users/sy_wei/Desktop/protome_enrichment/Months_follow_up_violin"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The first three rows of Table S2 contain the title and notes; row 4 is the
# actual header row.
df <- read_excel(
  path = input_file,
  sheet = input_sheet,
  skip = 3
)

if (!"Months_follow_up" %in% names(df)) {
  stop("Column 'Months_follow_up' was not found in Table S2.")
}

df_fm <- df %>%
  mutate(Months_follow_up = suppressWarnings(as.numeric(Months_follow_up))) %>%
  filter(!is.na(Months_follow_up))

if (nrow(df_fm) == 0) {
  stop("No valid numeric values were found in 'Months_follow_up'.")
}

q90 <- quantile(df_fm$Months_follow_up, probs = 0.90, na.rm = TRUE)
max_months <- max(df_fm$Months_follow_up, na.rm = TRUE)
upper_limit <- max(1, max_months * 1.25)
annotation_y <- max_months * 1.15

# Fix the random jitter positions so repeated runs produce the same figure.
set.seed(20260726)

p <- ggplot(df_fm, aes(x = "Samples pool", y = Months_follow_up)) +
  geom_violin(
    width = 0.9,
    trim = FALSE,
    fill = NA,
    color = "#FF6666",
    linewidth = 1.2
  ) +
  geom_boxplot(
    width = 0.18,
    fill = "white",
    color = "#FF6666",
    linewidth = 1.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.015,
    height = 0,
    size = 1.6,
    alpha = 0.8,
    color = "#FF6666"
  ) +
  annotate(
    "text",
    x = 1,
    y = annotation_y,
    label = paste0("Quantile90: ", round(q90, 2)),
    size = 5,
    hjust = 0.5
  ) +
  labs(
    x = NULL,
    y = "Months of follow-up"
  ) +
  scale_y_continuous(
    breaks = seq(0, ceiling(upper_limit / 6) * 6, by = 6),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(ylim = c(0, upper_limit)) +
  theme_classic(base_size = 14)

print(p)

ggsave(
  filename = file.path(output_dir, "Months_follow_up_violin.png"),
  plot = p,
  width = 5,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Months_follow_up_violin.pdf"),
  plot = p,
  width = 5,
  height = 6,
  units = "in",
  device = "pdf",
  bg = "white"
)

message("Valid samples: ", nrow(df_fm))
message("90th percentile: ", round(q90, 2), " months")
message("Plots saved in: ", output_dir)
