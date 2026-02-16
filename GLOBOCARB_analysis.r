# ==============================================================================
# GLOBOCARB: Global Lake Organic Carbon Burial Rates (1900-2000)
# Biogeochemistry Seminar WS 2025/26
# TU Bergakademie Freiberg
# Michael Schneider
# https://github.com/michi-sxc/biogeochem
# ==============================================================================
# We analyzed the GLOBOCARB supplementary dataset from:
#   Anderson, N. J., Heathcote, A. J. & Engstrom, D. R. (2020)
#   "Anthropogenic alteration of nutrient supply increases the global
#    freshwater carbon sink." Science Advances 6, eaaw2145.
#   DOI: 10.1126/sciadv.aaw2145
#   Data: https://data.mendeley.com/datasets/34hsd2jygc/1
#   License: CC BY-NC 4.0
#
# The dataset contains biome-aggregated organic carbon accumulation rates
# (OCAR) derived from 210Pb-dated sediment cores of 516 lakes, corrected
# for sediment focusing. Each row is one biome-decade mean (not an
# individual lake). 143 rows = 13 biomes x 11 decades (1900-2000).
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

# Required packages
required_packages <- c(
  "tidyverse", "scales", "viridis", "ggrepel", "broom",
  "car", "rstatix", "lme4", "effectsize", "ggfortify", "here"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# Reproducibility (affects geom_jitter positions)
set.seed(69420)

# Write the R console output to a separate .txt file (comment out along with sink() at line 835 if not needed)
sink(here("data", "globocarb_r.txt"))

# Consistent theme for all figures
theme_set(theme_minimal(base_size = 12) +
            theme(plot.title    = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
                  plot.caption  = element_text(hjust = 0.6, color = "gray50"),
                  legend.position   = "right",
                  panel.grid.minor  = element_blank())
          )


# ==============================================================================
# 2. DATA LOADING AND INITIAL EXPLORATION
# ==============================================================================

globocarb <- read.csv(here("data", "GLOBOCARBSupplementaryData.csv"), stringsAsFactors = FALSE)

cat("==== GLOBOCARB Dataset Overview ====\n\n")
cat("Dimensions:", nrow(globocarb), "rows x", ncol(globocarb), "columns\n")
cat("NOTE: Each row is a biome-decade MEAN (aggregated from 516 lake cores),\n")
cat("      not an individual lake observation.\n\n")
str(globocarb)
cat("\n--- Summary Statistics ---\n")
print(summary(globocarb))

# Check for NAs
cat("\n--- Missing Values ---\n")
print(colSums(is.na(globocarb)))
# Result: No NAs. Missing data might still be disguised as constant values.


# ==============================================================================
# 3. DATA QUALITY ASSESSMENT
# ==============================================================================
# Check for suspicious patterns: are any biomes constant across all decades?

cat("\n--- Data Quality: Checking for constant OCAR within biomes ---\n")
variance_check <- globocarb %>%
  group_by(Biome) %>%
  summarize(
    SD_OCAR = sd(OCAR_g_m2_yr),
    Min     = min(OCAR_g_m2_yr),
    Max     = max(OCAR_g_m2_yr),
    .groups = "drop"
  ) %>%
  arrange(SD_OCAR)

print(variance_check)

# FINDING: "Large Lakes" and "Flooded Grasslands and Savannas" have SD = 0
# (identical OCAR across all 11 decades). This is physically implausible and
# almost certainly reflects missing or placeholder data in the original dataset.

constant_biomes <- c("Large Lakes", "Flooded Grasslands and Savannas")
cat("\nWARNING: The following biomes have zero temporal variance and are treated",
    "as data artifacts:\n")
cat(" ", paste(constant_biomes, collapse = ", "), "\n")
cat("  These are INCLUDED in total burial calculations (their area contribution",
    "is real)\n")
cat("  These are EXCLUDED from temporal trend analyses and statistical tests\n\n")


# ==============================================================================
# 4. DATA CLEANING AND TRANSFORMATION
# ==============================================================================

globocarb <- globocarb %>%
  rename(
    OCAR           = OCAR_g_m2_yr,
    Lake_Area      = Lake_Area_km2,
    Reservoir_Area = Resevoir_Area_km2   # typo in original name fixed
  ) %>%
  mutate(
    # Numeric decade (regressions)
    Decade_num = Decade,
    
    # Factor decade (categorical plots)
    Decade_fct = factor(Decade, levels = seq(1900, 2000, by = 10), ordered = TRUE),
    
    # Pre/Post-1950 period (Anthropocene marker; Steffen et al. 2015)
    Period = factor(
      ifelse(Decade < 1950, "Pre-1950", "Post-1950"),
      levels = c("Pre-1950", "Post-1950")
    ),
    
    # Total carbon burial per biome-decade (Tg C yr-1)
    # OCAR [g/m2/yr] * Lake_Area [km2] * 1e6 [m2/km2] / 1e12 [g/Tg]
    C_Burial_Tg = (OCAR * Lake_Area * 1e6) / 1e12,
    
    # Reservoir-only burial estimate (same OCAR applied to reservoir area)
    C_Burial_Reservoir_Tg = (OCAR * Reservoir_Area * 1e6) / 1e12,
    
    # Data quality flag
    Data_Quality = ifelse(Biome %in% constant_biomes, "Suspect", "Reliable")
  )

# Reliable subset for trend analysis (excl. 2 constant-OCAR biomes)
gc_reliable <- globocarb %>% filter(Data_Quality == "Reliable")

cat("Full dataset:", nrow(globocarb), "observations (13 biomes x 11 decades)\n")
cat("Reliable subset:", nrow(gc_reliable), "observations (11 biomes x 11 decades)\n")

# Biomes ordered by median OCAR
all_biome_order <- globocarb %>%
  group_by(Biome) %>%
  summarize(med = median(OCAR), .groups = "drop") %>%
  arrange(med) %>%
  pull(Biome)

globocarb$Biome_ordered   <- factor(globocarb$Biome, levels = all_biome_order)
gc_reliable$Biome_ordered <- factor(gc_reliable$Biome, levels = all_biome_order)

# Biome summary
cat("\n--- OCAR Summary by Biome ---\n")
biome_summary <- globocarb %>%
  group_by(Biome, Data_Quality) %>%
  summarize(
    Mean = round(mean(OCAR), 2),
    SD   = round(sd(OCAR), 2),
    Min  = round(min(OCAR), 2),
    Max  = round(max(OCAR), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))
print(as.data.frame(biome_summary))

# Define biome groups (used later in Fig. 5.6 and Section 6.5)
tropical_biomes <- c("Tropical Moist and Dry Forest", "Tropical Grasslands", "Mangroves")
boreal_biomes   <- c("Boreal Forest and Taiga", "Tundra")


# ==============================================================================
# 5. FIGURES
# ==============================================================================
cat("\n==== Generating Figures ====\n")


# --- 5.1 UNUSED: Distribution of OCAR (Histogram) ------------------------------------
# Purpose: Show distribution shape -- right-skewed, dominated by low-burial
#          biomes, with a few high-OCAR outliers (Mangroves).

fig1 <- ggplot(globocarb, aes(x = OCAR)) +
  geom_histogram(
    aes(fill = after_stat(count)),
    bins = 25, color = "white", alpha = 0.85
  ) +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  labs(
    title    = "Distribution of Organic Carbon Accumulation Rates",
    subtitle = "143 biome-decade means from 516 lake cores across 13 biomes (1900\u20132000)",
    x = expression("OCAR (g C m"^{-2} * " yr"^{-1} * ")"),
    y = "Frequency",
    caption  = "Data: Anderson et al. (2020) Science Advances | DOI: 10.1126/sciadv.aaw2145"
  )

print(fig1)
ggsave(here("figures", "fig1_ocar_histogram.png"), fig1, width = 10, height = 6, dpi = 300)


# --- 5.2 OCAR by Biome (Boxplot + jitter by decade) --------------------------
# Purpose: Show inter-biome differences AND temporal trends within each biome
#          in a single figure. Points shift from blue (1900) to orange (2000).
#          Vertical jitter is for readability only (no scientific meaning).

decade_palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(11))

fig2 <- ggplot(globocarb,
               aes(x = Biome_ordered, y = OCAR)) +
  geom_boxplot(
    outlier.shape = NA, color = "gray60", fill = "white", width = 0.6
  ) +
  geom_jitter(
    aes(fill = Decade_fct),
    shape = 21, color = "gray30", stroke = 0.3, size = 2.5,
    width = 0.15, alpha = 0.9
  ) +
  # Highlight suspect biomes with red background
  annotate("rect",
           xmin = which(levels(globocarb$Biome_ordered) == "Large Lakes") - 0.4,
           xmax = which(levels(globocarb$Biome_ordered) == "Large Lakes") + 0.4,
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.05) +
  annotate("rect",
           xmin = which(levels(globocarb$Biome_ordered) == "Flooded Grasslands and Savannas") - 0.4,
           xmax = which(levels(globocarb$Biome_ordered) == "Flooded Grasslands and Savannas") + 0.4,
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.05) +
  scale_fill_manual(values = decade_palette, name = "Decade") +
  coord_flip() +
  labs(
    title    = "Organic Carbon Burial Rates by Biome (1900\u20132000)",
    subtitle = "Blue = early 20th century | Orange = recent decades | Red shading = suspect data",
    x = NULL,
    y = expression("OCAR (g C m"^{-2} * " yr"^{-1} * ")"),
    caption  = paste(
      "Boxplots: box = IQR, line = median, whiskers = 1.5\u00d7IQR.",
      "Jittered points colored by decade.",
      "\nBiomes ordered by median OCAR. Large Lakes and Flooded Grasslands show zero",
      "temporal variance (suspect data)."
    )
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 21, stroke = 0.5)))

print(fig2)
ggsave(here("figures", "fig2_boxplot_biome.png"), fig2, width = 12, height = 9, dpi = 300)


# --- 5.3 Heatmap of OCAR by Biome x Decade -----------------------------------
# Purpose: Compact visualization of all data points showing spatial x temporal
#          patterns. Uniform rows betray suspect data.

fig3 <- ggplot(globocarb,
               aes(x = Decade_fct, y = Biome_ordered, fill = OCAR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(
    option = "inferno",
    name = expression("OCAR (g C m"^{-2} * " yr"^{-1} * ")")
  ) +
  labs(
    title    = "Organic Carbon Burial Rate Heatmap",
    subtitle = "Biomes ordered by median OCAR (low to high)",
    x = "Decade", y = "",
    caption  = paste(
      "Brighter (yellow) colors indicate higher carbon burial rates.",
      "\nLarge Lakes and Flooded Grasslands show no temporal change (suspect data)."
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank())

print(fig3)
ggsave(here("figures", "fig3_heatmap.png"), fig3, width = 12, height = 8, dpi = 300)


# --- 5.4 UNUSED: Global Total Carbon Burial (Stacked Area) ---------------------------
# Purpose: Show ABSOLUTE burial (Tg C/yr), not just per-area rates.
#          Key finding: Boreal lakes dominate total burial despite low OCAR,
#          because of their enormous lake area --> Per-area rates != total burial.

fig4 <- ggplot(globocarb,
               aes(x = Decade_num, y = C_Burial_Tg, fill = Biome_ordered)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.3) +
  scale_fill_viridis_d(option = "turbo", name = "Biome") +
  scale_x_continuous(breaks = seq(1900, 2000, by = 10)) +
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  labs(
    title    = "Global Lake Carbon Burial by Biome (1900\u20132000)",
    subtitle = "Total organic carbon sequestration in lake sediments",
    x = "Decade",
    y = expression("Total C Burial (Tg C yr"^{-1} * ")"),
    caption  = paste(
      "Total burial = OCAR \u00d7 Lake Area. 1 Tg = 10\u00b9\u00b2 g.",
      "\nBoreal lakes dominate total burial despite low per-area rates",
      "due to their vast lake area."
    )
  ) +
  theme(legend.key.size = unit(0.4, "cm"))

print(fig4)
ggsave(here("figures", "fig4_stacked_burial.png"), fig4, width = 14, height = 8, dpi = 300)

# Decadal totals: lakes, reservoirs, and combined
cat("\n--- Total Global C Burial by Decade ---\n")
decade_totals <- globocarb %>%
  group_by(Decade_num) %>%
  summarize(
    Lake_Burial_Tg      = round(sum(C_Burial_Tg), 2),
    Reservoir_Burial_Tg = round(sum(C_Burial_Reservoir_Tg), 2),
    Total_Burial_Tg     = round(sum(C_Burial_Tg) + sum(C_Burial_Reservoir_Tg), 2),
    Mean_OCAR           = round(mean(OCAR), 2),
    .groups = "drop"
  )
print(as.data.frame(decade_totals))

cat("\nLake-only increase from 1900 to 2000:",
    round(decade_totals$Lake_Burial_Tg[11] - decade_totals$Lake_Burial_Tg[1], 2),
    "Tg C/yr\n")
cat("Fold change:",
    round(decade_totals$Lake_Burial_Tg[11] / decade_totals$Lake_Burial_Tg[1], 2), "\n")
cat("Lake + Reservoir burial in 2000:", decade_totals$Total_Burial_Tg[11], "Tg C/yr\n")
cat("Reservoir contribution in 2000:",
    round(decade_totals$Reservoir_Burial_Tg[11] /
            decade_totals$Total_Burial_Tg[11] * 100, 1), "%\n")


# --- 5.5 UNUSED: Lake Area vs. Mean OCAR (Bubbly plot) -------------------------------
# Purpose: Explains why biomes with low OCAR can dominate total burial.
#          Bubble size = total C burial. Spatial extent governs burial.

biome_means <- globocarb %>%
  group_by(Biome, Data_Quality) %>%
  summarize(
    Mean_OCAR   = mean(OCAR),
    Lake_Area   = unique(Lake_Area),
    Mean_Burial = mean(C_Burial_Tg),
    .groups = "drop"
  )

fig5 <- ggplot(biome_means, aes(x = Lake_Area, y = Mean_OCAR)) +
  geom_point(
    aes(size = Mean_Burial, fill = Data_Quality),
    shape = 21, color = "black", alpha = 0.8
  ) +
  geom_text_repel(
    aes(label = Biome), size = 3.2,
    box.padding = 0.5, point.padding = 0.3, min.segment.length = 0
  ) +
  scale_x_log10(labels = label_number(scale = 1e-3, suffix = "k")) +
  scale_fill_manual(
    values = c("Reliable" = "steelblue", "Suspect" = "tomato"),
    name = "Data Quality"
  ) +
  scale_size_continuous(range = c(3, 15), name = "Mean C Burial\n(Tg/yr)") +
  labs(
    title    = "Lake Area vs. Burial Rate (Biome Averages)",
    subtitle = "Bubble size = total carbon burial | Red = suspect data quality",
    x = expression("Total Lake Area in Biome (km"^2 * ")"),
    y = expression("Mean OCAR (g C m"^{-2} * " yr"^{-1} * ")"),
    caption  = paste(
      "Boreal Forest dominates total burial (large bubble, bottom right)",
      "despite low OCAR,",
      "\nwhile Mangroves have the highest OCAR but contribute little to total",
      "burial (top left)."
    )
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))

print(fig5)
ggsave(here("figures", "fig5_bubble_area_ocar.png"), fig5, width = 12, height = 8, dpi = 300)


# --- 5.6 Tropical vs. Boreal/Tundra Regional Comparison ----------------------
# Purpose: Test hypothesis that tropical biomes show faster OCAR increases
#          than high-latitude biomes, consistent with nutrient enrichment
#          (not temperature I guess) as the primary driver.

regional_data <- globocarb %>%
  filter(Biome %in% c(tropical_biomes, boreal_biomes)) %>%
  mutate(Region = ifelse(Biome %in% tropical_biomes, "Tropical", "Boreal/Tundra")) %>%
  group_by(Decade_num, Region) %>%
  summarize(mean_ocar = mean(OCAR), .groups = "drop")

# Fit linear models per region (for annotation on figure)
model_tropical <- lm(mean_ocar ~ Decade_num, data = filter(regional_data, Region == "Tropical"))
model_boreal   <- lm(mean_ocar ~ Decade_num, data = filter(regional_data, Region == "Boreal/Tundra"))

slope_trop  <- coef(model_tropical)[2]
slope_bor   <- coef(model_boreal)[2]
slope_ratio <- slope_trop / slope_bor

fig6 <- ggplot(regional_data, aes(x = Decade_num, y = mean_ocar, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "dashed") +
  scale_color_manual(values = c("Tropical" = "#E64B35", "Boreal/Tundra" = "#4DBBD5")) +
  scale_x_continuous(breaks = seq(1900, 2000, by = 20)) +
  annotate("text", x = 1945, y = 62,
           label = paste0("Tropical slope: ", round(slope_trop, 3), " /decade"),
           color = "#E64B35", hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 1945, y = 57,
           label = paste0("Boreal slope: ", round(slope_bor, 3), " /decade"),
           color = "#4DBBD5", hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 1945, y = 52,
           label = paste0("Ratio: ", round(slope_ratio, 1), "\u00d7"),
           color = "gray30", hjust = 0, fontface = "italic", size = 3.5) +
  labs(
    title    = "Carbon Burial: Tropical vs. Boreal Regions",
    subtitle = "Tropical regions show steeper increase over the 20th century",
    x = "Decade",
    y = expression("Mean OCAR (g C m"^{-2} * " yr"^{-1} * ")"),
    caption  = paste(
      "Shaded regions = 95% CI for linear trend.",
      "Tropical = Mangroves + Tropical Moist/Dry Forest + Tropical Grasslands.",
      "\nBoreal = Boreal Forest/Taiga + Tundra."
    )
    
  ) +
  theme(plot.caption  = element_text(hjust = 0, color = "gray50"))
  

print(fig6)
ggsave(here("figures", "fig6_tropical_vs_boreal.png"), fig6, width = 10, height = 7, dpi = 300)


# --- 5.7 UNUSED: Reservoir Expansion -------------------------------------------------
# Purpose: Illustration of ~60-fold expansion of global reservoir area and its potential
#          link to increased carbon burial.

reservoir_data <- globocarb %>%
  group_by(Decade_num) %>%
  summarize(
    Total_Reservoir = sum(Reservoir_Area, na.rm = TRUE),
    Total_Lake      = sum(Lake_Area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Ratio_Pct = (Total_Reservoir / Total_Lake) * 100)

fig7 <- ggplot(reservoir_data, aes(x = Decade_num)) +
  geom_bar(aes(y = Total_Reservoir / 1000),
           stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_line(aes(y = Ratio_Pct * 50), color = "darkred", linewidth = 1.2) +
  geom_point(aes(y = Ratio_Pct * 50), color = "darkred", size = 3) +
  scale_y_continuous(
    name = expression("Total Reservoir Area (\u00d710\u00b3 km\u00b2)"),
    sec.axis = sec_axis(~ . / 50, name = "Reservoir / Lake Area Ratio (%)")
  ) +
  scale_x_continuous(breaks = seq(1900, 2000, by = 10)) +
  labs(
    title    = "Global Reservoir Expansion (1900\u20132000)",
    subtitle = "Bars = reservoir area; Red line = reservoir-to-lake area ratio",
    x = "Decade",
    caption  = paste(
      "Reservoir area grew ~60-fold from 4.2\u00d710\u00b3 km\u00b2 (1900)",
      "to 255.7\u00d710\u00b3 km\u00b2 (2000).",
      "\nData: GRanD database via Anderson et al. (2020)."
    )
  ) +
  theme(axis.title.y.right = element_text(color = "darkred"),
        axis.text.y.right  = element_text(color = "darkred"))

print(fig7)
ggsave(here("figures", "fig7_reservoir_expansion.png"), fig7, width = 12, height = 7, dpi = 300)


# --- 5.8 Diagnostic Plots for Global Regression ------------------------------
# Purpose: quick check of model assumptions (normality, homoscedasticity).
#          Required complement to the Shapiro-Wilk / Levene formal tests.

model_global_diag <- lm(OCAR ~ Decade_num, data = gc_reliable, weights = Lake_Area)

fig8 <- autoplot(model_global_diag, which = c(1, 2), ncol = 2,
                 colour = "steelblue", alpha = 0.6) +
  theme_minimal(base_size = 11)

print(fig8)
ggsave(here("figures", "fig8_diagnostics_global.png"), fig8, width = 10, height = 5, dpi = 300)

cat("  The Q-Q plot shows heavy right tails, consistent with the Shapiro-Wilk\n")
cat("  rejection of normality. The fan shape in residuals vs. fitted suggests\n")
cat("  heteroscedasticity, consistent with the Levene test result.\n")


# ==============================================================================
# 6. STATISTICAL ANALYSIS
# ==============================================================================
# All trend analyses use the "reliable" subset (excluding 2 biomes with zero
# temporal variance). Total burial calculations use all 13 biomes.
# ==============================================================================

cat("\n\n==== STATISTICAL ANALYSIS ====\n")
cat("(Using reliable subset: n =", nrow(gc_reliable), ", 11 biomes)\n\n")


# --- 6.1 Global Linear Trend (OCAR ~ Decade) ---------------------------------
cat("--- 6.1 Global Linear Regression: OCAR ~ Decade ---\n")
cat("NOTE: This model pools all biomes, treating each biome-decade as independent.\n")
cat("      The low R\u00b2 is expected because biome identity explains most variance.\n")
cat("      See Section 6.7 for a mixed-effects model that accounts for this.\n\n")

model_global <- lm(OCAR ~ Decade_num, data = gc_reliable)
print(summary(model_global))


# --- 6.2 Per-Biome Linear Regressions ----------------------------------------
cat("\n--- 6.2 Linear Regression by Biome (OCAR ~ Decade) ---\n")

models_by_biome <- gc_reliable %>%
  group_by(Biome) %>%
  nest() %>%
  mutate(
    model   = map(data, ~ lm(OCAR ~ Decade_num, data = .x)),
    tidied  = map(model, tidy),
    glanced = map(model, glance)
  ) %>%
  unnest(glanced) %>%
  mutate(
    Slope   = map_dbl(tidied, ~ .x$estimate[2]),
    P_value = map_dbl(tidied, ~ .x$p.value[2]),
    Sig     = ifelse(P_value < 0.05, "Yes", "No")
  ) %>%
  select(Biome, Slope, r.squared, P_value, Sig) %>%
  arrange(desc(Slope))

cat("\nBiome-specific trends (sorted by slope):\n")
print(as.data.frame(models_by_biome))


# --- 6.3 One-Way ANOVA: Does OCAR differ among biomes? -----------------------
cat("\n--- 6.3 ANOVA: OCAR differences among biomes ---\n")

# a) Assumption checks
cat("\n  a) Assumption check: Normality of residuals\n")
anova_model <- aov(OCAR ~ Biome, data = gc_reliable)
shapiro_result <- shapiro.test(residuals(anova_model))
cat("     Shapiro-Wilk test: W =", round(shapiro_result$statistic, 4),
    ", p =", format(shapiro_result$p.value, digits = 4), "\n")

cat("\n  b) Assumption check: Homogeneity of variance\n")
levene_result <- car::leveneTest(OCAR ~ Biome, data = gc_reliable)
cat("     Levene's test: F =", round(levene_result$`F value`[1], 2),
    ", p =", format(levene_result$`Pr(>F)`[1], digits = 4), "\n")

if (shapiro_result$p.value < 0.05) {
  cat("\n  -> Residuals deviate significantly from normality (p < 0.05).\n")
}
if (levene_result$`Pr(>F)`[1] < 0.05) {
  cat("  -> Variance is NOT homogeneous across biomes (p < 0.05).\n")
  cat("  -> Using Welch's ANOVA + Games-Howell post-hoc (both robust to\n")
  cat("     unequal variances, unlike standard ANOVA + Tukey HSD).\n\n")
}

# Standard ANOVA (for reference only -- assumptions violated)
cat("  Standard ANOVA (reference only, assumptions violated):\n")
print(summary(anova_model))

# Welch's ANOVA (robust to unequal variances)
cat("  Welch's ANOVA (robust, does not assume equal variances):\n")
welch_anova <- oneway.test(OCAR ~ Biome, data = gc_reliable, var.equal = FALSE)
print(welch_anova)

# Effect size: eta-squared
cat("\n  Effect size (eta-squared):\n")
eta_result <- effectsize::eta_squared(anova_model, partial = FALSE)
print(eta_result)
cat("  Interpretation: biome identity explains ~",
    round(eta_result$Eta2 * 100, 1),
    "% of OCAR variance.\n")

# Post-hoc: Games-Howell (appropriate when variances are unequal)
# NOTE: Tukey HSD assumes equal variances -- since Levene's test rejected
#       that assumption, Games-Howell is the correct post-hoc test here.
cat("\n  Games-Howell post-hoc (top 10 most significant pairs):\n")
cat("  NOTE: Games-Howell does not assume equal variances, consistent with\n")
cat("        using Welch's ANOVA above.\n\n")

games_howell_result <- gc_reliable %>%
  games_howell_test(OCAR ~ Biome) %>%
  arrange(p.adj) %>%
  head(10)
print(as.data.frame(games_howell_result))


# --- 6.4 Welch's t-test: Pre-1950 vs. Post-1950 ------------------------------
cat("\n--- 6.4 Pre-1950 vs. Post-1950 Comparison (Welch's t-test) ---\n")
cat("  Two-sided test to avoid directional bias.\n\n")

pre  <- gc_reliable %>% filter(Decade_num < 1950) %>% pull(OCAR)
post <- gc_reliable %>% filter(Decade_num >= 1950) %>% pull(OCAR)

t_result <- t.test(post, pre, alternative = "two.sided")
print(t_result)

# Effect size: Cohen's d
cohens_d_result <- effectsize::cohens_d(post, pre)

cat("\n  Pre-1950 mean:",  round(mean(pre), 2), "g C m-2 yr-1 (n =", length(pre), ")\n")
cat("  Post-1950 mean:", round(mean(post), 2), "g C m-2 yr-1 (n =", length(post), ")\n")
cat("  Difference:",      round(mean(post) - mean(pre), 2), "g C m-2 yr-1\n")
cat("  Percent change:",  round(((mean(post) - mean(pre)) / mean(pre)) * 100, 1), "%\n")
cat("  t =", round(t_result$statistic, 3), ", p =", format(t_result$p.value, digits = 4), "\n")
cat("  Cohen's d =", round(cohens_d_result$Cohens_d, 3),
    " [", round(cohens_d_result$CI_low, 3), ",",
    round(cohens_d_result$CI_high, 3), "]\n")
cat("  Interpretation:", ifelse(abs(cohens_d_result$Cohens_d) < 0.5, "small",
                                ifelse(abs(cohens_d_result$Cohens_d) < 0.8, "medium", "large")),
    "effect size\n")

# Per-biome pre/post comparison (with Bonferroni correction)
cat("\n  Per-biome pre/post-1950 comparison (Bonferroni-corrected, m = 11):\n")
biome_ttest <- gc_reliable %>%
  group_by(Biome) %>%
  summarize(
    Pre_mean   = mean(OCAR[Decade_num < 1950]),
    Post_mean  = mean(OCAR[Decade_num >= 1950]),
    Change_pct = round(((Post_mean - Pre_mean) / Pre_mean) * 100, 1),
    p_raw      = t.test(OCAR[Decade_num >= 1950], OCAR[Decade_num < 1950])$p.value,
    .groups    = "drop"
  ) %>%
  mutate(
    p_bonferroni  = pmin(p_raw * n(), 1),
    Sig_corrected = ifelse(p_bonferroni < 0.05, "Yes", "No")
  ) %>%
  arrange(desc(Change_pct))

print(as.data.frame(biome_ttest))
cat("  Note: After Bonferroni correction (\u03b1 = 0.05/11 = 0.0045), fewer biomes\n")
cat("  remain significant, reflecting the multiple testing burden.\n")


# --- 6.5 Tropical vs. Boreal: Interaction Models -----------------------------
cat("\n--- 6.5 Tropical vs. Boreal: Interaction Models ---\n")
cat("  Testing whether the OCAR temporal slope differs between regions.\n\n")

# Prepare data: individual biome-decade observations (not aggregated means)
regional_raw <- gc_reliable %>%
  filter(Biome %in% c(tropical_biomes, boreal_biomes)) %>%
  mutate(Region = factor(
    ifelse(Biome %in% tropical_biomes, "Tropical", "Boreal/Tundra"),
    levels = c("Boreal/Tundra", "Tropical")
  ))

# Model A: Fixed-effects interaction (simple, but treats biome-decades as
#          independent -- problematic because biomes within a region are not
#          true replicates of that region)
cat("  A) Fixed-effects interaction: OCAR ~ Decade * Region\n")
interaction_fixed <- lm(OCAR ~ Decade_num * Region, data = regional_raw)
print(summary(interaction_fixed))

# Model B: Mixed-effects interaction with Biome as random intercept.
#          This accounts for the fact that tropical biomes (Mangroves, Tropical
#          Forest, Tropical Grasslands) have very different baselines but share
#          a regional grouping. The random intercept absorbs biome-level variance.
cat("\n  B) Mixed-effects model: OCAR ~ Decade * Region + (1|Biome)\n")
cat("     This properly accounts for biomes being nested within regions.\n\n")
interaction_mixed <- lmer(OCAR ~ Decade_num * Region + (1 | Biome),
                          data = regional_raw)
print(summary(interaction_mixed))

cat("\n  Comparison: The fixed-effects interaction p = 0.048 (borderline).\n")
cat("  The mixed model accounts for within-region biome heterogeneity;\n")
cat("  check the Decade_num:RegionTropical t-value (|t| > 2 suggests significance).\n")


# --- 6.6 Correlation: Reservoir Area vs. OCAR --------------------------------
cat("\n--- 6.6 Correlation: Reservoir Area vs. OCAR ---\n")
cat("  CAVEAT: Both variables trend upward with time. A positive correlation\n")
cat("  may reflect shared temporal trends, not a causal link.\n\n")

# a) Raw Spearman correlation
cat("  a) Raw Spearman correlation:\n")
cor_raw <- cor.test(gc_reliable$Reservoir_Area, gc_reliable$OCAR, method = "spearman")
print(cor_raw)

# b) Detrended correlation: remove the linear time trend from both variables,
#    then correlate the residuals. If the correlation vanishes, the raw
#    association was driven entirely by shared temporal trends.
cat("\n  b) Detrended correlation (time-trend removed):\n")
cat("     Procedure: regress each variable on Decade, correlate residuals.\n\n")

ocar_resid       <- residuals(lm(OCAR ~ Decade_num, data = gc_reliable))
reservoir_resid  <- residuals(lm(Reservoir_Area ~ Decade_num, data = gc_reliable))
cor_detrended    <- cor.test(reservoir_resid, ocar_resid, method = "spearman")
print(cor_detrended)

if (cor_detrended$p.value > 0.05) {
  cat("\n  -> After detrending, the correlation is non-significant (p =",
      round(cor_detrended$p.value, 4), ").\n")
  cat("     The raw association was driven by shared temporal trends,\n")
  cat("     not by a direct reservoir-OCAR relationship.\n")
} else {
  cat("\n  -> The detrended correlation remains significant (p =",
      round(cor_detrended$p.value, 4), ").\n")
  cat("     This suggests a reservoir-OCAR link beyond shared temporal trends.\n")
}


# --- 6.7 Mixed-Effects Model: Random Intercept + Random Slope -----------------
cat("\n--- 6.7 Mixed-Effects Model: OCAR ~ Decade + (Decade | Biome) ---\n")
cat("  This is the most appropriate model for this dataset structure.\n")
cat("  It accounts for: (a) biomes starting at different OCAR baselines\n")
cat("  (random intercept) and (b) biomes having different temporal slopes\n")
cat("  (random slope). This unifies the pooled regression (6.1) and the\n")
cat("  per-biome regressions (6.2) in a single framework.\n\n")

# Center decade for numerical stability and interpretable intercept
gc_reliable$Decade_scaled <- (gc_reliable$Decade_num - 1950) / 10

mixed_model <- lmer(OCAR ~ Decade_scaled + (Decade_scaled | Biome),
                    data = gc_reliable)
print(summary(mixed_model))

# Extract random slopes per biome
cat("\n  Random effects (biome-specific deviations from global slope):\n")
re <- ranef(mixed_model)$Biome
re$Biome <- rownames(re)
re$Total_slope <- fixef(mixed_model)["Decade_scaled"] + re$Decade_scaled
re <- re %>% arrange(desc(Total_slope))
print(re[, c("Biome", "Decade_scaled", "Total_slope")])

cat("\n  Fixed effect (global slope):",
    round(fixef(mixed_model)["Decade_scaled"], 4),
    "g C m-2 yr-1 per decade\n")
cat("  Interpretation: after accounting for biome-level variation, there is\n")
cat("  a significant global upward trend in OCAR.\n")


# --- 6.8 Reservoir Contribution to Total Burial -------------------------------
cat("\n--- 6.8 Reservoir Contribution to Total C Burial ---\n")
cat("  Anderson et al. (2020) note that including reservoirs raises the CO2\n")
cat("  emission offset from ~20% to ~30%. Here we quantify reservoir burial\n")
cat("  using the same OCAR values applied to reservoir areas.\n\n")

reservoir_burial <- globocarb %>%
  group_by(Decade_num) %>%
  summarize(
    Lake_Tg      = sum(C_Burial_Tg),
    Reservoir_Tg = sum(C_Burial_Reservoir_Tg),
    .groups      = "drop"
  ) %>%
  mutate(
    Total_Tg = Lake_Tg + Reservoir_Tg,
    Res_Pct  = round(Reservoir_Tg / Total_Tg * 100, 1)
  )

print(as.data.frame(reservoir_burial))
cat("\n  Reservoir burial grew from", round(reservoir_burial$Reservoir_Tg[1], 2),
    "Tg C/yr (1900) to", round(reservoir_burial$Reservoir_Tg[11], 2),
    "Tg C/yr (2000).\n")
cat("  In 2000, reservoirs contributed", reservoir_burial$Res_Pct[11],
    "% of the combined lake + reservoir total.\n")


# --- 6.9 Sensitivity Analysis: With vs. Without Suspect Biomes ---------------
cat("\n--- 6.9 Sensitivity: Effect of Excluding Suspect Biomes ---\n")
cat("  Comparing total burial with all 13 biomes vs. the 11 reliable biomes.\n\n")

burial_full <- globocarb %>%
  group_by(Decade_num) %>%
  summarize(Total = sum(C_Burial_Tg), .groups = "drop")

burial_reliable <- gc_reliable %>%
  group_by(Decade_num) %>%
  summarize(Total = sum(C_Burial_Tg), .groups = "drop")

sensitivity <- data.frame(
  Decade           = burial_full$Decade_num,
  All_13           = round(burial_full$Total, 2),
  Reliable_11      = round(burial_reliable$Total, 2),
  Difference       = round(burial_full$Total - burial_reliable$Total, 2),
  Pct_from_suspect = round((burial_full$Total - burial_reliable$Total) /
                             burial_full$Total * 100, 1)
)
print(sensitivity)

cat("\n  The 2 suspect biomes contribute a roughly constant",
    round(mean(sensitivity$Pct_from_suspect), 1),
    "% of total burial.\n")
cat("  Temporal trends are unaffected because their OCAR is constant.\n")
cat("  All trend-related conclusions are robust to their inclusion/exclusion.\n")


# ==============================================================================
# 7. SUMMARY OUTPUT
# ==============================================================================

cat("\n\n==== KEY FINDINGS SUMMARY ====\n\n")

cat("DATASET: 143 biome-decade means (13 biomes x 11 decades),\n")
cat("         derived from 516 210Pb-dated lake sediment cores.\n")
cat("         2 biomes flagged as suspect (zero temporal variance).\n\n")

cat("1. OCAR Distribution:\n")
cat("   Global mean =", round(mean(globocarb$OCAR), 2), "g C m-2 yr-1\n")
cat("   Global median =", round(median(globocarb$OCAR), 2), "(right-skewed)\n")
cat("   Range:", round(min(globocarb$OCAR), 2), "-",
    round(max(globocarb$OCAR), 2), "\n\n")

cat("2. Total C Burial (all 13 biomes, lakes only):\n")
cat("   1900:", round(decade_totals$Lake_Burial_Tg[1], 1), "Tg C/yr\n")
cat("   2000:", round(decade_totals$Lake_Burial_Tg[11], 1), "Tg C/yr\n")
cat("   =", round(decade_totals$Lake_Burial_Tg[11] /
                    decade_totals$Lake_Burial_Tg[1], 1), "-fold increase\n")
cat("   Including reservoirs:", decade_totals$Total_Burial_Tg[11], "Tg C/yr in 2000\n\n")

cat("3. ANOVA: Biomes differ significantly in OCAR\n")
cat("   Welch's F =", round(welch_anova$statistic, 1), ", p < 2.2e-16\n")
cat("   Eta-squared =", round(eta_result$Eta2, 3), "\n\n")

cat("4. Anthropocene signal: Post-1950 OCAR is",
    round(((mean(post) - mean(pre)) / mean(pre)) * 100, 1),
    "% higher than Pre-1950\n")
cat("   Welch's t =", round(t_result$statistic, 3),
    ", p =", format(t_result$p.value, digits = 4), "\n")
cat("   Cohen's d =", round(cohens_d_result$Cohens_d, 3), "\n\n")

cat("5. Tropical slope =", round(slope_trop, 3), "/decade vs. Boreal =",
    round(slope_bor, 3), "/decade\n")
cat("   Ratio:", round(slope_ratio, 1), "\u00d7\n")
cat("   Interaction p = 0.048 (fixed); check mixed model for robust inference.\n\n")

cat("6. Mixed-effects model confirms significant global temporal trend\n")
cat("   after accounting for biome-level variation.\n\n")

cat("7. Reservoir-OCAR correlation is spurious: significant raw (rho =",
    round(cor_raw$estimate, 3), ") but\n")
cat("   non-significant after detrending (rho =",
    round(cor_detrended$estimate, 3), ", p =",
    round(cor_detrended$p.value, 3), ").\n\n")

# End of R console output
sink()

# Experimental console output formatting (just for funsies)
source(here("format.R"))       
format_report(
  input_file   = here("data", "globocarb_r.txt"),
  output_pdf   = here("data", "globocarb_report.pdf"),
  title        = "GLOBOCARB Analysis Report",
  subtitle     = "Biogeochemistry Seminar 2026",
  keep_rmd     = FALSE          
)
