# GLOBOCARB: Global Lake Organic Carbon Burial Rates (1900–2000)

> **Biogeochemistry Seminar WS 2025/26 · TU Bergakademie Freiberg**  
> Michael Schneider · [`michi-sxc/biogeochem`](https://github.com/michi-sxc/biogeochem)

Analysis of the GLOBOCARB supplementary dataset from:

> Anderson, N. J., Heathcote, A. J. & Engstrom, D. R. (2020). *Anthropogenic alteration of nutrient supply increases the global freshwater carbon sink.* **Science Advances** 6, eaaw2145. [DOI: 10.1126/sciadv.aaw2145](https://doi.org/10.1126/sciadv.aaw2145)  
> Data: [Mendeley Data](https://data.mendeley.com/datasets/34hsd2jygc/1) · License: CC BY-NC 4.0

---

## Table of Contents

1. [Dataset Overview](#1-dataset-overview)
2. [Data Quality Assessment](#2-data-quality-assessment)
3. [OCAR by Biome](#3-ocar-by-biome)
4. [Total Global C Burial by Decade](#4-total-global-c-burial-by-decade)
5. [Statistical Analysis](#5-statistical-analysis)
   - [6.1 Global Linear Regression](#61-global-linear-regression-ocar--decade)
   - [6.2 Per-Biome Regressions](#62-linear-regression-by-biome)
   - [6.3 ANOVA: Biome Differences](#63-anova-ocar-differences-among-biomes)
   - [6.4 Pre- vs. Post-1950 Comparison](#64-pre-1950-vs-post-1950-comparison)
   - [6.5 Tropical vs. Boreal Interaction](#65-tropical-vs-boreal-interaction-models)
   - [6.6 Reservoir–OCAR Correlation](#66-correlation-reservoir-area-vs-ocar)
   - [6.7 Mixed-Effects Model](#67-mixed-effects-model-ocar--decade--decade--biome)
   - [6.8 Reservoir Contribution](#68-reservoir-contribution-to-total-c-burial)
   - [6.9 Sensitivity Analysis](#69-sensitivity-effect-of-excluding-suspect-biomes)
6. [Key Findings Summary](#6-key-findings-summary)

---

## 1. Dataset Overview

| Property | Value |
|---|---|
| Dimensions | 143 rows × 5 columns |
| Unit of observation | Biome-decade **mean** (not individual lake) |
| Source lakes | 516 ²¹⁰Pb-dated sediment cores |
| Coverage | 13 biomes × 11 decades (1900–2000) |

> **Note:** Each row is a biome-decade mean aggregated from 516 lake cores — not an individual lake observation.

### Variables

| Column | Type | Description |
|---|---|---|
| `Decade` | int | Decade start year (1900–2000) |
| `Biome` | chr | Biome name |
| `OCAR_g_m2_yr` | num | Organic carbon accumulation rate (g C m⁻² yr⁻¹) |
| `Lake_Area_km2` | num | Total lake area in biome (km²) |
| `Resevoir_Area_km2` | int | Total reservoir area in biome (km²) *(typo in original)* |

### Summary Statistics

| Statistic | Decade | OCAR (g C m⁻² yr⁻¹) | Lake Area (km²) | Reservoir Area (km²) |
|---|---|---|---|---|
| Min | 1900 | 1.864 | 31,083 | 0 |
| 1st Qu. | 1920 | 9.625 | 92,499 | 10 |
| Median | 1950 | 17.846 | 211,708 | 1,176 |
| **Mean** | **1950** | **28.845** | **376,421** | **7,430** |
| 3rd Qu. | 1980 | 38.152 | 394,236 | 7,586 |
| Max | 2000 | 145.666 | 1,644,866 | 54,197 |

 **Missing values:** None across all 5 columns.

---

## 2. Data Quality Assessment

Biomes checked for suspicious constant OCAR values across all decades:

| Biome | SD (OCAR) | Min | Max |
|---|---|---|---|
| Flooded Grasslands and Savannas | 0.00 | 31.00 | 31.00 |
| Large Lakes | 0.00 | 10.43 | 10.43 |
| Tundra | 1.21 | 1.86 | 5.70 |
| Temperate Conifer Forest | 1.77 | 3.76 | 9.47 |
| Boreal Forest and Taiga | 3.77 | 6.12 | 18.5 |
| Montane Grasslands | 4.05 | 10.9 | 24.0 |
| Deserts and Xeric Scrub | 4.94 | 1.91 | 18.5 |
| Temperate Mixed Forest | 10.1 | 17.7 | 49.4 |
| Mediterranean Forest | 14.2 | 4.37 | 39.1 |
| Temperate Grasslands | 18.3 | 24.7 | 87.8 |
| Tropical Grasslands | 19.5 | 8.64 | 71.7 |
| Tropical Moist and Dry Forest | 21.6 | 22.9 | 98.4 |
| Mangroves | 33.2 | 52.1 | 146 |

> **WARNING:** **Large Lakes** and **Flooded Grasslands and Savannas** have zero temporal variance - physically implausible, likely missing or placeholder data.
> - **Included** in total burial calculations (their lake area contribution is real)
> - **Excluded** from all temporal trend analyses and statistical tests

| Subset | Observations | Biomes |
|---|---|---|
| Full dataset | 143 | 13 |
| Reliable subset | 121 | 11 |

---

## 3. OCAR by Biome

Biomes ordered by mean OCAR (descending):

| Rank | Biome | Data Quality | Mean | SD | Min | Max |
|---|---|---|---|---|---|---|
| 1 | Mangroves |  Reliable | 92.90 | 33.25 | 52.13 | 145.67 |
| 2 | Temperate Grasslands |  Reliable | 56.18 | 18.28 | 24.73 | 87.83 |
| 3 | Tropical Moist and Dry Forest |  Reliable | 50.36 | 21.57 | 22.86 | 98.43 |
| 4 | Temperate Mixed Forest |  Reliable | 37.69 | 10.12 | 17.73 | 49.41 |
| 5 | Flooded Grasslands and Savannas | ⚠️ Suspect | 31.00 | 0.00 | 31.00 | 31.00 |
| 6 | Tropical Grasslands |  Reliable | 28.99 | 19.47 | 8.64 | 71.66 |
| 7 | Mediterranean Forest |  Reliable | 18.92 | 14.15 | 4.37 | 39.09 |
| 8 | Montane Grasslands |  Reliable | 16.10 | 4.05 | 10.92 | 23.95 |
| 9 | Boreal Forest and Taiga |  Reliable | 12.17 | 3.77 | 6.12 | 18.46 |
| 10 | Large Lakes | ⚠️ Suspect | 10.43 | 0.00 | 10.43 | 10.43 |
| 11 | Deserts and Xeric Scrub |  Reliable | 9.39 | 4.94 | 1.91 | 18.55 |
| 12 | Temperate Conifer Forest |  Reliable | 7.29 | 1.77 | 3.76 | 9.47 |
| 13 | Tundra |  Reliable | 3.57 | 1.21 | 1.86 | 5.70 |

*All values in g C m⁻² yr⁻¹*

---

## 4. Total Global C Burial by Decade

| Decade | Lake Burial (Tg/yr) | Reservoir Burial (Tg/yr) | **Total (Tg/yr)** | Mean OCAR |
|---|---|---|---|---|
| 1900 | 49.02 | 0.05 | **49.07** | 15.42 |
| 1910 | 58.26 | 0.09 | **58.36** | 18.81 |
| 1920 | 73.70 | 0.24 | **73.94** | 20.83 |
| 1930 | 74.73 | 0.44 | **75.17** | 24.23 |
| 1940 | 80.24 | 0.75 | **80.99** | 29.99 |
| 1950 | 85.56 | 1.24 | **86.80** | 29.66 |
| 1960 | 95.07 | 3.55 | **98.63** | 31.46 |
| 1970 | 110.54 | 6.18 | **116.72** | 37.89 |
| 1980 | 108.60 | 8.98 | **117.58** | 35.20 |
| 1990 | 103.70 | 9.32 | **113.03** | 34.73 |
| 2000 | 121.07 | 12.18 | **133.25** | 39.09 |

**Key metrics:**
- Lake-only increase (1900→2000): **+72.05 Tg C/yr** (×2.47 fold change)
- Total (lake + reservoir) in 2000: **133.25 Tg C/yr**
- Reservoir share in 2000: **9.1%**

> The Q-Q plot shows heavy right tails consistent with Shapiro-Wilk rejection of normality. The fan shape in residuals vs. fitted suggests heteroscedasticity, consistent with the Levene test result.

---

## 5. Statistical Analysis

*All trend analyses use the reliable subset (n = 121, 11 biomes). Total burial calculations use all 13 biomes.*

---

### 6.1 Global Linear Regression: OCAR ~ Decade

> **Note:** This model pools all biomes, treating each biome-decade as independent. The low R² is expected — biome identity explains most variance. See [§6.7](#67-mixed-effects-model-ocar--decade--decade--biome) for the appropriate mixed-effects model.
```
Call:
lm(formula = OCAR ~ Decade_num, data = gc_reliable)

Residuals:
   Min     1Q Median     3Q    Max 
-40.90 -20.11  -8.99  11.03 109.89 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)   
(Intercept) -501.6844   161.4887  -3.107  0.00237 **
Decade_num     0.2728     0.0828   3.295  0.00130 **

Residual standard error: 28.8 on 119 df
Multiple R-squared: 0.0836 | Adjusted R-squared: 0.0759
F-statistic: 10.86 on 1 and 119 DF | p-value: 0.001298
```

---

### 6.2 Linear Regression by Biome

Biome-specific trends, sorted by slope (g C m⁻² yr⁻¹ per decade):

| Rank | Biome | Slope | R² | p-value | Significant? |
|---|---|---|---|---|---|
| 1 | Mangroves | 0.606 | 0.365 | 0.049 | Yes |
| 2 | Tropical Moist and Dry Forest | 0.599 | 0.849 | < 0.001 | Yes |
| 3 | Temperate Grasslands | 0.538 | 0.953 | < 0.001 | Yes |
| 4 | Mediterranean Forest | 0.391 | 0.840 | < 0.001 | Yes |
| 5 | Tropical Grasslands | 0.348 | 0.351 | 0.055 | No |
| 6 | Temperate Mixed Forest | 0.297 | 0.945 | < 0.001 | Yes |
| 7 | Boreal Forest and Taiga | 0.093 | 0.676 | 0.002 | Yes |
| 8 | Montane Grasslands | 0.062 | 0.255 | 0.113 | No |
| 9 | Deserts and Xeric Scrub | 0.057 | 0.144 | 0.249 | No |
| 10 | Tundra | 0.026 | 0.522 | 0.012 | Yes |
| 11 | Temperate Conifer Forest | −0.015 | 0.083 | 0.391 | No |

---

### 6.3 ANOVA: OCAR Differences Among Biomes

#### Assumption Checks

| Test | Statistic | p-value | Result |
|---|---|---|---|
| Shapiro-Wilk (normality) | W = 0.9209 | p = 2.49 × 10⁻⁶ |  Normality violated |
| Levene's test (homogeneity) | F = 7.88 | p = 1.89 × 10⁻⁹ |  Equal variances violated |

→ Using **Welch's ANOVA + Games-Howell** post-hoc (both robust to unequal variances).

#### Standard ANOVA *(reference only; assumptions violated)*
```
             Df Sum Sq Mean Sq F value Pr(>F)    
Biome        10  81269    8127   33.78 <2e-16 ***
Residuals   110  26464     241                   
```

#### Welch's ANOVA
```
One-way analysis of means (not assuming equal variances)

F = 42.618 | num df = 10.000 | denom df = 42.532 | p-value < 2.2e-16
```

#### Effect Size

| Parameter | η² | 95% CI |
|---|---|---|
| Biome | **0.75** | [0.68, 1.00] |

> Biome identity explains **~75.4%** of OCAR variance.

#### Games-Howell Post-Hoc (Top 10 Most Significant Pairs)

*Does not assume equal variances; more consistent with Welch's ANOVA.*

| Comparison | Estimate | 95% CI | p (adj.) | Sig. |
|---|---|---|---|---|
| Tundra vs. Temperate Mixed Forest | −34.13 | [−46.48, −21.77] | 1.48 × 10⁻⁵ | **** |
| Tundra vs. Montane Grasslands | −12.53 | [−17.50, −7.55] | 1.72 × 10⁻⁵ | **** |
| Temperate Mixed Forest vs. Deserts & Xeric Scrub | 28.30 | [15.50, 41.11] | 2.52 × 10⁻⁵ | **** |
| Temperate Mixed Forest vs. Temperate Conifer Forest | 30.40 | [18.03, 42.77] | 3.74 × 10⁻⁵ | **** |
| Tundra vs. Temperate Grasslands | −52.61 | [−74.92, −30.30] | 7.00 × 10⁻⁵ | **** |
| Temperate Mixed Forest vs. Boreal Forest & Taiga | 25.53 | [12.98, 38.08] | 1.11 × 10⁻⁴ | *** |
| Temperate Grasslands vs. Deserts & Xeric Scrub | 46.79 | [24.34, 69.24] | 1.30 × 10⁻⁴ | *** |
| Temperate Grasslands vs. Temperate Conifer Forest | 48.89 | [26.58, 71.20] | 1.30 × 10⁻⁴ | *** |
| Tundra vs. Mangroves | −89.33 | [−129.89, −48.77] | 1.32 × 10⁻⁴ | *** |
| Temperate Conifer Forest vs. Mangroves | −85.61 | [−126.17, −45.04] | 1.90 × 10⁻⁴ | *** |

---

### 6.4 Pre-1950 vs. Post-1950 Comparison

*Two-sided Welch's t-test to avoid directional bias.*
```
Welch Two Sample t-test

t = 2.9348 | df = 116.56 | p-value = 0.004021
95% CI: [4.925, 25.368]
```

| Period | n | Mean OCAR (g C m⁻² yr⁻¹) |
|---|---|---|
| Pre-1950 | 55 | 22.06 |
| Post-1950 | 66 | 37.21 |
| **Difference** | — | **+15.15 (+68.7%)** |

**Effect size:** Cohen's d = **0.52** [0.155, 0.883] → *medium effect*

#### Per-Biome Pre/Post-1950 Comparison *(Bonferroni-corrected, m = 11)*

| Biome | Pre-1950 Mean | Post-1950 Mean | Change (%) | p (raw) | p (Bonf.) | Sig.? |
|---|---|---|---|---|---|---|
| Mediterranean Forest | 8.33 | 27.74 | +233.1% | 0.0114 | 0.126 | No |
| Tropical Grasslands | 14.67 | 40.93 | +178.9% | 0.0161 | 0.177 | No |
| Tropical Moist and Dry Forest | 32.97 | 64.85 | +96.7% | 0.0070 | 0.077 | No |
| Temperate Grasslands | 40.42 | 69.31 | +71.5% | 0.0018 | **0.019** | Yes |
| Tundra | 2.71 | 4.29 | +58.3% | 0.0191 | 0.210 | No |
| Temperate Mixed Forest | 28.61 | 45.26 | +58.2% | 0.0034 | **0.038** | Yes |
| Boreal Forest and Taiga | 9.52 | 14.38 | +51.1% | 0.0218 | 0.240 | No |
| Mangroves | 76.42 | 106.63 | +39.5% | 0.161 | 1.000 | No |
| Montane Grasslands | 13.27 | 18.45 | +39.1% | 0.0213 | 0.234 | No |
| Deserts and Xeric Scrub | 8.16 | 10.42 | +27.8% | 0.503 | 1.000 | No |
| Temperate Conifer Forest | 7.60 | 7.03 | −7.5% | 0.605 | 1.000 | No |

> After Bonferroni correction (α = 0.05/11 = 0.0045), fewer biomes remain significant, probably because of multiple testing burden.

---

### 6.5 Tropical vs. Boreal Interaction Models

*Testing whether the OCAR temporal slope differs between regions.*

#### Model A: Fixed-Effects Interaction (`OCAR ~ Decade * Region`)
```
Residuals:
    Min      1Q  Median      3Q     Max 
-61.447  -8.903  -3.609   5.216  80.855 

Coefficients:
                            Estimate Std. Error t value Pr(>|t|)  
(Intercept)               -108.880    341.620   -0.319   0.7512  
Decade_num                   0.060      0.175    0.342   0.7339  
RegionTropical            -842.881    441.030   -1.911   0.0616 .
Decade_num:RegionTropical    0.458      0.226    2.024   0.0482 *

Residual standard error: 25.98 on 51 df
Multiple R-squared: 0.5455 | Adjusted R-squared: 0.5188
F-statistic: 20.41 on 3 and 51 DF | p-value: 8.009e-09
```

#### Model B: Mixed-Effects Interaction (`OCAR ~ Decade * Region + (1|Biome)`)

*Accounts for biomes nested within regions via random intercept.*
```
REML criterion at convergence: 458.3

Random effects:
 Groups   Name        Variance  Std.Dev.
 Biome    (Intercept)  697.5    26.41   
 Residual              223.7    14.96   
Number of obs: 55 | Groups (Biome): 5

Fixed effects:
                            Estimate Std. Error t value
(Intercept)               -108.880    197.555   -0.551
Decade_num                   0.060      0.101    0.594
RegionTropical            -842.881    255.043   -3.305
Decade_num:RegionTropical    0.458      0.130    3.515
```

> The fixed-effects interaction p = 0.048 (borderline). The mixed model accounts for within-region biome heterogeneity; |t| = 3.515 for the interaction term (|t| > 2 suggests significance).

---

### 6.6 Correlation: Reservoir Area vs. OCAR

> **Caveat:** Both variables trend upward with time. A positive raw correlation may reflect shared temporal trends, not a causal link.

#### a) Raw Spearman Correlation
```
rho = 0.1673 | S = 245833 | p-value = 0.0666
```

*Not significant at α = 0.05.*

#### b) Detrended Correlation *(linear time trend removed from both variables)*
```
rho = 0.2275 | S = 228070 | p-value = 0.0121
```

→ The detrended correlation **remains significant** (p = 0.012), could be indictating a reservoir–OCAR link beyond shared temporal trends.

---

### 6.7 Mixed-Effects Model: OCAR ~ Decade + (Decade | Biome)

> **Most appropriate model for this dataset.** Accounts for (a) biomes starting at different OCAR baselines (random intercept) and (b) biomes having different temporal slopes (random slope). Unifies the pooled regression (§6.1) and per-biome regressions (§6.2) in a single framework.
```
Formula: OCAR ~ Decade_scaled + (Decade_scaled | Biome)

REML criterion at convergence: 957.2

Random effects:
 Groups   Name          Variance  Std.Dev.  Corr
 Biome    (Intercept)   728.63    26.993        
          Decade_scaled   4.73     2.175    0.97
 Residual               112.43    10.603       
Number of obs: 121 | Groups (Biome): 11

Fixed effects:
              Estimate Std. Error t value
(Intercept)    30.323      8.196   3.700
Decade_scaled   2.728      0.723   3.773
```

**Global slope: +2.728 g C m⁻² yr⁻¹ per decade** — significant after accounting for biome-level variation.

#### Biome-Specific Random Slopes

| Biome | Deviation from Global Slope | Total Slope |
|---|---|---|
| Mangroves | +4.455 | **7.183** |
| Temperate Grasslands | +2.157 | **4.886** |
| Tropical Moist and Dry Forest | +1.971 | **4.699** |
| Temperate Mixed Forest | +0.486 | **3.215** |
| Tropical Grasslands | +0.108 | **2.837** |
| Mediterranean Forest | −0.369 | **2.359** |
| Montane Grasslands | −1.348 | **1.380** |
| Boreal Forest and Taiga | −1.497 | **1.231** |
| Deserts and Xeric Scrub | −1.750 | **0.978** |
| Temperate Conifer Forest | −2.050 | **0.678** |
| Tundra | −2.163 | **0.565** |

*All slopes in g C m⁻² yr⁻¹ per decade.*

---

### 6.8 Reservoir Contribution to Total C Burial

*Same OCAR values applied to reservoir areas; cf. Anderson et al. (2020), who note that including reservoirs raises the CO₂ emission offset from ~20% to ~30%.*

| Decade | Lake (Tg/yr) | Reservoir (Tg/yr) | Total (Tg/yr) | Reservoir Share |
|---|---|---|---|---|
| 1900 | 49.02 | 0.054 | 49.07 | 0.1% |
| 1910 | 58.26 | 0.094 | 58.36 | 0.2% |
| 1920 | 73.70 | 0.235 | 73.94 | 0.3% |
| 1930 | 74.73 | 0.440 | 75.17 | 0.6% |
| 1940 | 80.24 | 0.747 | 80.99 | 0.9% |
| 1950 | 85.56 | 1.237 | 86.80 | 1.4% |
| 1960 | 95.07 | 3.554 | 98.63 | 3.6% |
| 1970 | 110.54 | 6.182 | 116.72 | 5.3% |
| 1980 | 108.60 | 8.976 | 117.58 | 7.6% |
| 1990 | 103.70 | 9.322 | 113.03 | 8.2% |
| 2000 | 121.07 | 12.176 | 133.25 | **9.1%** |

Reservoir burial grew from **0.05 Tg C/yr** (1900) to **12.18 Tg C/yr** (2000) — a ~226-fold increase.

---

### 6.9 Sensitivity: Effect of Excluding Suspect Biomes

| Decade | All 13 Biomes (Tg/yr) | Reliable 11 (Tg/yr) | Difference (Tg/yr) | Suspect Share |
|---|---|---|---|---|
| 1900 | 49.02 | 37.46 | 11.56 | 23.6% |
| 1910 | 58.26 | 46.71 | 11.56 | 19.8% |
| 1920 | 73.70 | 62.15 | 11.56 | 15.7% |
| 1930 | 74.73 | 63.18 | 11.56 | 15.5% |
| 1940 | 80.24 | 68.69 | 11.56 | 14.4% |
| 1950 | 85.56 | 74.00 | 11.56 | 13.5% |
| 1960 | 95.07 | 83.52 | 11.56 | 12.2% |
| 1970 | 110.54 | 98.99 | 11.56 | 10.5% |
| 1980 | 108.60 | 97.05 | 11.56 | 10.6% |
| 1990 | 103.70 | 92.15 | 11.56 | 11.1% |
| 2000 | 121.07 | 109.52 | 11.56 | 9.5% |

The 2 suspect biomes contribute a **roughly constant ~14.2%** of total burial. Because their OCAR is constant, temporal trends are entirely unaffected. All trend-related conclusions are robust to their inclusion or exclusion.

---

## 6. Key Findings Summary

**Dataset:** 143 biome-decade means (13 biomes × 11 decades) derived from 516 ²¹⁰Pb-dated lake sediment cores. 2 biomes flagged as suspect (zero temporal variance).

| # | Finding | Key Statistics |
|---|---|---|
| **1** | **OCAR Distribution** right-skewed, dominated by low-burial biomes | Mean = 28.84, Median = 17.85, Range: 1.86–145.67 g C m⁻² yr⁻¹ |
| **2** | **Total C Burial** increased 2.5-fold over the 20th century | 49 Tg/yr (1900) → 121 Tg/yr (2000); 133 Tg/yr incl. reservoirs |
| **3** | **Biomes differ significantly** in OCAR | Welch's F = 42.6, p < 2.2 × 10⁻¹⁶, η² = 0.754 |
| **4** | **Anthropocene signal:** Post-1950 OCAR is ~69% higher than Pre-1950 | t = 2.935, p = 0.004, Cohen's d = 0.52 (medium effect) |
| **5** | **Tropical vs. Boreal divergence** tropical slopes are ~8.6× steeper | Tropical: 0.518/decade vs. Boreal: 0.06/decade; interaction p = 0.048 |
| **6** | **Mixed-effects model** confirms significant global upward trend | Global slope: +2.728 g C m⁻² yr⁻¹ per decade |
| **7** | **Reservoir–OCAR correlation** persists after detrending | Raw ρ = 0.167 (p = 0.067); detrended ρ = 0.228 (p = 0.012) |

---

*Generated from R console output and partially formatted by Gemini 3. Source data: [Anderson et al. (2020)](https://doi.org/10.1126/sciadv.aaw2145) via [Mendeley Data](https://data.mendeley.com/datasets/34hsd2jygc/1).*
