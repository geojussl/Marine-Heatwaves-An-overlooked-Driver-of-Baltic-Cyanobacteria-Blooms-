# EnviroChange: statistical analysis of MHW and HAB
Justin Lingg-Laham and Jonah van den Bos

# Introduction

EnviroChange: Marine Heatwaves – An overlooked driver of Baltic
cyanobacteria blooms? Part 3: Statistical analysis of marine heatwaves
and cyanobacterial bloom intensity

This R script works on the post-processing, MHW event detection, and
statistical analysis linking marine heatwaves (MHWs) to cyanobacterial
bloom intensity (HABs) in the Baltic Sea.

The primary goal is to assess whether MHWs act as a trigger or as an
amplifying factor for HABs during summer months (June–August) between
2015 and 2025.

Worksteps: 1. Load and filter phycocyanin bloom metrics (monthly
aggregates) (CSV from Code: Part. 1) 2. Load daily SST and
climatological threshold time series (CSV from Code: Part. 2) 3. Detect
marine heatwave events using the heatwaveR implementation (Hobday et
al. 2016). 4. Compare bloom intensity between months with and without
MHW occurrence (Boxplot) 5. Calculate correlation analyses between MHW
duration and bloom intensity (Scatterplot) 6. Generate figure for
Presentation and Understanding.

Key assumptions: - Key assumptions from Part 1 and 2 - Analysis is
restricted to summer months (June–August). - Marine heatwaves are
detected from surface SST only. - Focus is on bloom intensity.

# Preparation

To ensure Reproducibility, it is crucial to understand that when the
document is rendered, the working directory is automatically set to the
folder in which it is located. Therefore, it is necessary to maintain
the attached folder structure! The dataimport takes place via relative
paths to be independent from lokal storage pathways from different
Computersystems

``` r
getwd()
```

    [1] "C:/Users/justi/OneDrive/Desktop/UGM Semester 1/environmental_systems/EnviroChange/korrelation"

## Libraries

The following Libraries are used to provide access to additional
functions:

``` r
library(dplyr)      # Package to manipulate Data (filter, summarise, etc. )
library(tidyverse)  # Meta-package with more Functions (ggplot2 as example)
library(lubridate)  # Package to work with dates 
library(data.table) # Package to work with big datasets
library(heatwaveR)  # Package for work with MHW
library(ggplot2)    # Package to visualize data
```

# 1. Phycocyanin

Load monthly phycocyanin bloom metrics from Part. 1

``` r
phycocyanin <- read.csv("./data/S2_PC_monthly_JJA_2015_2025_roi.csv", sep = ",", header = TRUE)
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
head(phycocyanin, 5)
```

      system.index   PC_max   PC_mean PC_median    PC_p90 bloom_count_approx
    1            0       NA        NA        NA        NA         0.00000000
    2            1 15.07162 6.8207444 6.8207444 6.8207444         0.11680663
    3            2 67.79510 0.9853975 0.9667403 1.2713874         0.04805542
    4            3 28.42011 9.8468171 9.8468171 9.8888413         1.35470263
    5            4 49.20474 0.9970491 0.9970491 0.9970491         0.01085846
      bloom_fraction        end image_count month      start threshold_used year
    1             NA 2015-06-30           0     6 2015-06-01             10 2015
    2    0.116806627 2015-07-31           1     7 2015-07-01             10 2015
    3    0.003003463 2015-08-31          16     8 2015-08-01             10 2015
    4    0.451567545 2016-06-30           3     6 2016-06-01             10 2016
    5    0.010858461 2016-07-31           1     7 2016-07-01             10 2016
                                        .geo
    1 {"type":"MultiPoint","coordinates":[]}
    2 {"type":"MultiPoint","coordinates":[]}
    3 {"type":"MultiPoint","coordinates":[]}
    4 {"type":"MultiPoint","coordinates":[]}
    5 {"type":"MultiPoint","coordinates":[]}

## Phycocyanin-Filter

apply quality filtering to retain reliable bloom observations

``` r
pc_filt <- phycocyanin %>%
  filter(
    image_count >= 3,   # Keep only data with more then three Images per monthly Scene
    bloom_fraction > 0, # exclude months without detected bloom-signal
    !is.na(PC_mean),    # remove records with missing PC-statistics
    !is.na(PC_p90)
  )
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
head(pc_filt, 5)
```

      system.index    PC_max   PC_mean PC_median    PC_p90 bloom_count_approx
    1            2  67.79510 0.9853975 0.9667403 1.2713874         0.04805542
    2            3  28.42011 9.8468171 9.8468171 9.8888413         1.35470263
    3            5  59.38373 0.2007151 0.2007151 0.2245405         0.03368644
    4            6 166.84210 2.1381653 2.0870251 2.7431478         0.58236025
    5            7 176.35062 1.9805666 1.7509694 3.5922519         1.48825511
      bloom_fraction        end image_count month      start threshold_used year
    1    0.003003463 2015-08-31          16     8 2015-08-01             10 2015
    2    0.451567545 2016-06-30           3     6 2016-06-01             10 2016
    3    0.005614407 2016-08-31           6     8 2016-08-01             10 2016
    4    0.011199236 2017-06-30          52     6 2017-06-01             10 2017
    5    0.016911990 2017-07-31          88     7 2017-07-01             10 2017
                                        .geo
    1 {"type":"MultiPoint","coordinates":[]}
    2 {"type":"MultiPoint","coordinates":[]}
    3 {"type":"MultiPoint","coordinates":[]}
    4 {"type":"MultiPoint","coordinates":[]}
    5 {"type":"MultiPoint","coordinates":[]}

# 2. Marine Heatwave

Load daily SST and climatological threshold time series from Part. 2

``` r
mhw_daily <- read.csv("./data/Daily_SST_p90_JJA_2015_2025_roi.csv", sep = ",", header = TRUE)
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
head(mhw_daily, 5)
```

            date      sst threshold   region
    1 2015-06-01 10.57029  12.76826 roi_test
    2 2015-06-02 10.76332  13.76968 roi_test
    3 2015-06-03 10.79149  14.04454 roi_test
    4 2015-06-04 10.83213  14.28565 roi_test
    5 2015-06-05 11.28999  14.62623 roi_test

# 3. MHW Detection (Hobday et al. 2016)

Marine Heatwave-detection using the Package: HeatwaveR

``` r
dat <- tibble(
  t    = as.Date(mhw_daily$date),       # Daily timestamps
  temp = mhw_daily$sst,                 # observed SST in °C
  seas = mhw_daily$threshold            # climatologic threshold (90th percentile)
)
dat <- dat %>%
  rename(thresh = seas) # rename treshold (to match heatwaveR naming convention)

dat <- dat %>%
  mutate(seas = thresh) #Duplicate threshold as seasonal climatology (required by                                detect_event)
colnames(dat) # check on Column-names
```

    [1] "t"      "temp"   "thresh" "seas"  

``` r
## Detect MHW events after Hobday et al. 2016 
## Definition: ≥ 5 consecutive days with SST exceeding the threshold
events <- detect_event(dat, minDuration = 5) 

event_table <- events$event # extract detected events
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
head(event_table, 5)
```

      event_no index_start index_peak index_end duration date_start  date_peak
    1        1          93         96        97        5 2016-06-01 2016-06-04
    2        2         120        123       124        5 2016-06-28 2016-07-01
    3        3         277        277       296       20 2018-06-01 2018-06-01
    4        4         322        339       348       27 2018-07-16 2018-08-02
    5        5         382        393       399       18 2019-06-14 2019-06-25
        date_end intensity_mean intensity_max intensity_var intensity_cumulative
    1 2016-06-05         0.6398        1.0469        0.2981               3.1990
    2 2016-07-02         0.1576        0.3406        0.1283               0.7878
    3 2018-06-20         1.0894        2.0381        0.6706              21.7880
    4 2018-08-11         1.1363        2.6735        0.7137              30.6810
    5 2019-07-01         0.4874        1.3478        0.3755               8.7732
      intensity_mean_relThresh intensity_max_relThresh intensity_var_relThresh
    1                   0.6398                  1.0469                  0.2981
    2                   0.1576                  0.3406                  0.1283
    3                   1.0894                  2.0381                  0.6706
    4                   1.1363                  2.6735                  0.7137
    5                   0.4874                  1.3478                  0.3755
      intensity_cumulative_relThresh intensity_mean_abs intensity_max_abs
    1                         3.1990            14.9287           15.6731
    2                         0.7878            17.8799           18.0956
    3                        21.7880            16.1414           17.4302
    4                        30.6810            21.7679           23.8210
    5                         8.7732            16.9523           18.4002
      intensity_var_abs intensity_cumulative_abs rate_onset rate_decline
    1            0.5624                  74.6433     0.4583       0.7031
    2            0.2253                  89.3995     0.1259       0.2941
    3            0.7424                 322.8283     3.5876       0.1047
    4            1.3507                 587.7332     0.1627       0.2779
    5            1.0440                 305.1415     0.1215       0.2534

``` r
## Derive annual marine heatwave characteristics
yearly <- event_table %>%
  mutate(year = year(date_start)) %>%
  group_by(year) %>%
  summarise(
    n_events = n(),                            # Number of MHW events per year
    mean_duration = mean(duration),            # Mean event duration (days)
    max_intensity = max(intensity_max),        # Max event intensity
    cum_intensity = sum(intensity_cumulative)  # cumulative annual intensity
  )

knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
head(yearly, 5) 
```

    # A tibble: 5 × 5
       year n_events mean_duration max_intensity cum_intensity
      <int>    <int>         <dbl>         <dbl>         <dbl>
    1  2016        2           5            1.05          3.99
    2  2018        2          23.5          2.67         52.5 
    3  2019        2          12            1.35         12.9 
    4  2020        2           7.5          2.31          9.06
    5  2021        2          22.5          2.65         68.0 

# 4. Compare bloom intensity between months with and without MHW occurrence

The Goal is to test whether cyanobacterial bloom intensity differs
between summer months affected by marine heatwaves and those without
MHWs

``` r
# Aggregate marine heatwave characteristics to monthly scale
# (Each month is classified by total MHW duration and intensity)
mhw_monthly <- event_table %>%
  mutate(
    year  = year(date_start),  # Extract year from event start date
    month = month(date_start)  # Extract month from event start date
  ) %>%
  group_by(year, month) %>%
  summarise(
    mhw_days = sum(duration),  # Total number of MHW days per month
    cum_intensity = sum(intensity_cumulative), # Cumulative MHW intensity per month
    .groups = "drop"
  ) %>%
  # classify months based on presence of MHW conditions
  mutate(
    mhw_present = ifelse(mhw_days > 0, "MHW month", "Without MHW")
  )

# Combination: monthly MHW with phycocyanin bloom data
pc_mhw <- pc_filt %>%
  left_join(mhw_monthly, by = c("year", "month")) %>%
  # lable months without deteced MHW 
  mutate(
    mhw_present = ifelse(is.na(mhw_present), "Without MHW", mhw_present)
  )

# ## Set factor order for plotting (left: No MHW, right: MHW)
pc_mhw$mhw_present <- factor(
  pc_mhw$mhw_present,
  levels = c("Without MHW", "MHW month")
)

# Visualization: Boxplot comparison of phycocyanin index between MHW and non-MHW months with ggplot2

boxplot <- ggplot(pc_mhw, aes(x = mhw_present, y = PC_mean))+
  geom_boxplot(
    fill = "#d0e0e1", color = "black", outlier.shape = NA, width = 0.5) +
  
  geom_jitter( # jittered points to avoid overlaping Points
    aes(color = mhw_present), width = 0.15, size = 2.5, alpha = 0.75) +
  
  scale_color_manual(values = c( "Without MHW"   = "#346a96", "MHW month" = "#355e73")) +
  
  labs(x = "", y = "Phycocyanin-Index [monthly mean]", title = "Phycocyanin-Index during summermonths (jun., jul., aug.) in the Baltic Sea \nwith and without marine heatwaves") +
  
  theme_minimal(base_size = 18) +
  
  theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold"), axis.text.x = element_text(face = "bold"))
print(boxplot)
```

![](korrelartion_mhw_pc_files/figure-commonmark/unnamed-chunk-7-1.png)

``` r
ggsave("figures/PC_Summer_MHW_500dpi.png", plot = boxplot, width=30, height=20, units="cm", dpi=500, bg="white")
```

# 5. Correlation: MHW duration and bloom intensity

The goal is to test whether the intensity of cyanobacterial blooms
increases with the duration of marine heatwave conditions within a
month.

For correlation the Spearman rank correlation is used.

``` r
# Monatliche MHW Tage zählen
mhw_monthly_count <- dat %>%
  mutate(is_mhw = temp > thresh) %>%                  # Binary MHW condition per day
  group_by(year = year(t), month = month(t)) %>%      # Aggregate to monthly scale
  summarise(mhw_days = sum(is_mhw), .groups = "drop")# Total number of MHW days per month

# merge monthly MHW duration and filtered phycocyanin bloom metrics
analysis_final <- pc_filt %>%
  left_join(mhw_monthly_count, by = c("year", "month"))

# Correlation: This test evaluates whether longer MHW durations are associated with higher cyanobacterial bloom intensity (PC 90th percentile)
cor.test(analysis_final$mhw_days, analysis_final$PC_p90, method = "spearman")
```


        Spearman's rank correlation rho

    data:  analysis_final$mhw_days and analysis_final$PC_p90
    S = 2374.6, p-value = 0.008497
    alternative hypothesis: true rho is not equal to 0
    sample estimates:
          rho 
    0.4717251 

## Visualize Correlation

``` r
p_final <- ggplot(analysis_final, aes(x = mhw_days, y = PC_p90)) +
  annotate("text", x = 0, y = max(analysis_final$PC_p90, na.rm=TRUE) * 0.78, label = "Summermonths without MHW", angle = 90, color = "#4A6A8A", fontface = "italic", size = 4.5)  +

  geom_smooth(method = "lm", color = "#355e73", fill = "#d0e0e1", alpha = 0.4, linewidth = 1) + # Trendline
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "#346a96", alpha = 0.5)+
  
  geom_point(color = "#355e73", size = 3, alpha = 0.6)+

  labs(title = "Correlation between marine heatwave duration and cyanobacterial bloom Intensity in the Baltic Sea \nduring summermonths (jun., jul., aug.) in research period (2015-2025)", 
    subtitle = "blue bar serves as baseline comparison",
    x = "Duration of MHW conditions [Days per month]",
    y = "Peak phycocyanin concentration\n[Phycocyanin 90th percentile]") +
  
  scale_x_continuous(breaks = seq(0, max(analysis_final$mhw_days, na.rm=TRUE), by = 5)) +

  theme_minimal(base_size = 14) +
  
  theme(plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(color = "grey40", face = "italic"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 18)),
    axis.title.x = element_text(margin = margin(t = 18))
  )+
  annotate("label", x = max(analysis_final$mhw_days, na.rm =TRUE),y = min(analysis_final$PC_p90, na.rm = TRUE), label = "p < 0.05\nn=30", hjust = 1.1, vjust = -0.5, fill = "white", alpha = 0.9)

print(p_final)
```

![](korrelartion_mhw_pc_files/figure-commonmark/unnamed-chunk-9-1.png)

``` r
ggsave("figures/MHW_vs_AlgalIntensity.png", plot = p_final, width=40, height=20, units="cm", dpi=500, bg="white")
```

# 6. Generate figure for Presentation and Understanding

``` r
library(ggplot2)

concept_data <- data.frame(
  mhw_days = c(0, 5, 10, 20, 30),
  bloom_intensity = c(3, 3.2, 4, 6, 9)
)

ggplot(concept_data, aes(x = mhw_days, y = bloom_intensity)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 8,
    ymin = -Inf, ymax = Inf,
    fill = "#d0e0e1", alpha = 0.6
  ) +
  annotate(
    "rect",
    xmin = 8, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "#346a96", alpha = 0.15
  ) +
  geom_line(
    color = "#355e73",
    linewidth = 1.2
  ) +
  geom_point(
    size = 4,
    color = "#355e73"
  ) +
  annotate(
    "text",
    x = 3, y = 9,
    label = "Bloom initiation\n(no or short MHW)",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 20, y = 9,
    label = "Bloom amplification\n(prolonged MHW)",
    size = 5,
    fontface = "bold",
    color = "#355e73"
  ) +
  
  annotate(
    "text",
    x = 15, y = 2,
    label = "Trigger ≠ Amplifier",
    size = 5,
    fontface = "italic"
  ) +
  
  labs(
    title = "Conceptual model: Marine heatwaves amplify cyanobacterial blooms",
    x = "Marine heatwave duration",
    y = "Bloom intensity (phycocyanin)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
```

![](korrelartion_mhw_pc_files/figure-commonmark/unnamed-chunk-10-1.png)
