# Code for creating Figure 2
# Code by: Shizhou Ma
# Version Date: October 2025

# This script combines CH4 and CO2 sequestration model outputs and creates
# Figure 2 showing historical and future trajectories across scenarios.
# Required packages to run the code:
# tidyverse V2.0.0, zoo V1.8-14, ggpubr V0.6.0


#Fill in your own file path
setwd_path <- "/Users/shizhouma/Desktop/Nature Submission/Ma et al. R code&data/Ma et al. Input Dataset/"
# CH4 input files
ch4paleo  <- read.csv(paste0(setwd_path, "CH4randomforestpredictpre1984.csv"))
ch4       <- read.csv(paste0(setwd_path, "CH4randomforestpredictcontemporary.csv"))
ch4ssp126 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp126.csv"))
ch4ssp245 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp245.csv"))
ch4ssp370 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp370.csv"))
ch4ssp585 <- read.csv(paste0(setwd_path, "CH4randomforestpredictssp585.csv"))

# CO2 input files
cseqpaleo  <- read.csv(paste0(setwd_path, "CO2randomforestpredictpre1984.csv"))
cseq       <- read.csv(paste0(setwd_path, "CO2randomforestpredictcontemporary.csv"))
cseqssp126 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp126.csv"))
cseqssp245 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp245.csv"))
cseqssp370 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp370.csv"))
cseqssp585 <- read.csv(paste0(setwd_path, "CO2randomforestpredictssp585.csv"))

# Inner join for paleo data (ch4 and cseq)

paleo <- inner_join(ch4paleo, cseqpaleo, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "historical")

# Inner join for historical data (ch4 and cseq)
historical <- inner_join(ch4, cseq, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "historical")

# Inner join for ssp126 data (ch4ssp126 and cseqssp126)
ssp126 <- inner_join(ch4ssp126, cseqssp126, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp126")

# Inner join for ssp245 data (ch4ssp245 and cseqssp245)
ssp245 <- inner_join(ch4ssp245, cseqssp245, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp245")

# Inner join for ssp370 data (ch4ssp370 and cseqssp370)
ssp370 <- inner_join(ch4ssp370, cseqssp370, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp370")

# Inner join for ssp585 data (ch4ssp585 and cseqssp585)
ssp585 <- inner_join(ch4ssp585, cseqssp585, by = c("Wetland.ID", "Year")) %>%
  select(-any_of(c("temporal.x", "temporal.y"))) %>%
  mutate(temporal = "ssp585")

# Combine all joined datasets into one data frame
inputforcbudget <- bind_rows(paleo, historical, 
                             ssp126,
                             ssp245, ssp370, ssp585) %>%
  arrange(Wetland.ID)

sitespred <- read.csv(paste0(setwd_path, "cmodelpredictionsinput.csv"))



#get the AREA info for each wetland
AREA <- sitespred %>% 
  select(Wetland.ID, AREA) %>% 
  distinct()
AREAtotal <- sum(AREA$AREA)
inputforcbudget <- inputforcbudget %>%
  left_join(AREA, by = "Wetland.ID")


# Create annual summary
annual_summary <- inputforcbudget %>%
  group_by(Year, temporal) %>%
  summarize(
    # Sum the fluxes themselves
    total_annualCH4  = sum(mean_annualch4*AREA,  na.rm = TRUE)/ AREAtotal,
    total_cseqmean   = sum(mean_cseq*AREA,   na.rm = TRUE)/ AREAtotal,
    
    # Combine standard errors using Root-Sum-of-Squares (RSS)
    total_annualSE   = sqrt(sum((sd_annualch4*AREA)^2,  na.rm = TRUE))/ AREAtotal,
    total_cseqSE     = sqrt(sum((sd_cseq*AREA)^2,    na.rm = TRUE))/ AREAtotal)


# Compute Regression Parameters by Scenario
# For CH₄:
reg_results_ch4 <- annual_summary %>%
  group_by(temporal) %>%
  summarize(
    ch4_slope = summary(lm(total_annualCH4 ~ Year, data = cur_data()))$coefficients["Year", "Estimate"],
    ch4_p     = summary(lm(total_annualCH4 ~ Year, data = cur_data()))$coefficients["Year", "Pr(>|t|)"]
  )

# For Cseq:
reg_results_cseq <- annual_summary %>%
  group_by(temporal) %>%
  summarize(
    cseq_slope = summary(lm(total_cseqmean ~ Year, data = cur_data()))$coefficients["Year", "Estimate"],
    cseq_p     = summary(lm(total_cseqmean ~ Year, data = cur_data()))$coefficients["Year", "Pr(>|t|)"]
  )

# Create Custom Legend Labels with Regression Info
ch4_labels <- reg_results_ch4 %>%
  mutate(label = paste0(temporal, 
                        " (Slope: ", round(ch4_slope, 2),
                        ", p: ", ifelse(ch4_p < 0.01, "< 0.01", formatC(ch4_p, format = "e", digits = 2)),
                        ")")) %>%
  select(temporal, label)

cseq_labels <- reg_results_cseq %>%
  mutate(label = paste0(temporal, 
                        " (Slope: ", round(cseq_slope, 2),
                        ", p: ", ifelse(cseq_p < 0.01, "< 0.01", formatC(cseq_p, format = "e", digits = 2)),
                        ")")) %>%
  select(temporal, label)

# Define a named vector for colors.
my_colors <- c("historical" = "blue", 
               "ssp126"     = "darkgreen", 
               "ssp245"     = "orange", 
               "ssp370"     = "purple", 
               "ssp585"     = "red")

# Reorder the custom labels to match the order in my_colors
ch4_labels <- ch4_labels[match(names(my_colors), ch4_labels$temporal), ]
cseq_labels <- cseq_labels[match(names(my_colors), cseq_labels$temporal), ]


library(zoo)

# Define labels:
scenario_labels <- c(
  historical = "Historical",
  ssp126     = "SSP1-2.6",
  ssp245     = "SSP2-4.5",
  ssp370     = "SSP3-7.0",
  ssp585     = "SSP5-8.5"
)

# include regression info:
ch4_labels <- reg_results_ch4 %>%
  mutate(label = paste0(
    scenario_labels[temporal],
    " (Slope: ", round(ch4_slope, 2),
    ", p: ", ifelse(ch4_p < 0.01, "< 0.01", formatC(ch4_p, "e", 2)),
    ")"
  )) %>%
  select(temporal, label)

annual_summary_ma <- annual_summary %>%
  # nest by scenario
  group_by(temporal) %>%
  group_modify(~{
    df <- .x %>% arrange(Year)
    if (.y$temporal == "historical") {
      # just 10-yr MA on historical itself
      df %>% mutate(ma_CH4 = rollmean(total_annualCH4, 10, 
                                      fill = NA, align = "right"))
    } else {
      # take last 9 yrs of historical
      hist_tail <- annual_summary %>%
        filter(temporal == "historical") %>%
        arrange(Year) %>%
        tail(9)
      # prepend, compute rollmean, then drop the prepended part
      combined <- bind_rows(hist_tail, df) %>% arrange(Year)
      ma_all   <- rollmean(combined$total_annualCH4, 10, 
                           fill = NA, align = "right")
      df %>% mutate(ma_CH4 = ma_all[-(1:nrow(hist_tail))])
    }
  }) %>%
  ungroup()

# overwrite annual_summary
annual_summary <- annual_summary_ma

# plot for CH4 (Figure 2a)
p_ch4 <- ggplot(annual_summary, aes(x = Year, y = total_annualCH4, color = temporal)) +
  # ribbon
  geom_ribbon(aes(ymin = total_annualCH4 - total_annualSE,
                  ymax = total_annualCH4 + total_annualSE,
                  fill = temporal),
              alpha = 0.2, color = NA) +
  # 2) raw data solid lines underneath
  geom_line(size = 1, alpha = 0.4) +
  # 3) white “halo” dashed line
  geom_line(aes(y = ma_CH4),
            linetype = "solid",
            size     = 5,
            color    = "white",
            alpha    = 0.6) +
  # 4) colored dashed MA on top
  geom_line(aes(y = ma_CH4),
            linetype = "dashed",
            size     = 1.3,
            alpha    = 1) +
  geom_vline(xintercept = 2020,
             linetype    = "dashed",
             size        = 1.2,
             color       = "black") +
  scale_color_manual(values = my_colors, labels = ch4_labels$label) +
  scale_fill_manual(values = my_colors, labels = ch4_labels$label) +
  scale_x_continuous(breaks = seq(1900, 2090, 10),
                     labels = seq(1900, 2090, 10),
                     limits = c(1900, 2090)) +
  scale_y_continuous(breaks = scales::pretty_breaks(5),
                     labels = scales::comma) +
  labs(x     = "Year",
       y     = expression(CH[4]~" (kg CH"[4]~ha^{-1}~")"),
       title = "",
       color = "Scenario",
       fill  = "Scenario") +
  theme_classic2(base_size = 14) +
  theme(
    axis.title   = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18)
  )

p_ch4




#plot for CO2 (Figure 2b)----
annual_summary_ma_cseq <- annual_summary %>%
  group_by(temporal) %>%
  group_modify(~{
    df <- .x %>% arrange(Year)
    if (.y$temporal == "historical") {
      df %>% mutate(ma_cseq = rollmean(total_cseqmean, 10, fill = NA, align = "right"))
    } else {
      # grab the last 9 yrs of historical
      hist_tail <- annual_summary %>%
        filter(temporal == "historical") %>%
        arrange(Year) %>%
        tail(9)
      combined <- bind_rows(hist_tail, df) %>% arrange(Year)
      ma_all <- rollmean(combined$total_cseqmean, 10, fill = NA, align = "right")
      df %>% mutate(ma_cseq = ma_all[-seq_len(nrow(hist_tail))])
    }
  }) %>%
  ungroup()

cseq_labels <- reg_results_cseq %>%
  mutate(label = paste0(
    scenario_labels[temporal], 
    " (Slope: ", round(cseq_slope, 2),
    ", p: ", ifelse(cseq_p < 0.01, "< 0.01",
                    formatC(cseq_p, format="e", digits=2)),
    ")"
  )) %>%
  select(temporal, label) %>%
  # reorder to match your my_colors order:
  slice(match(names(my_colors), temporal))

# 2) Plot with ribbon, faint solid lines, white halo dashes, then colored dashes on top
p_cseq <- ggplot(annual_summary, aes(x = Year, y = total_cseqmean, color = temporal, fill = temporal)) +
  # SE ribbon
  geom_ribbon(aes(
    ymin = total_cseqmean - total_cseqSE,
    ymax = total_cseqmean + total_cseqSE
  ), alpha = 0.2, color = NA) +
  # Raw data solid lines, lightly transparent
  geom_line(size = 1, alpha = 0.4) +
  # 1) White “halo” dashed layer
  geom_line(
    data     = annual_summary_ma_cseq,
    aes(y     = ma_cseq, group = temporal),
    linetype = "dashed",
    size     = 5,
    color    = "white",
    alpha    = 0.6
  ) +
  # 2) Colored dashed MA on top
  geom_line(
    data     = annual_summary_ma_cseq,
    aes(y     = ma_cseq, color = temporal),
    linetype = "dashed",
    size     = 1.2,
    alpha    = 1
  ) +
  geom_vline(xintercept = 2020,
             linetype    = "dashed",
             size        = 1.2,
             color       = "black") +
  # scales & labels
  scale_color_manual(values = my_colors, labels = cseq_labels$label) +
  scale_fill_manual(values  = my_colors, labels = cseq_labels$label) +
  scale_x_continuous(
    breaks = seq(1900, 2090, by = 10),
    labels = seq(1900, 2090, by = 10),
    limits = c(1900, 2090)
  ) +
  # scale_y_continuous(
  #   breaks = scales::pretty_breaks(n = 5),
  #   labels = scales::comma
  # ) +
  scale_y_continuous(
    limits = c(2600, 2800),
    breaks = seq(2600, 2800, length.out = 5),
    labels = scales::comma
  )+
  labs(
    x     = "Year",
    y     = expression(CO[2]~" (kg CO"[2]~ha^{-1}~")"),
    title = "",
    color = "Scenario",
    fill  = "Scenario"
  ) +
  theme_classic2(base_size = 15) +
  theme(
    axis.title   = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18)
  )

p_cseq



