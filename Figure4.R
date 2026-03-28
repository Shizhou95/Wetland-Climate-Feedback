# Code for creating Figure 4
# Code by: Shizhou Ma
# Version Date: October 2025

# This script creates Figure 4 using the `GHGpertubmerged` dataset generated
# from the radiative forcing analysis script.
# It generates:
#   1. Long-term net instantaneous radiative forcing trajectories
#   2. A zoomed-in view for inset plotting
# Note:
#   - Run the radiative forcing analysis script first to generate/load the
#     `GHGpertubmerged` dataset.
# Required packages to run the code:
# dplyr V1.1.4, ggplot2 V3.5.2

GHGpertubmerged <- GHGpertubmerged %>%
  filter(time >= -3000, time <= 2090) %>% 
  mutate(totalbase_color = "paleo")

# Define colour palette
data_colors <- c(
  paleo   = "black",
  hist    = "#A67C7A",
  ssp126  = "#BCCBE5",
  ssp245  = "#B1D0A9",
  ssp370  = "#C1BCBF",
  ssp585  = "#FF5733"
)

# Define the labels in the legend
scenario_labels <- c(
  paleo   = "Paleo",
  hist    = "Historical",
  ssp126  = "SSP1-2.6",
  ssp245  = "SSP2-4.5",
  ssp370  = "SSP3-7.0",
  ssp585  = "SSP5-8.5"
)

# Scale factor 1e12 for pW; 1e15 for fW
sf <- 1e12

ggplot(GHGpertubmerged, aes(x = time, group = scenario)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = net_RF * sf, color = scenario), size = 0.6) +
  geom_ribbon(aes(
    ymin = (net_RF - net_sd) * sf,
    ymax = (net_RF + net_sd) * sf,
    fill = scenario
  ), alpha = 0.2) +
  geom_line(
    data = filter(GHGpertubmerged, time >= 1901 & time <= 2091),
    aes(y = totalbase_cum * sf, color = totalbase_color),
    linetype = "solid", size = 1
  ) +
  geom_ribbon(
    data = filter(GHGpertubmerged, time >= 1901 & time <= 2090),
    aes(
      ymin = (totalbase_cum - totalbase_sd) * sf,
      ymax = (totalbase_cum + totalbase_sd) * sf,
      fill = totalbase_color
    ),
    alpha = 0.1
  ) +
  geom_line(
    data = filter(GHGpertubmerged, time >= 1901, scenario == "ssp126"),
    aes(y = totalbase_cum * sf, color = totalbase_color),
    linetype = "solid", size = 1
  ) +
  geom_vline(xintercept = -1600, linetype = "dashed", size = 1) +
  
  # Reference your named palettes here:
  scale_color_manual(
    name   = "Scenario",
    values = data_colors,
    labels = scenario_labels
  ) +
  scale_fill_manual(
    name   = "Scenario",
    values = data_colors,
    labels = scenario_labels
  ) +
  
  scale_x_continuous(
    breaks = pretty(GHGpertubmerged$time),
    labels = function(x) {
      ifelse(x < 0,
             paste0(abs(x), " BCE"),
             ifelse(x > 0,
                    paste0(x, " CE"),
                    "0"))
    }
  ) +
  labs(
    x = "Time",
    y = expression(
      "Net instantaneous radiative forcing (pW·m"^{"-2"}*")"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"   
  )


#zoom in version(inset figures)

GHGpertubmerged_limited <- GHGpertubmerged %>%
  filter(time >= 1901, time <= 2090)

GHGpertubmerged_limited <- GHGpertubmerged_limited %>%
  mutate(totalbase_color = "paleo")
ggplot(GHGpertubmerged_limited, aes(x = time, group = scenario)) +
  # Net RF lines & ribbons, scaled
  geom_line(aes(y = net_RF * sf, color = scenario), size = 1) +
  geom_ribbon(aes(
    ymin = (net_RF - net_sd) * sf,
    ymax = (net_RF + net_sd) * sf,
    fill = scenario
  ), alpha = 0.2) +
  
  # Dashed totalbase_cum for hist & ssp126, scaled
  geom_line(
    data = filter(GHGpertubmerged, time >= 1901, time <= 2090, scenario %in% c("hist", "ssp126")),
    aes(y = totalbase_cum * sf, color = totalbase_color),
    linetype = "dashed", size = 1
  ) +
  geom_ribbon(
    data = filter(GHGpertubmerged, time >= 1901, time <= 2090, scenario %in% c("hist", "ssp126")),
    aes(
      ymin = (totalbase_cum - totalbase_sd) * sf,
      ymax = (totalbase_cum + totalbase_sd) * sf,
      fill = totalbase_color
    ),
    alpha = 0.2
  ) +
  
  # Color & fill scales (no legend title)
  scale_color_manual(values = c(
    paleo   = "black",
    hist    = "#A67C7A",
    ssp126  = "#BCCBE5",
    ssp245  = "#B1D0A9",
    ssp370  = "#C1BCBF",
    ssp585  = "#FF5733"
  )) +
  scale_fill_manual(values = c(
    paleo   = "black",
    hist    = "#A67C7A",
    ssp126  = "#BCCBE5",
    ssp245  = "#B1D0A9",
    ssp370  = "#C1BCBF",
    ssp585  = "#FF5733"
  )) +
  
  # Custom x-axis: BCE / CE
  scale_x_continuous(
    breaks = pretty(GHGpertubmerged_limited$time),
    labels = function(x) {
      ifelse(x < 0,
             paste0(abs(x), " BCE"),
             ifelse(x > 0,
                    paste0(x, " CE"),
                    "0"))
    }
  ) +
  
  # Labels with ×10^15 on y-axis
  labs(
    x = "Time",
    y = expression(
      "Net instantaneous radiative forcing (pW·m"^{"-2"}*")"
    )
  ) +
  
  # Larger text and remove legend
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

