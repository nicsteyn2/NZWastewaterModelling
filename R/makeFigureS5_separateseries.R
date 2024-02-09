

rm(list=ls())

library(tidyverse)
library(cowplot)

source("R/supportLoadOutputs.R")

label_full = "Both series"
label_cases = "Reported cases only"
label_ww = "Wastewater only"

# Set parameters
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")

# Load data
df_full = read.csv(paste0("outputs/fullposterior/final.csv")) %>%
  mutate(date = as.Date(date), Data = label_full) %>%
  filter(date >= st_date, date <= en_date)

df_cases = read.csv(paste0("outputs/fullposterior/casesonly.csv")) %>%
  mutate(date = as.Date(date), Data = label_cases) %>%
  filter(date >= st_date, date <= en_date)

df_ww = read.csv(paste0("outputs/fullposterior/wastewateronly.csv")) %>%
  mutate(date = as.Date(date), Data = label_ww) %>%
  filter(date >= st_date, date <= en_date)

df = rbind(df_full, df_cases, df_ww)

# Plot Rt
plt_cases = ggplot(df %>% filter(Data != label_ww, variable=="Rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Data, alpha=Data)) +
  geom_line(aes(x=date,y=mean,color=Data)) +
  geom_abline(intercept=1, slope=0, linetype="dashed", color="black") +
  theme_bw() +
  ylab("Reproduction number") +
  xlab("") +
  theme(strip.background = element_rect(fill="transparent")) +
  scale_fill_manual(values=c("#3C6E47", "goldenrod")) +
  scale_color_manual(values=c("#3C6E47", "goldenrod")) + 
  scale_alpha_manual(values=c(0.5, 0.3))

# Plot where Data = Label is dark green and Data = cases only is orange
plt_wastewater = ggplot(df %>% filter(Data != label_cases, variable=="Rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Data, alpha=Data)) +
  geom_line(aes(x=date,y=mean,color=Data)) +
  geom_abline(intercept=1, slope=0, linetype="dashed", color="black") +
  theme_bw() +
  ylab("Reproduction number") +
  xlab("Date") +
  theme(strip.background = element_rect(fill="transparent")) +
  scale_fill_manual(values=c("#3C6E47", "darkblue")) +
  scale_color_manual(values=c("#3C6E47", "darkblue")) + 
  scale_alpha_manual(values=c(0.5, 0.2))
plt_wastewater

plt = plot_grid(plt_cases, plt_wastewater, ncol=1, align="v")
plt
ggsave("outputs/figures/figureS5_separateseries.png", plt, width=20, height=12, units="cm", dpi=300)
