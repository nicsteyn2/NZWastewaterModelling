
# This code is also used to produce Figure SX

rm(list=ls())

library(tidyverse)
library(cowplot)

source("R/supportLoadOutputs.R")

# Set parameters
filename = "final"
filedir = "outputs/fullposterior/"
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")
CARt_basetime = 141 # This is the time-point we compare CARt to to get our relative CARt

# Load data
df = read.csv(paste0(filedir, filename, ".csv")) %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= st_date, date <= en_date)

df_raw = read.csv(paste0(filedir, filename, "_rawdata.csv")) %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= st_date, date <= en_date)

# Process raw data into correct format for plotting
df_rawcases = df_raw %>% select(date, nCases) %>% rename(observed=nCases) %>% mutate(variable="Cases", variable_name="Reported cases")
df_rawww = df_raw %>% select(date, nGenomeCopies) %>% rename(observed=nGenomeCopies) %>% mutate(variable="Wastewater", variable_name="Wastewater")
df_raw = rbind(df_rawcases, df_rawww) %>%
  filter(!is.nan(observed))
rm(df_rawcases, df_rawww)

# Relabel variables to make things pretty
df$variable_name[df$variable=="CARt"] = "(b) Relative case ascertainment rate"
df$variable_name[df$variable=="It"] = "Estimated infection incidence"
df$variable_name[df$variable=="Rt"] = "(a) Instantaneous reproduction number"
df$variable_name[df$variable=="Wastewater"] = "(c) Wastewater"
df$variable_name[df$variable=="Cases"] = "(d) Reported cases"
df$variable_name = factor(df$variable_name,
                          levels=c("(a) Instantaneous reproduction number", "(b) Relative case ascertainment rate", "Estimated infection incidence", "(c) Wastewater", "(d) Reported cases"))

df_raw$variable_name[df_raw$variable=="Cases"] = "(d) Reported cases"
df_raw$variable_name[df_raw$variable=="Wastewater"] = "(c) Wastewater"


# Add on confidence:
df = df %>% mutate(confident_in_results = ifelse(date <= as.Date("2022-04-01"), "NotConfident", "Confident"))

# Scale CAR
CAR0 = df$mean[df$variable=="CARt" & df$t==CARt_basetime]
df$mean[df$variable=="CARt"] = df$mean[df$variable=="CARt"]/CAR0
df$lower[df$variable=="CARt"] = df$lower[df$variable=="CARt"]/CAR0
df$upper[df$variable=="CARt"] = df$upper[df$variable=="CARt"]/CAR0

# Produce hidden states plot
plt_hiddenstates = ggplot(df %>% filter(variable %in% c("CARt", "Rt")), aes(x=date, y=mean, fill=variable)) +
  geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
  facet_wrap("variable_name", scales="free_y", ncol=2) +
  theme_bw() +
  ylab("") +
  xlab("") +
  scale_fill_manual(values=c("CARt"="#9ACD32", "Rt"="#3C6E47")) +
  theme(legend.position="none", strip.background = element_rect(fill="transparent"))

plt_observeddata = ggplot(df %>% filter(variable %in% c("Cases", "Wastewater")), aes(x=date, y=mean)) +
  geom_line(aes(color = variable), alpha = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = variable), alpha=0.5) +
  geom_ribbon(aes(ymin=predictive_lower, ymax=predictive_upper, fill = variable), alpha=0.2) +
  geom_point(data=df_raw, aes(x=date, y=observed), color="black", size=0.1) +
  geom_line(data=df_raw, aes(x=date, y=observed), color="black", lwd=0.1) +
  facet_wrap("variable_name", scales="free_y", ncol=2) +
  scale_color_manual(values = c("Cases" = "goldenrod", "Wastewater" = "darkblue")) +
  scale_fill_manual(values = c("Cases" = "goldenrod", "Wastewater" = "darkblue")) +
  theme_bw() +
  ylab("") +
  xlab("Date") +
  theme(legend.position="none", strip.background = element_rect(fill="transparent"))

# Combine and save
plt = plot_grid(plt_hiddenstates, plt_observeddata, ncol=1, align="v")
plt

ggsave("outputs/results/figure3_results.png", plt, dpi=300, width=25, height=15, units="cm")

