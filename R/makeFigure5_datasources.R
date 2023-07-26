
rm(list=ls())

library(tidyverse)
library(cowplot)

source("R/supportLoadOutputs.R")

setwd("~/Documents/NZ Wastewater/NZWastewaterModelling")

# Set parameters
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")

# Load model results
df1 = read.csv("outputs/results/modeloutputs/uniform_casesonly.csv") %>% mutate(Data = "Reported cases")
df2 = read.csv("outputs/results/modeloutputs/uniform_wastewateronly.csv") %>% mutate(Data = "Wastewater")
df3 = read.csv("outputs/results/modeloutputs/uniform.csv") %>% mutate(Data = "Full model")
df = rbind(df1, df2, df3)
rm(df1, df2, df3)

df = df %>%
  mutate(date = as.Date(date), Data = factor(Data, levels=c("Full model", "Reported cases", "Wastewater"))) %>%
  filter(date >= st_date, date <= en_date)

df_plt = df %>% filter(variable=="Rt")




# Duplicate "Full model" rows
df_full_model_duplicate <- df_plt %>%
  filter(Data == "Full model") %>%
  mutate(facet = "Wastewater and Full model")

df_plt <- df_plt %>% mutate(facet = case_when(
  Data == "Reported cases" ~ "Reported cases and Full model",
  Data == "Wastewater" ~ "Wastewater and Full model",
  TRUE ~ as.character(Data)
))

df_plt <- bind_rows(df_plt, df_full_model_duplicate)
df_plt$facet[df_plt$facet=="Full model"] = "Reported cases and Full model"





# Plot
plt = ggplot(df_plt %>% filter(facet %in% c("Reported cases and Full model", "Wastewater and Full model")), 
                      aes(x=date, ymin=lower, ymax=upper, fill=Data, color=Data)) +
  geom_ribbon(data = . %>% filter(Data=="Wastewater"), alpha=0.15, lty=3) +
  geom_ribbon(data = . %>% filter(Data=="Reported cases"), alpha=0.4, lty=3) +
  geom_ribbon(data = . %>% filter(Data=="Full model"), alpha=0.4, lty=3) +
  geom_line(aes(y=mean)) +
  geom_hline(yintercept=1, lty=2) +
  facet_wrap(~facet, ncol=1) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.x=element_blank()) +
  ylab("Instantaneous reproduction number") +
  xlab("Date") +
  scale_color_manual(values = c("Full model" = "#3C6E47", "Reported cases" = "goldenrod", "Wastewater" = "darkblue")) +
  scale_fill_manual(values = c("Full model" = "#3C6E47", "Reported cases" = "goldenrod", "Wastewater" = "darkblue"))

plt
ggsave("outputs/results/figure5_datasource.png", plt, dpi=300, width=25, height=13, units="cm")


