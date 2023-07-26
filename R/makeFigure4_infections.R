
rm(list=ls())

library(tidyverse)
library(cowplot)

source("R/supportLoadOutputs.R")

filename = "final"
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")

# Load model results
df1 = read.csv(paste0("outputs/fullposterior/", filename, "_alpha2e9", ".csv")) %>%  mutate(alpha = 2e9)
df2 = read.csv(paste0("outputs/fullposterior/", filename, ".csv")) %>%  mutate(alpha = 3e9)
df3 = read.csv(paste0("outputs/fullposterior/", filename, "_alpha4e9", ".csv")) %>%  mutate(alpha = 4e9)

df_border = read.csv("data/borderWorkerInfections.csv") %>%
  filter(!is.nan(dailyInfectionsPer100K)) %>%
  mutate(incidence=dailyInfectionsPer100K * 7 * 51.5,
         date = as.Date(substr(date, 1, 10), "%d/%m/%Y"),
         car = NaN) %>%
  arrange(date) %>%
  mutate(cumulative=cumsum(incidence))

# Print initial CARs (for interest) and then tidy
CAR0_a = df1 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_b = df2 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_c = df3 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_a = CAR0_a[[1]]
CAR0_b = CAR0_b[[1]]
CAR0_c = CAR0_c[[1]]

# Update with relative CARt
df1_relcar = df1 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_a, lower=lower/CAR0_a, upper=upper/CAR0_a, variable="CARtRelative")
df2_relcar = df2 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_b, lower=lower/CAR0_b, upper=upper/CAR0_b, variable="CARtRelative")
df3_relcar = df3 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_c, lower=lower/CAR0_c, upper=upper/CAR0_c, variable="CARtRelative")

df = rbind(df1, df2, df3, df1_relcar, df2_relcar, df3_relcar) %>%
  mutate(date = as.Date(date), alpha=factor(alpha)) %>%
  filter(date >= st_date, date <= en_date)
rm(df1, df2, df3, df1_relcar, df2_relcar, df3_relcar)


# Make plot
df_plt = df %>%
  filter(variable %in% c("It", "CI", "Rt", "CARt", "CARtRelative")) %>%
  mutate(popline = ifelse(variable=="CI", 5.15e6, NA),
         Rline = ifelse(variable=="Rt", 1, NA),
         variable = case_when(variable=="It" ~ "(a) Infection incidence",
                              variable=="CI" ~ "(b) Cumulative infections",
                              variable=="CARt" ~ "(c) Absolute case ascertainment rate",
                              variable=="CARtRelative" ~ "(d) Relative case ascertainment rate",
                              variable=="Rt" ~ "(e) Instantaneous reproduction number"),
         variable = factor(variable, levels = c("(a) Infection incidence",
                                                "(b) Cumulative infections",
                                                "(c) Absolute case ascertainment rate",
                                                "(d) Relative case ascertainment rate",
                                                "(e) Instantaneous reproduction number")))

df_border_plt = rbind(
  df_border %>% select(date, incidence) %>% rename(value=incidence) %>% mutate(variable="It", value=value/7),
  df_border %>% select(date, cumulative) %>% rename(value=cumulative) %>% mutate(variable="CI"),
  df_border %>% select(date, car) %>% rename(value=car) %>% mutate(variable="CARt")) %>% 
  mutate(variable = case_when(variable=="It" ~ "(a) Infection incidence",
                              variable=="CI" ~ "(b) Cumulative infections",
                              variable=="CARt" ~ "(c) Case ascertainment rate"),
         variable = factor(variable, levels=c("(a) Infection incidence", "(b) Cumulative infections", "(c) Case ascertainment rate")))



plt = ggplot(df_plt) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=alpha), alpha=0.3) +
  geom_line(aes(x=date, y=mean, color=alpha)) +
  geom_point(data=df_border_plt %>% filter(variable %in% c("(a) Infection incidence", "(b) Cumulative infections")), aes(x=date,y=value), size=0.5) +
  geom_line(data=df_border_plt %>% filter(variable=="(b) Cumulative infections"), aes(x=date,y=value)) +
  geom_hline(aes(yintercept=popline), linetype=3) +
  geom_hline(aes(yintercept=Rline), linetype=3) +
  facet_wrap("variable", scales="free_y", ncol=1) +
  theme_bw() +
  theme(legend.position="right", strip.background = element_rect(fill="transparent")) +
  # scale_color_manual(values=c("#FF0000", "#EE82EE", "#800080")) +
  # scale_fill_manual(values=c("#FF0000", "#EE82EE", "#800080")) +
  scale_color_manual(values=c("#EE82EE", "#FF0000", "#800080")) +
  scale_fill_manual(values=c("#EE82EE", "#FF0000", "#800080")) +
  
  xlab("Date") +
  ylab("") +
  labs(color = expression(alpha), fill=expression(alpha))
plt
ggsave("outputs/results/figure4_infections.png", plt, dpi=300, width=25, height=25, units="cm")
