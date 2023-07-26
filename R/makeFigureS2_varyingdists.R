
rm(list=ls())

library(tidyverse)
library(cowplot)

source("R/supportLoadOutputs.R")

# Set parameters
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")

mainlabel = "reporting"
labela = "Earlier reporting"
labelb = "Standard distribution"
labelc = "Later reporting"

# Load data
df_a = read.csv(paste0("outputs/fullposterior/final_earlyreporting.csv")) %>%
  mutate(date = as.Date(date), Shift = labela) %>%
  filter(date >= st_date, date <= en_date)

df_b = read.csv(paste0("outputs/fullposterior/final.csv")) %>%
  mutate(date = as.Date(date), Shift = labelb) %>%
  filter(date >= st_date, date <= en_date)

df_c = read.csv(paste0("outputs/fullposterior/final_latereporting.csv")) %>%
  mutate(date = as.Date(date), Shift = labelc) %>%
  filter(date >= st_date, date <= en_date)

# If everything is setup properly the rest should just run

df = rbind(df_a, df_b, df_c) %>%
  filter(variable %in% c("Rt", "CARt")) %>%
  mutate(Shift = factor(Shift, levels = c(labela, labelb, labelc)))

df = rbind(df, df %>% filter(variable=="CARt") %>% mutate(variable="RelCARt"))

df$variable_name[df$variable=="Rt"] = "(a) Instantaneous reproduction number"
df$variable_name[df$variable=="CARt"] = "(b) Case ascertainment rate"
df$variable_name[df$variable=="RelCARt"] = "(c) Relative case ascertainment rate"

CAR0 = df$mean[df$variable=="CARt" & df$t==141 & df$Shift == labela]
df$mean[df$variable=="RelCARt" & df$Shift == labela] = df$mean[df$variable=="CARt" & df$Shift == labela]/CAR0
df$lower[df$variable=="RelCARt" & df$Shift == labela] = df$lower[df$variable=="CARt" & df$Shift == labela]/CAR0
df$upper[df$variable=="RelCARt" & df$Shift == labela] = df$upper[df$variable=="CARt" & df$Shift == labela]/CAR0

CAR0 = df$mean[df$variable=="CARt" & df$t==141 & df$Shift == labelb]
df$mean[df$variable=="RelCARt" & df$Shift == labelb] = df$mean[df$variable=="CARt" & df$Shift == labelb]/CAR0
df$lower[df$variable=="RelCARt" & df$Shift == labelb] = df$lower[df$variable=="CARt" & df$Shift == labelb]/CAR0
df$upper[df$variable=="RelCARt" & df$Shift == labelb] = df$upper[df$variable=="CARt" & df$Shift == labelb]/CAR0

CAR0 = df$mean[df$variable=="CARt" & df$t==141 & df$Shift == labelc]
df$mean[df$variable=="RelCARt" & df$Shift == labelc] = df$mean[df$variable=="CARt" & df$Shift == labelc]/CAR0
df$lower[df$variable=="RelCARt" & df$Shift == labelc] = df$lower[df$variable=="CARt" & df$Shift == labelc]/CAR0
df$upper[df$variable=="RelCARt" & df$Shift == labelc] = df$upper[df$variable=="CARt" & df$Shift == labelc]/CAR0


# Produce hidden states plot
plt_hiddenstates = ggplot(df, aes(x=date, y=mean, fill=Shift)) +
  geom_line(aes(color=Shift)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
  facet_wrap("variable_name", scales="free_y", ncol=1) +
  theme_bw() +
  ylab("") +
  xlab("") +
  theme(strip.background = element_rect(fill="transparent"))

plt_hiddenstates
ggsave(paste0("outputs/results/figureS2_", mainlabel, ".png"), plt_hiddenstates, dpi=300, width=20, height=15, units="cm")
