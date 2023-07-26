
rm(list=ls())

library(tidyverse)
library(cowplot)
library(ggforce)

label = "final"

modeldata = read.csv(paste0("outputs/fullposterior/", label, "_rawdata.csv")) %>%
  mutate(date=as.Date(date)) %>%
  filter(date >= as.Date("2022-01-01"), date <= as.Date("2023-03-31"))

mindate = min(modeldata$date)
maxdate = max(modeldata$date)

rawcases = read.csv("data/cases.csv") %>%
  mutate(date=as.Date(date, "%Y-%m-%d")) %>%
  filter(date >= mindate, date <= maxdate)

plt_cases = ggplot(modeldata) +
  geom_col(data=rawcases, aes(x=date, y=nCases), color="goldenrod") +
  geom_line(aes(x=date, y=nCases)) +
  ylab("Reported cases") +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(mindate, maxdate)

plt_ww = ggplot(modeldata %>% filter(!is.na(nGenomeCopies))) +
  geom_line(aes(x=date, y=nGenomeCopies), color="darkblue", lwd=0.25) +
  geom_point(aes(x=date, y=nGenomeCopies), color="darkblue") +
  ylab("Genome copies/person/day") +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(mindate, maxdate)
  

plt_ww_size = ggplot(modeldata %>% filter(!is.na(nGenomeCopies))) +
  geom_col(aes(x=date, y=nPopInCatchment/5.15e6)) +
  xlab("Date") +
  ylim(c(0,1)) +
  ylab("Prop. covered") +
  theme_bw() +
  xlim(mindate, maxdate)
  
plt = plot_grid(plt_cases, plt_ww, plt_ww_size,
                ncol=1,
                rel_heights = c(1,1,0.6),
                align="v")
plt
ggsave("outputs/results/figure1_data.png", plt, dpi=300, width=25, height=15, units="cm")


