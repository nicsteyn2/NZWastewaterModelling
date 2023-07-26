
rm(list=ls())

library(tidyverse)
library(cowplot)

infToRep = read.csv("outputs/results/misc/infToRepDist.csv") %>% data.frame()
infToShed = read.csv("outputs/results/misc/infToShedDist.csv") %>% data.frame() 
genTime = read.csv("outputs/results/misc/genTimeDist.csv") %>% data.frame()

onsToRep = read.csv("data/dist_esr_onset.csv") %>% data.frame() %>% rename(t=Days) %>% mutate(p=Count/sum(Count))
onsToShed = read.csv("data/dist_ww_shedding.csv") %>% data.frame() %>% rename(t=DaysFromSymptom, p=Shedding)
infToOns = data.frame(p=dweibull(seq(0,14),1.5, 3.2), t=seq(0,14)) %>% mutate(p=p/sum(p))


xstart = -3.6
titlefontsize = 12
axisfontsize = 9
axisticksize = 7

plt1 = ggplot(infToOns, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from infection") +
  ylab("Probability of symptom onset") +
  ggtitle("Incubation period (infection-to-onset)") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )

plt2 = ggplot(onsToShed, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from onset") +
  ylab("Proportion of total viral load shed") +
  ggtitle("Onset-to-shedding distribution") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )

plt3 = ggplot(onsToRep, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from onset") +
  ylab("Probability of reporting") +
  ggtitle("Onset-to-reporting distribution") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )

plt4 = ggplot(genTime, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from infection") +
  ylab("Probability of infection") +
  ggtitle("Generation time distribution") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )

plt5 = ggplot(infToShed, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from infection") +
  ylab("Proportion of total viral load shed") +
  ggtitle("Infection-to-shedding distribution") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )

plt6 = ggplot(infToRep, aes(x=t,y=p)) +
  geom_col() +
  xlab("Days from infection") +
  ylab("Probability of reporting") +
  ggtitle("Infection-to-reporting distribution") +
  theme_bw() +
  xlim(xstart, 20) +
  theme(
    plot.title = element_text(size = titlefontsize), 
    axis.title = element_text(size = axisfontsize), 
    axis.text = element_text(size = axisticksize)
  )


plt = plot_grid(plt1, plt2, plt3, plt4, plt5, plt6, ncol=3, align="v")
plt
ggsave("outputs/results/figureS1_distributions.png", plt, dpi=300, width=25, height=13, units="cm")

