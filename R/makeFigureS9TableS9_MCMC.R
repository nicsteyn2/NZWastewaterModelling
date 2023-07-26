
rm(list=ls())

library(tidyverse)
library(GGally)
library(coda)
library(jsonlite)
library(cowplot)

source("R/supportLoadOutputs.R")

# Specify inputs
PERIOD = 5
MODEL_NAME = "final"
N_CHAINS = 8
MIN_ITER = 101
THIN = 1

# Load data and make stuff
fname = paste0(MODEL_NAME, "/", MODEL_NAME, "_period", PERIOD)
df_samples = loadCleanSamples(N_CHAINS, fname, MIN_ITER, THIN) %>% mutate(chain=as.factor(chain)) %>% rename(sigCAR=sigRho)
df_samples_long = df_samples %>% pivot_longer(cols=c(sigCAR, sigR, kc, kw), names_to="param") %>% mutate(param=factor(param, levels=c("sigCAR", "sigR", "kc", "kw")))

mcmc_list <- mcmc.list(
  lapply(split(df_samples, df_samples$chain), function(chain_df) {
    chain_mcmc = mcmc(as.matrix(chain_df[c("sigCAR", "sigR", "kc", "kw")]))
    return(chain_mcmc)
  })
)

# Produce and save trace plot
trace_plot = ggplot(df_samples_long) +
  geom_line(aes(x=iter, y=value, color=chain)) +
  facet_grid("param", scales="free_y") +
  xlab("Iteration") +
  ylab("Parameter value") +
  theme_minimal()
trace_plot
ggsave(paste0("outputs/results/figureMCMC_period", PERIOD, "_trace.png"), trace_plot, width=25, height=10, dpi=300, units="cm")

# Produce and save pairs plot
pairs_plot = ggpairs(df_samples,
                     columns = c("sigCAR", "sigR", "kc", "kw"),
                     aes(color = factor(chain)),
                     upper=list(continuous = wrap("cor", size=3)),
                     lower = list(continuous = wrap("points", size = 1))) +
  theme_minimal()
pairs_plot
ggsave(paste0("outputs/results/figureMCMC_period", PERIOD, "_pairs.png"), pairs_plot, width=25, height=15, dpi=300, units="cm")


# Produce diagnostics table
samplesize = df_samples_long %>% group_by(param) %>% summarise(n=n()) 
ess = effectiveSize(mcmc_list)
gd = gelman.diag(mcmc_list, transform=TRUE)
gd

tblout = data.frame(sample_size = samplesize$n,
                    ess = round(ess),
                    gd=paste0(round(gd$psrf[,"Point est."], digits=2), " (", round(gd$psrf[,"Upper C.I."], digits=2), ")")) %>%
  t()
write.csv(tblout, paste0("outputs/results/tableMCMC_period", PERIOD,".csv"))

