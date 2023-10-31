
rm(list=ls())

library(tidyverse)
library(GGally)
library(coda)
library(jsonlite)

source("R/supportLoadOutputs.R")

# Specify inputs
PERIOD = 7
MODEL_NAME = "updates"
N_CHAINS = 8
MIN_ITER = 0
THIN = 1

# # Set filenames (use this code when we're not working with periods)
# fname = paste0(MODEL_NAME, "/", MODEL_NAME)
# fname_out = paste0("outputs/pmmh.nosync/", MODEL_NAME, "/", MODEL_NAME, "_diagnostics/")

# Set filenames (use this code when we are working with periods)
fname = paste0(MODEL_NAME, "/", MODEL_NAME, "_period", PERIOD)
fname_out = paste0("outputs/pmmh.nosync/", MODEL_NAME, "/", MODEL_NAME, "_diagnostics/", "period", PERIOD)

# Specify results to produce
PAIRS_PLOTS = TRUE
TRACE_PLOTS = TRUE
ESS = TRUE
ACCEPTANCEPROB = TRUE
GELMAN_DIAG = TRUE
PARAMSUMMARY = TRUE # TO WRITE FULL LATER, THIS IS A LITTLE DIFFERENT AS IT LOOKS AT ALL PERIODS. CODE IN MAKEPARAMSUMMARYTABLE.R
MCMC_HYPERPARAM_SUMMARY = TRUE

# PAIRS_PLOTS = FALSE
# TRACE_PLOTS = FALSE
# ESS = FALSE
# ACCEPTANCEPROB = FALSE
# GELMAN_DIAG = FALSE
# PARAMSUMMARY = TRUE # TO WRITE FULL LATER, THIS IS A LITTLE DIFFERENT AS IT LOOKS AT ALL PERIODS. CODE IN MAKEPARAMSUMMARYTABLE.R
# MCMC_HYPERPARAM_SUMMARY = TRUE


# Load data
df_samples = loadCleanSamples(N_CHAINS, fname, MIN_ITER, THIN)
df_accprob = loadAcceptanceProbs(N_CHAINS, fname, MIN_ITER)
mmhopts = fromJSON(paste0("outputs/pmmh.nosync/", fname, "_mmhopts.json"))


mcmc_list <- mcmc.list(
  lapply(split(df_samples, df_samples$chain), function(chain_df) {
    chain_mcmc = mcmc(as.matrix(chain_df[c("sigRho", "sigR", "kc", "kw")]))
    return(chain_mcmc)
  })
)



produceResults = function() {
  
  print(paste0("Total samples after cleaining: ", nrow(df_samples), ", or ", length(unique(df_samples$iter)), " per chain."))
  
  # Make pairs plot
  if (PAIRS_PLOTS) {
    pairs_plot = ggpairs(df_samples, columns = c("sigRho", "sigR", "kc", "kw"), aes(color = factor(chain))) +
      theme_bw() +
      ggtitle(paste0("Period ", PERIOD))
    pairs_plot
    ggsave(paste0(fname_out, "_pairs.png"), pairs_plot, width=30, height=20, units="cm")
    print(pairs_plot)
    
  }
  
  
  # Make trace plots
  if (TRACE_PLOTS) {
    
    df_samples_nonthin = loadCleanSamples(N_CHAINS, fname, MIN_ITER, 1)
    df_samples_long = df_samples_nonthin %>% pivot_longer(cols=c(sigRho, sigR, kc, kw), names_to="param") #%>% filter(iter %in% keepiters)
    
    plt_trace = ggplot(df_samples_long) +
      geom_line(aes(x=iter, y=value)) +
      facet_grid(param ~ chain, scales="free_y") +
      theme_bw() +
      ggtitle(paste0("Trace plots: period ", PERIOD))
    plt_trace
    ggsave(paste0(fname_out, "_trace.png"), plt_trace, width=30, height=20, units="cm")
    print(plt_trace)
    
  }
  
  # Print ESS
  if (ESS) {
    ess <- effectiveSize(mcmc_list)
    print("-----------------------------------------------------------------")
    print("Effective sample sizes: ")
    print("-----------------------------------------------------------------")
    print(ess)
    print("-----------------------------------------------------------------")
    cat("\n\n\n")
  }
  
  # Print acceptance prob
  if (ACCEPTANCEPROB) {
    df_accprob = loadAcceptanceProbs(N_CHAINS, fname, MIN_ITER)
    print("-----------------------------------------------------------------")
    print("Acceptance probabilities:")
    print("-----------------------------------------------------------------")
    print(df_accprob)
    print("-----------------------------------------------------------------")
    cat("\n\n\n")
  }
  
  # Print gelman diagnostic
  if (GELMAN_DIAG) {
    gd = gelman.diag(mcmc_list, transform=TRUE)
    print("-----------------------------------------------------------------")
    print("Gelman-Rubin convergence diagnostics:")
    print("-----------------------------------------------------------------")
    print(gd)
    print("-----------------------------------------------------------------")
    cat("\n\n\n")
  }
  
  # Print parameter summary
  if (PARAMSUMMARY) {
    print("-----------------------------------------------------------------")
    print("Parameter summary:")
    print("-----------------------------------------------------------------")
    print(summary(mcmc_list))
    print("-----------------------------------------------------------------")
    cat("\n\n\n")
  }
  
  # Print MCMC hyperparameter summaries
  if (MCMC_HYPERPARAM_SUMMARY) {
    s = summary(mcmc_list)$statistics
    hyperparam_table = data.frame(proposal_variances = mmhopts$proposalVariances, row.names = mmhopts$paramNames,
                                  posterior_stddev = s[,2]) %>%
      mutate(theoretical_prop_var = (2.38 * posterior_stddev) / 4,
             relative_prop_var = proposal_variances/theoretical_prop_var) %>%
      select(-posterior_stddev)
    print("-----------------------------------------------------------------")
    print("MCMC hyperparameter summary:")
    print("-----------------------------------------------------------------")
    print(t(hyperparam_table))
  }
  
}


dir.create(paste0("outputs/pmmh.nosync/", MODEL_NAME, "/", MODEL_NAME, "_diagnostics/"))
file.create(paste0(getwd(), "/", fname_out, ".txt"))

sink(paste0(fname_out, ".txt"), split = TRUE)
produceResults()
sink()

