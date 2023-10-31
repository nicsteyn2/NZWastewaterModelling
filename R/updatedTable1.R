library(tidyverse)
library(coda)

source("R/supportLoadOutputs.R")


getMCMCTable_AllPeriodsUpdates = function(MODEL_NAME, nChains, miniter, thin) {
  
  df_all = data.frame()
  for (period in seq(1,length(MODEL_NAME))) {
    
    print(period)
    
    if (length(MODEL_NAME) == 1) {
      fname = paste0(MODEL_NAME, "/", MODEL_NAME, "_period", period)
    } else {
      fname = paste0(MODEL_NAME[period], "/", MODEL_NAME[period], "_period", period)
    }
    
    # Load data
    tmp = loadMCMCOutputs(nChains, fname)
    df_samples = tmp[[1]]
    rm(tmp)
    
    # Calculate the min-max iter we can use
    minmaxiter_table = df_samples %>% group_by(chain) %>% summarise(maxiter=max(iter))
    minmaxiter = min(minmaxiter_table$maxiter)
    
    # Filter
    keepiters = seq(miniter, minmaxiter)
    keepiters = keepiters[seq(1, length(keepiters), by=thin)]
    df_in = df_samples %>% filter(iter %in% keepiters)
    
    # Make MCMC list
    mcmc_list <- mcmc.list(
      lapply(split(df_in, df_in$chain), function(chain_df) {
        mcmc(as.matrix(chain_df[c("sigRho", "sigR", "kc", "kw")]))
      })
    )
    
    # Calculate qunatiles
    df_tmp = getMCMCTable_OneList(mcmc_list) %>% mutate(period=period)
    
    # Add onto big boy
    df_all = rbind(df_all, df_tmp)
    
  }
  
  return(df_all)
  
}



label = c("final", "final", "final", "final", "final", "updates", "updates")
outlabel = "updates"


summaries = getMCMCTable_AllPeriodsUpdates(label, 8, 500, 1)

summaries$mean[summaries$var=="kw"] = summaries$mean[summaries$var=="kw"] * 1e6
summaries$lower[summaries$var=="kw"] = summaries$lower[summaries$var=="kw"] * 1e6
summaries$upper[summaries$var=="kw"] = summaries$upper[summaries$var=="kw"] * 1e6

df_summarytable = summaries
df_summarytable[] = lapply(df_summarytable, function(x) if(is.numeric(x)) signif(x, 2) else x)
df_summarytable = df_summarytable %>%
  mutate(value = paste0(mean, " (", lower, ", ", upper, ")")) %>%
  select(var, value, period) %>%
  spread(key=var, value=value)


write.csv(df_summarytable, paste0("outputs/updates/parametertable.csv"))
