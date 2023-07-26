
# Load raw MCMC outputs

loadMCMCOutputs = function(nChains, fname) {
  
  # Setup empty dfs
  df_samples = data.frame()
  df_allests = data.frame()
  df_loglik = data.frame()
  df_accprob = data.frame()
  
  for (ii in seq(1,nChains)) {
    
    # Set filenames
    fname_samples = paste0("outputs/pmmh.nosync/", fname, "_chainno", ii, "_retainedparams.csv")
    fname_allests = paste0("outputs/pmmh.nosync/", fname, "_chainno", ii, "_allparamlogliks.csv")
    fname_logliks = paste0("outputs/pmmh.nosync/", fname, "_chainno", ii, "_logliks.csv")
    fname_accprobs = paste0("outputs/pmmh.nosync/", fname, "_chainno", ii, "_acceptanceprobs.csv")
    
    # Load data
    df_samples_in = read.csv(fname_samples)
    df_allests_in = read.csv(fname_allests)
    df_loglik_in = read.csv(fname_logliks, col.names =c("loglik"))
    df_accprob_in = read.csv(fname_accprobs, col.names =c("accprob"))
    
    # Append info
    df_samples_in = df_samples_in %>% mutate(iter=seq(1,nrow(df_samples_in)), chain=ii) %>% filter(kw!=0)
    df_allests_in = df_allests_in %>% mutate(iter=seq(1,nrow(df_allests_in)), chain=ii) %>% filter(kw!=0)
    df_loglik_in = df_loglik_in %>% mutate(iter=seq(2,nrow(df_loglik_in)+1), chain=ii)
    df_accprob_in = df_accprob_in %>% mutate(iter=seq(2,nrow(df_accprob_in)+1), chain=ii) %>% filter(!is.na(accprob))
    
    # Stack onto main
    df_samples = rbind(df_samples, df_samples_in)
    df_allests = rbind(df_allests, df_allests_in)
    df_loglik = rbind(df_loglik, df_loglik_in)
    df_accprob = rbind(df_accprob, df_accprob_in)
    
  }
  
  return(list(df_samples, df_allests, df_loglik, df_accprob))
  
}


loadCleanSamples = function(nChains, fname, min_iter, thin) {
  
  # Load data
  tmp = loadMCMCOutputs(nChains, fname)
  df_samples = tmp[[1]]
  
  # Get min-max iteration (so all chains have the same number of iterations)
  maxiter = df_samples %>% group_by(chain) %>% summarise(maxiter=max(iter))
  minmaxiter = min(maxiter$maxiter)
  
  # Choose iters to keep
  keep_iters = seq(min_iter, minmaxiter, by=thin)
  df_samples = df_samples %>% filter(iter %in% keep_iters)
  
  return(df_samples)
  
}


loadAcceptanceProbs = function(nChains, fname, min_iter) {
  
  # Load data
  tmp = loadMCMCOutputs(nChains, fname)
  df_accprob = tmp[[4]]
  
  if (typeof(df_accprob$accprob) == "character") {
    df_accprob = df_accprob %>% mutate(accprob = as.numeric(accprob))
  }
  
  df_accprob = df_accprob %>%
    filter(iter>=min_iter) %>%
    group_by(chain) %>%
    summarise(mean_acc=mean(accprob, na.rm=TRUE))
  return(df_accprob)
  
}



getMCMCTable_AllPeriods = function(MODEL_NAME, nChains, miniter, thin) {
  
  df_all = data.frame()
  for (period in seq(1,5)) {
    
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


getMCMCTable_OneList = function(mcmc_list) {
  
  sum = summary(mcmc_list)
  stats_table = data.frame(sum$statistics)
  quant_table = data.frame(sum$quantiles)
  varnames = rownames(stats_table)
  
  df_out = data.frame(var=varnames,
                      mean=stats_table[varnames, "Mean"],
                      median=quant_table[varnames, "X50."],
                      lower=quant_table[varnames, "X2.5."],
                      upper=quant_table[varnames, "X97.5."])
  
  return(df_out)  
  
}

