
# rm(list=ls())

library(tidyverse)
library(coda)

source("R/supportLoadOutputs.R")

# label = "normal001v1"
label = "final"
outlabel = "final"


summaries = getMCMCTable_AllPeriods(label, 8, 500, 1)

summaries$mean[summaries$var=="kw"] = summaries$mean[summaries$var=="kw"] * 1e6
summaries$lower[summaries$var=="kw"] = summaries$lower[summaries$var=="kw"] * 1e6
summaries$upper[summaries$var=="kw"] = summaries$upper[summaries$var=="kw"] * 1e6

df_summarytable = summaries
df_summarytable[] = lapply(df_summarytable, function(x) if(is.numeric(x)) signif(x, 2) else x)
df_summarytable = df_summarytable %>%
  mutate(value = paste0(mean, " (", lower, ", ", upper, ")")) %>%
  select(var, value, period) %>%
  spread(key=var, value=value)

# df_summarytable = summaries
# df_summarytable = df_summarytable %>%
#   select(var, mean, period) %>%
#   spread(key=var, value=mean)

write.csv(df_summarytable, paste0("outputs/results/table1_params_", outlabel, ".csv"))


