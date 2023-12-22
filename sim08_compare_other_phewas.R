library(dplyr)

### files available from figshare
load("data/sim/mcmc_dat_eff/sim_data.RData")
load("data/sim/mcmc_dat_eff/sim_effect.RData")

#######################################################################
#### Running either Bonferroni or Benjamini Hochberg methods
get_metrics <- function(df_res, ref_df, method = "BH") {

  if (all.equal(df_res$truth_grp, rownames(ref_df)) != TRUE)
    ref_df <- ref_df[ order(match(rownames(ref_df), df_res$truth_grp)), ]
  print(all.equal(df_res$truth_grp, rownames(ref_df)))
  df_res[, colnames(ref_df)] <-
    ref_df[match(df_res$truth_grp, rownames(ref_df)), ]
  df_res <- df_res %>%
    mutate(method = method) %>%
    mutate(p.val.adj = ifelse(method == "BH", p.val.adj.bh, p.val.adj.bon)) %>%
    mutate(true = ifelse(grepl("Hc", truth), "Hc", "notHc")) %>%
    mutate(observed = ifelse(p.val.adj < 0.05, "Hc", "notHc")) %>%
    mutate(for_counts = case_when(
      true == "Hc" & observed == "Hc" ~ "TP",
      true == "Hc" & observed == "notHc" ~ "FN",
      true == "notHc" & observed == "notHc" ~ "TN",
      true == "notHc" & observed == "Hc" ~ "FP"
    ))
  counts_df <- as.data.frame(t(unclass(table(df_res$for_counts)))) %>%
    mutate(TPR = TP / (TP + FN)) %>%
    mutate(TNR = TN / (TN + FP)) %>%
    mutate(FDR = FP / (FP + TP))
  print(counts_df)
  return(list(ret_all = df_res, ret_fdr = counts_df))

}


load("data/sim/mcmc_results/mcmc_res_eff_all_df.RData")
load("data/sim/mcmc_results/processed_sim.RData")
load("data/sim/mcmc_dat_eff/sim_data_pmin_withp.RData")
df <- cbind(df, p.val.adj.bh = p.adjust(df[, "pval"], method = "BH")
           , p.val.adj.bon = p.adjust(df[, "pval"], method = "bonferroni"))
df_list_bh <- lapply(processed_res$ret_all, function(x) {
  get_metrics(x, df, method = "BH")
})

df_list_bon <- lapply(processed_res$ret_all, function(x) {
  get_metrics(x, df, method = "bonferroni")
})

save(df_list_bh, df_list_bon,
  file = "data/sim/mcmc_results/processed_sim_withFDR_bh_bon.RData")
#######################################################################