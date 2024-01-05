library(coda)
library(bayesplot)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggtext)
library(ggh4x)
library(gridExtra)
library(patchwork)
library(ggthemes)
library(dplyr)
library(cophescan)

#######################################################################
## mcmc diagnostic plots from multiple chains
plot_mcmc_diagnostics <- function(pref_path, suff_path = NULL, title_i = "",
    num_chains = 4, only_trace = FALSE, start = 101) {
    old <- bayesplot_theme_set(theme_minimal())
    color_scheme_set("mix-red-teal-yellow")
    if (!is.null(pref_path)){
      ### results path
      lst <- paste0(pref_path, title_i, "_chain_", 1:num_chains, ".RData")
      } else if (!is.null(suff_path)) {
      lst <- paste0(suff_path, "chain_", 1:num_chains, "_", title_i, ".RData")
    }
    mcmc_cophe <- list()

    for (chain in seq_along(lst)){
        load(lst[chain])
        if (!is.null(suff_path))
        cophe.hier.rcpp <- result
        chain_nm <- paste0("chain_", chain)
        mcmc_cophe[[chain_nm]] <- cbind(
        (cophe.hier.rcpp$ll[start:length(cophe.hier.rcpp$ll)]),
        t(cophe.hier.rcpp$parameters[, start:ncol(cophe.hier.rcpp$parameters)]))
        colnames(mcmc_cophe[[chain_nm]])[1] <- "ll"
        mcmc_cophe[[chain_nm]] <- mcmc( mcmc_cophe[[chain_nm]])
    }
    cm <- mcmc.list(mcmc_cophe)
    g1 <- list()
    g1[["trace"]] <- bayesplot::mcmc_trace(cm,
    facet_args <- list(ncol = 1, labeller = label_parsed)) +
    ggtitle(title_i) +
    theme_minimal() %+replace% theme(text = element_text(size = 13))

    color_scheme_set("blue")

    if (!only_trace) {
        g1[["hist"]] <- bayesplot::mcmc_hist(cm,
        facet_args = list(ncol = 1, strip.position = "left"))  +
        theme_minimal()
        color_scheme_set("red")
        g1[["int"]] <- bayesplot::mcmc_intervals(cm,
        pars = rownames(cophe.hier.rcpp$parameters))  +
        theme_minimal()
        color_scheme_set("teal")
        g1[["pairs"]] <- bayesplot::mcmc_pairs(cm,
            pars = rownames(cophe.hier.rcpp$parameters))
    }
    return(g1)
}

### function for collation of results from different cophescan approacges
return_processed_simres <- function(res_df_list, truth_vals) {
  ret <- vector("list", 6)
  ret_all <- vector("list", 6)
  names(ret) <- ret_names <-  c(
    "fixed_single", "fixed_susie",
    "hier_single", "hier_susie",
    "hier_single_rg", "hier_susie_rg"
  )
  ret_names_all <- names(ret_all) <- paste0(ret_names, "_full")
  res_df_names <- c("sing_fixd", "sus_fixd", 
                    "sing_mcmc_norg", "sus_mcmc_norg",
                    "sing_mcmc_rg", "sus_mcmc_rg")
  out_labels <- c(
    "Fxd.Pr|ABF", "Fxd.Pr|SuSIE.BF",
    "Hier.Pr|ABF|No.rg", "Hier.Pr|SuSIE.BF|No.rg",
    "Hier.Pr|ABF|With.rg", "Hier.Pr|SuSIE.BF|With.rg"
  )

  for (idx in 1:6){
    df <- res_df_list[[res_df_names[idx]]]
    ### prepare data to work with cran version 1.3.2 of cophescan
    ### (Not required if all processing done with the same version)
    df[, c("PP.Hn", "PP.Ha", "PP.Hc")] <- df[, c("Hn", "Ha", "Hc")]
    df$querytrait <- df$truth_grp

    ### Predict labels using cophescan
    df_call <- cophe.hyp.predict(df,
      grouping.vars = c("querysnp", "querytrait"))
    cont_mat <- as.data.frame.matrix(table(
      df_call$truth, df_call$cophe.hyp.call))[truth_vals, c("Hn", "Ha", "Hc")]

    ret[[ret_names[idx]]] <- data.frame(
      type = out_labels[idx], truth = truth_vals,
      total_truth = rowSums(cont_mat), OHn = cont_mat$Hn,
      OHa = cont_mat$Ha, OHc = cont_mat$Hc)

    df_call[, c("PP.Hn", "PP.Ha", "PP.Hc")] <- NULL
    colnames(df_call)[which(colnames(df_call) == "cophe.hyp.call")] <- "call"

    ret_all[[ret_names_all[idx]]] <- df_call
  }

  processed_result <- list(ret = ret, ret_all = ret_all)
  return(processed_result)
}

### calculate prediction metrics
calc_metrics <- function(df_pred, grp_nm) {
  TPR <- df_pred %>%
    filter(true == "Hc") %>%
    mutate(TP = OHc_sum) %>%
    mutate(FN = OHn_sum + OHa_sum) %>%
    mutate(TPR = TP / (TP + FN)) %>%
    mutate(merge = gsub("Hc_", "", grp))
  TNR <- df_pred %>%
    filter(true == "notHc") %>%
    mutate(TN = OHn_sum + OHa_sum) %>%
    mutate(FP = OHc_sum) %>%
    mutate(TNR = TN / (TN + FP)) %>%
    mutate(merge = gsub("notHc_", "", grp))
  metrics_df <- merge(TPR[, c(grp_nm, "type", "TP", "FN", "TPR", "merge")],
    TNR[, c("TN", "FP", "TNR", "merge")],
    by = "merge") %>%
    mutate(FDR = FP / (FP + TP))
  return(metrics_df)
}
#######################################################################
