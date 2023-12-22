library(data.table)
library(dplyr)
library(ggpubr)
source("utils.R")

#######################################################################
path <- c(
  "data/sim/mcmc_res_eff_300k/", "data/sim/mcmc_res_eff_300k/",
  "data/sim/mcmc_res_eff_1000k/", "data/sim/mcmc_res_eff_1000k/")
iter <- c(300, 300, 1000, 1000)
path <- paste0("data/sim/mcmc_res_eff_", iter, "k/")
suff <-  c("res_sing_cophe_mcmc_norg_result", "res_sus_cophe_mcmc_norg_result",
           "res_sing_cophe_mcmc_rg_result", "res_sus_cophe_mcmc_rg_result")

#######################################################################
## save all dfs of results with avg.posterior and avg.pik
cop_res_sim <- list()
load("data/sim/mcmc_dat_eff/sim_data.RData")
for (idx in seq_along(suff)){
  load(file.path(path[idx], paste0("chain_1_", suff[idx], ".RData")))
  if (grepl("sing", suff[idx])) {
    res_df <- simdata_sing
  } else {
    res_df <- simdata_sus
  }
  if (idx == 1) {
    cop_res_sim[["sing_fixd"]] <- simdata_sing
    cop_res_sim[["sus_fixd"]] <- simdata_sus
  }
  nm <- gsub("res_|_cophe|_result", "", suff[idx])
  out_df <- res_df
  out_df[,  c("Hn", "Ha", "Hc")] <- result$avg.posterior
  out_df[,  c("pn", "pa", "pc")] <- result$avg.pik
  cop_res_sim[[nm]] <- out_df
}
save(cop_res_sim, file = "data/sim/mcmc_results/mcmc_res_eff_all_df.RData")

#######################################################################
## save all dfs of results with avg.posterior and avg.pik of varying Hc
hc_ls <- c(50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000)
hc_res <- list()
for (hc in hc_ls){
  cop_res_sim_hc <- list()
  load(paste0("data/sim/mcmc_data_hc_vary/sim_data_hcsub_", hc, ".RData"))
  for (idx in seq_along(suff)){
    print(file.path(path[idx], paste0(
        "chain_1_hc_", hc, "_", suff[idx], ".RData")))
    load(file.path(path[idx], paste0(
        "chain_1_hc_",  hc, "_", suff[idx], ".RData")))
    if (grepl("sing", suff[idx])) {
      res_df <- simdata_sing_sub
    } else {
      res_df <- simdata_sus_sub
    }
    if (idx == 1) {
      cop_res_sim_hc[["sing_fixd"]] <- simdata_sing_sub
      cop_res_sim_hc[["sus_fixd"]] <- simdata_sus_sub
    }
    nm <- gsub("res_|_cophe|_result", "", suff[idx])
    out_df <- res_df
    out_df[,  c("Hn", "Ha", "Hc")] <- result$avg.posterior
    out_df[,  c("pn", "pa", "pc")] <- result$avg.pik
    cop_res_sim_hc[[nm]] <- out_df
  }
  hc_res[[as.character(hc)]] <- cop_res_sim_hc
}
save(hc_res, file = "data/sim/mcmc_results/mcmc_res_eff_all_df_hc.RData")

#######################################################################
#### save observed and true counts of different hypotheses
load("data/sim/mcmc_results/mcmc_res_eff_all_df.RData")
processed_res <- return_processed_simres(cop_res_sim)
save(processed_res, file = "data/sim/mcmc_results/processed_sim.RData")

#######################################################################
#### save observed and true counts of different hypotheses for varying Hc
## takes a long time
load("data/sim/mcmc_results/mcmc_res_eff_all_df_hc.RData")
hc_sub_list <- c(5000, 4000, 3000, 2000, 1000, 500, 200, 100, 50)

for (hc_sub in hc_sub_list) {
  processed_res_hc <- return_processed_simres(hc_res[[as.character(hc_sub)]],
    truth_vals = c("Hn", "Ha", "Hc")
  )
  save(processed_res_hc,
    file = paste0(
      "data/sim/mcmc_results/processed_sim_hc_",
      hc_sub, ".RData"))
  rm(processed_res_hc)
}

#######################################################################

##### Prepare data for plotting

### Prepare data for figures
## Data for figure 2a
## (Barplots of observed hypotheses)
load("data/sim/mcmc_results/mcmc_res_eff_all_df.RData")
load(file = "data/sim/mcmc_results/processed_sim.RData")
load("data/sim/mcmc_dat_eff/sim_data.RData")

df <- rbind(processed_res$ret$fixed_single, processed_res$ret$fixed_susie
           , processed_res$ret$hier_single, processed_res$ret$hier_susie,
           processed_res$ret$hier_single_rg,  processed_res$ret$hier_susie_rg)

k <- sapply(seq_len(nrow(df)), function(x) {
  ((df[x, c(4:6)]) /
    (sum((df[x, c(4:6)])))) * 100})

df <- cbind(df, t(k))
colnames(df)[7:9] <- c("Hn", "Ha", "Hc")
df <- df[, c(1:3, 7:9)]
df <- tidyr::gather(df, Hypothesis, value,  -truth, -type, -total_truth)

df$Hypothesis <- factor(df$Hypothesis,
  levels = c("Hn", "Ha", "Hc"))
df$truth <- paste0("True ", df$truth)
df$truth <- factor(df$truth,
  c("True Hn", "True Ha", "True Ha2", "True Hc", "True Hc2"))

## prepare labels for plotting
p_df <- df
p_df$`BF|Covariate` <- gsub("Fxd.Pr[|]|Hier.Pr[|]", "", p_df$type)
p_df$`BF|Covariate` <- gsub(".rg", ".r<sub>g</sub>", p_df$`BF|Covariate`)
p_df$`BF|Covariate`[which(p_df$`BF|Covariate` == "ABF")] <-
  "ABF|No.r<sub>g</sub>"
p_df$`BF|Covariate`[which(p_df$`BF|Covariate` == "SuSIE.BF")] <-
  "SuSIE.BF|No.r<sub>g</sub>"
p_df$`BF|Covariate` <- gsub(
  '"\\bSuSIE.BF\\b"',
  "SuSIE.BF|No.r<sub>g</sub>", p_df$`BF|Covariate`
)
p_df$Priors <- gsub(
  "[|]ABF|[|]No.rg|[|]Hier.Pr|[|]With.rg|[|]|SuSIE.BF",
  "", p_df$type)
p_df$Priors <- ifelse(p_df$Priors == "Fxd.Pr", "Fixed", "Hierarchical")
p_df$Priors <- factor(p_df$Priors, levels = c("Fixed", "Hierarchical"))
p_df$rg <- p_df$`BF|Covariate`
p_df$rg[grep("No", p_df$`BF|Covariate`)] <-  "No.rg"
p_df$rg[grep("With", p_df$`BF|Covariate`)] <-  "With.rg"
p_df$BF[grep("ABF", p_df$`BF|Covariate`)] <- "ABF"
p_df$BF[grep("SuSIE", p_df$`BF|Covariate`)] <- "SuSIE"
p_df$value <- unlist(p_df$value)
p_df$value_label <- round((p_df$value), 2)
p_df$value_label <- paste0(p_df$value_label, "%")
p_df$value_label[which(gsub("2|True ", "", p_df$truth) !=
  p_df$Hypothesis)] <- NA
fig2_df <- p_df
save(fig2_df, file = "data/sim/fig_data/fig2.RData")

## Data for figure 2b
## compare CoPheScan and conventional PheWAS
load("data/sim/mcmc_dat_eff/sim_data_pmin_withp.RData")
df <- cbind(df,
  p.val.adj.bh = p.adjust(df[, "pval"], method = "BH"),
  p.val.adj.bon = p.adjust(df[, "pval"], method = "bonferroni"))
simdata_sing[, colnames(df)] <-
  df[match(simdata_sing$truth_grp, rownames(df)), ]
simdata_sing <- simdata_sing %>%
  mutate(observed_bh = ifelse(p.val.adj.bh < 0.05, "Hc", "notHc")) %>%
  mutate(observed_bon = ifelse(p.val.adj.bon < 0.05, "Hc", "notHc"))
fig2b_bh <- as.data.frame(unclass(table(simdata_sing$truth,
  simdata_sing$observed_bh))) %>%
  mutate(type = "BH") %>%
  tibble::rownames_to_column("truth") %>%
  mutate(truth = paste0("True ", truth))
fig2b_bon <- as.data.frame(unclass(table(simdata_sing$truth,
  simdata_sing$observed_bon))) %>%
  mutate(type = "Bonf") %>%
  tibble::rownames_to_column("truth") %>%
  mutate(truth = paste0("True ", truth))

fig2b <- rbind(fig2b_bh, fig2b_bon)
fig2b <- fig2b %>%
  mutate(total = Hc + notHc) %>%
  mutate(Assoc = (Hc / total) * 100) %>%
  mutate(notAssoc = (notHc / total) * 100)
fig2b_melt <- reshape2::melt(fig2b[, c("type", "truth", "Assoc", "notAssoc")])
fig2b_melt <- fig2b_melt %>%
  mutate(variable = factor(variable, levels = c("notAssoc", "Assoc"))) %>%
  mutate(truth = factor(truth, levels = c("True Hn", "True Ha", "True Ha2",
                                          "True Hc", "True Hc2"))) %>%
  mutate(value_label = round((value), 2)) %>%
  mutate(value_label = paste0(value_label, "%"))
fig2b_melt$value_label[which(gsub("2|True ", "", fig2b_melt$truth) !=
                               fig2b_melt$Hypothesis)] <- NA
fig2b_melt$value_label[(grepl(c("True Ha|True Ha2|True Hn"), fig2b_melt$truth) &
  fig2b_melt$variable == "Assoc")] <- NA
fig2b_melt$value_label[(grepl(c("True Hc|True Hc2"), fig2b_melt$truth) &
  fig2b_melt$variable == "notAssoc")] <- NA

save(fig2b_melt, file = "data/sim/fig_data/fig2b.RData")

#######################################################################
## Data for figure 3a

g <- merge(processed_res$ret_all$hier_susie_full,
  processed_res$ret_all$hier_susie_rg_full,
  all.x = TRUE, by = "truth_grp", sort = FALSE)
g <- merge(g, simdata_sing[, c("truth_grp", "rg")],
  all.x = TRUE,
  by = "truth_grp", all.y = TRUE, sort = FALSE)
df <- data.frame(x = g$Hc.x, y = g$Hc.y, rg = g$rg, truth = g$truth.x)
df$truth <- paste0("True ", df$truth)
df$truth <- factor(df$truth, c("True Hn", "True Ha", "True Ha2",
                               "True Hc", "True Hc2"))
fig3a_df <- df
fig3a_df$hc <- fig3a_df$x < 0.6 & fig3a_df$y >= 0.6
fig3a_df$nohc <- fig3a_df$x < 0.6

table(fig3a_df$truth, fig3a_df$hc)
table(fig3a_df$truth, fig3a_df$alwhc)
table(fig3a_df$truth, fig3a_df$nohc)
save(fig3a_df, file = "data/sim/fig_data/fig3a.RData")

#######################################################################
