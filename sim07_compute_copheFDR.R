source("utils.R")
load("data/sim/mcmc_results/mcmc_res_eff_all_df.RData")

#######################################################################

processed_res <- list()

for (Hc.cutoff in seq(0, 1, by = 0.1)){
  print(Hc.cutoff)
  processed_res[[as.character(Hc.cutoff)]] <-
      return_processed_simres(cop_res_sim, Hc.cutoff = Hc.cutoff)
}
  processed_res[[as.character(0.55)]] <-
      return_processed_simres(cop_res_sim, Hc.cutoff = 0.55)

save(processed_res,
    file = "data/sim/mcmc_results/mcmc_processed_power_gr.RData")

for (nm in names(processed_res)) {
    proc <- processed_res[[nm]]
    save(proc, file =
        paste0("data/sim/mcmc_results/mcmc_processed_power_gr_",
                nm, ".RData"))
}

df_list <- lapply(names(processed_res[[1]]$ret),
    function(o) {
        lapply(processed_res, "[[", "ret") %>%
            lapply("[[", o) %>%
            do.call("rbind", .) %>%
            mutate(Hc.cutoff = sub(".[^.]+$", "", rownames(.)))})
names(df_list) <- names(processed_res[[1]]$ret)

save(df_list,
    file = "data/sim/mcmc_results/mcmc_processed_power_gr_df.RData")


df <- do.call("rbind", df_list)
df <- df %>%
  mutate(true = ifelse(grepl("Hc", truth), "Hc", "notHc")) %>%
  mutate(grp = paste0(true, "_", Hc.cutoff, "_", type))

a <- df %>%
  group_by(grp) %>%
  summarise(across(where(is.numeric), list(sum = sum)))
a[, c("true", "Hc.cutoff", "type")] <-
    do.call(rbind, stringr::str_split(a$grp, "_"))

df_metrics <- calc_metrics(a, grp_nm = "Hc.cutoff")

### plot TPR and FDR for different Hc.cutoff
df_metrics %>% ggplot(aes(Hc.cutoff, TPR)) +
    geom_point(size = 3, color = "grey33") +
    facet_wrap(~type, nrow = 2) +
    theme_bw() %+replace% theme(text = element_text(size = 12)) +
    scale_y_continuous(limits = c(0, 1))

df_metrics %>% ggplot(aes(Hc.cutoff, FDR)) +
    geom_point(size = 3, color = "grey33") +
    facet_wrap(~type, nrow = 2) +
    theme_bw() %+replace% theme(text = element_text(size = 12)) +
    scale_y_continuous(limits = c(0, 1))


#######################################################################
#### FDR for Varying Hc
hc_sub_list <- c(5000, 4000, 3000, 2000, 1000, 500, 200, 100, 50)
hc_perc_list <- round(hc_sub_list / (88048 + 4700 + hc_sub_list), 4) * 100
names(hc_perc_list) <- as.character(hc_sub_list)
hc_d <- list()
for (hc_sub in hc_sub_list){
  load(file = paste0("data/sim/mcmc_results/processed_sim_hc_",
                     hc_sub, ".RData"))
  df <- do.call("rbind", processed_res_hc$ret)
  df <- df %>%
    mutate(true = ifelse(grepl("Hc", truth), "Hc", "notHc")) %>%
    mutate(hcgrp = paste0(true, "_", hc_sub, "_", type))
  a <- df %>%
    group_by(hcgrp) %>%
    summarise(across(where(is.numeric), list(sum = sum)))
  a[, c("true", "hc_sub", "type")] <-
      do.call(rbind, stringr::str_split(a$hcgrp, "_"))
  hc_perc <- hc_perc_list[[as.character(hc_sub)]]
  hc_d[[as.character(hc_sub)]] <- calc_metrics(a, grp_nm = "hc_sub") %>%
    mutate(true.num.Hc = hc_perc)
}
processed_reslist

#### save data for figure 3b
p_df <- do.call("rbind", hc_d)
p_df$`BF|Covariate` <- gsub("Fxd.Pr[|]|Hier.Pr[|]", "", p_df$type)
p_df$`BF|Covariate` <- gsub(".rg", ".r<sub>g</sub>", p_df$`BF|Covariate`)
p_df$`BF|Covariate`[which(p_df$`BF|Covariate` == "ABF")] <-
    "ABF|No.r<sub>g</sub>"
p_df$`BF|Covariate`[which(p_df$`BF|Covariate` == "SuSIE.BF")] <-
    "SuSIE.BF|No.r<sub>g</sub>"
p_df$`BF|Covariate` <- gsub('"\\bSuSIE.BF\\b"', "SuSIE.BF|No.r<sub>g</sub>",
    p_df$`BF|Covariate`)
p_df$Priors <-
    gsub("[|]ABF|[|]No.rg|[|]Hier.Pr|[|]With.rg|[|]|SuSIE.BF", "", p_df$type)
p_df$Priors <- ifelse(p_df$Priors == "Fxd.Pr", "Fixed", "Hierarchical")
p_df$true.perc.Hc <- gsub(".*[[]", "", p_df$true.num.Hc)
p_df$true.perc.Hc <- as.numeric(gsub("[%].*", "", p_df$true.perc.Hc))

fig3b_df_fdr <- p_df
save(fig4b_df_fdr, file = "data/sim/fig_data/fig3b_fdr.RData")
#######################################################################
