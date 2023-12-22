### Diagnostics of the hierarchical models, check for convergence of chains
library(dplyr)
source("utils.R")

load("data/sim/mcmc_dat_eff/sim_data.RData")

#######################################################################
### Supplementary Figure 3

nits <- c(300000, 300000, 1000000, 1000000)
path <- paste0("data/sim/mcmc_results/mcmc_res_eff_", nits / 1000, "k/")
suff <-  c("res_sing_cophe_mcmc_norg_result", "res_sus_cophe_mcmc_norg_result",
            "res_sing_cophe_mcmc_rg_result", "res_sus_cophe_mcmc_rg_result")
threeb_lab <- c(
    expression(paste("ABF|No.r")[g]), expression(paste("SuSIE.BF|No.r")[g]),
    expression(paste("ABF|With.r")[g]), expression(paste("SuSIE.BF|With.r")[g])
)
gg_trace <- list()
for (idx in seq_along(suff)) {
    g1 <- plot_mcmc_diagnostics(
        pref_path = NULL, suff_path = path[idx], title_i = suff[idx],
        num_chains = 4, only_trace = TRUE, start = 101
    )
    gg_trace[[idx]] <- g1$trace + ggtitle(threeb_lab[idx])
}

ggarrange(
    plotlist = gg_trace, nrow = 2, ncol = 2,
    common.legend = TRUE, legend = "right"
)

ggsave("figures/final_figures/suppfig3.pdf",
        bg = "white", width = 11.6, height = 11, units = "in")
ggsave("figures/final_figures/suppfig3.png",
    bg = "white", width = 11.6, height = 11, units = "in"
)

#######################################################################
