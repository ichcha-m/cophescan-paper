cmdargs <- commandArgs(trailingOnly = TRUE)
# type: either 0 or 1
rg_inp <- (as.numeric(cmdargs[1]))
# chain 1
chain <- as.numeric(cmdargs[2])
# type: either single or susie
type <- cmdargs[3]

library(cophescan)
library(data.table)

#######################################################################

load("data/sim/mcmc_dat_eff/sim_data.RData")
if (rg_inp %in% c(1, 0)) {
  RG <- as.logical(rg_inp)
} else {
  stop("RG_inp should be either 0 or 1")
}

if (type == "single"){
  df <- simdata_sing_sub
  lbf_mat_out <- "res_sing_cophe_mcmc"
} else {
  df <- simdata_s_sub
  lbf_mat_out <- "res_sus_cophe_mcmc"
}

#######################################################################

nsnps <- df$nsnps

if (RG){
  covar_vec <- df$rg
  lbf_mat_out <- paste0(lbf_mat_out, "_rg")
} else{
  covar_vec <- rep(1, nrow(df))
  lbf_mat_out <- paste0(lbf_mat_out, "_norg")
}

## nits <- 300000
## thin <- 30

nits <- 1000000
thin <- 100

result <- run_metrop_priors(df,
    covar = RG, nits = nits, thin = thin,
    covar_vec = covar_vec, avg_pik = TRUE, avg_posterior = TRUE
)

### self: original folder mcmc_res_eff_pmin
# out_folder <- paste0("/mcmc_res_eff/", nits / 1000, "k/")
# ifelse(!dir.exists(out_folder), dir.create(out_folder), "Exists!")
save(result, nits, thin, file = paste0(
    "data/sim/mcmc_results/mcmc_res_eff",
    nits / 1000, "k/", "chain_", chain, "_", lbf_mat_out, "_result.RData"
))

#######################################################################
