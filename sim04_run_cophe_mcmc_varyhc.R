##### Running the MCMC model on thearying the number of Hc traits 

cmdargs <- commandArgs(trailingOnly = TRUE)
# type: either 0 or 1
rg_inp <- (as.numeric(cmdargs[1]))
# chain 1
chain <- as.numeric(cmdargs[2])
# hc_sub: one of these: 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000
hc_sub <- as.numeric(cmdargs[3])
# type: either single or susie
type <- cmdargs[4]
#######################################################################
library(cophescan)
library(data.table)

load(paste0(
  "data/sim/mcmc_dat_eff_hc_vary/sim_data_hcsub_",
  hc_sub, ".RData"
))

if (rg_inp %in% c(1, 0)) {
  RG <- as.logical(rg_inp)
} else {
  stop("RG_inp should be either 0 or 1")
}


print(hc_sub)

if (type == "single"){
  df <- simdata_sing_sub
  lbf_mat_out <- "res_sing_cophe_mcmc"
} else {
  df <- simdata_s_sub
  lbf_mat_out <- "res_sus_cophe_mcmc"
}

print(table(df$truth))
nsnps <- df$nsnps

if (RG) {
  covar_vec <- df$rg
  lbf_mat_out <- paste0(lbf_mat_out, "_rg")
} else {
  covar_vec <- rep(1, nrow(df))
  lbf_mat_out <- paste0(lbf_mat_out, "_norg")
}

## nits <- 300000
## thin <- 30

nits <- 1000000
thin <- 100

print(chain)
result <- run_metrop_priors(df,
    covar = RG, nits = nits, thin = thin,
    covar_vec = covar_vec, avg_pik = TRUE, avg_posterior = TRUE
)

save(result, nits, thin, file = paste0(
  "data/sim/mcmc_results/mcmc_res_eff_",
  nits / 1000, "k/", "chain_", chain, "_hc_", hc_sub, "_",
  lbf_mat_out, "_result.RData"
))

#######################################################################
