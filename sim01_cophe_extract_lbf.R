args <- commandArgs(trailingOnly = TRUE)
run <- as.numeric(args[1])

##########################################################################
library(cophescan)
library(dplyr)
library(magrittr)

### summary statistics files available in figshare
load("data/sim/simulated_cophescan/files_for_sim.RData")

### function to run cophe.susie and catch any errors thrown by susie and
### run cophe.single instead . Note: If there are no signals =with SuSIE
### cophescan switches to the ABF method aauomatically
ret_cophe_susie <- function(d, LD, querytrait) {
  d$LD <- LD
  d$N <- 20000
  cophe_susie_results <- tryCatch(
    expr = {
        cophe.susie(d, d$cv, querytrait)
    },
    error = function(e) {
        cophe.single(d, d$cv, querytrait)
    }
)
  return(cophe_susie_results)
}

## function to extract bayes factors
get_bf <- function(dataset, run) {

  print(paste("run:", run))
  (load(dataset))

  ## str(datasets[[1]])
  allsnps <- colnames(datasets[[1]]$beta)

  if (length(allsnps) < 1000) {
      return(NULL)
  }

  single_m <- list()
  susie_m <- list()
  beta_m <- list()
  vbeta_m <- list()
  cv <- list()
  snps <- list()

  for (j in seq_along(datasets)) {
    # reorder so first snp (column) is causal for Hc in each set
    cophe_dat <- lapply(seq_len(nrow(datasets[[j]]$beta)), function(x) {
        list(
            beta = datasets[[j]]$beta[x, ],
            varbeta = datasets[[j]]$vbeta[x, ], 
            snp = colnames(datasets[[j]]$beta), type = "cc",
            cv = datasets[[j]]$cv$Hc[1], 
            allcv = bind_rows(lapply(
                datasets[[j]]$cv,
                function(x) paste(x, collapse = ",")
            ))
        )
    })
    names(cophe_dat) <- paste(rownames(datasets[[j]]$beta),
        j, run, sep = "_")

    ## Run cophe.single
    single_m[[j]] <- multitrait.simplify(lapply(
        seq_along(cophe_dat),
        function(x) {
            cophe.single(cophe_dat[[x]],
                querysnpid = datasets[[j]]$cv$Hc[1],
                querytrait = names(cophe_dat)[x]
            )
        }
    ))
    ## Run cophe.susie
    susie_m[[j]] <- multitrait.simplify(lapply(
        seq_along(cophe_dat),
        function(x) {
            ret_cophe_susie(cophe_dat[[x]], LD,
                querytrait = names(cophe_dat)[x]
            )
        }
    ))

    beta_m[[j]] <- lapply(cophe_dat, function(x) x$beta)
    vbeta_m[[j]] <- lapply(cophe_dat, function(x) x$varbeta)
    cv[[j]] <- bind_rows(lapply(cophe_dat, function(x) x$allcv))
    snps[[j]] <- lapply(cophe_dat, function(x) x$snp)
  }

  single_df <- bind_rows(single_m)
  susie_df <- bind_rows(susie_m)
  beta_df <- bind_cols(beta_m)
  vbeta_df <- bind_cols(vbeta_m)
  cv_df <- bind_rows(cv)
  snp_df <- bind_cols(snps)

  single_df <- cbind(single_df,
      truth = sub("_.*", "", single_df$querytrait),
      truth_grp = single_df$querytrait
  )

  susie_df <- cbind(susie_df,
      truth = sub("_.*", "", susie_df$querytrait),
      truth_grp = susie_df$querytrait
  )

  cv_df$truth_grp <- single_df$truth_grp
  rownames(cv_df) <- single_df$truth_grp

  result <- list(
      simdata_sing = single_df, simdata_sus = susie_df, cv_df = cv_df,
      snp_df = snp_df, beta_df = beta_df, vbeta_df = vbeta_df
  )
}

simdata <- make_bf(file.path("data/sim/simulated_cophescan/", files_run[run]), run)
flnm <- files_run[run]
## Note version cophescan creates a new column called querytrait
save(simdata, alpha, beta, gamma, flnm, run,
    file = paste0("data/sim/mcmc_dat_eff/sim_data_", run, ".RData")
)

##########################################################################
