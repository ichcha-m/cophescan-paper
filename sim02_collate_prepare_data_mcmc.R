library(dplyr)
library(cophescan)
library(magrittr)

#########################################################################
### collate the data generated in sim01_cophe_extract_lbf.R
## this file is available in figshare
load("data/sim/simulated_cophescan/files_for_sim.RData")

collate_bf <- function(dataset, run) {
  print(paste("run:", run))
  (load(dataset))
  result <- list(
      single_df = simdata$simdata_sing,
      susie_df = simdata$simdata_sus,
      beta_df = simdata$beta_df, cv_df = simdata$cv_df,
      vbeta_df = simdata$vbeta_df, snp_df = simdata$snp_df
  )
  return(result)
}


simdata_both <- lapply(seq_along(files_run), function(x) {
    collate_bf(dataset =
            paste0("sim_data_", x, ".RData"), x
    )
})

simdata <- lapply(simdata_both, "[[", "single_df") %>%
    do.call("rbind", .) %>%
    subset(., nsnps == 1000)
simdata_sus <- lapply(simdata_both, "[[", "susie_df") %>%
    do.call("rbind", .) %>%
    subset(., nsnps == 1000)
simdata_cv <- lapply(simdata_both, "[[", "cv_df") %>% do.call("rbind", .)

simdata_beta <- as.data.frame(matrix())
simdata_vbeta <- data.frame(matrix())
simdata_snp <- data.frame(matrix())

for (x in seq_along(simdata_both)){
    if (!is.null(simdata_both[[x]]$cv_df)){
        simdata_beta <- cbind(simdata_beta, simdata_both[[x]]$beta_df)
        simdata_vbeta <- cbind(simdata_vbeta, simdata_both[[x]]$vbeta_df)
        simdata_snp <- cbind(simdata_snp, simdata_both[[x]]$snp_df)
    }
    if (x == 1) {
        simdata_beta[, 1] <- NULL
        simdata_vbeta[, 1] <- NULL
        simdata_snp[,  1] <- NULL
    }
}

dat$nsnps <- dat$lBF.Hc <- dat$lBF.Ha <- as.numeric(NA)
dat$truth_grp <- NA
message("gaps to start: ", sum(is.na(dat$nsnps)))

for (h in unique(dat$truth)) {
  # rows of dat that need filling
  wdat <- which(dat$truth == h & is.na(dat$nsnps))
  if (!length(wdat)) # no gaps left for this hypothesis
    next
  wsim <- which(simdata$truth == h) # rows of simdata to fill gaps
  if (length(wsim) > length(wdat)) # can only do length(wdat) gaps
    wsim <- head(wsim, length(wdat))
  dat[wdat, colnames(simdata)] <- simdata[wsim, , drop = FALSE] # fill
}

message("gaps remaining: ", sum(is.na(dat$nsnps)))

simdata_sus <- simdata_sus[which(simdata_sus$truth_grp %in%
    dat$truth_grp), ]
simdata_sus$rg   <- dat$rg[match(simdata_sus$truth_grp, dat$truth_grp)]
simdata_sing <- dat

simdata_beta_p <- simdata_beta[, simdata_sing$truth_grp]
simdata_snp_p <- simdata_snp[, simdata_sing$truth_grp]
simdata_vbeta_p <- simdata_vbeta[, simdata_sing$truth_grp]
simdata_cv_p <- simdata_cv[which(simdata_cv$truth_grp %in%
    simdata_sing$truth_grp), ]

print(table(dat$truth))
print(table(simdata_sus$truth))

colnames(simdata_sing) <- gsub("PP.", "", colnames(simdata_sing))
colnames(simdata_sus) <- gsub("PP.", "", colnames(simdata_sus))

save(files_run, simdata_sus, simdata_sing,
    file = "data/sim/mcmc_dat_eff/sim_data.RData")

save(simdata_beta_p, simdata_vbeta_p, simdata_cv_p, simdata_snp_p,
    file = "data/sim/mcmc_dat_eff/sim_effect.RData")

#########################################################################
### Simulate genetic correlation values

load("data/sim/mcmc_dat_eff/sim_data.RData")
set.seed(1)
simdata_sing$rg[grep("Hn|Ha", simdata_sing$truth)] <-
    pmin(rexp(length(grep("Hn|Ha", simdata_sing$truth)), 10))
set.seed(1)
simdata_sing$rg[grep("Hc", simdata_sing$truth)] <- 
    pmin(rexp(length(grep("Hc", simdata_sing$truth)), 5))
simdata_sing$rg[simdata_sing$rg > 1] <- 1

hist(simdata_sing$rg[grep("Hc", simdata_sing$truth)])
hist(simdata_sing$rg[grep("Hn|Ha", simdata_sing$truth)])

simdata_sus$rg = simdata_sing$rg[match(
    simdata_sus$truth_grp,
    simdata_sing$truth_grp
)]
# hist(simdata_sus$rg[grep("Hc", simdata_sus$truth)])
# hist(simdata_sus$rg[grep("Hn|Ha", simdata_sus$truth)])

### Note to self that originally stored as sim_data_pmin.RData but changed name
### to sim_data.RData in figshare
save(simdata_sing, simdata_sus, file = "data/sim/mcmc_dat_eff/sim_data.RData")
## -------------------------------------------------------------------

###########################################################################
### Prepare data for tests with varying proportion of Hc
load("data/sim/mcmc_dat_eff/sim_data.RData")

hc_sub_list <- c(5000, 4000, 3000, 2000, 1000, 500, 200, 100, 50)
for (hc_sub in hc_sub_list){
  simdata_sing_sub <- simdata_sing
  simdata_sing_sub <- simdata_sing_sub[which(simdata_sing_sub$truth %in%
      c("Hn", "Ha", "Hc")), ]
  idx <- which(simdata_sing_sub$truth == "Hc")
  set.seed(1)
  rmv <- sample(idx, length(idx) - hc_sub)
  simdata_sing_sub <- simdata_sing_sub[-rmv, ]
  print(table(simdata_sing_sub$truth))
  simdata_sus_sub <- simdata_sus[which(simdata_sus$truth_grp %in%
      simdata_sing_sub$truth_grp), ]
  table(simdata_sus_sub$truth)

  save(simdata_sing_sub, simdata_sus_sub,
      file = paste0("data/sim/mcmc_dat_eff_hc_vary/sim_data_hcsub_",
              hc_sub, ".RData"))
}
#############################################################################
### compute the p-values for the QV-QT pairs and save
## make sure lbf matrix matches the effect size matrices
simdata_cv_p <- simdata_cv_p[order(match(
  simdata_cv_p$truth_grp, simdata_sing$truth_grp)), ]
all.equal(colnames(simdata_snp_p), simdata_cv_p$truth_grp)
all.equal(simdata_sing$truth_grp, simdata_cv_p$trusimdatath_grp)
all.equal(simdata_sing$querysnp, simdata_cv_p$Hc)

### check if the datasets are the same
library(cophescan)
d <- list(beta = simdata_beta_p[, 10], varbeta = simdata_vbeta_p[, 10],
  snp = simdata_snp_p[, 10], N = 20000,
  cv = as.character(simdata_cv_p[10, "Hc"]), type = "cc")
c <- cophe.single(d, d$cv)

pval <- vector()
log10pval <- vector()
beta_val <- vector()
vbeta_val <- vector()
ri <- vector()

for (grp in simdata_cv_p$truth_grp){
  row_idx <- which(simdata_snp_p[, grp] %in%
    simdata_cv_p$Hc[simdata_cv_p$truth_grp == grp])
  ri[grp] <- row_idx
  beta_val[grp] <- simdata_beta_p[row_idx, grp]
  vbeta_val[grp] <- simdata_vbeta_p[row_idx, grp]
  log10pval[grp] <- -(pnorm(-abs(simdata_beta_p[row_idx, grp]) /
    sqrt(simdata_vbeta_p[row_idx, grp]), log.p = TRUE) + log(2)) / log(10)
  pval[grp] <-  10^(-log10pval[grp])
}

save(beta_val, vbeta_val, log10pval, pval, file = "sim_data_pmin_eff.RData")

df <- cbind(beta_val, vbeta_val, log10pval, pval)

simdata_sing[, colnames(df)] <- df[match(simdata_sing$truth_grp,
  rownames(df)), ]
simdata_sing[, colnames(df)] <- df[match(simdata_sing$truth_grp,
  rownames(df)), ]
simdata_sing$p.val.adj <- p.adjust(simdata_sing$pval, method = "BH")

simdata_sus[, colnames(df)] <- df[match(simdata_sus$truth_grp, rownames(df)), ]
simdata_sus[, colnames(df)] <- df[match(simdata_sus$truth_grp, rownames(df)), ]
simdata_sus$p.val.adj <- p.adjust(simdata_sus$pval, method = "BH")

save(simdata_sus, simdata_sing, ri, df,
  file = "data/sim/mcmc_dat_eff/sim_data_pmin_withp.RData")
#############################################################################