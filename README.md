## CoPheScan paper

Code for simulation experiments in  *[CoPheScan: phenome-wide association studies accounting for linkage disequilibrium](https://doi.org/10.1038/s41467-024-49990-8)* [1].

### Download the simulated data 
Available on figshare: [https://doi.org/10.6084/m9.figshare.24939408](https://doi.org/10.6084/m9.figshare.24939408). 

### Install CoPheScan

From CRAN (stable)
```r
install.packages("cophescan")
```

or from GitHub
```r
if(!require("remotes"))
   install.packages("remotes") # if necessary
remotes::install_github("ichcha-m/cophescan")
```

### Code for analyses of simulated data
The R scripts needed to extract data from the summary statistics, process and run CoPheScan is available in this repository:

 - sim01_cophe_extract_lbf.R 

    Summary statistics simulated using https://github.com/chr1swallace/cophescan-manuscript-sim-summary-data/ and available on figshare: [https://doi.org/10.6084/m9.figshare.24939408](https://doi.org/10.6084/m9.figshare.24939408). 

    Process these files and extract log Bayes factors (SuSIE and ABF) using CoPheScan. 

 - sim02_collate_prepare_data_mcmc.R
 
    Collate data from all files generated using sim01_cophe_extract_lbf.R as input to the hierarchical model. [The files files_for_sim.RData, sim_data.RData and sim_effect.RData (in sim_processed.tar.gz) are available on figshare [https://doi.org/10.6084/m9.figshare.24939408](https://doi.org/10.6084/m9.figshare.24939408).

 - sim03_run_cophe_mcmc_full.R
        
    Run the MCMC hierarchical model.

 - sim04_run_cophe_mcmc_varyhc.R

    Run MCMC hierarchical model on data with varying number of Hc traits.

 - sim05_run_diagnostics.R

    Plot diagnostics for the chains to check for convergence.

 - sim06_collate_results.R
    
    Collate MCMC results.

 - sim07_compute_copheFDR.R
    
    Calculating FDR at different Hc.cutoff.

 - sim08_compare_other_phewas.R
    
    Compare the results of CoPheScan with conventional PheWAS approaches.

 - sim09_figures_paper.R
    
    Figures of simulation results from the paper.

 - utils.R
    
    helper functions

### Reference
[1] Manipur, I., Reales, G., Sul, J.H., Shin, M.K., Longerich, S., Cortes, A. and Wallace, C., 2024. CoPheScan: phenome-wide association studies accounting for linkage disequilibrium. Nature Communications, 15(1), p.5862. [https://doi.org/10.1038/s41467-024-49990-8](https://doi.org/10.1038/s41467-024-49990-8)
