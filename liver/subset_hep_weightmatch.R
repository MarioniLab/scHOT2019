if (!require(DCARS)) {
  if (!require(usethis)) {
  install.packages("usethis",repos = "http://cran.us.r-project.org")
  }
  if (!require(devtools)) {
  install.packages("devtools",repos = "http://cran.us.r-project.org")
  }
  library(devtools)
  devtools::install_github("shazanfar/DCARS")
}
library(DCARS)
library(parallel)
library(minerva)
library(energy)

#pwd = getwd()
#write(pwd, file = "wd.txt")

load("output/all_wcors.RData")
hep_p_all <- readRDS("output/hep_p_all.Rds")
load("output/exprs_HVG.RData")


subset_values = rep(c(0.9,0.8,0.7,0.6,0.5), each = 10)
names(subset_values) <- make.unique(paste0("subset_", sprintf("%.2f", subset_values)), sep = "_")

subset_scenarios = sapply(subset_values, function(subsetval) {
  ind = sort(sample(ncol(hep_exprs_HVG), ceiling(subsetval*ncol(hep_exprs_HVG)), replace = FALSE))
  return(list(
    ind = ind,
    #data = hep_exprs_HVG[,ind],
    weightMatrix = DCARS::weightMatrix(length(ind), span = 0.25/subsetval))
  )
  }, simplify = FALSE)

parallel = TRUE
ncores = 10

p_all_only = FALSE

# calculate global values for each and take stratified sample of size 100 from 
# each and take the union of these
if (!file.exists("output/subset_global_weightmatch.RData")) {
global_spearman = lapply(subset_scenarios, function(scenario) {
  apply(pairs_hep_all,1,function(x){
  cor(hep_exprs_HVG[x[1],scenario$ind], hep_exprs_HVG[x[2],scenario$ind], method = "spearman")
})
  })
save(global_spearman, file = "output/subset_global_weightmatch.RData")
} else {
  load("output/subset_global_weightmatch.RData")
}


if (!file.exists("output/subset_output_weightmatch.RData")) {
subset_nsample = 20
subset_sampled_ind = unique(c(
  unlist(lapply(global_spearman, function(g) stratifiedSample(g, subset_nsample)))
))
print(length(subset_sampled_ind))

subset_perms = list() # save the permutations
subset_obs = list() # save the obs stats
subset_p_all = list() # save the results tables
subset_global = global_spearman
} else {
  load("output/subset_output_weightmatch.RData")
}

for (scenario_name in names(subset_scenarios)) {
  scenario = subset_scenarios[[scenario_name]]

  if (!p_all_only) { # this redoes p_all calculation for all

  if (scenario_name %in% names(subset_p_all)) next
  
  }  

  W_hep_scenario = scenario$"weightMatrix"
  exprs = hep_exprs_HVG[, scenario$ind]

  if (!p_all_only) {
  subset_obs[[scenario_name]] = DCARSacrossNetwork(
    exprs,
    edgelist = pairs_hep_all,
    W = W_hep_scenario,
    extractTestStatisticOnly = TRUE,
    weightedConcordanceFunction = weightedSpearman,
    weightedConcordanceFunctionW = "vector",
    verbose = TRUE
  )
}

  write(paste0(scenario_name, " obs calculated"), file = "subset_status.txt", append = TRUE)
  
  if (!p_all_only) {
  
  if (parallel) {
    permstats_raw = mclapply(as.list(subset_sampled_ind), function(h) {
      DCARSacrossNetwork(exprs,
                         edgelist = pairs_hep_all[h, , drop = FALSE],
                         W = W_hep_scenario,
                         weightedConcordanceFunction = weightedSpearman,
                         weightedConcordanceFunctionW = "vector",
                         verbose = TRUE,
                         extractPermutationTestStatistics = TRUE,
                         niter = 1000)
    }, mc.cores = ncores)
    permstats = lapply(permstats_raw, unlist)
    names(permstats) <- rownames(pairs_hep_all[subset_sampled_ind,])
  } else {
    permstats_raw = DCARSacrossNetwork(exprs,
                                       edgelist = pairs_hep_all[subset_sampled_ind,],
                                       W = W_hep_scenario,
                                       weightedConcordanceFunction = weightedSpearman,
                                       weightedConcordanceFunctionW = "vector",
                                       verbose = TRUE,
                                       extractPermutationTestStatistics = TRUE,
                                       niter = 1000)
    permstats = unlist(permstats_raw, recursive = FALSE)
    names(permstats) <- rownames(pairs_hep_all[subset_sampled_ind,])
  }
  subset_perms[[scenario_name]] <- permstats

  write(paste0(scenario_name, " perms calculated"), file = "subset_status.txt", append = TRUE)

} else {
  permstats <- subset_perms[[scenario_name]]
}

pairs_split = split(rownames(pairs_hep_all), rep(1:(2*ncores), length.out = nrow(pairs_hep_all)))

if (parallel) {

p_all_raw = mclapply(pairs_split, function(nms) {
  df = estimatePvaluesSpearman(subset_obs[[scenario_name]][nms], 
                                  subset_global[[scenario_name]], 
                                  permstats,
                                  usenperm = TRUE,
                                  nperm = 10000,
                                  plot = FALSE,
                                  maxDist = 0.1,
                                  verbose = TRUE)
  rownames(df) <- nms
  return(df)
  }, mc.cores = 2*ncores)
p_all = do.call(rbind, p_all_raw)
rownames(p_all) <- as.character(p_all$genepair)
p_all <- p_all[rownames(pairs_hep_all),]
} else {
 p_all = estimatePvaluesSpearman(subset_obs[[scenario_name]], 
                                  subset_global[[scenario_name]], 
                                  permstats,
                                  usenperm = TRUE,
                                  nperm = 10000,
                                  plot = FALSE,
                                  maxDist = 0.1,
                                  verbose = TRUE)
 rownames(p_all) <- names(subset_obs[[scenario_name]])
}
  p_all$fdr <- p.adjust(p_all$pval, method = "BH")
  subset_p_all[[scenario_name]] <- p_all
  
  write(paste0(scenario_name, " estimated pvals calculated"), file = "subset_status.txt", append = TRUE)

  # save and move to next one
  save(subset_perms,
       subset_obs,
       subset_p_all,
       subset_global,
       subset_sampled_ind,
       subset_scenarios,
       file = "output/subset_output_weightmatch.RData")

  write(paste0(scenario_name, " output saved"), file = "subset_status.txt", append = TRUE)

}

save(subset_obs,
     subset_p_all,
     subset_global,
     subset_sampled_ind,
     subset_scenarios,
     file = "output/subset_output_weightmatch_noperms.RData")