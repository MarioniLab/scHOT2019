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

# mine(rnorm(100), rnorm(100))$MIC
normMI = function(x,y,w = NULL) {
  mine(x, y)$MIC
}

#dcor(rnorm(100),rnorm(100))

# testing scenarios
test_statistics_scenarios = list(
  "Spearman" = list(type = "weighted",
                    hof = weightedSpearman_matrix,
                    hof_type = "matrix"),
  "Pearson" = list(type = "weighted",
                    hof = weightedPearson_matrix,
                    hof_type = "matrix"),
  "Spearman_block" = list(type = "block",
                    hof = weightedSpearman_matrix,
                    hof_type = "matrix"),
  "Pearson_block" = list(type = "block",
                    hof = weightedPearson_matrix,
                    hof_type = "matrix"),
  "MIC" = list(type = "block",
                    hof = normMI,
                    hof_type = "vector"),
  "BDC" = list(type = "block",
               hof = dcor,
               hof_type = "vector")
)

W_hep_thin = thin(W_hep, 30)

parallel = TRUE
ncores = 10

p_all_only = TRUE

# calculate global values for each and take stratified sample of size 100 from 
# each and take the union of these
if (!file.exists("output/test_statistics_global.RData")) {
global_spearman = apply(pairs_hep_all,1,function(x){
  cor(hep_exprs_HVG[x[1],], hep_exprs_HVG[x[2],], method = "spearman")
})
global_pearson = apply(pairs_hep_all,1,function(x){
  cor(hep_exprs_HVG[x[1],], hep_exprs_HVG[x[2],])
})
global_mic = apply(pairs_hep_all,1,function(x){
  normMI(hep_exprs_HVG[x[1],], hep_exprs_HVG[x[2],])
})
global_bdc = apply(pairs_hep_all,1,function(x){
  dcor(hep_exprs_HVG[x[1],], hep_exprs_HVG[x[2],])
})
save(global_spearman, global_pearson, global_mic, global_bdc,
  file = "output/test_statistics_global.RData")
} else {
  load("output/test_statistics_global.RData")
}


if (!file.exists("output/test_statistics_output.RData")) {
test_statistics_nsample = 350
test_statistics_sampled_ind = unique(c(
  stratifiedSample(global_spearman, test_statistics_nsample),
  stratifiedSample(global_pearson, test_statistics_nsample),
  stratifiedSample(global_mic, test_statistics_nsample),
  stratifiedSample(global_bdc, test_statistics_nsample)
))
length(test_statistics_sampled_ind)

test_statistics_perms = list() # save the permutations
test_statistics_obs = list() # save the obs stats
test_statistics_p_all = list() # save the results tables
test_statistics_global = list(
  "Spearman" = global_spearman,
  "Pearson" = global_pearson,
  "Spearman_block" = global_spearman,
  "Pearson_block" = global_pearson,
  "MIC" = global_mic,
  "BDC" = global_bdc
)
} else {
  load("output/test_statistics_output.RData")
}

for (scenario_name in names(test_statistics_scenarios)) {
  scenario = test_statistics_scenarios[[scenario_name]]

  if (!p_all_only) { # this redoes p_all calculation for all

  if (scenario_name %in% names(test_statistics_p_all)) next
  
  }

  if (scenario$type == "block") {
    W_hep_scenario = 1*(W_hep_thin > 0)
  } else {
    W_hep_scenario = W_hep_thin
  }
  

  if (!p_all_only) {
  test_statistics_obs[[scenario_name]] = DCARSacrossNetwork(
    hep_exprs_HVG,
    edgelist = pairs_hep_all,
    W = W_hep_scenario,
    extractTestStatisticOnly = TRUE,
    weightedConcordanceFunction = scenario$hof,
    weightedConcordanceFunctionW = scenario$hof_type,
    verbose = TRUE
  )
}

  write(paste0(scenario_name, " obs calculated"), file = "test_statistics_status.txt", append = TRUE)
  
  if (!p_all_only) {
  
  if (parallel) {
    permstats_raw = mclapply(as.list(test_statistics_sampled_ind), function(h) {
      DCARSacrossNetwork(hep_exprs_HVG,
                         edgelist = pairs_hep_all[h, , drop = FALSE],
                         W = W_hep_scenario,
                         weightedConcordanceFunction = scenario$hof,
                         weightedConcordanceFunctionW = scenario$hof_type,
                         verbose = TRUE,
                         extractPermutationTestStatistics = TRUE,
                         niter = 1000)
    }, mc.cores = ncores)
    permstats = lapply(permstats_raw, unlist)
    names(permstats) <- rownames(pairs_hep_all[test_statistics_sampled_ind,])
  } else {
    permstats_raw = DCARSacrossNetwork(hep_exprs_HVG,
                                       edgelist = pairs_hep_all[test_statistics_sampled_ind,],
                                       W = W_hep_scenario,
                                       weightedConcordanceFunction = scenario$hof,
                                       weightedConcordanceFunctionW = scenario$hof_type,
                                       verbose = TRUE,
                                       extractPermutationTestStatistics = TRUE,
                                       niter = 1000)
    permstats = unlist(permstats_raw, recursive = FALSE)
    names(permstats) <- rownames(pairs_hep_all[test_statistics_sampled_ind,])
  }
  test_statistics_perms[[scenario_name]] <- permstats

  write(paste0(scenario_name, " perms calculated"), file = "test_statistics_status.txt", append = TRUE)

} else {
  permstats <- test_statistics_perms[[scenario_name]]
}

pairs_split = split(rownames(pairs_hep_all), rep(1:(2*ncores), length.out = nrow(pairs_hep_all)))

if (parallel) {

p_all_raw = mclapply(pairs_split, function(nms) {
  df = estimatePvaluesSpearman(test_statistics_obs[[scenario_name]][nms], 
                                  test_statistics_global[[scenario_name]][nms], 
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
 p_all = estimatePvaluesSpearman(test_statistics_obs[[scenario_name]], 
                                  test_statistics_global[[scenario_name]], 
                                  permstats,
                                  usenperm = TRUE,
                                  nperm = 10000,
                                  plot = FALSE,
                                  maxDist = 0.1,
                                  verbose = TRUE)
 rownames(p_all) <- names(test_statistics_obs[[scenario_name]])
}
  p_all$fdr <- p.adjust(p_all$pval, method = "BH")
  test_statistics_p_all[[scenario_name]] <- p_all
  
  write(paste0(scenario_name, " estimated pvals calculated"), file = "test_statistics_status.txt", append = TRUE)

  # save and move to next one
  save(test_statistics_perms,
       test_statistics_obs,
       test_statistics_p_all,
       test_statistics_global,
       test_statistics_sampled_ind,
       W_hep_thin,
       file = "output/test_statistics_output.RData")

  write(paste0(scenario_name, " output saved"), file = "test_statistics_status.txt", append = TRUE)

}
