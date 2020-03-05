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
#library(minerva)
#library(energy)

#pwd = getwd()
#write(pwd, file = "wd.txt")

load("output/counts_coords.RData")
source("../functions.R")
#load("output/exprs_HVG.RData")
nonDEgenes = readRDS("output/nonDEgenes.Rds")
load("output/sampled_permstats.Rdata")

pair_genes = intersect(HVG, nonDEgenes)

pairs = t(combn(intersect(HVG, nonDEgenes),2))
rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")

counts_genes = counts[pair_genes,]
#span_sampled_ind = which(rownames(pairs) %in% rownames(pairs_sampled))

#dim(pairs)

span_values = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70)/100
names(span_values) <- paste0("span_", sprintf("%.2f", span_values))

span_scenarios = sapply(span_values, function(spanval) {
  list(weightMatrix = weightMatrix_nD(coords, span = spanval))
  }, simplify = FALSE)

parallel = TRUE
ncores = 10

p_all_only = FALSE

# calculate global values for each and take stratified sample of size 100 from 
# each and take the union of these
if (!file.exists("output/span_global.RData")) {
global_spearman = apply(pairs, 1, function(x) cor(counts[x[1],], counts[x[2],],
  method = "spearman"))
save(global_spearman,
  file = "output/span_global.RData")
} else {
  load("output/span_global.RData")
}


if (!file.exists("output/span_output.RData")) {
span_nsample = 500
span_sampled_ind = unique(c(
  stratifiedSample(global_spearman, span_nsample)
))
length(span_sampled_ind)

span_perms = list() # save the permutations
span_obs = list() # save the obs stats
span_p_all = list() # save the results tables
span_global = list(
  "Spearman" = global_spearman
)
} else {
  load("output/span_output.RData")
}

for (scenario_name in names(span_scenarios)) {
  scenario = span_scenarios[[scenario_name]]

  if (!p_all_only) { # this redoes p_all calculation for all

  if (scenario_name %in% names(span_p_all)) next
  
  }  

  W_scenario = scenario$"weightMatrix"

  if (!p_all_only) {
  span_obs[[scenario_name]] = DCARSacrossNetwork(
    counts_genes,
    edgelist = pairs,
    W = W_scenario,
    extractTestStatisticOnly = TRUE,
    weightedConcordanceFunction = weightedSpearman,
    weightedConcordanceFunctionW = "vector",
    verbose = TRUE
  )
}

  write(paste0(scenario_name, " obs calculated"), file = "span_status.txt", append = TRUE)
  
  if (!p_all_only) {
  
  if (parallel) {
    permstats_raw = mclapply(as.list(span_sampled_ind), function(h) {
      DCARSacrossNetwork(counts_genes,
                         edgelist = pairs[h, , drop = FALSE],
                         W = W_scenario,
                         weightedConcordanceFunction = weightedSpearman,
                         weightedConcordanceFunctionW = "vector",
                         verbose = TRUE,
                         extractPermutationTestStatistics = TRUE,
                         niter = 1000)
    }, mc.cores = ncores)
    permstats = lapply(permstats_raw, unlist)
    names(permstats) <- rownames(pairs[span_sampled_ind,])
  } else {
    permstats_raw = DCARSacrossNetwork(counts_genes,
                                       edgelist = pairs[span_sampled_ind,],
                                       W = W_scenario,
                                       weightedConcordanceFunction = weightedSpearman,
                                       weightedConcordanceFunctionW = "vector",
                                       verbose = TRUE,
                                       extractPermutationTestStatistics = TRUE,
                                       niter = 1000)
    permstats = unlist(permstats_raw, recursive = FALSE)
    names(permstats) <- rownames(pairs[span_sampled_ind,])
  }
  span_perms[[scenario_name]] <- permstats

  write(paste0(scenario_name, " perms calculated"), file = "span_status.txt", append = TRUE)

} else {
  permstats <- span_perms[[scenario_name]]
}

pairs_split = split(rownames(pairs), rep(1:(2*ncores), length.out = nrow(pairs)))

if (parallel) {

p_all_raw = mclapply(pairs_split, function(nms) {
  df = estimatePvaluesSpearman(span_obs[[scenario_name]][nms], 
                                  span_global[["Spearman"]], 
                                  permstats,
                                  usenperm = TRUE,
                                  nperm = 10000,
                                  plot = FALSE,
                                  maxDist = 2,
                                  verbose = TRUE)
  rownames(df) <- nms
  return(df)
  }, mc.cores = 2*ncores)
p_all = do.call(rbind, p_all_raw)
rownames(p_all) <- as.character(p_all$genepair)
p_all <- p_all[rownames(pairs),]
} else {
 p_all = estimatePvaluesSpearman(span_obs[[scenario_name]], 
                                  span_global[["Spearman"]], 
                                  permstats,
                                  usenperm = TRUE,
                                  nperm = 10000,
                                  plot = FALSE,
                                  maxDist = 2,
                                  verbose = TRUE)
 rownames(p_all) <- names(span_obs[[scenario_name]])
}
  p_all$fdr <- p.adjust(p_all$pval, method = "BH")
  span_p_all[[scenario_name]] <- p_all
  
  write(paste0(scenario_name, " estimated pvals calculated"), file = "span_status.txt", append = TRUE)

  # save and move to next one
  save(span_perms,
       span_obs,
       span_p_all,
       span_global,
       span_sampled_ind,
       file = "output/span_output.RData")

  write(paste0(scenario_name, " output saved"), file = "span_status.txt", append = TRUE)

}
