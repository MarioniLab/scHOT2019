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

#pwd = getwd()
#write(pwd, file = "wd.txt")

load("output/all_wcors.RData")
hep_p_all <- readRDS("output/hep_p_all.Rds")
load("output/exprs_HVG.RData")


if (!file.exists("output/hep_full_pvals.Rds")) {
  hep_full_pvals = rep(NA, nrow(pairs_hep_all))
  names(hep_full_pvals) <- rownames(pairs_hep_all)
} else {
  hep_full_pvals = readRDS("output/hep_full_pvals.Rds")
}

for (i in 1:nrow(pairs_hep_all)) {
  print(i)
  
  write(i, file = "status.txt", append = TRUE)
  
  if (!is.na(hep_full_pvals[i])) next
  
  if (i %% 20 == 0) {
    print("saving!")
    saveRDS(hep_full_pvals, "output/hep_full_pvals.Rds")
   # plot(-log10(hep_full_pvals), -log10(hep_p_all$pval)); abline(c(0,1))
  }
  
  hepgenes = as.character(pairs_hep_all[i,])
  
  obs = DCARS(hep_exprs_HVG,
                            xname = pairs_hep_all[i,1],
                            yname = pairs_hep_all[i,2],
                            # edgelist = pairs_hep_all[1:2, , drop = FALSE],
                            W = W_hep,
                            weightedConcordanceFunction = weightedSpearman,
                            weightedConcordanceFunctionW = "vector",
                            verbose = TRUE,
                            extractTestStatisticOnly = TRUE,
                            niter = 10000)
  
  hep_dat = hep_exprs_HVG[hepgenes,]
  
  mcPerms = unlist(mclapply(1:40, function(j) {
    DCARS(hep_dat,
                            xname = pairs_hep_all[i,1],
                            yname = pairs_hep_all[i,2],
                            # edgelist = pairs_hep_all[1:2, , drop = FALSE],
                            W = W_hep,
                            weightedConcordanceFunction = weightedSpearman,
                            weightedConcordanceFunctionW = "vector",
                            verbose = TRUE,
                            extractPermutationTestStatistics = TRUE,
                            niter = 250)
  }, mc.cores = 40))
  
  hep_full_pvals[i] = mean(mcPerms >= obs)
   
}

 print("saving!")
 saveRDS(hep_full_pvals, "output/hep_full_pvals.Rds")