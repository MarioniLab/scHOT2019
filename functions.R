weightMatrix_nD = function(x, span = 0.5) {
  
  # x is a matrix with rows corresponding to cells and columns corresponding to dimensions to calculate distance
  
  ncells = nrow(x)
 
  # calculate euclidean distance of points on 2D
  # coords = cbind(x,y)
  coords = as.matrix(x)
  d = as.matrix(dist(coords))
  
  # extract a weights vector per cell
  
  W_raw = sapply(seq_len(ncells), function(cell) {
    dvec = d[cell,]
    vals = rep(0, ncells)
    vals[order(dvec)[1:ceiling(span*ncells)]] = seq(1, 0, length.out = ceiling(span*ncells))
    return(vals)
  }, simplify = FALSE)
  
  W = do.call(rbind, W_raw)
  
  return(W)
  
}

cordist = function(x) {
  
  cor = cor(t(x), method = "spearman")
  
  return(as.dist((1-cor)/2))
  
}

weightedMeanMatrixStats = function(x,y,w) {
  require(matrixStats)
  weightedMean(x, w)
}

genesetGOtest = function(set, universe, termsList) {
  # set is the character vector of genes in geneset
  # universe is character vector of genes to include in universe
  # termsList is a list of the pathways
  
  termsListfiltered = lapply(termsList, function(x) return(intersect(universe, x)))
  keepForTesting = unlist(lapply(termsListfiltered, length)) >= 8 & unlist(lapply(termsListfiltered, length)) <= 500
  
  pval = sapply(1:length(termsList), function(i){
    
    if (!keepForTesting[i]) return(1)
    
    termGenes = termsListfiltered[[i]]
    return(fisher.test(table(universe %in% set, universe %in% termGenes), alt = "g")$p.value)
  })
  names(pval) <- names(termsListfiltered)
  return(pval)
}
