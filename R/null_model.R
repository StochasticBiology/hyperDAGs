#' Simulate from a null model of no character interactions, given a dataset
#'
#' @param m vector of data
#' @param acc.rate numeric (default 1). For phylogenetic data, the null model involves simulation on the given tree, with character acquisition rates proportional to each character's prevalence in the data. This param gives the constant of proportionality. Choose it so that the number of positive observations in the null simulation is comparable to that in the original data.
#' @param tree tree (optional). Tree linking the observations.
#' @return a named list, containing the solution object for this null model instance, and a dataframe with fit_properties of the solution
#' @export
simulate_null_model = function(m, acc.rate = 1,
                               tree = NULL) {
char.probs = colSums(m/nrow(m))
L = ncol(m)

if(is.null(tree)) {
  r.set = list()
  r.set$x = 0
  r.set[["ancnames"]] = c()
  r.set[["descnames"]] = c()
  for(i in 1:nrow(m)) {
    tmp = rbinom(L, 1, prob=char.probs)
    src = tmp*0
    r.set$x = r.set$x + sum(tmp)
    r.set[["ancnames"]] = c(r.set[["ancnames"]], paste0(src, collapse=""))
    r.set[["descnames"]] = c(r.set[["descnames"]], paste0(tmp, collapse=""))
  }
  tree = star_tree(nrow(m))
} else {
    r.set = simulate_accumulation(0, L, use.tree=tree,
                                      accumulation.rate = acc.rate,
                                      dynamics="heterogeneous",
                                      char.probs = char.probs)
}
    n1s = sum(unlist(r.set$x))
    soln = simplest_DAG(r.set[["ancnames"]], r.set[["descnames"]])
    sim.df = fit_properties(soln)
    sim.df$n1 = n1s
    return(list(
                soln=soln,
                sim.df=sim.df))
}

#simulate_null_model(mro.m, 0.08, mro.c$tree)
#simulate_null_model(tb.m, 110, src.data$tree)
