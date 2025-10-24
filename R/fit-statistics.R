

#' Compute excess branching score of a graph
#'
#' @param g graph to analyse
#' @return numeric excess branching score
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' branching_count(sol$best.graph)
#' @export
branching_count = function(g) {
  b = igraph::degree(g, mode="out")-1
  b[b<0] = 0
  return(sum(b))
}

#' Compute layer sum of a graph
#'
#' @param g graph to analyse
#' @return numeric layer sum
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' layer_sum(sol$best.graph)
#' @export
layer_sum = function(g) {
  b = igraph::degree(g, mode="out")-1
  # get the name, and excess branching, of all nodes with excess branching > 0
  b.names = igraph::V(g)$name[b > 0]
  b.outs = b[b>0]
  # count the 1s in each name to get the layer of each node
  count1s = sapply(b.names, stringr::str_count, "1")
  # return summed product of excess branching and layer
  # i.e. the sum of source layer over all edges that contribute to excess branching
  if(length(b.outs) == 0) { return(0) }
  return(sum(b.outs*count1s))
}

#' Summarise a set of properties of a fitted HyperDAGs model
#'
#' @param fit model structure to analyse
#' @param verbose Boolean, whether to output summary as a message (default false)
#' @return data frame of summary properties
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' fit_properties(sol)
#' @export
fit_properties = function(fit, verbose=FALSE) {
  thisL = stringr::str_length(fit$dataset$ancestors[1])
  df = data.frame(L = stringr::str_length(fit$dataset$ancestors[1]),
                  ntrans = nrow(fit$dataset),
                  nuniqtrans = nrow(unique(fit$dataset)),
                  nuniqstates = length(unique(c(fit$dataset$ancestors, fit$dataset$descendants))),
                  S = round(1-fit$best.bc/(nrow(fit$dataset)-1), digits=2),
                  Sprime = round(1-fit$best.bc/(nrow(unique(fit$dataset))-1), digits=2),
                  Sstar = round(1-fit$best.bc/(choose(thisL, floor(thisL/2))-1), digits=2),
                  modE = length(igraph::E(fit$best.graph)),
                  B = branching_count(fit$best.graph),
                  LS = layer_sum(fit$best.graph))
  if(verbose == TRUE) {
    str = paste0("L = ", df$L,
                 "; ntrans = ", df$ntrans,
                 "; ntotal = ", df$ntotal,
                 "; nuniq = ", df$nuniq,
                 "; S = ", df$S,
                 "; S' = ", df$Sprime,
                 "; S* = ", df$Sstar,
                 "; |E| = ", df$modE,
                 "; B = ", df$B,
                 "; LS = ", df$LS
    )
    message(str)
  }
  return(df)
}

#' Verify that a solution graph spans the required transitions
#'
#' @param g graph of fitted model
#' @param ancnames vector of ancestor states
#' @param descnames vector of descendant states
#' @return Boolean: are edges spanned?
#' @examples
#' ancnames = c("000", "010", "100")
#' descnames = c("110", "110", "101")
#' sol = simplest_DAG(ancnames, descnames)
#' transitions_spanned(sol$best.graph, ancnames, descnames)
#' @export
transitions_spanned = function(g, ancnames, descnames) {
  edges.ok = TRUE
  for(i in 1:length(ancnames)) {
    anc.name = ancnames[i]
    desc.name = descnames[i]
    anc.ref = which(igraph::V(g)$name==anc.name)
    desc.ref = which(igraph::V(g)$name==desc.name)
    if(!is.finite(igraph::distances(g, anc.ref, desc.ref, mode="out"))) {
      message(paste0("- Lost connection ", anc.name, " -> ", desc.name))
      edges.ok = FALSE
    }
  }
  return(edges.ok)
}

