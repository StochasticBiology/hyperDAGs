# Kostas' implementation of the "Gutin algorithm"
# https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 (page 3)
# https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5
# Iain's edits following this, finding a first try for a DAG that spans observation pairs

#' Binary to decimal helper function
#'
#' @param x string binary value to convert
#' @return numeric decimal value
#' @examples
#' BinToDec("1001")
#' @export
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

#' Decimal to binary helper function
#'
#' @param x numeric value to convert
#' @param len numeric length of required binary string
#' @return string binary value
#' @examples
#' DecToBin(9, 4)
#' @export
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

#' Decimal to (numeric) binary helper function
#'
#' @param x numeric value to convert
#' @param len numeric length of required binary string
#' @return numeric binary vector
#' @examples
#' DecToBinV(9, 4)
#' @export
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

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
#' @param verbose Boolean, whether to output summary as a message
#' @return data frame of summary properties
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' fit_properties(sol)
#' @export
fit_properties = function(fit, verbose=TRUE) {
  thisL = stringr::str_length(fit$dataset$ancestors[1])
  df = data.frame(L = stringr::str_length(fit$dataset$ancestors[1]),
                  ntrans = stringr::str_length(fit$dataset$ancestors[1]),
                  ntotal = nrow(fit$dataset),
                  nuniq = length(unique(c(fit$dataset$ancestors, fit$dataset$descendants))),
                  S = round(1-fit$best.bc/nrow(fit$dataset), digits=2),
                  Sprime = round(1-fit$best.bc/nrow(unique(fit$dataset)), digits=2),
                  Sstar = round(1-fit$best.bc/choose(thisL, floor(thisL/2)), digits=2),
                  modE = length(igraph::E(fit$best.graph)),
                  B = branching_count(fit$best.graph),
                  LS = layer_sum(fit$best.graph))
  if(verbose == TRUE) {
  str = paste0("L = ", stringr::str_length(fit$dataset$ancestors[1]),
               "; ntrans = ", nrow(unique(fit$dataset)),
               "; ntotal = ", nrow(fit$dataset),
               "; nuniq = ", length(unique(c(fit$dataset$ancestors, fit$dataset$descendants))),
               "; S = ", round(1-fit$best.bc/nrow(fit$dataset), digits=2),
               "; S' = ", round(1-fit$best.bc/nrow(unique(fit$dataset)), digits=2),
               "; S* = ", round(1-fit$best.bc/choose(thisL, floor(thisL/2)), digits=2),
               "; |E| = ", length(igraph::E(fit$best.graph)),
               "; B = ", branching_count(fit$best.graph),
               "; LS = ", layer_sum(fit$best.graph)
  )
  message(str)
  }
  return(df)
}

#' Write a solution to a file, choosing a particular set of single steps for each multistep change (used to compare HyperHMM output)
#'
#' @param trans XXX
#' @param L XXX
#' @param fname string filename
#' @return XXX
#' @examples
#' #sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' #write_single_steps(sol)
#' print(999)
#' @export
write_single_steps = function(trans, L, fname) {
  trans.set = matrix(ncol=2)
  for(i in 1:nrow(trans)) {
    src = DecToBinV(trans[i,1], L)
    dest = DecToBinV(trans[i,2], L)
    changes = which(dest-src == 1)
    curr = src
    for(j in changes) {
      trans.set = rbind(trans.set, c(BinToDec(curr), BinToDec(curr)+2**(L-j) ))
      curr[j] = 1
    }
  }
  utils::write.table(unique(trans.set[2:nrow(trans.set),]), fname, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#' Difference helper function for vectorised calculations
#'
#' @param x numeric
#' @param y numeric
#' @return numeric difference
#' @examples
#' matdiff(5, 2)
#' @export
matdiff <- function(x, y) {
  return(x-y)
}

#' Check compatibility between two bitstrings, reported as a distance
#'
#' @param s1 character bit string 1
#' @param s2 character bit string 2
#' @return numeric difference
#' @examples
#' is_compat("101", "001")
#' @export
is_compat = function(s1, s2) {
  diffs = as.numeric(strsplit(s2, "")[[1]]) - as.numeric(strsplit(s1, "")[[1]])
  if(min(diffs) < 0) {
    return(Inf)
  } else if(sum(diffs) == 0) {
    return(Inf)
  } else {
    return(sum(diffs))
  }
}

#' Check compatibility between two numeric bitstrings, reported as a distance
#'
#' @param s1 numeric bit vector 1
#' @param s2 numeric bit vector 2
#' @return numeric difference
#' @examples
#' is_compat_num(c(1,0,1), c(0,0,1))
#' @export
is_compat_num = function(s1, s2) {
  diffs = s2-s1
  sumdiffs = sum(diffs)
  if(min(diffs) < 0) {
    return(Inf)
  } else if(sumdiffs == 0) {
    return(Inf)
  } else {
    return(sumdiffs)
  }
}

#' Compute the simplest spanning arborescence for a dataset
#'
#' @param ancnames character vector of states
#' @param descnames character vector of descendant states (optional)
#' @return list: names (state list), raw.graph (first solution), rewired.graph (minimised layer sum), raw.bc (branching count of first solution), rewired.bc (branching count of rewired solution), raw.ls (layer sum of first solution), rewired.ls (layer sum of rewired solution)
#' @examples
#' simplest_arborescence(c("1011", "1001", "0100"))
#' @export
simplest_arborescence = function(ancnames, descnames=NULL) {

  # for the arborescence picture we don't care about ancestor/descendant relationships; just take the union of all data
  names = c(ancnames, descnames)

  message("Starting Algorithm 1")

  ## compute the distance matrix among all entries

  # sort data, this way the always the minimum nearest neighbor will be chosen
  len=stringr::str_length(names[1])
  distMat=stringdist::stringsimmatrix(names,method = "hamming")
  distMat=round(len-len*distMat)
  tdistMat = distMat
  distMat=distMat/distMat
  diag(distMat) <- len+1

  skeleton= data.frame(Anc=character(),
                       Desc=character())

  message(". Computing difference statistics")
  onecounts = sapply(names, stringr::str_count, "1")

  onecountsdiff = outer(onecounts, onecounts, Vectorize(matdiff))
  absonecountsdiff = abs(onecountsdiff)

  # Kostas first condition:
  set.1 = which(tdistMat > absonecountsdiff, arr.ind = TRUE)

  # Kostas third [our second] condition:
  set.2 = which(onecountsdiff > 0, arr.ind = TRUE)
  set.2 = set.2[set.2[,1] < set.2[,2],]

  # Kostas second [our third] condition:
  set.3 = which(onecountsdiff < 0, arr.ind = TRUE)
  set.3 = set.3[set.3[,1] < set.3[,2],]

  # condition 1 trumps conditions 2 and 3, so find indices unique to 2 and 3
  unique.2 = suppressMessages( as.matrix(dplyr::anti_join(as.data.frame(set.2), as.data.frame(set.1))) )
  unique.3 = suppressMessages( as.matrix(dplyr::anti_join(as.data.frame(set.3), as.data.frame(set.1))) )

  # apply distMat effects of these conditions
  distMat[set.1] = len+1
  distMat[unique.2] = 2
  distMat[unique.3] = 1

  # add the corresponding entries to the skeleton (ji for condition 2, ij for condition 3)
  skeleton = rbind(skeleton, data.frame(Anc=names[unique.2[,2]],
                                        Desc=names[unique.2[,1]]))
  skeleton = rbind(skeleton, data.frame(Anc=names[unique.3[,1]],
                                        Desc=names[unique.3[,2]]))

  # just to make sure we don't have self-loops
  skeleton=skeleton[which(skeleton$Anc != skeleton$Desc),]
  skeleton=unique(skeleton)

  message(". Creating working graphs")
  # create the two sets for the bipartite graph B
  V=unique(union(skeleton[,1],skeleton[,2]))
  Vprim=unique(skeleton[,2])

  # IGJ: paste is already vectorised
  Vprim = paste(Vprim, "P", sep="")

  # IGJ: vectorising this too
  B = data.frame(Anc=skeleton[,1], Desc=paste(skeleton[,2], "P", sep=""))

  graphB <- igraph::graph_from_edgelist(as.matrix(B),directed = F)

  types= c(rep(0,length(V)),rep(1,length(Vprim)))
  names(types) <- c(V,Vprim)

  # IGJ vectorising this
  edgesB = c(t(as.matrix(B)))

  message(". Working with graphs")
  # create the bipartite graph B
  bp=igraph::make_bipartite_graph(types,edgesB)
  igraph::is_bipartite(graphB)

  # find a maximum matching M in B.
  maxbp= igraph::max_bipartite_match(bp)

  # create M*
  Mstar=maxbp$matching
  Mstar=Mstar[ceiling(length(Mstar)/2+1):length(Mstar)]
  from=as.character(Mstar)
  to=names(Mstar)
  to=sub("P","",to)
  from[which(is.na(from))]= strrep("0",nchar(skeleton[1,1]))

  # final is the final MINLEAF arborescence
  final=cbind(from,to)
  graphB <- igraph::graph_from_edgelist(as.matrix(final),directed = T)
  tree.leaves=setdiff(final[,2],final[,1])

  ### rewire
  message(". Computing compatibilities for rewiring")

  num.names = lapply(strsplit(names, ""), as.numeric)

  compats = outer(num.names, num.names, Vectorize(is_compat_num))

  message(". Rewiring")

  new.final = final
  for(i in 1:nrow(new.final)) {
    anc = new.final[i,1]
    desc.refs = which(new.final[,1]==anc)
    outdeg = length(desc.refs)
    if(outdeg > 1) {
      for(j in 1:outdeg) {
        desc = new.final[desc.refs[j],2]
        desc.ref = which(names == desc)[1]
        desc.compats = compats[,desc.ref]
        best.new.ref = which(desc.compats == min(desc.compats) )[1]
        new.final[desc.refs[j],1] = names[best.new.ref]
      }
    }
  }

  message(". Wrapping up")

  new.graphB <- igraph::graph_from_edgelist(as.matrix(new.final),directed = T)

  new.final[new.final[,1]==paste0(rep("0", len), collapse=""),]

  rlist = list(len = len,
               names = names,
               raw.graph = graphB,
               rewired.graph = new.graphB,
               raw.bc = branching_count(graphB),
               rewired.bc = branching_count(new.graphB),
               raw.ls = layer_sum(graphB),
               rewired.ls = layer_sum(new.graphB))

  return(rlist)
}

#' Compute the simplest transition-spanning DAG for ancestor-descendant data
#'
#' @param ancnames character vector of states
#' @param descnames character vector of descendant states
#' @return list: names (state list), raw.graph (first solution), rewired.graph (minimised layer sum), raw.bc (branching count of first solution), rewired.bc (branching count of rewired solution), raw.ls (layer sum of first solution), rewired.ls (layer sum of rewired solution), dataset (data), best.graph (best solution), best.bc (branching count of best solution)
#' @examples
#' simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' @export
simplest_DAG = function(ancnames, descnames) {

  # first get the simplest spanning arborescence for the union of the data
  s.a = simplest_arborescence(ancnames, descnames)

  if(length(descnames) == 0) {
    rlist = s.a
    rlist$dataset = ancnames
    rlist$best.graph = rlist$rewired.graph
    rlist$best.bc = rlist$rewired.bc
    return(rlist)
  }

  message("Starting Algorithm 2")

  graphC = s.a$rewired.graph
  L = s.a$len

  message(". Adding edges")
  # loop through descendant nodes
  for(this.desc in unique(descnames)) {
    # if this descendant node doesn't exist in the graph, add it
    if(!(this.desc %in% igraph::V(graphC)$name)) {
      #      print(paste0("Adding ", this.desc))
      graphC = igraph::add_vertices(graphC, n = 1)
      igraph::V(graphC)$name[length(igraph::V(graphC)$name)] = this.desc
    }
    # get list of ancestral states for this descendant
    this.ancs = c()
    anc.refs = which(descnames == this.desc)
    for(i in anc.refs) {
      this.ancs = c(this.ancs, ancnames[i])
    }
    # count the 1s in each; we want to start connecting from the one with most 1s
    count1s = rep(0, length(this.ancs))
    for(i in 1:length(this.ancs)) {
      count1s[i] = stringr::str_count(this.ancs[i], "1")
    }
    count1order = order(count1s, decreasing=TRUE)
    # loop through ancestors for this descendant
    for(i in count1order) {
      this.this.ancs = this.ancs[i]
      # if we don't have a path from this ancestor, add an edge
      path = suppressWarnings( igraph::get.shortest.paths(graphC, from=this.this.ancs, to=this.desc) )
      if(length(path$vpath[[1]]) == 0) {
        from_id = which(igraph::V(graphC)$name == this.this.ancs)
        to_id = which(igraph::V(graphC)$name == this.desc)
        graphC = igraph::add_edges(graphC, c(from_id, to_id))
      } else {
        #        print(paste0(" Already have path"))
      }
    }
  }

  #######
  # now prune edges we don't need

  message(". Pruning edges")

  graphD = graphC
  paths = list()
  plengths = rep(0, length(ancnames))
  # for debugging, check that we have all the observations we need
  for(i in 1:length(ancnames)) {
    paths[[i]] = suppressWarnings( igraph::get.shortest.paths(graphD, from=ancnames[i], to=descnames[i]) )
    plengths[i] = length(paths[[i]]$vpath[[1]])
  }

  # loop through (dynamic) edge set. *not* a very efficient approach
  eref = 1
  while(eref < length(igraph::E(graphD))) {
    e = igraph::E(graphD)[eref]
    # prune this edge
    tmp = igraph::delete_edges(graphD, e)
    plengths = rep(0, length(ancnames))
    paths = list()
    # check if all our paths are still intact
    for(i in 1:length(ancnames)) {
      paths[[i]] = suppressWarnings( igraph::get.shortest.paths(tmp, from=ancnames[i], to=descnames[i]) )
      plengths[i] = length(paths[[i]]$vpath[[1]])
    }
    # if so, remove this edge
    if(min(plengths) > 0) {
      graphD = tmp
    } else {
      eref = eref+1
    }
  }

  paths = list()
  plengths = rep(0, length(ancnames))
  for(i in 1:length(ancnames)) {
    paths[[i]] = suppressWarnings( igraph::get.shortest.paths(graphD, from=ancnames[i], to=descnames[i]) )
    plengths[i] = length(paths[[i]]$vpath[[1]])
  }
  if(min(plengths) == 0) {
    message("Pruning's gone wrong!")
  }

  message(". Connecting to root")
  in_degrees <- igraph::degree(graphD, mode = "in")
  zeroes = which(in_degrees == 0)
  root.name = paste0(rep("0", stringr::str_length(ancnames[i])), collapse="")
  root.label = which(igraph::V(graphD)$name == root.name)
  if(length(root.label) == 0) {
    graphD = igraph::add_vertices(graphD, n = 1)
    igraph::V(graphD)$name[length(igraph::V(graphD)$name)] = root.name
    root.label = length(igraph::V(graphD)$name)
  }
  for(this.zero in zeroes) {
    if(this.zero != root.label) {
      graphD = igraph::add_edges(graphD, c(root.label, this.zero))
    }
  }

  message(". Wrapping up")
  graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")

  dataset = data.frame(ancestors = ancnames,
                       descendants = descnames )

  rlist = s.a
  rlist$dataset = dataset
  rlist$best.graph = graphD
  rlist$best.bc = branching_count(graphD)

  return(rlist)
}

#' Plot the unrewired arboresence (output from Gutin algorithm) and the rewired version increasing layer sum
#'
#' @param graphs list output of model fit
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_stage_1(sol)
#' @export
plot_stage_1 = function(graphs) {
  # extract graphs from algorithm output and assign nodes to layers
  graphB = graphs$raw.graph
  new.graphB = graphs$rewired.graph
  graphB.layers = sapply(igraph::V(graphB)$name, stringr::str_count, "1")
  new.graphB.layers = sapply(igraph::V(new.graphB)$name, stringr::str_count, "1")
  L = graphs$len

  # decide on label size (heuristic, for clarity)
  label.size = 3
  if(L > 5) {
    label.size = 2
  }
  if(L > 20) {
    label.size = 0
  }

  return(
    ggpubr::ggarrange(
      ggraph::ggraph(graphB, layout="sugiyama", layers=graphB.layers) +
        ggraph::geom_edge_link(color="#CCCCCC") +
        ggraph::geom_node_text(ggplot2::aes(label=name), angle=45, hjust=0, size=label.size) +
        ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
        ggplot2::scale_y_continuous(expand = c(0.1, 0.1)) +
        ggraph::theme_graph(),
      ggraph::ggraph(new.graphB, layout="sugiyama", layers=new.graphB.layers) +
        ggraph::geom_edge_link(color="#CCCCCC") +
        ggraph::geom_node_text(ggplot2::aes(label=name), angle=45, hjust=0, size=label.size) +
        ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
        ggplot2::scale_y_continuous(expand = c(0.1, 0.1)) +
        ggraph::theme_graph(),
      labels = c("A", "B")
    )
  )
}

#' Plot the simplest spanning arborescence and simplest transition-spanning DAG
#'
#' @param graphs list output of model fit
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_stage_2(sol)
#' @export
plot_stage_2 = function(graphs) {
  # extract graphs from algorithm output and assign layers
  graphB = graphs$rewired.graph
  graphD = graphs$best.graph
  graphB.layers = sapply(igraph::V(graphB)$name, stringr::str_count, "1")
  graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
  L = graphs$len

  # decide on label size (heuristic, for clarity)
  label.size = 3
  if(L > 5) {
    label.size = 2
  }
  if(L > 20) {
    label.size = 1
    igraph::V(graphB)$name = 1:length(igraph::V(graphB))
    igraph::V(graphD)$name = 1:length(igraph::V(graphD))
  }

  return(
    ggpubr::ggarrange( ggraph::ggraph(graphB, layout="sugiyama", layers=graphB.layers) +
                         ggraph::geom_edge_link(color="#CCCCCC") +
                         ggraph::geom_node_text(ggplot2::aes(label=name), size=label.size, angle=45, hjust=0) +
                         ggplot2::ggtitle(branching_count(graphB)) +
                         ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
                         ggplot2::theme_void(), #ggraph::theme_graph(),
                       ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
                         ggraph::geom_edge_link(color="#CCCCCC") +
                         ggraph::geom_node_text(ggplot2::aes(label=name), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) +
                         ggplot2::ggtitle(branching_count(graphD)) +
                         ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
                         ggplot2::theme_void(), #ggraph::theme_graph(),

                       nrow = 1, labels=c("A", "B")
    )
  )
}

#' Plot a solution graph
#'
#' @param graphD graph to plot
#' @param v.labels data frame of vertex labels (optional)
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_stage_p(sol$best.graph)
#' @export
plot_stage_p = function(graphD, v.labels = data.frame(Species=NA)) {
  # extract graphs from algorithm output and assign layers
  graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
  L = stringr::str_length(igraph::V(graphD)$name[1])

  # decide on label size (heuristic, for clarity)
  label.size = 3
  igraph::V(graphD)$plotname = igraph::V(graphD)$name
  if(L > 5) {
    label.size = 2
    igraph::V(graphD)$plotname = igraph::V(graphD)$name
  }
  if(L > 20) {
    label.size = 1
    igraph::V(graphD)$plotname = 1:length(igraph::V(graphD))
  }

  if(!is.na(v.labels$Species[1])) {
    igraph::V(graphD)$v.label = ""
    for(i in 1:nrow(v.labels)) {
      ref = which(igraph::V(graphD)$name == v.labels$label[i])
      igraph::V(graphD)$v.label[ref] = v.labels$Species[i]
    }
  }
  g = ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
    ggraph::geom_edge_link(color="#CCCCCC") +
    ggraph::geom_node_text(ggplot2::aes(label=plotname), size=label.size, angle=45, hjust=0) +
    ggplot2::ggtitle(paste0("B = ", branching_count(graphD), collapse="")) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_void() #ggraph::theme_graph()
  if(!is.na(v.labels$Species[1])) {
    g = g +
      ggraph::geom_node_text(ggplot2::aes(label=v.label), size=2, angle=45, hjust=0)  #, check_overlap = TRUE) +
  }
  return(g)
}

#' Plot a solution graph with edges weighted by downstream datapoints
#'
#' @param graphD graph to plot
#' @param labels vector of vertex labels (optional)
#' @param thresh numeric threshold below which edges are not plotted (default 5)
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_weights(sol$best.graph, thresh = 0)
#' @export
plot_weights = function(graphD, labels=c(""), thresh=5) {
  # extract graphs from algorithm output and assign layers
  graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
  L = stringr::str_length(igraph::V(graphD)$name[1])
  if(labels[1] == "") {
    labels = as.character(1:L)
  }
  graphD.size = igraph::neighborhood.size(graphD, L+1, mode="out")
  igraph::E(graphD)$thickness = as.numeric(graphD.size[igraph::ends(graphD, es = igraph::E(graphD), names = FALSE)[, 2]])
  this.ends = igraph::ends(graphD, es=igraph::E(graphD))
  srcs = strsplit(this.ends[,1], split="")
  dests = strsplit(this.ends[,2], split="")
  for(i in 1:nrow(this.ends)) {
    igraph::E(graphD)$label[i] = paste0("+", paste0(labels[which(srcs[[i]]!=dests[[i]])], collapse=","), collapse="")
  }
  # decide on label size (heuristic, for clarity)
  edge.label.size = 3
  label.size = 3
  igraph::V(graphD)$plotname = igraph::V(graphD)$name
  if(L > 5) {
    label.size = 2
    igraph::V(graphD)$plotname = igraph::V(graphD)$name
    edge.label.size = 3
  }
  if(L > 20) {
    label.size = 1
    igraph::V(graphD)$plotname = 1:length(igraph::V(graphD))
    edge.label.size = 2
  }

  full.plot = ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
    ggraph::geom_edge_link(color="#CCCCCC", label_size= edge.label.size, angle_calc = "across", label_colour = "grey", ggplot2::aes(label=label, alpha=log(as.numeric(thickness)+1))) +
    ggraph::geom_node_text(ggplot2::aes(label=plotname), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) +
    ggplot2::ggtitle(paste0("B = ", branching_count(graphD), collapse="")) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_void() #ggraph::theme_graph()

  graphE = igraph::delete_edges(graphD, which(as.numeric(igraph::E(graphD)$thickness) < thresh))
  degrees = igraph::degree(graphE)
  graphE = igraph::delete_vertices(graphE, which(degrees == 0))
  graphE.layers = sapply(igraph::V(graphE)$name, stringr::str_count, "1")

  label.size = 0
  thresh.plot = ggraph::ggraph(graphE, layout="sugiyama", layers=graphE.layers) +
    ggraph::geom_edge_link(color="#AAAAFF", label_size= edge.label.size, angle_calc = "across", label_colour = "black",
                           ggplot2::aes(label=label, edge_width = as.numeric(thickness), alpha=as.numeric(thickness))) +
    ggraph::geom_node_text(ggplot2::aes(label=name), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) +
    ggraph::scale_edge_width(limits=c(0,NA)) +
    ggplot2::theme_void() + #ggraph::theme_graph() +
    ggplot2::theme(legend.position="none")

  return(list(full.plot=full.plot, thresh.plot=thresh.plot,
              full.graph=graphD, thresh.graph=graphE))
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
