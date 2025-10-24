# Kostas' implementation of the "Gutin algorithm"
# https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 (page 3)
# https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5
# Iain's edits following this, finding a first try for a DAG that spans observation pairs

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
                                        Desc=names[unique.2[,1]],
                                        stringsAsFactors=FALSE))
  skeleton = rbind(skeleton, data.frame(Anc=names[unique.3[,1]],
                                        Desc=names[unique.3[,2]],
                                        stringsAsFactors=FALSE))

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

