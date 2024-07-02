# Kostas' implementation of the "Gutin algorithm"
# https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 (page 3)
# https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5
# Iain's edits following this, finding a first try for a DAG that spans observation pairs

library(igraph)
library(stringr)
library(DescTools)
library(dplyr)
library(stringdist)
library(ggraph)
library(ggplot2)
library(ggpubr)

# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

# decimal to binary function
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

# decimal to binary function, returning a numerical vector
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

# excess branching score for minimisation
branching.count = function(g) {
  b = degree(g, mode="out")-1
  b[b<0] = 0
  return(sum(b))
}

# layer sum
layer.sum = function(g) {
  b = degree(g, mode="out")-1
  b.names = V(g)$name[b > 0]
  b.outs = b[b>0]
  count1s = sapply(b.names, str_count, "1")
  return(sum(b.outs*count1s))
}

# Define the function f
matdiff <- function(x, y) {
  return(x-y)  
}

# check compatibility between two bitstrings, reported as a distance
is.compat = function(s1, s2) {
  diffs = as.numeric(strsplit(s2, "")[[1]]) - as.numeric(strsplit(s1, "")[[1]])
  if(min(diffs) < 0) {
    return(Inf) 
  } else if(sum(diffs) == 0) {
    return(Inf)
  } else {
    return(sum(diffs))
  }
}

# check compatibility between two bitstrings, reported as a distance
is.compat.num = function(s1, s2) {
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

# compute the simplest spanning arborescence for a dataset
simplest.arborescence = function(ancnames, descnames=NULL) {
  
  # for the arborescence picture we don't care about ancestor/descendant relationships; just take the union of all data
  names = c(ancnames, descnames)
  
  message("Starting Algorithm 1")
  
  ## compute the distance matrix among all entries
  
  # sort data, this way the always the minimum nearest neighbor will be chosen
  len=str_length(names[1])
  distMat=stringsimmatrix(names,method = "hamming")
  distMat=round(len-len*distMat)
  tdistMat = distMat
  distMat=distMat/distMat
  diag(distMat) <- len+1
  
  skeleton= data.frame(Anc=character(),
                       Desc=character())
  
  message(". Computing difference statistics")
  onecounts = sapply(names, str_count, "1")
  
  onecountsdiff = outer(onecounts, onecounts, Vectorize(matdiff)) 
  absonecountsdiff = abs(onecountsdiff)
  
  # Kostas first condition:
  # if ( as.numeric(StrDist(names[i],names[j],method = "hamming"))> abs((str_count(names[i], "1")-str_count(names[j], "1"))) ) {
  set.1 = which(tdistMat > absonecountsdiff, arr.ind = TRUE)
  
  # Kostas third [our second] condition:
  # (str_count(names[i], "1")-str_count(names[j], "1"))>0 [where i<j]
  set.2 = which(onecountsdiff > 0, arr.ind = TRUE)
  set.2 = set.2[set.2[,1] < set.2[,2],]
  
  # Kostas second [our third] condition:
  # (str_count(names[i], "1")-str_count(names[j], "1"))<0 ) [where i<j]
  set.3 = which(onecountsdiff < 0, arr.ind = TRUE)
  set.3 = set.3[set.3[,1] < set.3[,2],]
  
  # condition 1 trumps conditions 2 and 3, so find indices unique to 2 and 3
  unique.2 = as.matrix(anti_join(as.data.frame(set.2), as.data.frame(set.1)))
  unique.3 = as.matrix(anti_join(as.data.frame(set.3), as.data.frame(set.1)))
  
  # apply distMat effects of these conditions
  distMat[set.1] = len+1
  distMat[unique.2] = 2
  distMat[unique.3] = 1
  
  # add the corresponding entries to the skeleton (ji for condition 2, ij for condition 3)
  skeleton = rbind(skeleton, data.frame(Anc=names[unique.2[,2]],
                                        Desc=names[unique.2[,1]]))
  skeleton = rbind(skeleton, data.frame(Anc=names[unique.3[,1]],
                                        Desc=names[unique.3[,2]]))
  
  # reduced version of Kostas' conditionals
  #for (i in 1:(nrow(distMat)-1)) {
  #  for (j in (i+1):nrow(distMat)) {
  #    if ( tdistMat[i,j] > absonecountsdiff[i,j]) { #abs(onecounts[i]-onecounts[j]) ) {
  #      distMat[i,j]=len+1
  #      distMat[j,i]=len+1
  #    } else if ( onecountsdiff[i,j]>0 ) {
  #      distMat[i,j]=2
  #      skeleton[nrow(skeleton)+1,]=c(names[j],names[i])
  #    } else if ( onecountsdiff[i,j]<0 ) {
  #      distMat[i,j]=1
  #      skeleton[nrow(skeleton)+1,]=c(names[i],names[j])
  #    }
  #  }
  #}
  
  # just to make sure we don't have self-loops
  skeleton=skeleton[which(skeleton$Anc != skeleton$Desc),]
  skeleton=unique(skeleton)
  
  message(". Creating working graphs")
  # create the two sets for the bipartite graph B
  V=unique(union(skeleton[,1],skeleton[,2]))
  Vprim=unique(skeleton[,2])
  
  # IGJ: paste is already vectorised
  Vprim = paste(Vprim, "P", sep="")
  # replacing
  #for (i in 1:length(Vprim)) {
  #  Vprim[i]=paste(Vprim[i],"P",sep = "")
  #}
  
  # IGJ: vectorising this too
  B = data.frame(Anc=skeleton[,1], Desc=paste(skeleton[,2], "P", sep=""))
  # replacing
  #B= data.frame(Anc=character(),
  #              Desc=character())
  #for (i in 1:nrow(skeleton)) {
  #  B[nrow(B)+1,]=c(skeleton[i,1],paste(skeleton[i,2],"P",sep = ""))
  #}
  
  graphB <- graph_from_edgelist(as.matrix(B),directed = F)
  
  types= c(rep(0,length(V)),rep(1,length(Vprim)))
  names(types) <- c(V,Vprim)
  
  # IGJ vectorising this
  edgesB = c(t(as.matrix(B)))
  # replacing
  #edgesB=vector()
  #for (i in 1:nrow(B)) {
  #  edgesB=c(edgesB,B[i,1],B[i,2])
  #}
  
  message(". Working with graphs")
  # create the bipartite graph B
  bp=make_bipartite_graph(types,edgesB)
  is_bipartite(graphB)
  
  # find a maximum matching M in B.
  maxbp= max_bipartite_match(bp)
  
  # create M* 
  Mstar=maxbp$matching
  Mstar=Mstar[ceiling(length(Mstar)/2+1):length(Mstar)]
  from=as.character(Mstar)
  to=names(Mstar)
  to=sub("P","",to)
  from[which(is.na(from))]= strrep("0",nchar(skeleton[1,1]))
  
  # final is the final MINLEAF arborescence
  final=cbind(from,to)
  graphB <- graph_from_edgelist(as.matrix(final),directed = T)
  #write.table(final, file="gutinsAlgorithm-out.csv", row.names=FALSE, col.names=T,sep = ",",quote = F)
  tree.leaves=setdiff(final[,2],final[,1])
  #length(tree.leaves)
  
  ### IGJ code from here
  #ggraph(graphB) + geom_edge_link() + geom_node_text(aes(label=name))
  
  ### rewire
  message(". Computing compatibilities for rewiring")
  
 num.names = lapply(strsplit(names, ""), as.numeric)
  
  #Sys.time()
  #compats = outer(names, names, Vectorize(is.compat)) 
  #Sys.time()
  compats = outer(num.names, num.names, Vectorize(is.compat.num))
  #Sys.time()
    
  message(". Rewiring")
  
  new.final = final
  for(i in 1:nrow(new.final)) {
    anc = new.final[i,1]
    desc.refs = which(new.final[,1]==anc) 
    outdeg = length(desc.refs)
    if(outdeg > 1) {
      #      print(i)
      #      print(paste("Thinking about ", anc, " descendants ", new.final[desc.refs,2]))
      for(j in 1:outdeg) {
        desc = new.final[desc.refs[j],2]
        #        print(paste("  Thinking about ", desc))
        desc.ref = which(names == desc)[1]
        desc.compats = compats[,desc.ref]
        best.new.ref = which(desc.compats == min(desc.compats) )[1]
        #        print(paste("  best compats ", names[best.new.ref]))
        
        new.final[desc.refs[j],1] = names[best.new.ref]
      }
    }
  }
  
  message(". Wrapping up")
  
  new.graphB <- graph_from_edgelist(as.matrix(new.final),directed = T)
  
  #ggraph(new.graphB, layout="sugiyama", layers = new.graphB.layers) + geom_edge_link() + geom_node_text(aes(label=name))
  
  new.final[new.final[,1]==paste0(rep("0", len), collapse=""),]
  
  rlist = list(len = len,
               names = names,
               raw.graph = graphB,
               rewired.graph = new.graphB,
               raw.bc = branching.count(graphB),
               rewired.bc = branching.count(new.graphB),
               raw.ls = layer.sum(graphB),
               rewired.ls = layer.sum(new.graphB))
  
  return(rlist)
}

# compute the simplest transition-spanning DAG for ancestor-descendant data
simplest.DAG = function(ancnames, descnames) {
  
  # first get the simplest spanning arborescence for the union of the data
  s.a = simplest.arborescence(ancnames, descnames)
  
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
    #    print(paste0("Thinking about ", this.desc))
    
    # if this descendant node doesn't exist in the graph, add it
    if(!(this.desc %in% V(graphC)$name)) {
      #      print(paste0("Adding ", this.desc))
      graphC = add_vertices(graphC, n = 1)
      V(graphC)$name[length(V(graphC)$name)] = this.desc
    }
    # get list of ancestral states for this descendant
    this.ancs = c()
    anc.refs = which(descnames == this.desc)
    for(i in anc.refs) {
      this.ancs = c(this.ancs, ancnames[i])
    }
    #    print(paste0("Ancestors are ", this.ancs))
    # count the 1s in each; we want to start connecting from the one with most 1s
    count1s = rep(0, length(this.ancs))
    for(i in 1:length(this.ancs)) {
      count1s[i] = str_count(this.ancs[i], "1")
    }
    count1order = order(count1s, decreasing=TRUE)
    #    print(count1order)
    # loop through ancestors for this descendant
    for(i in count1order) {
      this.this.ancs = this.ancs[i]
      #      print(paste0(" Thinking about ", this.this.ancs))
      
      # if we don't have a path from this ancestor, add an edge
      path = get.shortest.paths(graphC, from=this.this.ancs, to=this.desc)
      if(length(path$vpath[[1]]) == 0) {
        from_id = which(V(graphC)$name == this.this.ancs)
        to_id = which(V(graphC)$name == this.desc)
        #        print(paste0(" Adding ", this.this.ancs, "-", this.desc, " = ", from_id, " ", to_id))
        graphC = add_edges(graphC, c(from_id, to_id))
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
    paths[[i]] = get.shortest.paths(graphD, from=ancnames[i], to=descnames[i])
    plengths[i] = length(paths[[i]]$vpath[[1]])
  }
  
  # loop through (dynamic) edge set. *not* a very efficient approach
  eref = 1
  while(eref < length(E(graphD))) {
    e = E(graphD)[eref]
    # prune this edge
    tmp = delete_edges(graphD, e)
    plengths = rep(0, length(ancnames))
    paths = list()
    # check if all our paths are still intact
    for(i in 1:length(ancnames)) {
      paths[[i]] = get.shortest.paths(tmp, from=ancnames[i], to=descnames[i])
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
    paths[[i]] = get.shortest.paths(graphD, from=ancnames[i], to=descnames[i])
    plengths[i] = length(paths[[i]]$vpath[[1]])
  }
  if(min(plengths) == 0) {
    message("Pruning's gone wrong!")
  }
  
  message(". Connecting to root")
  in_degrees <- degree(graphD, mode = "in")
  zeroes = which(in_degrees == 0)
  root.name = paste0(rep("0", str_length(ancnames[i])), collapse="")
  root.label = which(V(graphD)$name == root.name)
  if(length(root.label) == 0) {
    graphD = add_vertices(graphD, n = 1)
    V(graphD)$name[length(V(graphD)$name)] = root.name
    root.label = length(V(graphD)$name)
  }
  for(this.zero in zeroes) {
    if(this.zero != root.label) {
      graphD = add_edges(graphD, c(root.label, this.zero))
    }
  }
  
  message(". Wrapping up")
  graphD.layers = sapply(V(graphD)$name, str_count, "1")
  
  dataset = data.frame(ancestors = ancnames,
                       descendants = descnames )
  
  rlist = s.a
  rlist$dataset = dataset
  rlist$best.graph = graphD
  rlist$best.bc = branching.count(graphD)
  
  return(rlist)
}

# plot the unrewired arboresence (output from Gutin algorithm) and the rewired version increasing layer sum
plot.stage.1 = function(graphs) {
  # extract graphs from algorithm output and assign nodes to layers 
  graphB = graphs$raw.graph
  new.graphB = graphs$rewired.graph
  graphB.layers = sapply(V(graphB)$name, str_count, "1")
  new.graphB.layers = sapply(V(new.graphB)$name, str_count, "1")
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
    ggarrange(
      ggraph(graphB, layout="sugiyama", layers=graphB.layers) + 
        geom_edge_link(color="#CCCCCC") + 
        geom_node_text(aes(label=name), angle=45, hjust=0, size=label.size) + 
        scale_x_continuous(expand = c(0.1, 0.1)) +
        scale_y_continuous(expand = c(0.1, 0.1)) +
        theme_graph(),
      ggraph(new.graphB, layout="sugiyama", layers=new.graphB.layers) + 
        geom_edge_link(color="#CCCCCC") + 
        geom_node_text(aes(label=name), angle=45, hjust=0, size=label.size) + 
        scale_x_continuous(expand = c(0.1, 0.1)) +
        scale_y_continuous(expand = c(0.1, 0.1)) +
        theme_graph(),
      labels = c("A", "B")
    )
  )
}

# plot the simplest spanning arborescence and simplest transition-spanning DAG
plot.stage.2 = function(graphs) {
  # extract graphs from algorithm output and assign layers
  graphB = graphs$rewired.graph
  graphD = graphs$best.graph
  graphB.layers = sapply(V(graphB)$name, str_count, "1")
  graphD.layers = sapply(V(graphD)$name, str_count, "1")
  L = graphs$len
  
  # decide on label size (heuristic, for clarity)
  label.size = 3
  if(L > 5) {
    label.size = 2
  }
  if(L > 20) {
    label.size = 1
    V(graphB)$name = 1:length(V(graphB))
    V(graphD)$name = 1:length(V(graphD))
  }
  
  return(
    ggarrange( ggraph(graphB, layout="sugiyama", layers=graphB.layers) + 
                 geom_edge_link(color="#CCCCCC") + 
                 geom_node_text(aes(label=name), size=label.size, angle=45, hjust=0) + 
                 ggtitle(branching.count(graphB)) + scale_x_continuous(expand = c(0.1, 0.1)) +
                 theme_graph(),
               ggraph(graphD, layout="sugiyama", layers=graphD.layers) + 
                 geom_edge_link(color="#CCCCCC") + 
                 geom_node_text(aes(label=name), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) + 
                 ggtitle(branching.count(graphD)) + scale_x_continuous(expand = c(0.1, 0.1)) +
                 theme_graph(),
               #ggtexttable(dataset, theme=ttheme(base_size=12)), nrow=1, labels=c("A", "B", "C"), 
               nrow = 1, labels=c("A", "B")
               )
  )
}

# plot the simplest spanning arborescence and simplest transition-spanning DAG
plot.stage.p = function(graphD) {
  # extract graphs from algorithm output and assign layers
  graphD.layers = sapply(V(graphD)$name, str_count, "1")
  L = str_length(V(graphD)$name[1])
  
  # decide on label size (heuristic, for clarity)
  label.size = 3
  if(L > 5) {
    label.size = 2
  }
  if(L > 20) {
    label.size = 1
    V(graphD)$name = 1:length(V(graphD))
  }
  
  return(
    ggraph(graphD, layout="sugiyama", layers=graphD.layers) + 
      geom_edge_link(color="#CCCCCC") + 
      geom_node_text(aes(label=name), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) + 
      ggtitle(paste0("B = ", branching.count(graphD), collapse="")) + scale_x_continuous(expand = c(0.1, 0.1)) +
      theme_graph()
  )
}

# untested: plot the simplest graphs embedded in the full hypercube
plot.full.cube = function(graphs) {
  L = graphs$len
  graphB = graphs$rewired.graph
  graphC = graphs$best.graph
  
  # construct full hypercube for comparison
  pow2 = 2**((L-1):0)
  am = matrix(ncol=2)
  # produce list of decimal edges
  for(i in 1:(2**L-1)) {
    anc = DecToBinV(i-1, len=L)
    to.1 = which(anc == 0)
    for(j in 1:length(to.1)) {
      desc = i-1+pow2[to.1[j]]
      am = rbind(am, c(i-1, desc))
    }
  }
  
  # convert to graph with binary labels
  ambin = apply(am[2:nrow(am),], c(1,2), DecToBin, len=L)
  graphO <- graph_from_edgelist(as.matrix(ambin),directed = T)
  graphO = graph.union(graphO, graphB, graphC)
  graphO.layers = sapply(V(graphO)$name, str_count, "1")
  
  # figure out which edges in the complete hypercube are those that we found in graph B/C
  rowsB <- apply(get.edgelist(graphB), 1, paste, collapse = ",")
  rowsC <- apply(get.edgelist(graphC), 1, paste, collapse = ",")
  rowsO <- apply(get.edgelist(graphO), 1, paste, collapse = ",")
  common_rows_B <- intersect(rowsB, rowsO)
  common_rows_C <- intersect(rowsC, rowsO)
  indexes_B <- which(rowsO %in% common_rows_B)
  indexes_C <- which(rowsO %in% common_rows_C)
  
  # set a variable for those edges included in graphC
  E(graphO)$skeleton_B = E(graphO)$skeleton_C = 0
  E(graphO)[indexes_B]$skeleton_B = 1
  E(graphO)[indexes_C]$skeleton_C = 1
  
  return(
    ggarrange(
      ggraph(graphO) + 
        geom_edge_link(aes(edge_color=factor(skeleton_B), edge_alpha=factor(skeleton_B))) + 
        scale_edge_alpha_manual(values=c("0"=1/L, "1"=1)) + 
        #  scale_edge_color_manual(values=c("0"="#EEEEEE", "1"="black")) + 
        theme_graph() + 
        theme(legend.position="none"),
      ggraph(graphO) + 
        geom_edge_link(aes(edge_color=factor(skeleton_C), edge_alpha=factor(skeleton_C))) + 
        scale_edge_alpha_manual(values=c("0"=1/L, "1"=1)) + 
        #  scale_edge_color_manual(values=c("0"="#EEEEEE", "1"="black")) + 
        theme_graph() + 
        theme(legend.position="none")
    )
  )
}

transitions.spanned = function(g, ancnames, descnames) {
  edges.ok = TRUE
  for(i in 1:length(ancnames)) {
    anc.name = ancnames[i]
    desc.name = descnames[i]
    anc.ref = which(V(g)$name==anc.name)
    desc.ref = which(V(g)$name==desc.name)
    if(!is.finite(distances(g, anc.ref, desc.ref, mode="out"))) {
      message(paste0("- Lost connection ", anc.name, " -> ", desc.name))
      edges.ok = FALSE
    }
  }
  return(edges.ok)
}
