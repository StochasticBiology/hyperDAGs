#' Plot the unrewired arboresence (output from Gutin algorithm) and the rewired version increasing layer sum
#'
#' @param soln list output of model fit
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_stage_1(sol)
#' @export
plot_stage_1 = function(soln) {
  # extract graphs from algorithm output and assign nodes to layers
  graphB = soln$raw.graph
  new.graphB = soln$rewired.graph
  graphB.layers = sapply(igraph::V(graphB)$name, stringr::str_count, "1")
  new.graphB.layers = sapply(igraph::V(new.graphB)$name, stringr::str_count, "1")
  L = soln$len

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
                         ggplot2::labs(caption = paste0("B = ", branching_count(graphB), collapse="")) +
                         ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
                         ggplot2::theme_void(), #ggraph::theme_graph(),
                       ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
                         ggraph::geom_edge_link(color="#CCCCCC") +
                         ggraph::geom_node_text(ggplot2::aes(label=name), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) +
                         ggplot2::labs(caption = paste0("B = ", branching_count(graphD), collapse="")) +
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
    ggplot2::labs(caption = paste0("B = ", branching_count(graphD), collapse="")) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_void() #ggraph::theme_graph()
  if(!is.na(v.labels$Species[1])) {
    g = g +
      ggraph::geom_node_text(ggplot2::aes(label=v.label), size=2, angle=45, hjust=0)  #, check_overlap = TRUE) +
  }
  return(g)
}

#' Plot a solution graph with control over styling
#'
#' @param graphD graph to plot
#' @param v.labels data frame of vertex labels (optional)
#' @param label.size numeric node label size (default chooses according to length)
#' @param point.size numeric node point size (default 0)
#' @param edge.alpha numeric edge transparency (default 1)
#' @param label.style character node plot style ("full" binary strings, "points" just points, or "labels" individual codes for nodes, default chooses according to length)
#' @return ggarrange object containing plots
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_stage_gen(sol$best.graph)
#' @export
plot_stage_gen = function(graphD,
                          v.labels = data.frame(Species=NA),
                          label.size = 0,
                          point.size = 0,
                          edge.alpha = 1,
                          label.style = 0) {
  # extract graphs from algorithm output and assign layers
  graphD.layers = sapply(igraph::V(graphD)$name, stringr::str_count, "1")
  L = stringr::str_length(igraph::V(graphD)$name[1])

  # decide on label size (heuristic, for clarity)
  igraph::V(graphD)$plotname = igraph::V(graphD)$name
  if(L > 5) {
    if(label.size == 0) {
      label.size = 2
    }
    if(label.style == 0) {
      igraph::V(graphD)$plotname = igraph::V(graphD)$name
    }
  }
  if(L > 20) {
    if(label.size == 0) {
      label.size = 1
    }
    if(label.style == 0) {
      igraph::V(graphD)$plotname = 1:length(igraph::V(graphD))
    }
  }

  if(label.style == "full") {
    igraph::V(graphD)$plotname = igraph::V(graphD)$name
  }
  if(label.style == "labels") {
    igraph::V(graphD)$plotname = 1:length(igraph::V(graphD))
  }
  if(label.style == "points") {
    igraph::V(graphD)$plotname = rep("",length(igraph::V(graphD)))
  }
  if(!is.na(v.labels$Species[1])) {
    igraph::V(graphD)$v.label = ""
    for(i in 1:nrow(v.labels)) {
      ref = which(igraph::V(graphD)$name == v.labels$label[i])
      igraph::V(graphD)$v.label[ref] = v.labels$Species[i]
    }
  }
  g = ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
    ggraph::geom_edge_link(color="#CCCCCC", alpha = edge.alpha) +
    ggraph::geom_node_text(ggplot2::aes(label=plotname), size=label.size, angle=45, hjust=0) +
    # ggplot2::ggtitle(paste0("B = ", branching_count(graphD), collapse="")) +
    #  ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_void() #ggraph::theme_graph()
  if(!is.na(v.labels$Species[1])) {
    g = g +
      ggraph::geom_node_text(ggplot2::aes(label=v.label), size=label.size, angle=45, hjust=0)  #, check_overlap = TRUE) +
  }
  if(label.style == "points" & point.size > 0) {
    g = g + ggraph::geom_node_point(size=point.size)
  }
  return(g)
}

#' Plot a solution graph with edges weighted by downstream datapoints
#'
#' @param graphD graph to plot
#' @param labels vector of vertex labels (optional)
#' @param thresh numeric threshold below which edges are not plotted (default 5)
#' @param edge.label.size numeric edge label size (default 3)
#' @param check.overlaps Boolean, check for label overlaps (default FALSE)
#' @param orient string, edge label orientation "across", "along" (default "across")
#' @return list object containing a collection of plots: full and thresholded plots, and full and thresholded graphs
#' @examples
#' sol = simplest_DAG(c("000", "010", "100"), c("110", "110", "101"))
#' plot_weights(sol$best.graph, thresh = 0)
#' @export
plot_weights = function(graphD, labels=c(""), thresh=5, edge.label.size=3, check.overlaps=FALSE, orient = "across") {
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
    igraph::E(graphD)$label[i] = paste0("+", paste0(labels[which(srcs[[i]]!=dests[[i]])], collapse="\n"), collapse="")
  }
  # decide on label size (heuristic, for clarity)
  #edge.label.size = 3
  label.size = 3
  igraph::V(graphD)$plotname = igraph::V(graphD)$name
  if(L > 5) {
    label.size = 2
    igraph::V(graphD)$plotname = igraph::V(graphD)$name
    #  edge.label.size = 3
  }
  if(L > 20) {
    label.size = 1
    igraph::V(graphD)$plotname = 1:length(igraph::V(graphD))
    #  edge.label.size = 2
  }

  full.plot = ggraph::ggraph(graphD, layout="sugiyama", layers=graphD.layers) +
    ggraph::geom_edge_link(color="#CCCCCC", label_size= edge.label.size, angle_calc = orient, label_colour = "grey", check_overlap = check.overlaps, ggplot2::aes(label=label, alpha=log(as.numeric(thickness)+1))) +
    ggraph::geom_node_text(ggplot2::aes(label=plotname), size=label.size, angle=45, hjust=0  , check_overlap = check.overlaps) +
    ggplot2::ggtitle(paste0("B = ", branching_count(graphD), collapse="")) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_void() #ggraph::theme_graph()

  graphE = igraph::delete_edges(graphD, which(as.numeric(igraph::E(graphD)$thickness) < thresh))
  degrees = igraph::degree(graphE)
  graphE = igraph::delete_vertices(graphE, which(degrees == 0))
  graphE.layers = sapply(igraph::V(graphE)$name, stringr::str_count, "1")

  label.size = 0
  thresh.plot = ggraph::ggraph(graphE, layout="sugiyama", layers=graphE.layers) +
    ggraph::geom_edge_link(color="#AAAAFF", label_size= edge.label.size, angle_calc = orient, label_colour = "black",check_overlap = check.overlaps,
                           ggplot2::aes(label=label, edge_width = as.numeric(thickness), alpha=as.numeric(thickness))) +
    ggraph::geom_node_text(ggplot2::aes(label=name), size=label.size, angle=45, hjust=0, check_overlap = check.overlaps) +
    ggraph::scale_edge_width(limits=c(0,NA)) +
    ggplot2::theme_void() + #ggraph::theme_graph() +
    ggplot2::theme(legend.position="none")

  return(list(full.plot=full.plot, thresh.plot=thresh.plot,
              full.graph=graphD, thresh.graph=graphE))
}
