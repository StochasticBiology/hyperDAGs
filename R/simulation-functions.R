
#' Produce the set of binary strings of length nstr with k 1's
#'
#' @param nstr numeric length of binary strings
#' @param k numeric number of 1s in the strings
#' @return character vector set of strings
#' @examples
#' binary_strings_with_k_ones(5, 3)
#' @export
binary_strings_with_k_ones <- function(nstr, k) {
  if (k > nstr || k < 0) return(character(0))
  idx <- combn(nstr, k, simplify = FALSE)
  sapply(idx, function(pos) {
    bits <- rep("0", nstr)
    bits[pos] <- "1"
    paste(bits, collapse = "")
  })
}

#' Produce a set of n random binary strings of length nstr
#'
#' @param n numeric number of strings to produce
#' @param nstr numeric length of binary strings
#' @return character vector set of strings
#' @examples
#' random_binary_strings(5, 3)
#' @export
random_binary_strings <- function(n, nstr) {
  rm = matrix(round(runif(n*nstr)), ncol=n)
  return(apply(rm, 1, paste0, collapse=""))
}

#' Produce a set of n character strings of length nstr reflecting two distinct pathways
#' "leftwards" and "rightwards" through the bitstring
#'
#' @param n numeric number of strings to produce
#' @param nstr numeric length of binary strings
#' @return character vector set of strings
#' @examples
#' bilinear_binary_strings(10, 4)
#' @export
bilinear_binary_strings <- function(n, nstr) {
  rm = matrix(0, nrow=nstr, ncol=n)
  one.count = 0
  for(i in 0:(floor(nstr/2) - 1)) {
    for(j in 0:one.count) {
      rm[2*i+1,j+1] = 1
      rm[2*i+2,n-j] = 1
    }
    one.count = one.count+1
    if(one.count >= n) { one.count = 0 }
  }
  return(apply(rm, 1, paste0, collapse=""))
}

#' Create a star phylogeny with n tips
#'
#' @param n numeric number of tips
#' @return tree object
#' @examples
#' star_tree(4)
#' @export
star_tree = function(n) {
  # optional: tip labels
  tip.labels <- paste0("t", 1:n)

  # create a star tree
  star.tree <- ape::stree(n = n, type = "star")
  star.tree$tip.label <- tip.labels
  return(star.tree)
}

#' Plot a phylogeny and associated set of bitstrings at the tips
#'
#' @param my.tree tree object
#' @param tip.data matrix or data frame containing one row for each tip
#' @return ggplot object
#' @examples
#' #my.tree = star_tree(4)
#' #tip.data = as.list(matrix(0, nrow=4, ncol=5))
#' #plot_tree_data(my.tree, tip.data)
#' print(999)
#' @export
plot_tree_data = function(my.tree, tip.data) {
  data.m = do.call(rbind, tip.data)[1:length(my.tree$tip.label),]
  rownames(data.m) = my.tree$tip.label
  colnames(data.m) = 1:ncol(data.m)
  g.core = ggtree::ggtree(my.tree)
  this.plot = ggtree::gheatmap(g.core, as.data.frame(data.m), low="white", high="#AAAAAA", colnames=FALSE) +
    theme(legend.position="none")

  return(this.plot)
}

#' Simulate dynamics for L characters on a given tree, or a simulated tree with n tips
#'
#' @param tree.size numeric number of tips on tree. If 0, we assume that a tree is being provided in use.tree
#' @param L numeric number of characters to simulate
#' @param tree.type character type of tree to simulate. "random" is only option for now
#' @param birth.rate numeric birth rate for birth-death process simulating random tree
#' @param death.rate numeric death rate for birth-death process simulating random tree
#' @param accumulation.rate numeric rate of accumulation of features
#' @param dynamics character type of accumulation dynamics to simulate. "random" (all characters random, max one change per branch); "linear" (single pathway dependence); "poisson" (all characters random, Poisson accumulation); "mixed" (mix of linear and random); "heterogeneous" (random, with different relative probabilities in char.probs). default "random"
#' @param use.tree NULL or tree object. if tree.size is 0, the tree provided here will be used instead of simulating a new one
#' @param char.probs numeric vector of relative probabilities of each character acquisition, for "heterogeneous" accumulation
#' @return a list: x (data at nodes and tips); my.tree (tree structure); ancnames (labels of ancestors for each transition); descnames (labels of descendants for each transition)
#' @examples
#' simulate_accumulation(32, 4)
#' @export
# simulate dynamics for L "iid" characters on a tree with n tips
# dynamics can be random, linear, or mixed
# tree type can be random or star
# other parameters take default values but can be specified
simulate_accumulation = function(tree.size,
                                    L,
                                    tree.type = "random",
                                    birth.rate = 1,
                                    death.rate = 0.5,
                                    accumulation.rate = 1,
                                    dynamics = "random",
                                    use.tree = NULL,
                                    char.probs=NULL) {

  # create random phylogeny with tree.size nodes from birth-death process parameterised as above
  if(tree.size == 0) {
    my.tree = use.tree
  } else {
    my.tree = ape::rphylo(tree.size, birth=birth.rate, death=death.rate)
    my.tree$node.label = as.character(1:my.tree$Nnode)
    tree.labels = c(my.tree$tip.label, my.tree$node.label)
  }

  # generate state for all nodes traversing tree breadth-first
  # and setting the state of the child nodes according to
  # accumulation rate and branch length
  my.root = phangorn::getRoot(my.tree)
  to.do = c(my.root)
  # initialise state list
  x = list()
  ancnames = c()
  descnames = c()
  x[[my.root]] = rep(0,L)
  # while we still have vertices to simulate
  while(length(to.do) > 0) {
    # initialise a new to-do list for next iteration
    new.to.do = c()
    # loop through each node in current to-do list
    for(i in to.do) {
      this.outgoing.edges = which(my.tree$edge[,1] == i)
      # loop over this node's children
      for(j in this.outgoing.edges) {
        this.child = my.tree$edge[j,2]
        this.branch.length = my.tree$edge.length[j]
        # construct state for this child based on its parent
        x[[this.child]] = x[[i]]
        ref = which(x[[this.child]] == 0)
        if(length(ref) > 0) {
          if(dynamics == "linear") {
            # find leftmost zero in current state, and change with some probability
            # recall dynamics here are 00000 -> 10000 -> 11000 -> 11100 -> 11110 -> 11111
            if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref[1]] = 1 }
          }
          if(dynamics == "random") {
            if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][sample(ref, 1)] = 1}
          }
          if(dynamics == "heterogeneous") {
            char.hazards = accumulation.rate*this.branch.length*char.probs[ref]
            char.hazards[char.hazards > 1] = 1
            these.gains = rbinom(length(ref), 1, char.hazards)
            x[[this.child]][these.gains==1] = 1
          }
          if(dynamics == "poisson") {
            n.accum = rpois(1, accumulation.rate*this.branch.length)
            x[[this.child]][sample(ref, n.accum, replace=TRUE)] = 1
          }
          if(dynamics == "mixed") {
            if(runif(1) < 0.1) {
              if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][sample(ref, 1)] = 1}
            } else {
              if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref[1]] = 1 }
            }
          }
          anc.str = paste0(x[[i]], collapse="")
          desc.str = paste0(x[[this.child]], collapse="")
          #cat(anc.str, "->", desc.str, "\n")
          ancnames = c(ancnames, anc.str)
          descnames = c(descnames, desc.str)
          # add this child to to state list, and to next iteration's to-do
          new.to.do = c(new.to.do, this.child)
        }
      }
      # update to-do list
      to.do = new.to.do
    }
  }
  return(list(x=x,
              my.tree=my.tree,
              ancnames=ancnames,
              descnames=descnames))
}
