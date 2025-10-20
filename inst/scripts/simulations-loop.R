library(ape)
library(phangorn)
library(parallel)
library(hypertrapsct)
library(ggpubr)

set.seed(1)

binary_strings_with_k_ones <- function(n, k) {
  if (k > n || k < 0) return(character(0))
  idx <- combn(n, k, simplify = FALSE)
  sapply(idx, function(pos) {
    bits <- rep("0", n)
    bits[pos] <- "1"
    paste(bits, collapse = "")
  })
}

L = 10
dyn.set = c("linear", "random", "max.spread")
tree.size = 128

birth.rate = 1
death.rate = 0.5

solns = list()
ancnames = descnames = list()
fits = data.frame()

for(L in c(3, 5, 7, 9)) {
  for(tree.size in c(32, 64, 128)) {
    for(dynamics in dyn.set) {
      if(dynamics == "random") {
        n.seed = 5
      } else {
        n.seed = 1
      }
      for(seed in 1:n.seed) {
        set.seed(seed)
        expt.label = paste0(L, ".", tree.size, ".", dynamics, ".", seed, collapse="")
        message(expt.label)
        ancnames[[expt.label]] = c()
        descnames[[expt.label]] = c()

        if(dynamics == "max.spread") {
          descnames[[expt.label]] = rep(binary_strings_with_k_ones(L, floor(L/2)), length.out = tree.size-1)
          ancnames[[expt.label]] = rep(binary_strings_with_k_ones(L, 0), length.out = tree.size-1)
        } else {
          # accumulation rate for features (and loss rate, for reversible setup)
          accumulation.rate = 1
          # create random phylogeny with tree.size nodes from birth-death process parameterised as above
          my.tree = rphylo(tree.size, birth=birth.rate, death=death.rate)
          my.tree$node.label = as.character(1:my.tree$Nnode)
          tree.labels = c(my.tree$tip.label, my.tree$node.label)

          # generate state for all nodes traversing tree breadth-first
          # and setting the state of the child nodes according to
          # accumulation rate and branch length
          my.root = getRoot(my.tree)
          to.do = c(my.root)
          # initialise state list
          x = list()
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
                }

                anc.str = paste0(x[[i]], collapse="")
                desc.str = paste0(x[[this.child]], collapse="")
                #cat(anc.str, "->", desc.str, "\n")
                ancnames[[expt.label]] = c(ancnames[[expt.label]], anc.str)
                descnames[[expt.label]] = c(descnames[[expt.label]], desc.str)
                # add this child to to state list, and to next iteration's to-do
                new.to.do = c(new.to.do, this.child)
              }
            }
            # update to-do list
            to.do = new.to.do
          }
        }
        solns[[expt.label]] = simplest_DAG(ancnames[[expt.label]], descnames[[expt.label]])
        fits = rbind(fits, cbind(data.frame(label=expt.label),
                                 fit_properties(solns[[expt.label]], verbose=FALSE)))
      }
    }
  }
}

sol = simplest_DAG(ancnames[["linear"]], descnames[["linear"]])
plot_stage_p(sol$best.graph)
sol = simplest_DAG(ancnames[["random"]], descnames[["random"]])
plot_stage_p(sol$best.graph)
sol = simplest_DAG(ancnames[["max.spread"]], descnames[["max.spread"]])
plot_stage_p(sol$best.graph)
fit_properties(sol)
ggarrange(plot_stage_1(solns[[1]]), plot_stage_1(solns[[2]]), nrow=2)
