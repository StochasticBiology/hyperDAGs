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

random_binary_strings <- function(n, nstr) {
  rm = matrix(round(runif(n*nstr)), ncol=n)
  return(apply(rm, 1, paste0, collapse=""))
}

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

star_tree = function(n) {
  # optional: tip labels
tip.labels <- paste0("t", 1:n)

# create a star tree
star.tree <- ape::stree(n = n, type = "star")
star.tree$tip.label <- tip.labels
return(star.tree)
}

plot_tree_data = function(my.tree, tip.data) {
  data.m = do.call(rbind, tip.data)[1:length(my.tree$tip.label),]
  rownames(data.m) = my.tree$tip.label
  colnames(data.m) = 1:ncol(data.m)
  g.core = ggtree::ggtree(my.tree)
  this.plot = ggtree::gheatmap(g.core, as.data.frame(data.m), low="white", high="#AAAAAA", colnames=FALSE) +
    theme(legend.position="none")

  return(this.plot)
}

L = 10
dyn.set = c("linear", "random", "mixed", "bilinear", "max.spread", "spread")
tree.size = 128
#dyn.set = "bilinear"
birth.rate = 1
death.rate = 0.5

solns = list()
ancnames = descnames = list()
x.set = tree.set = list()
fits.raw = data.frame()
data.plots = list()

for(L in c(3, 5, 7, 9, 20)) {
  for(tree.size in c(32, 64, 128)) {
    for(dynamics in dyn.set) {
      if(dynamics == "random" | dynamics == "spread") {
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

        if(dynamics == "spread" | dynamics == "max.spread" | dynamics == "bilinear") {
          my.tree = star_tree(tree.size)
          descnames[[expt.label]] = rep(binary_strings_with_k_ones(L, floor(L/2)), length.out = tree.size-1)
          ancnames[[expt.label]] = rep(binary_strings_with_k_ones(L, 0), length.out = tree.size-1)
          if(dynamics == "spread") {
            descnames[[expt.label]][1:floor((tree.size-1)/2)] = random_binary_strings(L, floor((tree.size-1)/2))
          }
          if(dynamics == "bilinear") {
            descnames[[expt.label]] = bilinear_binary_strings(L, tree.size-1)
          }
          x = strsplit(descnames[[expt.label]], split="")
          x[[length(x)+1]] = strsplit(ancnames[[expt.label]][1], split="")[[1]]
          x = lapply(x, as.numeric)
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
                  if(dynamics == "mixed") {
                    if(runif(1) < 0.1) {
                      if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][sample(ref, 1)] = 1}
                    } else {
                      if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref[1]] = 1 }
                    }
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
        fits.raw = rbind(fits.raw, cbind(data.frame(label=expt.label),
                                 fit_properties(solns[[expt.label]], verbose=FALSE)))
        x.set[[expt.label]] = x
        tree.set[[expt.label]] = my.tree
      }
    }
  }
}

save(fits, file="sims-fits.Rdata")
save(solns, file="sims-solns.Rdata")
save(x.set, file="sims-x-set.Rdata")
save(tree.set, file="sims-tree-set.Rdata")


#### pull example for simple demo

name = names(x.set)[grep("5.32.random.1", names(x.set))]
this.x = x.set[[name]]
this.tree = tree.set[[name]]
plot.1 = plot_tree_data(this.tree, this.x)
this.soln = solns[[name]]

this.soln$best.ls = layer_sum(this.soln$best.graph)
this.soln$Sprime = fit_properties(this.soln)$Sprime
plot.2 = plot_stage_gen(this.soln$raw.graph,
                                         label.size = 4) +
  coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift)) +
  labs(caption = paste0("B = ", this.soln$raw.bc, ", LS = ", this.soln$raw.ls, collapse="")) +
  theme(plot.caption = element_text(size = 14, hjust = 0.5))

plot.3 = plot_stage_gen(this.soln$rewired.graph,
                        label.size = 4) +
  coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift)) +
  labs(caption = paste0("B = ", this.soln$rewired.bc, ", LS = ", this.soln$rewired.ls, collapse="")) +
  theme(plot.caption = element_text(size = 14, hjust = 0.5))

plot.4 = plot_stage_gen(this.soln$best.graph,
                        label.size = 4) +
  coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift)) +
  labs(caption = paste0("B = ", this.soln$best.bc, ", S' = ", this.soln$Sprime, " (LS = ", this.soln$best.ls, ")", collapse="")) +
  theme(plot.caption = element_text(size = 14, hjust = 0.5))

simple.plot = ggarrange(plot.1, plot.2, plot.3, plot.4,
          labels=c("A", "B", "C", "D"))

sf = 3
png("sim-simple-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(simple.plot)
dev.off()

#### pull examples for demonstrations of different dynamics

names.demo.plot = names(x.set)[grep("7.32.*1", names(x.set))]

margin.shift = 25
demo.soln.plots = demo.data.plots = list()
for(name in names.demo.plot) {
  this.x = x.set[[name]]
  this.tree = tree.set[[name]]
  demo.data.plots[[name]] = plot_tree_data(this.tree, this.x)
  this.soln = solns[[name]]
  demo.soln.plots[[name]] = plot_stage_gen(this.soln$best.graph,
                                      label.size = 4) +
    coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift))
}

sim.plot = ggarrange(
  ggarrange(plotlist = demo.data.plots, nrow=1),
  ggarrange(plotlist = demo.soln.plots, nrow=1), nrow = 2
)

png("sim-plot-data.png", width=1300*sf, height=600*sf, res=72*sf)
print(sim.plot)
dev.off()

### statistics under different dynamics

fits = fits.raw
fits$type = "random"
fits$type[grep("linear", fits$label)] = "linear"
fits$type[grep("bilinear", fits$label)] = "bilinear"
fits$type[grep("spread", fits$label)] = "spread"
fits$type[grep("mixed", fits$label)] = "mixed"
fits$type[grep("max.spread", fits$label)] = "max.spread"
fits$tree.size =as.numeric(sapply(strsplit(fits$label, "\\."), `[`, 2))


ggarrange(
  ggplot(fits, aes(x=type, y=S, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)),
  ggplot(fits, aes(x=type, y=Sprime, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)),
  ggplot(fits, aes(x=type, y=Sstar, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)),
  nrow=3)

sprime.plot = ggplot(fits, aes(x=type, y=Sprime, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Dynamics", y = "S'", color = "L:", shape = "n:")

sf = 3
png("sprime-plot.png", width=400*sf, height=200*sf, res=72*sf)
print(sprime.plot)
dev.off()

margin.shift = 25
ggarrange(
  plot_stage_gen(solns[["5.128.mixed.1"]]$raw.graph, label.size = 4) +
    coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift)),
  plot_stage_gen(solns[["5.128.mixed.1"]]$rewired.graph, label.size = 4) +
    coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift)),
  plot_stage_gen(solns[["5.128.mixed.1"]]$best.graph, label.size = 4) +
  coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift))
)

ggarrange(
plot_stage_gen(solns[["20.128.mixed.1"]]$raw.graph,
               label.size = 4, label.style = "labels"),
plot_stage_gen(solns[["20.128.mixed.1"]]$rewired.graph,
               label.size = 4, label.style = "labels"),
plot_stage_gen(solns[["20.128.mixed.1"]]$best.graph,
               label.size = 4, label.style = "labels"))

