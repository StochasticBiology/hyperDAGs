library(ape)
library(phangorn)
library(parallel)
library(hypertrapsct)
library(ggpubr)

set.seed(1)

L = 10
dyn.set = c("linear", "random", "mixed", "bilinear", "max.spread", "spread")
star.phylo = c("spread", "max.spread", "bilinear")
tree.size = 128
birth.rate = 1
death.rate = 0.5

solns = list()
ancnames = descnames = list()
x.set = tree.set = list()
fits.raw = data.frame()
data.plots = list()

run.sims = FALSE

if(run.sims == TRUE) {
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

          # the dynamics we'll "simulate" on a star phlogeny
          if(dynamics %in% star.phylo) {
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
            # otherwise, we'll simulate a random phylogeny with these dynamics
            r.set = simulate_accumulation(tree.size, L, dynamics=dynamics)
            x = r.set[["x"]]
            my.tree = r.set[["my.tree"]]
            ancnames[[expt.label]] = r.set[["ancnames"]]
            descnames[[expt.label]] = r.set[["descnames"]]
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

  save(fits, file="sims-fits.Rdata")
  save(solns, file="sims-solns.Rdata")
  save(x.set, file="sims-x-set.Rdata")
  save(tree.set, file="sims-tree-set.Rdata")
} else {
  sims.fits.file = system.file("extdata", "sims-fits.Rdata", package="hyperdags")
  load(sims.fits.file)
  sims.solns.file = system.file("extdata", "sims-solns.Rdata", package="hyperdags")
  load(sims.solns.file)
  sim.x.set.file = system.file("extdata", "sims-x-set.Rdata", package="hyperdags")
  load(sim.x.set.file)
  sim.tree.set.file = system.file("extdata", "sims-tree-set.Rdata", package="hyperdags")
  load(sim.tree.set.file)
}

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
### statistics under different dynamics

fits = fits.raw
fits$type = "random"
fits$type[grep("linear", fits$label)] = "linear"
fits$type[grep("bilinear", fits$label)] = "bilinear"
fits$type[grep("spread", fits$label)] = "spread"
fits$type[grep("mixed", fits$label)] = "mixed"
fits$type[grep("max.spread", fits$label)] = "max.spread"
fits$tree.size =as.numeric(sapply(strsplit(fits$label, "\\."), `[`, 2))

level.order = c("linear", "bilinear", "mixed", "random", "spread", "max.spread")

#names.demo.plot = names(x.set)[grep("7.32.*1", names(x.set))]
names.demo.plot = paste0("7.32.", level.order, ".1")

margin.shift = 25
demo.soln.plots = demo.data.plots = list()
for(name in names.demo.plot) {
  this.x = x.set[[name]]
  this.tree = tree.set[[name]]
  demo.data.plots[[name]] = plot_tree_data(this.tree, this.x)
  this.soln = solns[[name]]
  if(grepl("linear", name)) {
    demo.soln.plots[[name]] = plot_stage_gen(this.soln$best.graph,
                                             label.size = 4) +
      coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift))
  } else {
    demo.soln.plots[[name]] = plot_stage_gen(this.soln$best.graph,
                                             label.style = "points",
                                             label.size = 2) +
      coord_cartesian(clip = "off") + theme(plot.margin = margin(margin.shift, margin.shift, margin.shift, margin.shift))
  }
}

sim.plot = ggarrange(
  ggarrange(plotlist = demo.data.plots, nrow=1),
  ggarrange(plotlist = demo.soln.plots, nrow=1), nrow = 2
)

png("sim-plot-data.png", width=1000*sf, height=600*sf, res=72*sf)
print(sim.plot)
dev.off()

#########

png("sim-plot-data-alls.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange(
  ggplot(fits, aes(x=type, y=S, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)) ,#+ theme(legend.position="none"),
  ggplot(fits, aes(x=type, y=Sprime, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)),
  ggplot(fits, aes(x=type, y=Sstar, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5)),# + theme(legend.position="none"),
  nrow=3)
dev.off()

sprime.plot = ggplot(fits, aes(x=factor(type, levels=level.order), y=Sprime, color=factor(L), shape=factor(tree.size))) + geom_point(position = position_dodge(width = 0.5), size=4) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(x = "Dynamics", y = "S'", color = "L:", shape = "n:")

sf = 3
png("sprime-plot.png", width=400*sf, height=200*sf, res=72*sf)
print(sprime.plot)
dev.off()

both.plot = ggarrange(
  ggarrange(plotlist = demo.data.plots, nrow=1, labels=c("A", "B", "C", "D", "E", "F")),
  ggarrange(plotlist = demo.soln.plots, nrow=1),
  sprime.plot+theme_minimal(base_size=20)+
    theme(legend.position="bottom"), labels=c("", "", "G"),
  heights=c(1,1,1.5), nrow = 3
)

sf = 3
png("sprime-both-plot.png", width=1000*sf, height=700*sf, res=72*sf)
print(both.plot)
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

