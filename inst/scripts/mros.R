library(hypertrapsct)
library(phytools)
library(ggpubr)

### reproduce MRO paper figure
mro.df = read.csv(system.file("extdata/mro-barcodes-2025-2.csv", package="hyperdags"))
mro.df[,2:ncol(mro.df)] = 1-mro.df[,2:ncol(mro.df)]
mro.tree = read.tree(system.file("extdata/mro-ncbi-tree-2025.nwk", package="hyperdags"))
mro.tree$tip.label = gsub("_", " ", mro.tree$tip.label)
mro.c = curate.tree(mro.tree, mro.df)
#mro.src = apply(mro.c$srcs, 1, paste0, collapse="")
#mro.dest = apply(mro.c$dests, 1, paste0, collapse="")
mro.dest = apply(unique(rbind(mro.c$srcs, mro.c$dests)), 1, paste0, collapse="")
mro.src = rep("000000000", length(mro.dest))

### more rational?
mro.df = read.csv(system.file("extdata/mro-barcodes-2025-1.csv", package="hyperdags"))
mro.df[,2:ncol(mro.df)] = 1-mro.df[,2:ncol(mro.df)]
mro.tree = read.tree(system.file("extdata/mro-ncbi-tree-2025.nwk", package="hyperdags"))
mro.tree$tip.label = gsub("_", " ", mro.tree$tip.label)
mro.c = curate.tree(mro.tree, mro.df)
mro.src = apply(mro.c$srcs, 1, paste0, collapse="")
mro.dest = apply(mro.c$dests, 1, paste0, collapse="")
#mro.dest = apply(unique(rbind(mro.c$srcs, mro.c$dests)), 1, paste0, collapse="")
#mro.src = rep("000000000", length(mro.dest))


mro.soln.real = simplest_DAG(mro.src, mro.dest)
fit_properties(mro.soln.real)
plot_stage_gen(mro.soln.real$best.graph)

####### heterogeneous null model

mro.char.probs = colSums(mro.df[,2:ncol(mro.df)])/nrow(mro.df)

mro.n1s = c()
mro.r.sim.df = data.frame()
for(i in 1:10) {
mro.r.set = simulate_accumulation(0, 9, use.tree=mro.c$tree,
                              accumulation.rate = 0.08,
                              dynamics="heterogeneous",
                              char.probs = mro.char.probs)
mro.n1s = c(mro.n1s, sum(unlist(mro.r.set$x)))
mro.soln = simplest_DAG(mro.r.set[["ancnames"]], mro.r.set[["descnames"]])
mro.r.sim.df = rbind(mro.r.sim.df, fit_properties(mro.soln))
}
sum(mro.df[,2:ncol(mro.df)])
mro.n1s

mro.plot = ggarrange(plotHypercube.curated.tree(mro.c),
          plot_tree_data(mro.r.set$my.tree, mro.r.set$x),
          plot_stage_gen(mro.soln.real$best.graph, label.style = "points", label.size=1) +
            labs(caption=paste0("B = 10, S' = 0.77", collapse="")) +
            theme(plot.caption = element_text(hjust = 0.5)),
          plot_stage_gen(mro.soln$best.graph, label.style="points", label.size=1) +
            labs(caption=paste0("B = 33 ± 4, S' = 0.57 ± 0.04", collapse="")) +
            theme(plot.caption = element_text(hjust = 0.5)),
          labels = c("A", "B"))

sf = 2
png("mro-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(mro.plot)
dev.off()

fit_properties(mro.soln.real)
mro.r.sim.df

########## homogeneous

mro.0.n1s = c()
mro.0.r.sim.df = data.frame()
for(i in 1:10) {
  mro.0.r.set = simulate_accumulation(0, 9, use.tree=mro.c$tree,
                                    accumulation.rate = 0.4,
                                    dynamics="poisson")
  mro.0.n1s = c(mro.0.n1s, sum(unlist(mro.0.r.set$x)))
  mro.0.soln = simplest_DAG(mro.0.r.set[["ancnames"]], mro.0.r.set[["descnames"]])
  mro.0.r.sim.df = rbind(mro.0.r.sim.df, fit_properties(mro.0.soln))
}
sum(mro.df[,2:ncol(mro.df)])
mro.0.n1s

ggarrange(plotHypercube.curated.tree(mro.c),
          plot_tree_data(mro.0.r.set$my.tree, mro.0.r.set$x),
          plot_stage_gen(mro.soln.real$best.graph),
          plot_stage_gen(mro.0.soln$best.graph))

fit_properties(mro.soln)
mro.0.r.sim.df

cat("Solution S' = ", fit_properties(mro.soln.real)$Sprime, "\n",
    " compared to ", round(mean(mro.0.r.sim.df$Sprime), digits=2),
    " +- ", round(sd(mro.0.r.sim.df$Sprime), digits=2), " homogeneous\n",
    " and ", round(mean(mro.r.sim.df$Sprime), digits=2),
    " +- ", round(sd(mro.r.sim.df$Sprime), digits=2), " heterogeneous\n")

