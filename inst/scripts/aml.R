library(hypertrapsct)
library(igraph)
library(ggtree)
library(ggraph)
library(stringr)

load("~/Dropbox/Documents/2024_Projects/HyperEvol/Workspaces/cancer-mar-8.RData")
plotHypercube.motifs(parallelised.runs[[4]], 
                     label.size = 3,
                     featurenames = AML[[4]],
                     label.scheme = "sparse") + theme(legend.position="none")
g.cancer.graph2 = plotHypercube.sampledgraph2(parallelised.runs[[4]], use.arc = FALSE, featurenames = AML[[4]], 
                                              edge.label.size=3, edge.label.angle = "along", node.labels=FALSE,
                                              no.times=TRUE, small.times=FALSE, thresh=0.008, truncate=5,
                                              use.timediffs = FALSE, edge.check.overlap = FALSE) +
  theme(legend.position="none") + coord_flip() + scale_y_reverse()
g.cancer.graph2
