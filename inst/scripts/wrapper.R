library(hyperdags)
library(reshape2)
library(ggrepel)
library(ggplotify)
library(ggpubr)

# functions for reading data and plotting
srcfile <- system.file("scripts", "get_data.R", package = "hyperdags")
source(srcfile)

#### basic example for demonstration

# simple illustrative dataset
ancnames  = c("10000","01000","00000","10010","01011","10110")
descnames = c("11110","11110","01011","11110","01111","11111")
# get STSD
s.dag = simplest_DAG(ancnames, descnames)
# plot and summarise: A is SSA, B is STSD, numbers are branching counts
plot_stage_2(s.dag)
fit_properties(s.dag)
transitions_spanned(s.dag$best.graph, ancnames, descnames)

#### real experiments
# runtime ~2h on modern Mac

sf = 3

expt.index = 0
expt.out = list()

# mtDNA case study takes some minutes; ptDNA will take more
for(expt in c("inline", "TBsimp", "TB", "CGH", "cancer", "mtDNA", "ptDNA")) {
  if(expt == "mtDNA") { dset = get_mtDNA_data() }
  if(expt == "ptDNA") { dset = get_ptDNA_data() }
  if(expt == "cancer") { dset = get_cancer_data() }
  if(expt == "inline") { dset = get_inline_data() }
  if(expt == "TB") { dset = get_TB_data() }
  if(expt == "TBsimp") { dset = get_TB_data(simple=TRUE) }
  if(expt == "CGH") { dset = get_CGH_data() }

  L = dset[["L"]]
  ancnames = dset[["ancnames"]]
  descnames = dset[["descnames"]]

  dset.trans = unique(data.frame(Ancestor=ancnames, Descendant=descnames))
  states = unique(c(ancnames, descnames))
  dset.states = data.frame(States=states)

  # for the simplified TB experiment, store data for plotting later
  if(expt == "TBsimp") { tbsimp.dset.states = dset.states; tbsimp.dset.trans = dset.trans }

  # *** CENTRAL LINE ***: apply the algorithm to find our solution
  s.dag = simplest_DAG(ancnames, descnames)

  # verify that our solution works
  if(length(descnames) != 0) {
    spanned = transitions_spanned(s.dag$best.graph, ancnames, descnames)
    if(spanned == TRUE) {
      message("-- transitions verified")
    } else {
      message("-- transitions *NOT* verified")
    }
  }

  # output visualisations of the solution and steps along the way
  sf = 3
  png(paste0("stage-1-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot_stage_1(s.dag))
  dev.off()
  png(paste0("stage-2-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot_stage_2(s.dag))
  dev.off()

  expt.out[[expt]] = s.dag
}

### the following section produces graphical output for the article

# some options for graphical styling in output
scale.x = scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))
scale.y = scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
title.style = theme(plot.title = element_text(
  size = 14,     # Change font size
  family = "Arial",
  face = "plain", # Change font face to bold
  hjust = 0.5,   # Center align the title
  vjust = 0 #1.5    # Adjust vertical alignment
))


  tb.simp.data =  ggarrange(ggtexttable(tbsimp.dset.states, rows=NULL),
                            ggarrange(ggtexttable(dset.trans[1:(nrow(tbsimp.dset.trans)/3),], rows=NULL),
                                      ggtexttable(dset.trans[(nrow(tbsimp.dset.trans)/3+1):(2*nrow(tbsimp.dset.trans)/3),], rows=NULL),
                                      ggtexttable(dset.trans[(2*nrow(tbsimp.dset.trans)/3+1):(nrow(tbsimp.dset.trans)),], rows=NULL),
                                      nrow=1),
                            nrow = 1, labels=c("A", "B"), widths=c(0.3,1))
  png(paste0("stage-0-", expt, ".png", collapse=""), width=600*sf, height=480*sf, res=72*sf)
  print(tb.simp.data)
  dev.off()

# simplified TB plot
png("fig-tbsimp.png", width=850*sf, height=700*sf, res=72*sf)
print(ggarrange(
  tb.simp.data,
  ggarrange(plot_stage_p(expt.out[[2]]$raw.graph) + scale.y + title.style,
                plot_stage_p(expt.out[[2]]$rewired.graph) + scale.y + title.style,
                plot_stage_p(expt.out[[2]]$best.graph) + scale.y + title.style,
                nrow = 1,
                labels=c("C", "D", "E")),
  nrow=2, heights=c(1,0.66))
)
dev.off()

# full TB plot
png("fig-tb.png", width=1200*sf, height=500*sf, res=72*sf)
w = plot_weights(expt.out[[3]]$best.graph,
                 labels = c("INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP",
                            "MOX", "OFL", "PRO"))
print(ggarrange(plot_stage_p(expt.out[[3]]$best.graph) + scale.y + title.style,
                w$thresh.plot,
                labels=c("C", "D")))
dev.off()

# cancer progression plot. to match previous research we rotate the visualisations here
png("fig-cancer.png", width=1000*sf, height=500*sf, res=72*sf)
w = plot_weights(expt.out[[5]]$best.graph,
                 labels = c("FLT3", "NPM1", "WT1", "DNMT3A", "KRAS", "NRAS", "RUNX1",
                            "IDH1", "IDH2", "PTPN11", "SRSF2", "ASXL1", "BCOR", "STAG2",
                            "TP53", "U2AF1", "SF3B1", "TET2", "CSF3R", "JAK2", "GATA2",
                            "EZH2", "PPM1D", "SETBP1", "KIT", "CBL", "PHF6", "MYC", "ETV6",
                            "MPL", "SMC3"))
print(ggarrange(                 plot_stage_p(expt.out[[5]]$best.graph) + coord_flip() + scale.y + scale_y_reverse() + title.style+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 w$thresh.plot + scale_y_reverse() + coord_flip()+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 w$thresh.plot+coord_flip() + scale_y_reverse() + theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 nrow=1, labels = c("A", "B", "C"))
)
dev.off()

# organelle DNA plot

# highlight some species of interest
mt.orgs = c("homo sapiens", "arabidopsis thaliana", "reclinomonas americana",
         "saccharomyces cerevisiae", "plasmodium falciparum", "physcomitrium patens",
         "triticum aestivum")
mt.raw = read.csv("mt-barcodes-manual.csv")
L = ncol(mt.raw)-1
mt.names = apply(1-mt.raw[,2:(L+1)], 1, paste0, collapse="")
mt.labs = data.frame()
for(i in 1:length(mt.orgs)) {
  mt.labs = rbind(mt.labs, data.frame(Species = mt.orgs[i], label=mt.names[which(mt.raw$Species == mt.orgs[i])]))
}

# same for plastids
pt.orgs = c("arabidopsis thaliana", "physcomitrium patens", "hydnora visseri", "parasitaxus usta",
            "polysiphonia elongata", "gracilaria changii", "cladosiphon okamuranus")
pt.raw = read.csv("pt-barcodes-manual.csv")
L = ncol(pt.raw)-1
pt.names = apply(1-pt.raw[,2:(L+1)], 1, paste0, collapse="")
pt.labs = data.frame()
for(i in 1:length(pt.orgs)) {
  pt.labs = rbind(pt.labs, data.frame(Species = pt.orgs[i], label=pt.names[which(pt.raw$Species == pt.orgs[i])]))
}

# plot the solution, with labelled species highlights
png("fig-odna.png", width=1000*sf, height=1000*sf, res=72*sf)
print(ggarrange(plot_stage_p(expt.out[[6]]$best.graph, v.labels=mt.labs) + title.style,
                plot_stage_p(expt.out[[7]]$best.graph, v.labels=pt.labs) + title.style,
                nrow = 2,
                labels=c("A", "B")))
dev.off()

# additional organelle DNA and cancer plots, reflecting the orderings of individual features

oDNA.types = c("MT", "PT", "AML")
oDNA.g = list()
oDNA.g2 = list()

# loop through these case studies
for(oDNA.expt in 1:3) {
  # pull the solution for this case study
  oDNA = oDNA.types[oDNA.expt]
  if(oDNA == "MT") {
    wg = expt.out[[6]]$best.graph
  } else if(oDNA == "PT") {
    wg = expt.out[[7]]$best.graph
  } else if(oDNA == "AML") {
    wg = expt.out[[5]]$best.graph
  }
  # get the node names in this solution
  these.nodes = ends(wg, E(wg)[1])
  src = strsplit(these.nodes[1,1], "")[[1]]
  L = length(src)
  if(oDNA == "MT") {
    genes = colnames(mt.raw)[2:(L+1)]
  } else if(oDNA == "PT") {
    genes = colnames(pt.raw)[2:(L+1)]
  } else if(oDNA == "AML") {
    genes = c("FLT3", "NPM1", "WT1", "DNMT3A", "KRAS", "NRAS", "RUNX1",
                       "IDH1", "IDH2", "PTPN11", "SRSF2", "ASXL1", "BCOR", "STAG2",
                       "TP53", "U2AF1", "SF3B1", "TET2", "CSF3R", "JAK2", "GATA2",
                       "EZH2", "PPM1D", "SETBP1", "KIT", "CBL", "PHF6", "MYC", "ETV6",
                       "MPL", "SMC3")
  }

  # populate a matrix that will record features present when each feature is gained, and vice versa
  mat.dep = matrix(0, nrow=L, ncol=L)

  # go through the edges in the solution, recording how many features are present when each feature is gained
  for(i in 1:length(E(wg))) {
    these.nodes = ends(wg, E(wg)[i])
    src = strsplit(these.nodes[1,1], "")[[1]]
    dest = strsplit(these.nodes[1,2], "")[[1]]
    changes = which(src != dest)
    presents = which(src == "1")
    absents = which(dest == "0")
    mat.dep[presents, absents] = mat.dep[presents, absents]+1
  }

  # label this matrix and summarise in a dataframe for plotting
  colnames(mat.dep) = genes
  rownames(mat.dep) = genes
  gene.stats = data.frame(gene=genes, cols=colSums(mat.dep), rows=rowSums(mat.dep))

  # produce the corresponding plots
  oDNA.g[[oDNA.expt]] = ggplot(gene.stats, aes(x=cols, y=rows, label=gene)) +
    geom_point() + geom_text_repel(max.overlaps = 30) + scale_x_sqrt() + scale_y_sqrt() +
    labs(x="Other losses preceding loss", y="Other losses following loss")

  oDNA.g2[[oDNA.expt]] = as.ggplot(pheatmap(mat.dep,
                                  treeheight_row = 0, treeheight_col = 0,
                                  fontsize_row = 6, fontsize_col = 6))
}

png("oDNA-points.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(oDNA.g[[1]] + theme_minimal(), oDNA.g[[2]] + theme_minimal(),
          labels=c("A", "B"))
dev.off()

# output these plots
png("oDNA-heatmaps.png", width=800*sf, height=1200*sf, res=72*sf)
ggarrange(oDNA.g2[[1]] + theme_void(), oDNA.g2[[2]] + theme_void(),
          labels=c("A", "B"), nrow=2, heights=c(1, 2.5))
dev.off()

png("AML-points.png", width=400*sf, height=300*sf, res=72*sf)
oDNA.g[[3]] + theme_minimal() +
  labs(x = "Other gains preceding gain", y = "Other gains following gain")
dev.off()
