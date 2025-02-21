source("paired-algorithm.R")

library(reshape2)
library(ggrepel)
library(ggplotify)

#### basic example for demonstration

# simple illustrative dataset
ancnames  = c("10000","01000","00000","10010","01011","10110")
descnames = c("11110","11110","01011","11110","01111","11111")
# get STSD
s.dag = simplest.DAG(ancnames, descnames)
# plot and summarise: A is SSA, B is STSD, numbers are branching counts
plot.stage.2(s.dag)
fit.properties(s.dag)

#### real experiments
# runtime ~2h on modern Mac

sf = 3
# if there are multiple python3 installs on the machine, use this to set the required path
# local.python = "python3"
run.python = FALSE    # don't run the Python version by default; change to do so
local.python = "/opt/homebrew/Caskroom/miniconda/base/bin/python3"
expt.index = 0
expt.out = list()

# mtDNA case study takes some minutes; ptDNA will take more
for(expt in c("inline", "TBsimp", "TB", "CGH", "cancer", "mtDNA", "ptDNA")) {
  expt.index = expt.index + 1
  if(expt == "mtDNA") {
    dfraw = read.csv("mt-trans-manual.csv")
    dfraw[,3:ncol(dfraw)] = 1-dfraw[,3:ncol(dfraw)]
    write.table(dfraw, "mt-trans-manual-gains.csv", sep=",", row.names=FALSE, col.names=TRUE)
    L = (ncol(dfraw)-2)/2
    ancnames = apply(dfraw[,3:(2+L)], 1, paste0, collapse="")
    descnames = apply(dfraw[,(2+L+1):(2+2*L)], 1, paste0, collapse="")
  } else if(expt == "ptDNA") {
    dfraw = read.csv("pt-trans-manual.csv")
    dfraw[,3:ncol(dfraw)] = 1-dfraw[,3:ncol(dfraw)]
    write.table(dfraw, "pt-trans-manual-gains.csv", sep=",", row.names=FALSE, col.names=TRUE)
    L = (ncol(dfraw)-2)/2
    ancnames = apply(dfraw[,3:(2+L)], 1, paste0, collapse="")
    descnames = apply(dfraw[,(2+L+1):(2+2*L)], 1, paste0, collapse="")
  } else if(expt == "cancer") {
    df1 = read.csv("cancer-srcs.csv", header=TRUE)
    df2 = read.csv("cancer-dests.csv", header=TRUE)
    rownames(df1) = NULL
    rownames(df2) = NULL
    ancnames = apply(df1, 1, paste0, collapse="")
    descnames = apply(df2, 1, paste0, collapse="")
  } else if(expt == "inline") {
    L = 5
    df = matrix(c( 1,0,0,0,0,
                   0,1,0,0,0,
                   0,0,0,0,0,
                   1,0,0,1,0,
                   0,1,0,1,1,
                   1,0,1,1,0), ncol=5, byrow=TRUE)
    # store the binary name of each entry as a character
    
    ancnames=apply(df, 1, paste0, collapse = '')
    
    # model sets of descendant nodes
    
    dfdesc = matrix(c(1,1,1,1,0,
                      1,1,1,1,0,
                      0,1,0,1,1,
                      1,1,1,1,0,
                      0,1,1,1,1,
                      1,1,1,1,1), ncol=5, byrow = TRUE)
    
    # get string labels
    descnames=apply(dfdesc, 1, paste0, collapse = '')
  } else if(expt == "TB") {
    tbdf = read.table("tb_drug.txt", colClasses = "character")
    L = str_length(tbdf[1,1])
    ancnames = tbdf[,1]
    descnames = tbdf[,2]
  } else if(expt == "TBsimp") {
    tbdf = read.table("tb_drug.txt", colClasses = "character")
    ancnames = tbdf[,1]
    descnames = tbdf[,2]
    ancnames = substr(ancnames, 1, 5)
    descnames = substr(descnames, 1, 5)
    L = 5
  } else if(expt == "CGH") {
    tbdf = read.table("ovarian_cgh_header.csv", header=TRUE, colClasses = "character", sep = " ")
    zero.site = paste0(rep("0", ncol(tbdf)), collapse="")
    L = nrow(tbdf)
    ancnames = rep(zero.site, nrow(tbdf))
    descnames = apply(tbdf, 1, paste0, collapse="")
  } 
  
  dset.trans = unique(data.frame(Ancestor=ancnames, Descendant=descnames))
  states = unique(c(ancnames, descnames))
  dset.states = data.frame(States=states)
  
  if(expt == "TBsimp") {
    df = 3
    tb.simp.data =  ggarrange(ggtexttable(dset.states, rows=NULL),
                              ggarrange(ggtexttable(dset.trans[1:(nrow(dset.trans)/3),], rows=NULL), 
                                        ggtexttable(dset.trans[(nrow(dset.trans)/3+1):(2*nrow(dset.trans)/3),], rows=NULL),
                                        ggtexttable(dset.trans[(2*nrow(dset.trans)/3+1):(nrow(dset.trans)),], rows=NULL),
                                        nrow=1),
                              nrow = 1, labels=c("A", "B"), widths=c(0.3,1))
    png(paste0("stage-0-", expt, ".png", collapse=""), width=600*sf, height=480*sf, res=72*sf)
    print(tb.simp.data)
    dev.off()
  }
  
  s.dag = simplest.DAG(ancnames, descnames)
  write.table(apply(as_edgelist(s.dag$best.graph), c(1,2), BinToDec), paste0(expt, "-table.txt", collapse=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.single.steps(apply(as_edgelist(s.dag$best.graph), c(1,2), BinToDec), L, paste0(expt, "-table-ss.txt", collapse=""))
  if(length(descnames) != 0) {
    spanned = transitions.spanned(s.dag$best.graph, ancnames, descnames)
    if(spanned == TRUE) {
      message("-- transitions verified") 
    } else {
      message("-- transitions *NOT* verified")
    }
  }
  
  sf = 3
  png(paste0("stage-1-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.1(s.dag))
  dev.off()
  png(paste0("stage-2-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.2(s.dag))
  dev.off()
  
  if(run.python) {
  message("Python code:")
  write.table(data.frame(anc=ancnames,desc=descnames), paste0(expt, "-trans-basic.txt", collapse=""), quote = FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
  write.table(data.frame(anc=ancnames,desc=descnames), "_in", quote = FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
  system(paste0(local.python, " disjointp.py", collapse=""))
  p.df = as.matrix(read.table("_out", header=FALSE, sep = " ", colClasses = rep("character",2)))
  p.g = graph_from_edgelist(p.df)
  spanned = transitions.spanned(s.dag$best.graph, ancnames, descnames)
  if(spanned == TRUE) {
    message("-- transitions verified") 
  } else {
    message("-- transitions *NOT* verified")
  }
  png(paste0("stage-p-", expt, ".png", collapse=""), width=1100*sf, height=300*sf, res=72*sf)
  print(ggarrange(plot.stage.2(s.dag), plot.stage.p(p.g), 
                  nrow=1, widths=c(2,1), labels=c("", "C")))
  dev.off()
  }
  
  expt.out[[expt.index]] = s.dag
}

scale.x = scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))
scale.y = scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
title.style = theme(plot.title = element_text(
  size = 14,     # Change font size
  family = "Arial",
  face = "plain", # Change font face to bold
  hjust = 0.5,   # Center align the title
  vjust = 0 #1.5    # Adjust vertical alignment
)) 

png("fig-tbsimp.png", width=850*sf, height=700*sf, res=72*sf)
print(ggarrange(
  tb.simp.data, 
  ggarrange(plot.stage.p(expt.out[[2]]$raw.graph) + scale.y + title.style, 
                plot.stage.p(expt.out[[2]]$rewired.graph) + scale.y + title.style, 
                plot.stage.p(expt.out[[2]]$best.graph) + scale.y + title.style, 
                nrow = 1,
                labels=c("C", "D", "E")),
  nrow=2, heights=c(1,0.66))
)
dev.off()

png("fig-tb.png", width=1200*sf, height=500*sf, res=72*sf)
w = plot.weights(expt.out[[3]]$best.graph, 
                 labels = c("INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", 
                            "MOX", "OFL", "PRO"))
print(ggarrange(plot.stage.p(expt.out[[3]]$best.graph) + scale.y + title.style,
                w$thresh.plot,
                labels=c("C", "D")))
dev.off()

png("fig-cancer.png", width=1000*sf, height=500*sf, res=72*sf)
w = plot.weights(expt.out[[5]]$best.graph, 
                 labels = c("FLT3", "NPM1", "WT1", "DNMT3A", "KRAS", "NRAS", "RUNX1",
                            "IDH1", "IDH2", "PTPN11", "SRSF2", "ASXL1", "BCOR", "STAG2",
                            "TP53", "U2AF1", "SF3B1", "TET2", "CSF3R", "JAK2", "GATA2",
                            "EZH2", "PPM1D", "SETBP1", "KIT", "CBL", "PHF6", "MYC", "ETV6",
                            "MPL", "SMC3"))
print(ggarrange(                 plot.stage.p(expt.out[[5]]$best.graph) + coord_flip() + scale.y + scale_y_reverse() + title.style+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 w$thresh.plot + scale_y_reverse() + coord_flip()+theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 w$thresh.plot+coord_flip() + scale_y_reverse() + theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1)),
                                 nrow=1, labels = c("A", "B", "C"))
)
dev.off()


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
#plot.stage.p(expt.out[[6]]$best.graph, v.labels=mt.labs)

pt.orgs = c("arabidopsis thaliana", "physcomitrium patens", "hydnora visseri", "parasitaxus usta",
            "polysiphonia elongata", "gracilaria changii", "cladosiphon okamuranus")
pt.raw = read.csv("pt-barcodes-manual.csv")
L = ncol(pt.raw)-1
pt.names = apply(1-pt.raw[,2:(L+1)], 1, paste0, collapse="")
pt.labs = data.frame()
for(i in 1:length(pt.orgs)) {
  pt.labs = rbind(pt.labs, data.frame(Species = pt.orgs[i], label=pt.names[which(pt.raw$Species == pt.orgs[i])]))
}
#plot.stage.p(expt.out[[7]]$best.graph, v.labels=pt.labs)

png("fig-odna.png", width=1000*sf, height=1000*sf, res=72*sf)
print(ggarrange(plot.stage.p(expt.out[[6]]$best.graph, v.labels=mt.labs) + title.style, 
                plot.stage.p(expt.out[[7]]$best.graph, v.labels=pt.labs) + title.style, 
                nrow = 2,
                labels=c("A", "B")))
dev.off()

oDNA.types = c("MT", "PT")
oDNA.g = list()
oDNA.g2 = list()

for(oDNA.expt in 1:2) {
  oDNA = oDNA.types[oDNA.expt]
  if(oDNA == "MT") {
    wg = expt.out[[6]]$best.graph
  } else{
    wg = expt.out[[7]]$best.graph
  }
  these.nodes = ends(wg, E(wg)[1])
  src = strsplit(these.nodes[1,1], "")[[1]]
  L = length(src)
  if(oDNA == "MT") {
    genes = colnames(mt.raw)[2:(L+1)]
  } else {
    genes = colnames(pt.raw)[2:(L+1)]
  }
  
  mat.dep = matrix(0, nrow=L, ncol=L)
  
  for(i in 1:length(E(wg))) {
    these.nodes = ends(wg, E(wg)[i])
    src = strsplit(these.nodes[1,1], "")[[1]]
    dest = strsplit(these.nodes[1,2], "")[[1]]
    changes = which(src != dest)
    presents = which(src == "1")
    absents = which(dest == "0")
    mat.dep[presents, absents] = mat.dep[presents, absents]+1
  }
  
  colnames(mat.dep) = genes
  rownames(mat.dep) = genes
  gene.stats = data.frame(gene=genes, cols=colSums(mat.dep), rows=rowSums(mat.dep))
  
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

png("oDNA-heatmaps.png", width=800*sf, height=1200*sf, res=72*sf)
ggarrange(oDNA.g2[[1]] + theme_void(), oDNA.g2[[2]] + theme_void(), 
          labels=c("A", "B"), nrow=2, heights=c(1, 2.5))
dev.off()
