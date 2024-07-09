source("paired-algorithm.R")

# runtime ~2h on modern Mac

sf = 3
# if there are multiple python3 installs on the machine, use this to set the required path
# local.python = "python3"
local.python = "/opt/homebrew/Caskroom/miniconda/base/bin/python3"
expt.index = 0
expt.out = list()

data.properties = function(fit) {
  str = paste0("L = ", str_length(fit$dataset$ancestors[1]),
               "; ntrans = ", nrow(unique(fit$dataset)),
               "; ntotal = ", nrow(fit$dataset),
               "; nuniq = ", length(unique(c(fit$dataset$ancestors, fit$dataset$descendants))),
               "; S = ", round(1-fit$best.bc/nrow(fit$dataset), digits=2),
               "; S' = ", round(1-fit$best.bc/nrow(unique(fit$dataset)), digits=2),
               "; |E| = ", length(E(fit$best.graph)),
               "; B = ", branching.count(fit$best.graph),
               "; LS = ", layer.sum(fit$best.graph)
               )
  message(str)
}

write.single.steps = function(trans, L, fname) {
  trans.set = matrix(ncol=2)
  for(i in 1:nrow(trans)) {
    src = DecToBinV(trans[i,1], L)
    dest = DecToBinV(trans[i,2], L)
    changes = which(dest-src == 1)
    curr = src
    for(j in changes) {
     trans.set = rbind(trans.set, c(BinToDec(curr), BinToDec(curr)+2**(L-j) ))
      curr[j] = 1
    }
  }
  write.table(unique(trans.set[2:nrow(trans.set),]), fname, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

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
  
  expt.out[[expt.index]] = s.dag
}

title.style = theme(plot.title = element_text(
  size = 14,     # Change font size
  family = "Arial",
  face = "plain", # Change font face to bold
  hjust = 0.5,   # Center align the title
  vjust = 0 #1.5    # Adjust vertical alignment
))

png("fig-tbsimp.png", width=1000*sf, height=300*sf, res=72*sf)
print(ggarrange(plot.stage.p(expt.out[[2]]$raw.graph) + title.style, 
                plot.stage.p(expt.out[[2]]$rewired.graph) + title.style, 
                plot.stage.p(expt.out[[2]]$best.graph) + title.style, 
                nrow = 1,
                labels=c("A", "B", "C")))
dev.off()

png("fig-tb.png", width=1200*sf, height=500*sf, res=72*sf)
w = plot.weights(expt.out[[3]]$best.graph, 
                 labels = c("INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", 
                            "MOX", "OFL", "PRO"))
print(ggarrange(plot.stage.p(expt.out[[3]]$best.graph) + title.style,
                w$thresh.plot,
                labels=c("A", "B")))
dev.off()

png("fig-cancer.png", width=600*sf, height=600*sf, res=72*sf)
w = plot.weights(expt.out[[5]]$best.graph, 
                 labels = c("FLT3", "NPM1", "WT1", "DNMT3A", "KRAS", "NRAS", "RUNX1",
                            "IDH1", "IDH2", "PTPN11", "SRSF2", "ASXL1", "BCOR", "STAG2",
                            "TP53", "U2AF1", "SF3B1", "TET2", "CSF3R", "JAK2", "GATA2",
                            "EZH2", "PPM1D", "SETBP1", "KIT", "CBL", "PHF6", "MYC", "ETV6",
                            "MPL", "SMC3"))
print(ggarrange(plot.stage.p(expt.out[[5]]$best.graph) + title.style,
                w$thresh.plot,
                labels=c("A", "B"), nrow=2))
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
