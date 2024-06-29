source("paired-algorithm.R")

sf = 3
# if there are multiple python3 installs on the machine, use this to set the required path
# local.python = "python3"
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
    ancnames = tbdf[,1]
    descnames = tbdf[,2]
  } else if(expt == "TBsimp") {
    tbdf = read.table("tb_drug.txt", colClasses = "character")
    ancnames = tbdf[,1]
    descnames = tbdf[,2]
    ancnames = substr(ancnames, 1, 5)
    descnames = substr(descnames, 1, 5)
  } else if(expt == "CGH") {
    tbdf = read.table("ovarian_cgh_header.csv", header=TRUE, colClasses = "character", sep = " ")
    zero.site = paste0(rep("0", ncol(tbdf)), collapse="")
    ancnames = rep(zero.site, nrow(tbdf))
    descnames = apply(tbdf, 1, paste0, collapse="")
  } 
  
  s.dag = simplest.DAG(ancnames, descnames)
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

png("fig-tb.png", width=600*sf, height=500*sf, res=72*sf)
print(plot.stage.p(expt.out[[3]]$best.graph) + title.style)
dev.off()

png("fig-cancer.png", width=600*sf, height=300*sf, res=72*sf)
print(plot.stage.p(expt.out[[5]]$best.graph) + title.style)
dev.off()

png("fig-odna.png", width=1000*sf, height=1000*sf, res=72*sf)
print(ggarrange(plot.stage.p(expt.out[[6]]$best.graph) + title.style, 
                plot.stage.p(expt.out[[7]]$best.graph) + title.style, 
                nrow = 2,
                labels=c("A", "B")))
dev.off()
