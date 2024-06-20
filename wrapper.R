source("paired-algorithm.R")

sf = 3

# mtDNA case study takes some minutes; ptDNA will take more
for(expt in c("inline", "TBsimp", "TB", "CGH", "cancer", "mtDNA", "ptDNA")) {
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
    ancnames = apply(tbdf, 1, paste0, collapse="")
    descnames = NULL
  } 
  
  s.dag = simplest.DAG(ancnames, descnames)

  sf = 3
  png(paste0("stage-1-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.1(s.dag))
  dev.off()
  png(paste0("stage-2-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.2(s.dag))
  dev.off()
}

V(graphD)$name2 = 1:length(V(graphD))
png("debug.png", width=2000, height=2000)
ggraph(graphD, layout="sugiyama", layers=graphD.layers) + 
  geom_edge_link(color="#CCCCCC") + 
  geom_node_text(aes(label=name2), size=label.size, angle=45, hjust=0) + #, check_overlap = TRUE) + 
  ggtitle(branching.count(graphD)) + scale_x_continuous(expand = c(0.1, 0.1)) +
  theme_graph()
dev.off()

test.name = "000000000000000000000001000000000000000000000000000000000000000000"
which(descnames == test.name)
