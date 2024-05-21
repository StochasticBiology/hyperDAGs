source("paired-algorithm-2.R")

sf = 2

for(expt in c("inline", "TBsimp", "TB", "CGH", "cancer")) {
  
  if(expt == "mtDNA") {
    # inPut="test-tb"
    inPut="test-mtdna-full"
    # inPut="test-1"
    
    # read file with strings, add the all-zero string, and remove duplicates
    dfraw = read.table(paste(inPut,".csv",sep = ""),header = F,as.is = T,colClasses = "numeric",sep=",")
    dfraw[nrow(dfraw)+1,]=rep(0,ncol(dfraw))
    df=unique(dfraw)
    L = ncol(dfraw)
    # store the binary name of each entry as a character
    names=vector()
    for (i in 1:nrow(df)) {
      names[i]=paste(as.character(df[i,]),collapse = '')
    }
    
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
  png(paste0("stage-1-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.1(s.dag))
  dev.off()
  png(paste0("stage-2-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(plot.stage.2(s.dag))
  dev.off()
}


