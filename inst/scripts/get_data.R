require(ggpubr)
require(stringr)

get_mtDNA_data = function() {
  datafile <- system.file("extdata", "mt-trans-manual.csv", package = "hyperdags")
  dfraw <- read.csv(datafile)
  dfraw[,3:ncol(dfraw)] = 1-dfraw[,3:ncol(dfraw)]
  L = (ncol(dfraw)-2)/2
  ancnames = apply(dfraw[,3:(2+L)], 1, paste0, collapse="")
  descnames = apply(dfraw[,(2+L+1):(2+2*L)], 1, paste0, collapse="")
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

get_ptDNA_data = function() {
  datafile <- system.file("extdata", "pt-trans-manual.csv", package = "hyperdags")
  dfraw = read.csv(datafile)
  dfraw[,3:ncol(dfraw)] = 1-dfraw[,3:ncol(dfraw)]
  L = (ncol(dfraw)-2)/2
  ancnames = apply(dfraw[,3:(2+L)], 1, paste0, collapse="")
  descnames = apply(dfraw[,(2+L+1):(2+2*L)], 1, paste0, collapse="")
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

get_cancer_data = function() {
  datafile1 <- system.file("extdata", "cancer-srcs.csv", package = "hyperdags")
  datafile2 <- system.file("extdata", "cancer-dests.csv", package = "hyperdags")
  df1 = read.csv(datafile1, header=TRUE)
  df2 = read.csv(datafile2, header=TRUE)
  rownames(df1) = NULL
  rownames(df2) = NULL
  ancnames = apply(df1, 1, paste0, collapse="")
  descnames = apply(df2, 1, paste0, collapse="")
  L = ncol(df1)
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

get_inline_data = function() {
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
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

get_TB_data = function(simple = FALSE) {
  datafile <- system.file("extdata", "tb_drug.txt", package = "hyperdags")
  tbdf = read.table(datafile, colClasses = "character")
  L = str_length(tbdf[1,1])
  ancnames = tbdf[,1]
  descnames = tbdf[,2]
  if(simple == TRUE) {
    ancnames = substr(ancnames, 1, 5)
    descnames = substr(descnames, 1, 5)
  }
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

get_CGH_data = function() {
  datafile <- system.file("extdata", "ovarian_cgh_header.csv", package = "hyperdags")
  tbdf = read.table(datafile, header=TRUE, colClasses = "character", sep = " ")
  zero.site = paste0(rep("0", ncol(tbdf)), collapse="")
  L = nrow(tbdf)
  ancnames = rep(zero.site, nrow(tbdf))
  descnames = apply(tbdf, 1, paste0, collapse="")
  return(list(L=L, ancnames=ancnames, descnames=descnames))
}

