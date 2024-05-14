# Kostas' implementation of the "Gutin algorithm"
# https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 (page 3)
# https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5

# there are lots of data-cleaning steps that could be optimised, avoiding loops and data conversions, but this can come later

library(igraph)
library(stringr)
library(DescTools)
library(dplyr)
library(stringdist)
library(ggraph)
library(ggplot2)

set.seed(1)

# thoughts so far:
# 1 -- current -- build arborescence for ancestors, attach descendants starting from bottom and moving up
# 2 -- build arborescence for all nodes, add edges where needed to capture interactions, prune afterwards
# 3 -- (in parallel) shift branch points as far down the tree as possible

# so far this can be "inline", "file", "TB" (1 above) or "TB2" (arborescence for all nodes)
expt = "inline"

if(expt == "file") {
  # inPut="test-tb"
  # inPut="test-mtdna-full"
  inPut="test-1"
  
  # read file with strings, add the all-zero string, and remove duplicates
  dfraw = read.table(paste(inPut,".csv",sep = ""),header = F,as.is = T,colClasses = "numeric",sep=",")
  dfraw[nrow(dfraw)+1,]=rep(0,ncol(dfraw))
  df=unique(dfraw)
  
  # store the binary name of each entry as a character
  names=vector()
  for (i in 1:nrow(df)) {
    names[i]=paste(as.character(df[i,]),collapse = '')
  }
  
} else if(expt == "inline") {
  df = matrix(c( 1,0,0,0,0,
                 0,1,0,0,0,
                 0,0,0,0,0,
                 1,0,0,1,0,
                 0,1,0,1,1,
                 1,0,1,1,0), ncol=5, byrow=TRUE)
  # store the binary name of each entry as a character
  names=vector()
  for (i in 1:nrow(df)) {
    names[i]=paste(as.character(df[i,]),collapse = '')
    
    # model sets of descendant nodes
    dfdesc = matrix(c(1,0,1,1,0,
                      0,1,1,1,1,
                      0,1,0,1,1,
                      1,1,1,1,0,
                      0,1,1,1,1,
                      1,1,1,1,1), ncol=5, byrow = TRUE)
    
    dfdesc = matrix(c(1,1,1,1,0,
                      1,1,1,1,0,
                      0,1,0,1,1,
                      1,1,1,1,0,
                      0,1,1,1,1,
                      1,1,1,1,1), ncol=5, byrow = TRUE)
    
    # get string labels
    descnames=vector()
    for (i in 1:nrow(dfdesc)) {
      descnames[i]=paste(as.character(dfdesc[i,]),collapse = '')
    }
    
  }
} else if(expt == "TB") {
  tbdf = read.table("tb_drug.txt", colClasses = "character")
  names = tbdf[,1]
  descnames = tbdf[,2]
} else if(expt == "TB2") {
  tbdf = read.table("tb_drug.txt", colClasses = "character")
  names = c(tbdf[,1],tbdf[,2])
}

# compute the distance matrix among all entries

# sort data, this way the always the minimum nearest neighbor will be chosen
len=str_length(names[1])
distMat=stringsimmatrix(names,method = "hamming")
distMat=round(len-len*distMat)
tdistMat = distMat
distMat=distMat/distMat
diag(distMat) <- len+1

skeleton= data.frame(Anc=character(),
                     Desc=character())

onecounts = sapply(names, str_count, "1")
# Define the function f
matdiff <- function(x, y) {
  return(x-y)  # Example function, replace with your own function
}

onecountsdiff = outer(onecounts, onecounts, Vectorize(matdiff)) 
absonecountsdiff = abs(onecountsdiff)

# Kostas first condition:
# if ( as.numeric(StrDist(names[i],names[j],method = "hamming"))> abs((str_count(names[i], "1")-str_count(names[j], "1"))) ) {
set.1 = which(tdistMat > absonecountsdiff, arr.ind = TRUE)

# Kostas third [our second] condition:
# (str_count(names[i], "1")-str_count(names[j], "1"))>0 [where i<j]
set.2 = which(onecountsdiff > 0, arr.ind = TRUE)
set.2 = set.2[set.2[,1] < set.2[,2],]

# Kostas second [our third] condition:
# (str_count(names[i], "1")-str_count(names[j], "1"))<0 ) [where i<j]
set.3 = which(onecountsdiff < 0, arr.ind = TRUE)
set.3 = set.3[set.3[,1] < set.3[,2],]

# condition 1 trumps conditions 2 and 3, so find indices unique to 2 and 3
unique.2 = as.matrix(anti_join(as.data.frame(set.2), as.data.frame(set.1)))
unique.3 = as.matrix(anti_join(as.data.frame(set.3), as.data.frame(set.1)))

# apply distMat effects of these conditions
distMat[set.1] = len+1
distMat[unique.2] = 2
distMat[unique.3] = 1

# add the corresponding entries to the skeleton (ji for condition 2, ij for condition 3)
skeleton = rbind(skeleton, data.frame(Anc=names[unique.2[,2]],
                                      Desc=names[unique.2[,1]]))
skeleton = rbind(skeleton, data.frame(Anc=names[unique.3[,1]],
                                      Desc=names[unique.3[,2]]))

# reduced version of Kostas' conditionals
#for (i in 1:(nrow(distMat)-1)) {
#  for (j in (i+1):nrow(distMat)) {
#    if ( tdistMat[i,j] > absonecountsdiff[i,j]) { #abs(onecounts[i]-onecounts[j]) ) {
#      distMat[i,j]=len+1
#      distMat[j,i]=len+1
#    } else if ( onecountsdiff[i,j]>0 ) {
#      distMat[i,j]=2
#      skeleton[nrow(skeleton)+1,]=c(names[j],names[i])
#    } else if ( onecountsdiff[i,j]<0 ) {
#      distMat[i,j]=1
#      skeleton[nrow(skeleton)+1,]=c(names[i],names[j])
#    }
#  }
#}

# just to make sure we don't have self-loops
skeleton=skeleton[which(skeleton$Anc != skeleton$Desc),]
skeleton=unique(skeleton)

# create the two sets for the bipartite graph B
V=unique(union(skeleton[,1],skeleton[,2]))
Vprim=unique(skeleton[,2])

# IGJ: paste is already vectorised
Vprim = paste(Vprim, "P", sep="")
# replacing
#for (i in 1:length(Vprim)) {
#  Vprim[i]=paste(Vprim[i],"P",sep = "")
#}

# IGJ: vectorising this too
B = data.frame(Anc=skeleton[,1], Desc=paste(skeleton[,2], "P", sep=""))
# replacing
#B= data.frame(Anc=character(),
#              Desc=character())
#for (i in 1:nrow(skeleton)) {
#  B[nrow(B)+1,]=c(skeleton[i,1],paste(skeleton[i,2],"P",sep = ""))
#}

graphB <- graph_from_edgelist(as.matrix(B),directed = F)

types= c(rep(0,length(V)),rep(1,length(Vprim)))
names(types) <- c(V,Vprim)

# IGJ vectorising this
edgesB = c(t(as.matrix(B)))
# replacing
#edgesB=vector()
#for (i in 1:nrow(B)) {
#  edgesB=c(edgesB,B[i,1],B[i,2])
#}

# create the bipartite graph B
bp=make_bipartite_graph(types,edgesB)
is_bipartite(graphB)

# find a maximum matching M in B.
maxbp= max_bipartite_match(bp)

# create M* 
Mstar=maxbp$matching
Mstar=Mstar[ceiling(length(Mstar)/2+1):length(Mstar)]
from=as.character(Mstar)
to=names(Mstar)
to=sub("P","",to)
from[which(is.na(from))]= strrep("0",nchar(skeleton[1,1]))

# final is the final MINLEAF arborescence
final=cbind(from,to)
graphB <- graph_from_edgelist(as.matrix(final),directed = T)
write.table(final, file="gutinsAlgorithm-out.csv", row.names=FALSE, col.names=T,sep = ",",quote = F)
tree.leaves=setdiff(final[,2],final[,1])
length(tree.leaves)

######### Iain's code from here

ggraph(graphB) + geom_edge_link() + geom_node_text(aes(label=name))

#### the following code is only applicable to the "inline" case study
if(expt != "inline") {
  stop("Use the inline case study for this bit!")
}

# model sets of descendant nodes
dfdesc = matrix(c(1,0,1,1,0,
                  0,1,1,1,1,
                  0,1,0,1,1,
                  1,1,1,1,0,
                  0,1,1,1,1,
                  1,1,1,1,1), ncol=5, byrow = TRUE)

dfdesc = matrix(c(1,1,1,1,0,
                  1,1,1,1,0,
                  0,1,0,1,1,
                  1,1,1,1,0,
                  0,1,1,1,1,
                  1,1,1,1,1), ncol=5, byrow = TRUE)

# get string labels
descnames=vector()
for (i in 1:nrow(dfdesc)) {
  descnames[i]=paste(as.character(dfdesc[i,]),collapse = '')
}

# ancestral states are those in the original data
ancs = names(V(graphB))

graphC = graphB

# loop through descendant nodes
for(this.desc in unique(descnames)) {
  print(paste0("Thinking about ", this.desc))
  
  # if this descendant node doesn't exist in the graph, add it
  if(!(this.desc %in% V(graphC)$name)) {
    print(paste0("Adding ", this.desc))
    graphC = add_vertices(graphC, n = 1)
    V(graphC)$name[length(V(graphC)$name)] = this.desc
  }
  # get list of ancestral states for this descendant
  this.ancs = c()
  anc.refs = which(descnames == this.desc)
  for(i in anc.refs) {
    this.ancs = c(this.ancs, paste0(as.character(df[i,]), collapse=""))
  }
  print(paste0("Ancestors are ", this.ancs))
  # count the 1s in each; we want to start connecting from the one with most 1s
  count1s = rep(0, length(this.ancs))
  for(i in 1:length(this.ancs)) {
    count1s[i] = str_count(this.ancs[i], "1")
  }
  count1order = order(count1s, decreasing=TRUE)
  print(count1order)
  # loop through ancestors for this descendant
  for(i in count1order) {
    this.this.ancs = this.ancs[i]
    print(paste0(" Thinking about ", this.this.ancs))
    
    # if we don't have a path from this ancestor, add an edge
    path = get.shortest.paths(graphC, from=this.this.ancs, to=this.desc)
    if(length(path$vpath[[1]]) == 0) {
      from_id = which(V(graphC)$name == this.this.ancs)
      to_id = which(V(graphC)$name == this.desc)
      print(paste0(" Adding ", this.this.ancs, "-", this.desc, " = ", from_id, " ", to_id))
      graphC = add_edges(graphC, c(from_id, to_id))
    } else {
      print(paste0(" Already have path"))
    }
  }
}

ggraph(graphC) + geom_edge_link() + geom_node_text(aes(label=name))

dataset = data.frame(ancestors = apply(df, 1, paste0, collapse=""),
                     descendants = apply(dfdesc, 1, paste0, collapse="") )

png("ex-output.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange( ggraph(graphB) + geom_edge_link() + geom_node_text(aes(label=name)) + scale_x_continuous(expand = c(0.1, 0.1)),
           ggraph(graphC) + geom_edge_link() + geom_node_text(aes(label=name)) + scale_x_continuous(expand = c(0.1, 0.1)),
           ggtexttable(dataset, theme=ttheme(base_size=12)), nrow=1, labels=c("A", "B", "C"))
dev.off()
