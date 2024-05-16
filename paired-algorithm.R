# Kostas' implementation of the "Gutin algorithm"
# https://link.springer.com/chapter/10.1007/978-3-540-68880-8_23 (page 3)
# https://link.springer.com/chapter/10.1007/978-3-319-94830-0_5
# Iain's edits following this, finding a first try for a DAG that spans observation pairs

library(igraph)
library(stringr)
library(DescTools)
library(dplyr)
library(stringdist)
library(ggraph)
library(ggplot2)
library(ggpubr)

# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

# decimal to binary function
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

# decimal to binary function, returning a numerical vector
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

set.seed(1)

# thoughts so far:
# 1 -- current -- build arborescence for ancestors, attach descendants starting from bottom and moving up (TB)
# 2 -- build arborescence for all nodes, add edges where needed to capture interactions, prune afterwards (TB2, without pruning!)
# 3 -- (in parallel) shift branch points as far down the tree as possible

# so far this can be "inline", "file", "TB" (1 above) or "TB2" (part of 2 above; arborescence for all nodes)
expt = "TB2"

# excess branching score for minimisation
branching.count = function(g) {
  b = degree(g, mode="out")-1
  b[b<0] = 0
  return(sum(b))
}

if(expt == "file") {
  # inPut="test-tb"
  # inPut="test-mtdna-full"
  inPut="test-1"
  
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
  
} else if(expt == "inline") {
  L = 5
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
  L = 10
  tbdf = read.table("tb_drug.txt", colClasses = "character")
  names = tbdf[,1]
  descnames = tbdf[,2]
} else if(expt == "TB2") {
  L = 10
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

### IGJ code from here
ggraph(graphB) + geom_edge_link() + geom_node_text(aes(label=name))

### rewire
is.compat = function(s1, s2) {
  diffs = as.numeric(strsplit(s2, "")[[1]]) - as.numeric(strsplit(s1, "")[[1]])
  if(min(diffs) < 0) {
    return(Inf) 
  } else if(sum(diffs) == 0) {
    return(Inf)
  } else {
    return(sum(diffs))
  }
}
is.compat("001", "011")

compats = outer(names, names, Vectorize(is.compat)) 

new.final = final
for(i in 1:nrow(new.final)) {
  anc = new.final[i,1]
  desc.refs = which(new.final[,1]==anc) 
  outdeg = length(desc.refs)
  if(outdeg > 1) {
    print(i)
    print(paste("Thinking about ", anc, " descendants ", new.final[desc.refs,2]))
    for(j in 1:outdeg) {
      desc = new.final[desc.refs[j],2]
      print(paste("  Thinking about ", desc))
      desc.ref = which(names == desc)[1]
      desc.compats = compats[,desc.ref]
      best.new.ref = which(desc.compats == min(desc.compats) )[1]
      print(paste("  best compats ", names[best.new.ref]))
      
      new.final[desc.refs[j],1] = names[best.new.ref]
    }
  }
}


new.graphB <- graph_from_edgelist(as.matrix(new.final),directed = T)
new.graphB.layers = sapply(V(new.graphB)$name, str_count, "1")

ggraph(new.graphB, layout="sugiyama", layers = new.graphB.layers) + geom_edge_link() + geom_node_text(aes(label=name))

graphB = new.graphB

#### the following code is only applicable to the "inline" case study
if(expt == "inline") {
  
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
  ancs = names(V(graphB))
} else if(expt == "TB") {
  ancs = names(V(graphB))
} else if(expt == "TB2") {
  ancs = names = tbdf[,1]
  descnames = tbdf[,2]
} else {
  stop("This case study isn't appropriate for the WIP descendant attachment!")
}

# ancestral states are those in the original data


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
    this.ancs = c(this.ancs, names[i])
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

graphC.layers = sapply(V(graphC)$name, str_count, "1")

ggraph(graphC, layout="sugiyama", layers = graphC.layers) + geom_edge_link() + geom_node_text(aes(label=name))

ggraph(graphC) + geom_edge_link() + geom_node_text(aes(label=name))

dataset = data.frame(ancestors = names,
                     descendants = descnames )

sf = 2
png(paste0("output-rewire-", expt, ".png", collapse=""), width=1600*sf, height=600*sf, res=72*sf)
ggarrange( ggraph(graphB) + geom_edge_link() + geom_node_text(aes(label=name), angle=45, hjust=0) + ggtitle(branching.count(graphB)) + scale_x_continuous(expand = c(0.1, 0.1)),
           ggraph(graphC) + geom_edge_link() + geom_node_text(aes(label=name), angle=45, hjust=0) + ggtitle(branching.count(graphC)) + scale_x_continuous(expand = c(0.1, 0.1)),
           ggtexttable(dataset, theme=ttheme(base_size=12)), nrow=1, labels=c("A", "B", "C"), label.y=0.1)
dev.off()

# construct full hypercube for comparison
pow2 = 2**((L-1):0)
am = matrix(ncol=2)
# produce list of decimal edges
for(i in 1:(2**L-1)) {
  anc = DecToBinV(i-1, len=L)
  to.1 = which(anc == 0)
  for(j in 1:length(to.1)) {
    desc = i-1+pow2[to.1[j]]
    am = rbind(am, c(i-1, desc))
  }
}

# convert to graph with binary labels
ambin = apply(am[2:nrow(am),], c(1,2), DecToBin, len=L)
graphO <- graph_from_edgelist(as.matrix(ambin),directed = T)
graphO = graph.union(graphO, graphB, graphC)
graphO.layers = sapply(V(graphO)$name, str_count, "1")

# figure out which edges in the complete hypercube are those that we found in graph B/C
rowsB <- apply(get.edgelist(graphB), 1, paste, collapse = ",")
rowsC <- apply(get.edgelist(graphC), 1, paste, collapse = ",")
rowsO <- apply(get.edgelist(graphO), 1, paste, collapse = ",")
common_rows_B <- intersect(rowsB, rowsO)
common_rows_C <- intersect(rowsC, rowsO)
indexes_B <- which(rowsO %in% common_rows_B)
indexes_C <- which(rowsO %in% common_rows_C)

# set a variable for those edges included in graphC
E(graphO)$skeleton_B = E(graphO)$skeleton_C = 0
E(graphO)[indexes_B]$skeleton_B = 1
E(graphO)[indexes_C]$skeleton_C = 1
ggarrange(
ggraph(graphO) + 
  geom_edge_link(aes(edge_color=factor(skeleton_B), edge_alpha=factor(skeleton_B))) + 
  scale_edge_alpha_manual(values=c("0"=1/L, "1"=1)) + 
#  scale_edge_color_manual(values=c("0"="#EEEEEE", "1"="black")) + 
  theme_graph() + 
  theme(legend.position="none"),
ggraph(graphO) + 
  geom_edge_link(aes(edge_color=factor(skeleton_C), edge_alpha=factor(skeleton_C))) + 
  scale_edge_alpha_manual(values=c("0"=1/L, "1"=1)) + 
#  scale_edge_color_manual(values=c("0"="#EEEEEE", "1"="black")) + 
  theme_graph() + 
  theme(legend.position="none")
)
# not quite right yet...