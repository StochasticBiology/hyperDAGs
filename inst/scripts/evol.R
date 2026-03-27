library(phytools)
library(hyperinf)
library(ggpubr)

# from https://academic.oup.com/sysbio/article/59/6/723/1711587?guestAccessKey=#supplementary-data

tree = read.tree(system.file("extdata/squamates.tree", package="hyperdags"))
tips = read.table(system.file("extdata/squamates-treetips.txt", package="hyperdags"), sep=" ")
data = read.table(system.file("extdata/squamates.txt", package="hyperdags"))

tips$V2 = gsub("[,;]", "", tips$V2)
genera = strsplit(tips$V2, "_")
gen.df = data.frame()
for(i in 1:nrow(tips)) {
  gen.df = rbind(gen.df, data.frame(tip=as.numeric(tips$V1[i]),
                                    genus=genera[[i]][length(genera[[i]])]))
}
for(i in 1:length(tree$tip.label)) {
  ref = which(gen.df$tip == tree$tip.label[i])
  tree$tip.label[i] = gen.df$genus[ref]
}
plotTree(tree)

data$V1 = gsub(";", "", data$V1)
data$V2 = as.numeric(gsub(",", "", data$V2))
data$V3 = as.numeric(gsub(",", "", data$V3))

data$a.1 = ifelse(data$V2 < 5, 1, 0)
data$a.2 = ifelse(data$V2 < 4, 1, 0)
data$a.3 = ifelse(data$V2 < 3, 1, 0)
data$a.4 = ifelse(data$V2 < 2, 1, 0)
data$a.5 = ifelse(data$V2 < 1, 1, 0)
data$b.1 = ifelse(data$V3 < 5, 1, 0)
data$b.2 = ifelse(data$V3 < 4, 1, 0)
data$b.3 = ifelse(data$V3 < 3, 1, 0)
data$b.4 = ifelse(data$V3 < 2, 1, 0)
data$b.5 = ifelse(data$V3 < 1, 1, 0)

df = as.data.frame(data[,c(1, 4:ncol(data))])
row.names(df) = NULL
colnames(df) = c("ID", "m<5", "m<4", "m<3", "m<2", "m<1", "p<5", "p<4", "p<3", "p<2", "p<1")
plot_hyperinf_data(df, tree) 

mfit = hyperinf(df, tree, method="hyperdags")
ggarrange(plot_hyperinf_data(df, tree),
          plot_hyperinf(mfit))

mfit.prop = fit_properties(mfit)
squamate.m = as.matrix(df[,2:ncol(df)])

null.df = data.frame()
for(i in 1:100) {
  tmp = simulate_null_model(m=squamate.m, acc.rate=3, tree=tree)
  null.df = rbind(null.df, tmp$sim.df)
}
cat(sum(squamate.m), "vs", mean(null.df$n1), "\n", 
    mfit.prop$Sprime, "vs ", mean(null.df$Sprime), "+-", sd(null.df$Sprime), "\n",
    length(which(null.df$Sprime > mfit.prop$Sprime))/nrow(null.df))

# from https://onlinelibrary.wiley.com/doi/epdf/10.1002/jez.b.21212

m2 = read.table(system.file("extdata/bones.txt", package="hyperdags"), sep = " ")
d2 = apply(m2[,3:ncol(m2)], c(1,2), as.numeric)
colnames(d2) = c("Capitate", "Hamate", "Triquetral", "Lunate", "Trapezium", "Trapezoid")
names.short = c("Cap", "Ham", "Tri", "Lun", "Trz", "Trd")
fit2 = hyperinf(d2)
plot_hyperinf(fit2)

fit2.dags = hyperinf(d2, method="hyperdags")
ggarrange(plot_hyperinf_data(d2), plot_hyperinf(fit2.dags))
plot_hyperinf(fit2.dags)
bfit.prop = fit_properties(fit2.dags)

null.b.df = data.frame()
for(i in 1:100) {
  tmp = simulate_null_model(m=d2, tree=NULL)
  null.b.df = rbind(null.b.df, tmp$sim.df)
}
cat(sum(d2), "vs", mean(null.b.df$n1), "\n", 
   "S: ", bfit.prop$S, "vs ", mean(null.b.df$S), "+-", sd(null.b.df$S), "\n",
   length(which(null.b.df$S > bfit.prop$S))/nrow(null.b.df), "\n",
   "S': ", bfit.prop$Sprime, "vs ", mean(null.b.df$Sprime), "+-", sd(null.b.df$Sprime), "\n",
    length(which(null.b.df$Sprime > bfit.prop$Sprime))/nrow(null.b.df))

#######

sf = 2
png("morpho-cases.png", width=600*sf, height=500*sf, res=72*sf)
ggarrange(plot_hyperinf_data(df, tree),
          plot_hyperinf(mfit),
          plot_hyperinf_data(d2), plot_hyperinf(fit2.dags, feature.names = names.short),
          labels = c("Ai", "ii", "Bi", "ii"))
dev.off()
