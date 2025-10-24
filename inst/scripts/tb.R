library(readxl)       # to read Excel datafile
library(dplyr)        # for data wrangling
library(phytools)     # to read external tree data
library(ggplot2)      # }
library(ggpubr)       # } for visualisation
library(ggraph)       # }
library(hypertrapsct)

# get profiles of drug resistance/susceptibility
# Supplementary Table 4 of Casali et al., https://www.nature.com/articles/ng.2878.s3
system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls")
o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")

# get phylogeny linking isolates
# Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt")
tree = read.tree("41588_2014_BFng2878_MOESM34_ESM.txt")

##########
##### Curate data

# extract isolate ID and resistance profiles to our ten drugs (discarding mutation info)
# remove any incomplete profiles and recast as binary strings
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL", "PRO")
df = o.df[,which(colnames(o.df) %in% col.interest)]
missing.rows = unique(which(df == ".", arr.ind=TRUE)[,1])
final.df = df[-missing.rows,]
final.df = final.df %>%
  mutate(across(-Isolate, ~ ifelse(. == "R", 1, 0)))
final.df = as.data.frame(final.df)

##########
##### Evolutionary accumulation modelling

# reconstruct ancestral states using parsimony picture and extract before-after transitions
src.data = curate.tree(tree, final.df)

ancnames.tb = apply(src.data$srcs, 1, paste0, collapse="")
descnames.tb = apply(src.data$dests, 1, paste0, collapse="")
tb.soln.real = simplest_DAG(ancnames.tb, descnames.tb)

####### heterogeneous null model

tb.char.probs = colSums(final.df[,2:ncol(final.df)])/nrow(final.df)

tb.n1s = c()
tb.r.sim.df = data.frame()
for(i in 1:10) {
  tb.r.set = simulate_accumulation(0, 10, use.tree=src.data$tree,
                                accumulation.rate = 110,
                                dynamics="heterogeneous",
                                char.probs = tb.char.probs)
  tb.n1s = c(tb.n1s, sum(unlist(tb.r.set$x)))
  tb.soln = simplest_DAG(tb.r.set[["ancnames"]], tb.r.set[["descnames"]])
  tb.r.sim.df = rbind(tb.r.sim.df, fit_properties(tb.soln))
}

sum(final.df[,2:ncol(final.df)])
tb.n1s

ggarrange(plotHypercube.curated.tree(src.data),
          plot_tree_data(tb.r.set$my.tree, tb.r.set$x),
          plot_stage_gen(tb.soln.real$best.graph),
          plot_stage_gen(tb.soln$best.graph))

fit_properties(tb.soln)
tb.r.sim.df

####### homogeneous null model

tb.0.n1s = c()
tb.0.r.sim.df = data.frame()
for(i in 1:10) {
  tb.0.r.set = simulate_accumulation(0, 10, use.tree=src.data$tree,
                                   accumulation.rate = 1300,
                                   dynamics="poisson")
  tb.0.n1s = c(tb.0.n1s, sum(unlist(tb.0.r.set$x)))
  tb.0.soln = simplest_DAG(tb.0.r.set[["ancnames"]], tb.0.r.set[["descnames"]])
  tb.0.r.sim.df = rbind(tb.0.r.sim.df, fit_properties(tb.0.soln))
}

sum(final.df[,2:ncol(final.df)])
tb.0.n1s

ggarrange(plotHypercube.curated.tree(src.data),
          plot_tree_data(tb.0.r.set$my.tree, tb.0.r.set$x),
          plot_stage_gen(tb.soln.real$best.graph),
          plot_stage_gen(tb.0.soln$best.graph))

fit_properties(tb.soln.real)
tb.0.r.sim.df

cat("Solution S' = ", fit_properties(tb.soln.real)$Sprime, "\n",
    " compared to ", round(mean(tb.0.r.sim.df$Sprime), digits=2),
    " +- ", round(sd(tb.0.r.sim.df$Sprime), digits=2), " homogeneous\n",
    " and ", round(mean(tb.r.sim.df$Sprime), digits=2),
    " +- ", round(sd(tb.r.sim.df$Sprime), digits=2), " heterogeneous\n")
