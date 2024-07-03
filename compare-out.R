library(reshape2)
library(ggplot2)
library(ggpubr)
tab.norm = as.matrix(read.table("mean_tb-normal.txt"))
tab.norm.one = as.matrix(read.table("mean_tb-normal-one.txt"))
tab.cons = as.matrix(read.table("mean_tb-constrained.txt"))
tab.cons.one = as.matrix(read.table("mean_tb-constrained-one.txt"))
df.norm = melt(tab.norm)
df.norm.one = melt(tab.norm.one)
df.cons = melt(tab.cons)
df.cons.one = melt(tab.cons.one)
ggarrange(ggplot(df.norm, aes(x=Var1,y=Var2,size=value)) + geom_point(),
          ggplot(df.norm.one, aes(x=Var1,y=Var2,size=value)) + geom_point(),
          ggplot(df.cons, aes(x=Var1,y=Var2,size=value)) + geom_point(),
          ggplot(df.cons.one, aes(x=Var1,y=Var2,size=value)) + geom_point())

tab.norm = t(as.matrix(read.table("mean_tb-full-normal.txt")))
tab.cons.one = t(as.matrix(read.table("mean_tb-full-constrained-one.txt")))
df.norm = melt(tab.norm)
df.cons.one = melt(tab.cons.one)
sf = 2
png("hyperhmm-hyperdag.png", width=600*sf, height=400*sf, res=72*sf)
ggarrange(ggplot(df.norm, aes(x=Var1,y=factor(Var2),size=value)) + geom_point() + labs(x="Time", y="Feature", size="P"),
          ggplot(df.cons.one, aes(x=Var1,y=factor(Var2),size=value)) + geom_point() + labs(x="Time", y="Feature", size="P"),
          labels=c("A", "B"))
dev.off()
