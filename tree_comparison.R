library(ape)

all_tax = read.tree(file = "")
add25 = read.tree(file = "")

targ = as.phylo(all_tax)
cur = as.phylo(add25)

#compare = 
all.equal.phylo(targ, cur)
dist.topo(targ, cur)
