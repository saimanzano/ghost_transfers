library(treeio)
library(ggtree)
library(stringr)
library(ggplot2)

tr <- read.nexus('../data/atps-stf2-sparse-braces.timetree.tree.taxa.combined')

ggtree(tr) +
  geom_tippoint(aes(colour = str_split(label, '_', simplify = TRUE)[, 1])) +
  geom_nodelab(aes(x = branch, label = node))

# Clades to remove (mitochondrial and chloroplastic)
# 1240
# 898

a <- phytools::getDescendants(tr, 1240)
a <- a[which(a <= Ntip(tr))]

b <- phytools::getDescendants(tr, 898)
b <- b[which(b <= Ntip(tr))]

tr <- drop.tip(tr, c(a, b))

phytools::tip

ggtree(tr) +
  geom_tippoint(aes(colour = str_split(label, '_', simplify = TRUE)[, 1])) +
  geom_nodelab(aes(x = branch, label = node))

# write.tree(tr, file = '../data/pruned_tree.nwk')

tr <- read.tree('../data/pruned_tree.nwk')

ggtree(tr) +
  geom_tippoint(aes(colour = str_split(label, '_', simplify = TRUE)[, 1])) +
  labs(colour = '')

