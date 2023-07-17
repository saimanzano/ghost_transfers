library(treeio)
library(ggtree)
library(stringr)
library(ggplot2)

tr <- read.nexus('../data/atps-stf2-sparse-braces.timetree.tree.taxa.combined')
source('../../../../software/developing/phygeno/gene_tree_plotting/src/plot_tree.R')

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

tr$tip.label <- str_replace_all(tr$tip.label, '\'', '')

ggtree(tr) +
  geom_tippoint(aes(colour = str_split(label, '_', simplify = TRUE)[, 2]), show.legend = FALSE) +
  # geom_nodelab(aes(x = branch, label = node))
  geom_nodelab(aes(x = branch, label = node), geom = 'label', alpha = 0.6)

# write.tree(tr, file = '../data/pruned_tree.nwk')

tr <- read.tree('../data/pruned_tree.nwk')
tr$tip.label <- str_replace_all(tr$tip.label, '\'', '')

p <- ggtree(tr) +
  geom_tippoint(aes(colour = str_split(label, '_', simplify = TRUE)[, 1])) +
  labs(colour = '') +
  geom_tiplab(aes(label = str_split(str_split(label, '_', simplify = TRUE)[, 2], '-', simplify = TRUE)[, 1]), size = 50) +
  geom_nodelab() +
  geom_nodelab(aes(x = branch, label = node), geom = 'label', alpha = 0.6, size = 150)
# p

pdf('../outputs/tree_big.pdf', width = 982, height = 1712 * 2)
p
dev.off()

clades <- c('Eukaryotes' = 1148,
            'Asgard' = 1147,
            'Asgard' = 1235,
            'Asgard' = 1238,
            'TACK' = 1251,
            'Euryarchaeota' = 1332,
            'DPANN' = 1476,
            'DPANN' = 1483,
            'DPANN' = 1557,
            'Alphaproteobacteria' = 1105,
            'Gammaproteobacteria' = 1080,
            'Betaproteobacteria' = 1098,
            'Epsilonproteobacteria' = 1074,
            'CPR like' = 1039,
            'Planctomycetes' = 1033,
            'Verrucomicrobia' = 1027,
            'Deltaproteobacteria' = 979,
            'Acidobacteria' = 1007,
            'Elusimicrobia' = 1013,
            'Bacteroidetes' = 941,
            'Chlorobi' = 939,
            'Cyanobacteria' = 821,
            'Chloroflexi' = 837,
            'Actinobacteria' = 849,
            'Firmiputes' = 879,
            'Fusobacteria' = 812,
            'Synergistetes' = 809,
            'Aquificae' = 807,
            'Thermotogae' = 800,
            'Deinococcus' = 803,
            'Other bacteria' = 797,
            'Other bacteria' = 792,
            'Candidatus' = 968,
            'Spirochaetes' = 1020,
            'Nitrospirae' = 1004)

p <- ggcollapsed(tr, clades, mode = 'max')

p + geom_tiplab() +
  xlim(0, 6000) +
  geom_nodelab(aes(label = node))
