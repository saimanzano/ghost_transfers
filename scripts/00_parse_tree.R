#!/usr/bin/env Rscript
#  Copyright Moisès Bernabeu, Saioa Manzano-Morales & Toni Gabaldón <saioa.manzano@bsc.es>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Loading libraries ----
library(treeio)
library(phytools)
library(stringr)

# Reading and identifying the chloroplastic and mitochondrial tips
# Clades to remove (mitochondrial and chloroplastic)
# 1240
# 898

tr <- read.nexus('../data/atps-stf2-sparse-braces.timetree.tree.taxa.combined')

a <- getDescendants(tr, 1240)
a <- a[which(a <= Ntip(tr))]

b <- getDescendants(tr, 898)
b <- b[which(b <= Ntip(tr))]

tr <- drop.tip(tr, c(a, b))

tr$tip.label <- str_replace_all(tr$tip.label, '\'', '')

write.tree(tr, file = '../data/pruned_tree.nwk')
