import ete3

treefile = '../data/pruned_tree.nwk'

t = ete3.PhyloTree(treefile)

euk_leaves = []
for leave in t.iter_leaves():
    leave.add_feature('clade', leave.name.split('_', 1)[0].replace("'", ''))
    leave.add_feature('phylum', leave.name.split('_', 1)[1].replace("'", ''))
    if leave.clade == 'Eukaryota':
        euk_leaves.append(leave.name)


euk_anc = t.get_common_ancestor(euk_leaves)
euk_ancestors = euk_anc.iter_ancestors ()

ancs = []
for anc in euk_ancestors:
    ancs.append(anc)

laeca_age = ancs[0].get_distance(euk_leaves[0])
euk_stem_length = euk_anc.dist
leca_age = laeca_age - euk_stem_length

print('node', 'birth', 'death', 'length', 'clades', sep='\t')
print('FECA-LECA', laeca_age, leca_age, euk_stem_length, 'Eukaryota', sep = '\t')

i = 0
for node in t.traverse():
    if not node.is_leaf():
        birth = node.get_distance(node.get_leaf_names()[0])
        length = node.dist
        death = birth - length

        clades = [x.clade for x in node.get_leaves()]

        print(i, birth, death, length, ';'.join(set(clades)), sep='\t')
    else:
        birth = node.dist
        length = node.dist
        print(i, birth, 0, length, ';'.join(set(clades)), sep='\t')
    
    i += 1
