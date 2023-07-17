import ete3

treefile = '../data/pruned_tree.nwk'

t = ete3.PhyloTree(treefile)

euk_leaves = []
for leaf in t.iter_leaves():
    leaf.add_feature('clade', leaf.name.split('_', 1)[0].replace("'", ''))
    leaf.add_feature('phylum', leaf.name.split('_', 1)[1].split('-', 1)[0].split('_', 1)[0].replace("'", ''))
    if '-' in leaf.name:
        leaf.add_feature('classes', leaf.name.split('_', 1)[1].split('-')[1].split('_', 1)[0].replace("'", ''))
    else:
        leaf.add_feature('classes', leaf.name.split('_', 1)[1].split('_')[1].split('_', 1)[0].replace("'", ''))
    if leaf.clade == 'Eukaryota':
        euk_leaves.append(leaf.name)


euk_anc = t.get_common_ancestor(euk_leaves)
euk_ancestors = euk_anc.iter_ancestors ()

ancs = []
for anc in euk_ancestors:
    ancs.append(anc)

laeca_age = ancs[0].get_distance(euk_leaves[0])
euk_stem_length = euk_anc.dist
leca_age = laeca_age - euk_stem_length

print('node', 'birth', 'death', 'length', 'clades', 'phylums', 'class', sep='\t')
print('FECA-LECA', laeca_age, leca_age, euk_stem_length, 'Eukaryota', 'Eukaryota', 'Eukaryota', sep = '\t')

i = 0
for node in t.traverse():
    if not node.is_leaf():
        birth = node.get_distance(node.get_leaf_names()[0])
        length = node.dist
        death = birth - length

        clades = ';'.join(set([x.clade for x in node.get_leaves()]))
        phylum = ';'.join(set([x.phylum for x in node.get_leaves()]))
        classes = ';'.join(set([x.classes for x in node.get_leaves()]))

        print(i, birth, death, length, clades, phylum, classes, sep='\t')
    else:
        birth = node.dist
        length = node.dist

        clades = ';'.join(set([x.clade for x in node.get_leaves()]))
        phylum = ';'.join(set([x.phylum for x in node.get_leaves()]))
        classes = ';'.join(set([x.classes for x in node.get_leaves()]))

        print(i, birth, 0, length, clades, phylum, classes, sep='\t')
    
    i += 1
