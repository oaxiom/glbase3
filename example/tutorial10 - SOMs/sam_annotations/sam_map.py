import itertools, math
from collections import defaultdict, Counter
from glbase3 import glload
import glbase3

from extended import derivation_type, gene_layer_name
from extended import rsem_tag_counts, rsem_tag_counts_average, total_tag_counts, bowtie2_tag_counts, bowtie2_tag_counts_average, bowtie2_tag_counts_percent_map
from extended import num_replicates, for_publication
from extended import sample_description, sample_to_gse, gse_to_reference, pmid_to_gse, cell_source
from extended import rnaseq_library_parameters, bad_samples

low_expressed_threshold = 10**1.6 # see te_counts/2.gc_norm.R
low_expressed_threshold_log2 = math.log(low_expressed_threshold, 2)
low_expressed_threshold_log10 = 1.6

#print low_expressed_threshold

desc_to_name = dict((v,k) for k, v in sample_description.items()) # get back the id name
# I artificially append this one, as filesystems may alter the name if the desc is used to save:
desc_to_name['Naive CD4+ T cells, anti-CD3-28 treated'] = 'naive cd4t 24h'

# remapped for description name
gene_layer = {}
for k in gene_layer_name:
    gene_layer[k] = []
    for n in gene_layer_name[k]:
        gene_layer[k].append(sample_description[n])


colour_guide = {"Embryonic": "orange",
    "Germ cells": "yellow",
    "Mesoderm": "green",
    "Blood mesoderm": "red",
    "Endoderm": "indigo",
    "Neurectoderm": "blue",
    "Surface ectoderm": "cyan",
    "Neural crest": "magenta"
    }

sample_colours = {} # coloured by germ layer
for k in sorted({x for v in gene_layer.values() for x in v}):
    sample_colours[k] = "grey" # no layer
    for gl in gene_layer:
        if k in gene_layer[gl]:
            sample_colours[k] = colour_guide[gl]

def get_cols(cond_names):
    cols = []
    for c in cond_names:
        if c in sample_colours:
            cols.append(sample_colours[c])
        else:
            print('{} not found in colors'.format(c))
            cols.append('grey')
    return cols

all_gses = [sample_to_gse[k] for k in sample_to_gse]
all_gses = set([x for sublist in all_gses for x in sublist])

def get_details(sample_name):
    # get the germ lineage:
    germ = "?"
    for k in gene_layer:
        if sample_name in gene_layer_name[k]:
            germ = k
            break

    # These key names are also the labels for Table S1:
    ret = {"name": sample_name,
        #"GEO Accession": ", ".join(sample_to_gse[sample_name]),
        #"PMID": ", ".join([pmid_to_gse[gse] for gse in sample_to_gse[sample_name]]),
        #"Reference(s)": "; ".join([gse_to_reference[gse] for gse in sample_to_gse[sample_name]]),
        "Germ layer": germ,
        "Description": sample_description[sample_name],
        'Cell type/tissue source': derivation_type[sample_name],
        #'Total Number of tags for each replicate': ', '.join([str(int(i)) for i in total_tag_counts[sample_name]]),
        #'Number of mapped tags for each replicate': ', '.join([str(int(i)) for i in bowtie2_tag_counts[sample_name]]),
        #'Percent mapped for each replicate':  ', '.join([str(int(i)) for i in bowtie2_tag_counts_percent_map[sample_name]]),
        #'Average number of mapped tags': int(bowtie2_tag_counts_average[sample_name]),
        #'Total number of replicates': num_replicates[sample_name],
        #'Number of studies': len(sample_to_gse[sample_name]),
        #'Pooled/Single-cell': cell_source[sample_name],
        #'Paired-end/Single-end': ", ".join([i[0] for i in rnaseq_library_parameters[sample_name]]),
        #'Read-lengths(s)': ", ".join([str(i[1]) for i in rnaseq_library_parameters[sample_name]]),
        #'Methods(s)': ", ".join([i[2] for i in rnaseq_library_parameters[sample_name]]),
        #'Sequence molecules(s)': ", ".join([i[3] for i in rnaseq_library_parameters[sample_name]]),
        #'Sequencing machine(s)': ", ".join([i[4] for i in rnaseq_library_parameters[sample_name]]),
        }

    return(ret)

def get_citation_list(list_of_gses):
    """
    Get a citation list e.g.

    GSE0000000 (Smith et al., 2012; PMID: 000000000), ...
    """
    gses = list(set(sum([sample_to_gse[k] for k in list_of_gses], [])))
    gses.sort()
    r = ", ".join(["%s (%s, PMID:%s)" % (gse, gse_to_reference[gse], pmid_to_gse[gse]) for gse in gses])
    return(r)

def remap_expn_sample_names(expn):
    """
    Super-specific method to remap the sample names.
    """
    names = expn.getConditionNames()
    new_names = []
    for n in names:
        if n in sample_description:
            new_names.append(sample_description[n])
        else:
            new_names.append('?%s' % n)
    expn.setConditionNames(new_names)

if __name__ == "__main__":
    # This should later be changed to mm10v79:
    expn = glload("../te_counts/genes_cpm_expression.glb") # Using the published one for now

    print(expn.getConditionNames())

    print("Testing coherence")
    print("\nSample names:")
    # Test that all desc names are actually in the expn and vice versa:
    expn_names = expn.getConditionNames()
    desc_names = list(sample_description.keys())
    ss1 = set(expn_names) - set(desc_names)
    ss2 = set(desc_names) - set(expn_names)
    if ss1:
        print(">>>Missing in sample_description:\n", '\n'.join(sorted(list(ss1))))
    print()
    if ss2:
        print(">>>Missing in exp:\n", '\n'.join(sorted(list(ss2))))
    print("\n>>>Checking for duplicate formal names in sample_description:")
    print(" failed\n".join([x for x, y in list(Counter(list(sample_description.values())).items()) if y > 1]))

    print("\n>>>Sample GEO:")
    for c in expn.getConditionNames():
        if c not in sample_to_gse:
            print("Sample '%s' has no GSE entry" % c)

    print("\n>>>Gene layer:")
    all_assigned_to_gene_layer = sorted({x for v in gene_layer.values() for x in v})
    remap_expn_sample_names(expn)
    for c in expn.getConditionNames():
        if True not in [c == i for i in all_assigned_to_gene_layer]:
            print("ERROR: Sample '%s' has not been assigned to a germ layer" % c)

    print("\n>>>Check each entry in the germ layers is actually in the expn data:")
    all_conds = expn.getConditionNames()
    for lin in gene_layer:
        for sam in gene_layer[lin]:
            if sam not in all_conds:
                print("ERROR: '%s' entry in germ layer, no actual sample in data" % sam)

    print("\n>>>Duplicates in germ layer:")
    for l in list(gene_layer.keys()):
        cc = [x for x, y in list(Counter(gene_layer[l]).items()) if y > 1]
        if cc:
            print(" failed\n".join(cc))

    print("\n>>>Number of samples in layer")
    for k in gene_layer:
        print(k, len(gene_layer[k]))

    print('\n', set(sum(list(sample_to_gse.values()), [])))
    print('>>>Number of unique GSE/study numbers:', len(set(sum(list(sample_to_gse.values()), []))))
    print('Make sure you check there is not "" entry in the GSEs')

    print("\n>>>GSE -> PMID")
    for gse in all_gses:
        res = [True, True]
        res_text = []
        if gse not in gse_to_reference:
            res[0] = False
            #res_text.append("no reference")

        if gse not in pmid_to_gse:
            res[1] = False
            res_text.append("no pmid")

        if False in res:
            print("'%s': '', " % (gse))

    import csv
    oh = open("sample_map.tsv", "w")
    writer = csv.writer(oh, dialect=csv.excel_tab)
    keys = list(sample_description.keys())
    keys.sort()
    key_order = ["Description",
        "Germ layer",
        'Cell type/tissue source',
        #'Pooled/Single-cell',
        #'Total Number of tags for each replicate',
        #'Number of mapped tags for each replicate',
        #'Total number of replicates',
        #'Average number of mapped tags',
        #'Number of studies',
        #'Sequencing machine(s)',
        #'Read-lengths(s)',
        #'Paired-end/Single-end',
        #'Sequence molecules(s)',
        #'Methods(s)',
        #"GEO Accession",
        #"PMID",
        #"Reference(s)",
        ]
    writer.writerow(key_order)
    newl = []
    for k in keys:
        d = get_details(k)
        writer.writerow([d[k] for k in key_order])
        newl.append(d)
    oh.close()

    print()
    print("\n>>>Citation:")
    gses = list(gse_to_reference.keys())
    gses.sort()
    r = ", ".join(["%s (%s:%s)" % (gse, gse_to_reference[gse], pmid_to_gse[gse]) for gse in gses])
    print(r)
    print()

    gl = glbase3.genelist()
    gl.load_list(newl)
    gl.save("sample_map.glb")

    print()
    print("Number of studies in domain:")
    res = {}
    for k in gene_layer_name:
        res[k] = []
        for ct in gene_layer_name[k]:
            res[k] += sample_to_gse[ct]

        res[k] = set(res[k])

    for k in sorted(res):
        print(k, len(res[k]))
