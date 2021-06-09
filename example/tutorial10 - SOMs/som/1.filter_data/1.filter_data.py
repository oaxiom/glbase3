"""

Basically a quick and dirty replacement for the old gene expression
classification system. The aim is to remove non-informative genes for things like PCA and
SOMs to get cleaner plots.

"""

import math, os, sys
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append("../../sam_annotations/")
import sam_map

expn = glload("../../te_counts/genes_ntc_expression.glb")
conds = expn.getConditionNames()
sam_map.remap_expn_sample_names(expn)
#expn.log(2, 0.1)

labs = ['Pou5f1', 'Sox2', 'Nanog', 'Cdx2', 'Esrrb', 'Dnmt3l', 'Utf1', # embryonic
    'Eomes', 'Hand1', 'Tfap2c', 'Tfap2a', # Other embryonic
    'Sox7', 'Hnf4a', 'Gata4', 'Gata6', 'Sox17', # endoderm
    'Foxa2', 'Myog',
    'Trp63', # Surface ectoderm
    'Myod1', # Muscle mesoderm
    'Nr2f1', 'Nr2f2', 'Msx1', 'Msx2', 'Tfap2a', # Neural crest
    'Tal1', 'Cebpe',# Blood mesoderm
    'Sox10','Olig2', 'Zfp533', 'Neurod1', 'Neurod6', 'Olig1', 'Myt1l', 'Pou3f1', 'Ascl1', # neurectoderm masters
    'Lin28a', 'Dmrtb1', # Germ cell master regulators
    #'Ep300', 'Stat3', 'Stat1', 'Junb', 'Jun', 'Klf4', 'Atf2', 'Mef2d', 'Max', 'Smad1', 'Gapdh', 'Actb', # C/T-independent
    #'Peg3', 'Zfp600', 'Zfp229', 'Mir17hg', '2610305D13Rik', '9830147E19Rik',
    ]

expn.draw_scatter_CV(filename='s1.scatter_CV.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-2, 2**17], vlines=[sam_map.low_expressed_threshold,],
    hlines=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])
expn.save('real_data_set_step1.glb')

expn = expn.filter_by_CV(2, 20)
expn.draw_scatter_CV(filename='s2.scatter_CV_afterCV.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-2, 2**17], vlines=[sam_map.low_expressed_threshold,],
    hlines=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])
expn.save('real_data_set_step2.glb')

expn = expn.filter_low_expressed(10**1.6, 7) # The magic number
expn.draw_scatter_CV(filename='s3.scatter_CV_afterLOW.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-2, 2**17], vlines=[sam_map.low_expressed_threshold,],
    hlines=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])
expn.save('real_data_set_step3.glb')

expn = expn.filter_high_expressed(2**7, 110)
expn.draw_scatter_CV(filename='s4.scatter_CV_afterHIGH.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-2, 2**17], vlines=[sam_map.low_expressed_threshold,],
    hlines=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

# The below feels too crude
#expn = expn.filter_by_mean_expression(2**1, 2**10)
#expn.draw_scatter_CV(filename='s5.scatter_CV_afterHIGHmean.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
#    ylims=[0,17], xlims=[2**-2, 2**17], vlines=[sam_map.low_expressed_threshold,],
#    hlines=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

expn.log(2,.1)
expn.save("real_data_set.glb")
expn[0:100].save("test_data_set.glb")

expn = expn.getColumns(['ensg', 'name'], strip_expn=True)
expn.saveTSV("real_data_set.tsv")
