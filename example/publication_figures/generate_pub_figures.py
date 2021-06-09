"""

This code produces all of the figures, as they appear in the glbase paper.

To make the final figures the svg files were imported into iDraw and arranged.

Some of the axis legends were re-typed as text into iDraw to improve clarity.

This will not work out of the box as it needs some raw data downloaded from UCSC, or GEO.
You can see the relevant tutorials for where and how to get that data.

"""

import sys, os
user_root = os.path.expanduser("~")
from glbase3 import *
format.minimal_bed.update({"strand": 5})

config.draw_mode = "pdf"

# -------- Figure 1
# panels are out of order so I don't re-load data multiple times

# -- panel A
fa = genelist(filename="../tutorial7/STAT3_motifs_sequence_100bp.fa", format=format.fasta, name="STAT3")
fa2 = genelist(filename="../tutorial7/random_background_100bp.fa", format=format.fasta)
stat3 = motif(name="stat3", sequence="ttcnnngaa")
stat3.scanMotifFrequency(genelist=[fa], filename="figure-1a.png", random_fastas=[fa2], xlabel="Base pairs around peak summit")

# -- panel D
# Tutorial 1 Histogram and heatmap, top10/bott10?
expn_format = {"force_tsv": True, "enst": 0, "name": 1, "tf_type": 2}
pec_expn = expression(filename="../tutorial1/macrophage-all-tfdbd.tsv", format=expn_format, expn="column[3:5]")
pec_expn = pec_expn.add_fc_key("fc", "pec_untre", "pec_il10")
pec_expn.mult(1e6)
pec_expn.log(2, 1)
pec_expn = pec_expn.removeDuplicates("name")
pec_expn.sort("fc")
top = pec_expn[0:20]
pec_expn.reverse()
bot = pec_expn[0:20]
expn = top + bot
expn.sort("fc")
expn.conditions = ["-IL-10", "+IL-10"] # make them nice
expn.heatmap(filename="figure-1d.pdf", bracket=[-2, 6], row_cluster=True, col_cluster=False)
expn.stats.print_stats()

# -- panel B
chip = genelist(filename="../shared_raw_data/pec_pil10_peaks.xls", format=format.macs_summit, name="STAT3 ChIP-seq in macrophages")
randoms = [genelist(filename="../shared_raw_data/random_peaks-%s.bed" % i, format=format.minimal_bed) for i in [1,2]]
mm9 = glload("../shared_raw_data/mm9_refGene.glb")
chip.genome_distribution(filename="figure-1b.pdf", genome=mm9, random_lists=randoms)

# -- panel C
print(expn.conditions)
pec_expn.conditions = ["PEC macrophages -IL-10 expression log2(TPM)", "PEC macrophages +IL-10 expression log2(TPM)"] # make labels nice
pec_expn._optimiseData() # after the expn.conditions name change I need to rebuild the indeces...
newl = []
for item in pec_expn:
    if item["conditions"][0] == 0 or item["conditions"][1] == 0: # I trim these to make the figure draw faster in the pdf
        pass
    else:
        newl.append(item)
pec_expn.load_list(newl)
pec_expn.log(2, 1)
pec_expn.scatter("PEC macrophages -IL-10 expression log2(TPM)", "PEC macrophages +IL-10 expression log2(TPM)", filename="figure-1c.png", spot_size=30)
1/0
# -- panel E
pp = flat_track(filename=os.path.join(user_root, "mm9/phastCons/placental.flat"), bin_format="f", name="Placental conservation")
bed = genelist(filename="../tutorial4/Sox2_Oct4_ol_w100_annot.bed", name="Sox2/Oct4 binding sites", format=format.minimal_bed)
pileup_data = pp.pileup(filename="figure-1e.png", genelists=bed, mask_zero=True, bandwidth=500,
    xlabel="Position around Sox2/Oct4 binding site", ylabel="Average phastCons conservation score")

# -- panel F
t = flat_track(filename="../tutorial4/p300_reads.flat")
pileup_data = t.heatmap(filename="figure-1f.pdf", imshow=True, genelist=bed[0:300], distance=2000, bins=100)
