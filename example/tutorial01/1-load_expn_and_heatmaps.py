"""

Code from tutorial 1

"""

import os
from glbase3 import *

config.draw_mode = 'pdf'

#docs_path = "../../docs/source/images/" # I want to save the image to the docs folder so it appears in the documentation.
#config.draw_size = "small" # I set this so the figure sizes are nice for the docs

expn_format = {"force_tsv": True, "enst": 0, "name": 1, "tf_type": 2}

expn = expression(filename="../shared_raw_data/macrophage-all-tfdbd.tsv.gz", format=expn_format, expn="column[3:5]", gzip=True)
print(expn)
print(expn.getConditionNames())

expn = expn.add_fc_key("fc", "pec_untre", "pec_il10")
print(expn)
print(expn.find("Stat3"))

expn.stats.print_stats()

expn.hist(filename="1.hist_no_log.pdf")

expn.mult(1e6)
expn.log(2)
expn.stats.print_stats()

expn.hist(filename="1.hist_after_log.pdf", range=[-5, 10])

expn = expn.removeDuplicates("name")

expn.pie(filename="1.tf_type_pie.pdf", key="tf_type", font_size=6)
expn.pie(filename="1.strand_pie.pdf", key="tf_type", font_size=6)

expn.heatmap(filename="1.heatmap.pdf")

expn.sort("fc")
print(expn)

top30 = expn[0:30]
top30.heatmap(filename="1.heatmap-top30.pdf")

expn.reverse()
bottom30 = expn[0:30]
bottom30.heatmap(filename="1.heatmap-bot30.pdf")

expn = expn.normaliseToCondition("pec_untre")
expn.sort("pec_il10")

expn.heatmap(filename="1.heatmap-normed.pdf", bracket=[0,1], row_cluster=False)

expn.reverse()
top30 = expn[0:30]
top30.heatmap(filename="1.heatmap-normed-top30.pdf", row_cluster=False, bracket=[0,1])

expn.reverse()
bottom30 = expn[0:30]
bottom30.heatmap(filename="1.heatmap-normed-bot30.pdf", row_cluster=False, bracket=[0,1])
