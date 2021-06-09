
"""

This is the code from the publication.

relative paths were removed in the text for clarity.

"""


import glbase3 as gl

peaks = gl.genelist(filename="../shared_raw_data/macs_list.xls.gz", format=gl.format.macs_summit, gzip=True)
bed = gl.genelist(filename="../shared_raw_data/Sox2_Oct4_ol_w100_annot.bed.gz", format=gl.format.bed, gzip=True)
fasta = gl.genelist(filename="../shared_raw_data/Fasta_file.fa.gz", format=gl.format.fasta, gzip=True)

print(peaks)

peaks = peaks[0:1] # list was truncated for clarity:
print(peaks.pointify()) # Take the middle point of the interval
print(peaks.expand("loc", 100)) # expand the left and right border by 100 bp

print(bed)

print(fasta)
