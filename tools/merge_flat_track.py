#! /usr/bin/python
"""

Part of glbase,

MErge flat tracks to a single flat

"""

import sys, os, csv, time, numpy
import gzip as opengzip
from collections.abc import Iterable
from .. import flat_track
from .. import config
from .. import location

def merge_flats(
    infilenames:list,
    outfilename:str,
    name:str,
    mode:str='average',
    skip_non_standard_chroms=False,
    **kargs):
    """
    **Purpose**
        Convert a list of genomic coordinates to a flat file (Actually an SQL
        db)

    **Arguments**
        infilenames
            a list of .flats to merge

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track

        mode (Optional, default='average')
            Way to merge the flats.

            average = take the mean of the tracks.
            sum = take the sum of all tracks.

        skip_non_standard_chroms (Optional, default=False)
            only use canonical chromsome names (i.e. chr1, chrX, chrM)
            and not scaffolds and unplaced (chr14_KIA... etc).

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert infilenames, "no filename specified"
    assert isinstance(infilenames, Iterable), 'infilenames is not an iterable'
    assert os.path.realpath(outfilename), "no save filename specified"

    n = 0
    m = 0
    total = 0

    out = flat_track(filename=outfilename, new=True, name=name, bin_format='f')

    merged_chrom_data = {}
    largest_lengths = {}

    s = time.time()
    for f in infilenames:
        f = flat_track(filename=f, new=False, name=name, bin_format='f')
        #print(largest_lengths)

        for chrom in f.get_all_chrom_names():
            config.log.info('Doing {}'.format(chrom))
            if skip_non_standard_chroms and '_' in lastchrom:
                continue

            chrom_array_to_add = numpy.array(f.get_array_chromosome(chrom))

            if chrom not in merged_chrom_data:
                merged_chrom_data[chrom] = chrom_array_to_add
                largest_lengths[chrom] = chrom_array_to_add.shape[0]
                continue

            largest_lengths[chrom] = max(largest_lengths[chrom], chrom_array_to_add.shape[0])
            # Check the sizes and pad if needed;
            if merged_chrom_data[chrom].shape[0] < largest_lengths[chrom]:
                merged_chrom_data[chrom] = numpy.pad(merged_chrom_data[chrom], (0, largest_lengths[chrom]-merged_chrom_data[chrom].shape[0]), constant_values=0)

            if chrom_array_to_add.shape[0] < largest_lengths[chrom]:
                chrom_array_to_add = numpy.pad(chrom_array_to_add, (0, largest_lengths[chrom]-chrom_array_to_add.shape[0]), constant_values=0)

            merged_chrom_data[chrom] += numpy.array(chrom_array_to_add)

        del f

    e = time.time()
    config.log.info("Took: {:.2f} seconds".format(e-s))

    # Write out merge
    for chrom in merged_chrom_data:
        chrom_array = numpy.array(merged_chrom_data[chrom])
        if mode == 'average':
            chrom_array /= len(infilenames)

        out.add_chromosome_array(chrom, chrom_array)

    out.finalise()
    config.log.info("Finalised library.")

    return True

