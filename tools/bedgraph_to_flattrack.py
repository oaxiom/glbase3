#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

"""

import sys, os, csv, time

import numpy
import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def bedgraph_to_flat(infilenames, outfilename, name, gzip=None, all_in_mem=False, **kargs):
    """
    **Purpose**
        Convert a bedGraph file to a flat file (Actually an h5 db)

        A bedGraph file should look like this:

        chr1	3000881	3000986	1
        chr1	3000986	3001081	2
        chr1	3001081	3001157	1
        chr1	3001157	3001186	2
        chr1	3001186	3001357	1
        chr1	3001535	3001735	1
        chr1	3002009	3002209	1
        chr1	3003894	3003917	1
        chr1	3003917	3004094	2


    **Arguments**
        infilename
            the name of the bedGraph file to read from.

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track

        gzip (Optional, default=False)
            The input file(s) is a gzip file.

        all_in_mem (Optional, default=False)
            If you have a lot of memory, set this to True and the flat creation is done
            all in memory and only committed to disk right at the end.

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert infilenames, "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"
    if not isinstance(infilenames, list):
        infilenames = [infilenames]

    bin_format = 'f' # Only one supported
    n = 0
    m = 0
    total = 0

    config.log.info("Started %s -> %s" % (infilenames, outfilename))

    open_mode = None
    open_mode = open if not gzip else opengzip.open
    s = time.time()
    chrom_arrays = {}
    chrom_sizes = {}
    num_bps = 0

    # Preparse to get chrom sizes;
    config.log.info('Preparse to get chromsome sizes')
    for f in infilenames:
        oh = open_mode(f, 'rt')

        for line in oh:
            if "#" in line:
                continue
            line = line.split()
            chrom = line[0]
            left = int(line[1])
            rite = int(line[2])
            if chrom not in chrom_sizes:
                chrom_sizes[chrom] = 0
            if rite > chrom_sizes[chrom]:
                chrom_sizes[chrom] = rite+1000
            num_bps += (rite - left) # bodge

    for chrom in chrom_sizes:
        config.log.info('{} = {:,} bp'.format(chrom, chrom_sizes[chrom]))
        chrom_arrays[chrom] = numpy.zeros((chrom_sizes[chrom]))

    config.log.info('{:,} bps (estimated) to do'.format(num_bps))

    config.log.info('Collect data')
    for f in infilenames:
        config.log.info("Started {}".format(f))
        oh = open_mode(f, 'rt')

        for line in oh:
            if "#" in line:
                continue

            line = line.strip().split()
            chrom = line[0]
            left = int(line[1])
            rite = int(line[2])
            score = float(line[3])

            chrom_arrays[chrom][left:rite] = score

            if n % 1e6 == 0:
                print('Processed: {:,} bp'.format(n))
            n += (rite - left)

    # Commit:
    flat = flat_track(filename=outfilename, new=True, name=name, bin_format=bin_format)
    for chrom in chrom_arrays:
        config.log.info('Committing: {}'.format(chrom))
        flat.add_chromosome_array(chromosome=chrom.replace('chr', ''), arr=chrom_arrays[chrom])
    del chrom_arrays

    e = time.time()

    config.log.info("Finalise library. Contains '{:,}' bps of data".format(n))
    flat.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return True

