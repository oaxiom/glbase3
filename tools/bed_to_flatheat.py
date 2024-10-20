#! /usr/bin/env python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

Uses the new flat_track2

"""

import sys, os, csv, time, numpy

from math import floor

import gzip as opengzip
from .. import flat_heat
from .. import config
from .. import location

__warning_lr_not_sorted = False

def is_pe_inner_loop(f, chr_sizes, infilename, gzip, tot_tag_count, ybins, ymax):
    # PE version, assumes l and r are frags
    n = 0
    for ch in sorted(chr_sizes.keys()):
        config.log.info(f'Extracting {ch}')
        
        # Make a 2D array:
        this_chrom = numpy.zeros((chr_sizes[ch] // 10, ybins)) 

        for file in infilename:
            config.log.info(f"Collecting from {file}")
            oh = open(file, "rt") if not gzip else opengzip.open(file, 'rt')
            for line in oh: # Skip glbase overhead:
                if "#" in line:
                    continue

                line = line.split('\t')
                if line[0] != ch:
                    continue

                # We assume a strict BED file
                l = int(line[1]) // 10 # in bins
                r = int(line[2]) // 10 # in bins
                w = int(line[2]) - int(line[1]) # Not divided by 10 into bins

                if w > ymax:
                    # issue warning?
                    continue
                    
                if w < 0:
                    if not __warning_lr_not_sorted:
                        config.log.waring('right coordinate is less than left')
                        __warning_lr_not_sorted = True
                    continue
                    
                ybin = floor((w / ymax) * ybins)# in bins
                
                for bp in range(l, r):
                    this_chrom[bp, ybin] += 1 # chrom data is per 10 bp.

                if n % 1e6 == 0 and n>0:
                    config.log.info("{0:,} tags parsed ({1:.1f}%)".format(n, n/tot_tag_count*100))
                n += 1 # need to do this
            oh.close()
        #config.log.info('Committing %s' % ch)
        f.add_chromosome_array(ch.strip(), numpy.array(this_chrom))
        del this_chrom

def bed_to_flatheat(
    infilename:str,
    outfilename:str,
    name:str,
    ymax:int,
    gzip:bool=False,
    ybins:int=50,
    ):
    """
    **Purpose**
        Convert a bed file containing reads into a flat
        db)

        A bed file should look like this:

        chr1    3000881 3000986
        chr1    3000986 3001081
        chr1    3001081 3001157
        chr1    3001157 3001186
        chr1    3001186 3001357
        chr1    3001535 3001735
        chr1    3002009 3002209
        chr1    3003894 3003917
        chr1    3003917 3004094

        Note, that the flatheats DO NOT use strand.
        
        flatheats must be from paired-end data

        Note, also, that normalisation is NOT required at generation, and can instead be
        performed at runtime.

    **Arguments**
        infilename
            the name of the bed file(s) to read from.

            You can send a list of filenames

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track


        ybins (int, optional, default=50)
            number of bins for the y-axis
        
        ymax (int, Required)
            size of the y-axis in base pairs.

        gzip (Optional, default=False)
            The input file(s) is a gzip file.

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert infilename, "no filename specified"
    assert outfilename, "no save filename specified"

    if not isinstance(infilename, list):
        infilename = [infilename]

    n = 0

    bin_format = 'i'

    f = flat_heat(filename=outfilename, new=True, name=name, ymax=ymax, ybins=ybins)

    config.log.info("Started %s -> %s" % (infilename, outfilename))

    s = time.time()
    config.log.info('Preparse BED(s)')
    chr_sizes = {}
    cleft = 0
    for file in infilename:
        config.log.info("Collecting %s" % (file,))

        oh = open(file, "rt") if not gzip else opengzip.open(file, 'rt')
        # Pre-parse to grab the maximum chromosome sizes:
        for line in oh: # Skip glbase overhead:
            if "#" in line:
                continue
            line = line.split('\t')
            # We assume a strict BED file,
            ch = line[0]
            r = int(line[2])+5000 # pad out any likely read extension

            if ch not in chr_sizes:
                chr_sizes[ch] = r
            if r > chr_sizes[ch]:
                chr_sizes[ch] = r
            if n % 1000000 == 0 and n>0:
                config.log.info("{:,}M reads parsed".format((n // 1e6),))
            n += 1 # need to do this
        oh.close()
    total_read_count = int(n)

    f.set_total_num_reads(total_read_count)

    config.log.info(f'Total read count {total_read_count:,}')

    config.log.info('Observed chromsomes sizes:')
    for ch in chr_sizes:
        config.log.info(f'Chromosome: {ch} = {chr_sizes[ch]:,} bp')

    config.log.info("Building library")

    is_pe_inner_loop(f, chr_sizes, infilename, gzip, total_read_count, ybins, ymax)
    f.meta_data['isPE'] = True
    f.meta_data['total_read_count'] = total_read_count
    f.finalise()

    e = time.time()
    config.log.info("Took: %.1f seconds" % (e-s))
    return True
