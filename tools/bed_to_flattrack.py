#! /usr/bin/env python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

Uses the new flat_track2

"""

import sys, os, csv, time, numpy

import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def is_pe_inner_loop(f, chr_sizes, infilename, gzip, tot_tag_count):
    # PE version, assumes l and r are frags
    n = 0
    for ch in sorted(chr_sizes.keys()):
        config.log.info('Extracting %s' % ch)
        this_chrom = [0] * (chr_sizes[ch]+1)

        for file in infilename:
            config.log.info("Collecting from %s" % (file,))
            oh = open(file, "rt") if not gzip else opengzip.open(file, 'rt')
            for line in oh: # Skip glbase overhead:
                if "#" in line:
                    continue

                line = line.split('\t')
                if line[0] != ch:
                    continue
                #print line

                # We assume a strict BED file
                l = int(line[1])
                r = int(line[2])

                for bp in range(l, r):
                    this_chrom[bp] += 1

                if n % 1000000 == 0 and n>0:
                    config.log.info("{0:,} tags parsed ({1:.1f}%)".format(n, n/tot_tag_count*100))
                n += 1 # need to do this
            oh.close()
        config.log.info('Committing %s' % ch)
        f.add_chromosome_array(ch.strip(), numpy.array(this_chrom))
        del this_chrom

def is_se_inner_loop(f, chr_sizes, infilename, gzip, read_extend, tot_tag_count):
    # SE version, assumes l and r need read extension based on strand
    n = 0
    for ch in sorted(chr_sizes.keys()):
        config.log.info('Extracting %s' % ch)
        this_chrom = [0] * (chr_sizes[ch]+1)

        for file in infilename:
            config.log.info("Collecting from %s" % (file,))
            oh = open(file, "rt") if not gzip else opengzip.open(file, 'rt')
            for line in oh: # Skip glbase overhead:
                if "#" in line:
                    continue

                line = line.split('\t')
                if line[0] != ch:
                    continue

                # We assume a strict BED file
                l = int(line[1])
                r = int(line[2])
                s = line[5].strip()
                if s == '+':
                    r += read_extend
                elif s == '-':
                    l -= read_extend
                else:
                    raise AssertionError('strand %s not found' % s)

                for bp in range(l, r):
                    this_chrom[bp] += 1

                if n % 1000000 == 0 and n>0:
                    config.log.info("{0:,} tags parsed ({1:.1f}%)".format(n, n/tot_tag_count*100))
                n += 1 # need to do this
            oh.close()
        config.log.info('Committing %s' % ch)
        f.add_chromosome_array(ch.strip(), numpy.array(this_chrom))
        del this_chrom

def bed_to_flat(infilename, outfilename, name, isPE, read_extend=None, strand=False, gzip=None, **kargs):
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

        Note, that the flats DO NOT use strand.

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

        isPE (Required)
            Does the BED contain PE fragments, or SE 5' ends orientated by strand?

        read_extend (Required if isPE=False)
            read_extension to perform if the reads are single ended.

        strand (Required to be True if isPE=False)
            If isPE is False then you must specify strand as True and a +/- strand
            must be in your column 6 of your BED (starting at column 1).

        gzip (Optional, default=False)
            The input file(s) is a gzip file.

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert infilename, "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"
    if isPE:
        read_extend = 0
        assert not strand, 'strand is ignored in isPE=True libraries'
    else:
        assert read_extend is not None, 'You must specify a read_extend if isPE=False'
        assert strand, 'if isPE=False, you must set strand=True and must have strang +/- in column 6 of the BED'

    if not isinstance(infilename, list):
        infilename = [infilename]

    n = 0

    bin_format = 'i'
    f = flat_track(filename=outfilename, new=True, name=name, bin_format=bin_format)

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
            r = int(line[2])+5000+read_extend # pad out any likely read extension

            if ch not in chr_sizes:
                chr_sizes[ch] = r
            if r > chr_sizes[ch]:
                chr_sizes[ch] = r
            if n % 1000000 == 0 and n>0:
                config.log.info("%s,000,000 reads parsed" % ((n // 1000000),))
            n += 1 # need to do this
        oh.close()
    total_read_count = int(n)
    f.set_total_num_reads(total_read_count)
    config.log.info('Total read count %s' % total_read_count)

    config.log.info('Observed chromsomes sizes:')
    for ch in chr_sizes:
        config.log.info('Chromosome: {0} = {1:,} bp'.format(ch, chr_sizes[ch]))

    if isPE:
        is_pe_inner_loop(f, chr_sizes, infilename, gzip, total_read_count)
        f.meta_data['isPE'] = True
    else:
        is_se_inner_loop(f, chr_sizes, infilename, gzip, read_extend, total_read_count)
        f.meta_data['isPE'] = False

    f.meta_data['total_read_count'] = total_read_count

    config.log.info("Finalising library")
    f.finalise()

    e = time.time()
    config.log.info("Took: %.1f seconds" % (e-s))
    return True
