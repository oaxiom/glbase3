#! /usr/bin/python
"""

Part of glbase,

"""

import sys, os, csv, time, numpy

import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def wig_to_flat(infilename, outfilename, name, gzip=None, **kargs):
    """
    **Purpose**
        Convert a variableStep wig file containing reads into a flat db)

        variableStep chrom=chr10
        3131448 0.00975
        3131449 0.00975
        3131450 0.00975
        3131451 0.00975
        3131452 0.00975
        3131453 0.00975
        3131454 0.00975
        3131455 0.00975
        3131456 0.00975

        Note, that the flats DO NOT use strand.
        
        Note, also, that normalisation is NOT required at generation, and can instead be 
        performed at runtime.

    **Arguments**
        infilename
            the name of the bed file(s) to read from. 
            
            You can send a lit of filenames

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track

        gzip (Optional, default=False)
            The input file(s) is a gzip file.
        
    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert infilename, "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"
    
    n = 0
    f = flat_track(filename=outfilename, new=True, name=name, bin_format='f') # always f

    config.log.info("Started %s -> %s" % (infilename, outfilename))

    s = time.time()
    
    config.log.info('Preparse BED(s)')
    chr_sizes = {}
    this_chrom = []
    this_chrom_pos = 0
    this_chrom_name = None
    
    cleft = 0
    # Can only have one file
    if not gzip:
        oh = open(infilename, "rU")
    else:
        oh = opengzip.GzipFile(infilename, 'rt')

    # Pre-parse to grab the maximum chromosome sizes:
    for line in oh: # Skip glbase overhead:
        if "#" in line:
            continue
        if 'variableStep' in line:
            # commit last chrom:
            if this_chrom_name:
                f.add_chromosome_array(this_chrom_name.strip(), this_chrom) 
                chr_sizes[this_chrom_name] = len(this_chrom)
            # get new chrom
            this_chrom_name = str(line.split(' ')[1].split('=')[1]).strip()
            config.log.info('Parsing chr=%s' % this_chrom_name)
            this_chrom_pos = 0
            this_chrom = []
            continue
                       
        line = line.split('\t')
        pos = int(line[0])
        sco = float(line[1])
        while pos >= len(this_chrom):
            this_chrom.append(0.0)
        this_chrom[pos] = sco

        if n % 1000000 == 0 and n>0:
            config.log.info("%s,000,000 bps parsed" % ((n // 1000000),))
        
        n += 1 # need to do this 
    oh.close()
    total_read_count = int(n)
    config.log.info('Total bp count %s' % total_read_count)
    
    config.log.info('Observed chromsomes sizes:')
    for ch in chr_sizes:
        config.log.info('Chromosome: %s = %s' % (ch, chr_sizes[ch]))

    config.log.info("Finalising library")
    f.finalise()

    e = time.time()
    config.log.info("Took: %.1f seconds" % (e-s))
    return(True)
    
