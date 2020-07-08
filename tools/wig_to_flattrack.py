#! /usr/bin/python
"""

Part of glbase,

"""

import sys, os, csv, time, numpy

import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def wig_to_flat(infilenames, outfilename, name, gzip=False, **kargs):
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
        infilenames
            a filename, or list of filenames (if you want to merge wiggles)

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
    assert infilenames, "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"
    if not isinstance(infilenames, list):
        infilenames = [infilenames]

    n = 0
    flat = flat_track(filename=outfilename, new=True, name=name, bin_format='f') # always f

    s = time.time()

    config.log.info('Preparse BED(s)')

    cleft = 0
    open_mode = None
    if not gzip:
        open_mode = open
    else:
        open_mode = opengzip.open

    chrom_name = None
    list_of_chroms = {}
    for f in infilenames:
        config.log.info("Started %s " % (f, ))
        oh = open_mode(f, 'rt')

        # Pre-parse to grab the maximum chromosome sizes:
        for line in oh: # Skip glbase overhead:
            if "#" in line:
                continue
            if 'variableStep' in line:
                # commit the size of the last chrom:
                if chrom_name and list_of_chroms[chrom_name] < int(pos):
                    list_of_chroms[chrom_name] = int(pos)+1 # store the chrom size

                chrom_name = str(line.split(' ')[1].split('=')[1]).strip()

                config.log.info('Found %s' % chrom_name)
                if chrom_name not in list_of_chroms:
                    list_of_chroms[chrom_name] = 0 # length of chrom;
                continue

            line = line.split('\t')
            pos = line[0]
        oh.close()
        # store the final chrom
        if chrom_name and list_of_chroms[chrom_name] < int(pos):
            list_of_chroms[chrom_name] = int(pos)+1 # store the chrom size

    config.log.info('Observed chromsomes sizes:')
    for ch in sorted(list_of_chroms):
        config.log.info('Chromosome: %s = %s' % (ch, list_of_chroms[ch]))

    # Go back through, once for each chromosome to make a numpy array for each chrom
    for chrom in list_of_chroms:
        config.log.info('Building array for %s' % chrom)
        arr = [0] * (list_of_chroms[chrom]+1000) # add a small pad for edge cases
        for f in infilenames:
            config.log.info("Started %s" % (f, ))
            oh = open_mode(f, 'rt')

            record = False
            for line in oh:
                if "#" in line:
                    continue
                if 'variableStep' in line:
                    # commit this chrom:
                    if record: # If I get here again, and record = True then the chrom is finished
                        break

                    # get new chrom
                    this_chrom_name = str(line.split(' ')[1].split('=')[1]).strip()
                    if this_chrom_name == chrom:
                        record = True
                        continue

                if record:
                    line = line.strip().split('\t')
                    pos = int(line[0])
                    sco = float(line[1])
                    arr[pos] += sco

                    n += 1
                    if n % 10000000 == 0:
                        config.log.info('Processed: {:,} bps'.format(n))
                        #break

            oh.close()

        # take the average
        arr = numpy.array(arr)
        arr /= len(infilenames)

        # Done this chrom, commit the array:
        config.log.info('Finished %s' % chrom)
        flat.add_chromosome_array(chrom, arr)

    config.log.info("Finalising library")
    flat.finalise()

    e = time.time()
    config.log.info("Took: %.1f seconds" % (e-s))
    return(True)

