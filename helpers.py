"""

helpers.py

a selection of functions to help various parts of the code.

Used in places like the CSV format language and in selecting criteria, etc...

"""

import cPickle, math, sys, os
import utils, config

from data import *
from location import location
from errors import BadBinaryFileFormatError

# ----------------------------------------------------------------------
# Some helper functions

def glload(filename):
    """
    **Purpose**
        Load a glbase binary file
        (Actually a Python pickle)
        
    **Arguments**
        filename (Required)
            the filename of the glbase binary file to load.
            
    **Returns**
        The glbase object previously saved as a binary file
    """
    assert os.path.exists(os.path.realpath(filename)), "File '%s' not found" % filename

    try:
        oh = open(os.path.realpath(filename), "rb")
        newl = cPickle.load(oh)
        oh.close()
    except cPickle.UnpicklingError:
        raise BadBinaryFileFormatError, filename
        
    # Recalculate the _optimiseData for old lists, and new features
    try:
        if newl.qkeyfind:
            pass
        if "loc" in newl[0].keys() or "tss_loc" in newl[0].keys(): # buckets are only present if a loc key is available.
            if newl.buckets: # added in 0.381, only in objects with tss_loc or loc key.
                pass
    except Exception:
        config.log.warning("Old glb format, will rebuild buckets and/or qkeyfind, consider resaving")
        newl._optimiseData()
        
    try:
        cons = len(newl._conditions) # expression-like object
        config.log.info("Loaded '%s' binary file with %s items, %s conditions" % (filename, len(newl), cons))
    except AttributeError:
        config.log.info("Loaded '%s' binary file with %s items" % (filename, len(newl)))   
    return(newl)

def change_drawing_mode(mode):
    """
    Change the output driver

    send me "png", "eps"|"ps" for your output needs
    
    This is a duplication of config.change_draw_mode() ...
    """
    assert mode in config.valid_draw_modes, "Draw output mode '%s' is not recognised" % mode

    config.draw_mode = mode

# ----------------------------------------------------------------------
# Criteria functions.

def fold2UpOrDown(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > 2:
                    return(True)
                elif normed_data < 0.5:
                    return(True)
                else:
                    return(False)
        return(False)

def fold2Down(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data < 0.5:
                    return(True)
                else:
                    return(False)
        return(False)

def fold2Up(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > 2:
                    return(True)
                else:
                    return(False)
        return(False)

def XDown(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    X = kargs["X"]
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data < X:
                    return(True)
                else:
                    return(False)
        return(False)

def XUp(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the Xdown data
    """
    X = kargs["X"]
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > X:
                    return(True)
                else:
                    return(False)
        return(False)

# For formatting the CSV loading.

def lst_find(lst, predicate): # I need a helper function to find the item
    return (i for i, j in enumerate(lst) if predicate(j)).next()

def cat_columns(c1, c2, sep=' '):
    # concatenate two columns together
    return('%s%s%s' % (c1, sep, c2))

def strandSorter(chr, left, right, strand):
    """
    A helper proc to extract the tss from a list of coords.
    """
    if strand in positive_strand_labels:
        return(location(chr=chr, left=left, right=left))
    elif strand in negative_strand_labels:
        return(location(chr=chr, left=right, right=right))
    return(None)

# various other helpers for normalisation etc..
