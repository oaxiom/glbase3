"""

Initialise glbase, import all the libraries, set up the environment etc.

Requires:
* numpy
* matplotlib
* scipy
* sklearn
* h5py
* networkx
"""

import sys, os

#-----------------------------------------------------------------------
# Load all of the global configuration options.
try:
    from . import config
    from .errors import LibraryNotFoundError
except:
    print("Error: Fatal - GLbase is not installed correctly, cannot find my own libraries")
    print("       Is the python 'sys.path' correct?")
    sys.exit() # no raise if I can't get errors, it's surely a fatal installation problem.

# ----------------------------------------------------------------------
# Test for availability of the core non-standard libs.
# These need to be available as the subsequent load/checking is weak/non-existent.

try:
    import numpy
    config.NUMPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - Numpy is not available or not installed")

try:
    import scipy
    config.SCIPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - Scipy is not available or not installed")

try:
    import matplotlib
    matplotlib.use("Agg") # cluster friendly!
    config.MATPLOTLIB_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - matplotlib not available or not installed")

try:
    import sklearn
    config.SKLEARN_AVAIL = True
except Exception:
    raise LibraryNotFoundError("Fatal - sklearn not available or not installed")

try:
    import h5py
    config.H5PY_AVAIL = True
except Exception:
    config.log.warning('Fatal - h5py not available or not installed')

try:
    import networkx
    config.NETWORKX_AVAIL = True
except Exception:
    config.log.warning('Fatal - networkx not available or not installed') # pass silently as networkx is optional.

try:
    import pygraphviz
    config.PYGRAPHVIZ_AVAIL = True
except Exception:
    pass # pass silently as pygraphviz is optional.

try:
    import graphviz # sometimes comes in the wrong namespace!
    config.PYGRAPHVIZ_AVAIL = True
except Exception:
    pass # pass silently as pygraphviz is optional.

try:
    import pydot
    config.PYDOT_AVAIL = True
except Exception:
    pass # pass silently as pygraphviz is optional.

try:
    import numexpr
    config.NUMEXPR_AVAIL = True
except Exception:
    pass # pass silently as numexpr is optional.

# ----------------------------------------------------------------------
# Now import the rest of my libraries - assumes here they are available.
# If I can get config and errors then these are probably available too.

from .helpers import glload, change_drawing_mode, fold2UpOrDown, fold2Down, fold2Up, XDown, XUp, lst_find, cat_columns, strandSorter
from .location import location
from .genelist import genelist
from .expression import expression
from .genome import genome
from .genome_sql import genome_sql # To replace genome
from .delayedlist import delayedlist
from .glglob import glglob
from .element import motif
from .track import track
from .flat_track import flat_track
from .progress import progressbar
from .pwm import pwm
from .pwms import pwms
from .ecrbase import ecrbase, tfbs_iter
from .region import region
from .expression import expression
from .logos import logo
from .draw import draw
from .format_container import fc
from .fastq import fastq
from .glgo import glgo
from .draw import adjust_text
from .hic import hic, merge_hiccys
from . import realtime
from . import gldata
from . import utils
from . import format
from . import cmaps

from .tools.seqToTrk import seqToTrk
from .tools.wigstep_to_flattrack import wigstep_to_flat
from .tools.gerp_to_flattrack import gerp_to_flat
from .tools.bedgraph_to_flattrack import bedgraph_to_flat
from .tools.bed_to_flattrack import bed_to_flat
from .tools.wig_to_flattrack import wig_to_flat
from .tools.rnaseq import rnaseqqc

def version():
    config.log.info("glbase - version: %s %s" % (config.version, config.DATE))
    config.log.info("The working directory is: '%s'" % (os.getcwd()))

config.set_log_level('info')

# export all of the libraries, methods and helpers.
__all__ = ["genelist", "fastq", "expression", "genome", "genome_sql", "track", "flat_track", "delayedlist",
            "glgo", "hic", # primary objects
            'config',
            #"rigidgrid", # Temporarily unavailable
            'merge_hiccys', # hic support
            "location",
            "pwm", "pwms", # PWM object support
            "flags",  "format",
            "utils", "glload", "seqToTrk", "logo",
            "glglob", "motif",  "wigstep_to_flat", "bedgraph_to_flat", 'bed_to_flat', 'wig_to_flat',
            "rnaseqqc", "gldata",
            "gerp_to_flat", "draw", "fc",
            "progressbar", "ecrbase", "region", "realtime",
            #"realtime2", # This is deprecated
            "tfbs_iter",
            "strandSorter",
            'adjust_text',
            "cmaps",
            "change_drawing_mode", "fold2UpOrDown", "fold2Down", 'fold2Up', 'XDown', 'XUp', 'lst_find', 'cat_columns', 'strandSorter'
            ]
            # in future I want to get rid of dir() and control what gets exported.
