"""

Initialise glbase, import all the libraries, set up the environment etc.

Requires:
* numpy
* matplotlib
* scipy
* sklearn
* h5py
"""

import sys, os, logging
from pkgutil import iter_modules

available_modules = list((name for loader, name, ispkg in iter_modules()))

#-----------------------------------------------------------------------
# Load all of the global configuration options.
try:
    from . import config
    from .errors import LibraryNotFoundError
except:
    print("Error: Fatal - glbase3 is not installed correctly, cannot find my own libraries")
    print("       Is the python 'sys.path' correct?")
    sys.exit() # no raise if I can't get errors, it's surely a fatal installation problem.

# ----------------------------------------------------------------------
# Test for availability of the core non-standard libs.
# These need to be available as the subsequent load/checking is weak/non-existent.

if 'numpy' in available_modules:
    config.NUMPY_AVAIL = True
else:
    raise LibraryNotFoundError("Fatal - Numpy is not available or not installed")

if 'scipy' in available_modules:
    config.SCIPY_AVAIL = True
else:
    raise LibraryNotFoundError("Fatal - Scipy is not available or not installed")

try:
    import matplotlib
    matplotlib.use("Agg") # cluster friendly!
    config.MATPLOTLIB_AVAIL = True
except ImportError:
    raise LibraryNotFoundError("Fatal - matplotlib not available or not installed")

if 'sklearn' in available_modules:
    config.SKLEARN_AVAIL = True
else:
    raise LibraryNotFoundError("Fatal - sklearn not available or not installed")

if 'h5py' in available_modules:
    config.H5PY_AVAIL = True
else:
    config.log.warning('h5py not available or not installed')

if 'networkx' in available_modules:
    config.NETWORKX_AVAIL = True
    # pass silently as pygraphviz is optional.

if 'pygraphviz' in available_modules:
    config.PYGRAPHVIZ_AVAIL = True
    # pass silently as pygraphviz is optional.

if 'graphviz' in available_modules: # sometimes comes in the wrong namespace!
    config.PYGRAPHVIZ_AVAIL = True
    # pass silently as pygraphviz is optional.

if 'pydot' in available_modules:
    config.PYDOT_AVAIL = True
    # pass silently as pydot is optional.

if 'statsmodels' in available_modules:
    config.STATSMODELS_AVAIL = True

#try:
#    import numexpr
#    config.NUMEXPR_AVAIL = True
#except Exception:
#    pass # pass silently as numexpr is optional.

if 'umap' in available_modules:
    config.UMAP_LEARN_AVAIL = True
    umap_log = logging.getLogger("umap")
    umap_log.setLevel(logging.CRITICAL) # silence debug output
else:
    pass # pass silently as umap is optional.

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
from .flat_track import flat_track
from .progress import progressbar
from .pwm import pwm
from .pwms import pwms
from .expression import expression
from .logos import logo
from .draw import draw
from .format_container import fc
from .fastq import fastq
from .draw import adjust_text
from .hic import hic, merge_hiccys
from .massspec import massspec
from .glgo import glgo
#from .ecrbase import ecrbase, tfbs_iter
#from .region import region
#from .intervaltree import intervaltree # Later integrate into genelist; expose here for now
from . import realtime
from . import utils
from . import format
from . import cmaps

from .tools.wigstep_to_flattrack import wigstep_to_flat
from .tools.merge_flat_track import merge_flats
from .tools.gerp_to_flattrack import gerp_to_flat
from .tools.bedgraph_to_flattrack import bedgraph_to_flat
from .tools.bed_to_flattrack import bed_to_flat
from .tools.wig_to_flattrack import wig_to_flat
from .tools.rnaseq import rnaseqqc

def version():
    config.log.info("glbase3 - version: {} {}".format(config.version, config.DATE))
    config.log.info("The working directory is: '{}'".format(os.getcwd()))

config.set_log_level('info')

# export all of the libraries, methods and helpers.
__all__ = [
    "genelist",
    "fastq",
    "expression",
    "genome",
    "genome_sql",
    #"track", # Deprecated. use flat_track
    "flat_track",
    "delayedlist",
    'massspec',
    "glglob",
    "hic",# primary objects
    'config',
    'merge_hiccys', # hic support
    "location",
    "pwm", "pwms", # PWM object support
    "format",
    "glload",
    # Tools/
    "wigstep_to_flat", "bedgraph_to_flat", 'bed_to_flat', 'wig_to_flat', "gerp_to_flat", 'merge_flats',
    "logo",
    "motif",
    "rnaseqqc",
    "glgo",
    "fc",
    #"rigidgrid", # Unavailable
    #"ecrbase",
    #"region",
    "realtime",
    #"realtime2", # This is deprecated
    #"tfbs_iter",
    'intervaltree',  # Useful utils to export
    "utils",
    'adjust_text',
    "change_drawing_mode",
    "progressbar",
    "draw",
    "cmaps", "strandSorter",# Miscellaneous
    "fold2UpOrDown", "fold2Down", 'fold2Up', 'XDown', 'XUp', 'lst_find', 'cat_columns', 'strandSorter'
    ]
