"""
config.py

config must be imported before any other glbase library.

"""

# config cannot import any other glbase module

import os, logging, subprocess

# -------------- Versioning data


try:
    oh = open(os.path.join(os.path.split(__file__)[0], "version.num"), "rU")
    VERSION = "1.%s" % oh.readline().strip().replace("+", "") # The hg hook will put a plus as I call just before committing.
    oh.close()
except Exception:
    VERSION = "version data not found"
DATE = ""

__version__ = VERSION # There's a lot of needless redundancy here...
version = VERSION # this is the valid one I want to use in future.

# -------------- General options

SILENT = False # set this to True to silence all glbase output. Only works at startup
DEBUG = True
do_logging = True

# flags for the availability of libraries
MATPLOTLIB_AVAIL = False # required
NUMPY_AVAIL = False # required
SCIPY_AVAIL = False # required
SKLEARN_AVAIL = False # required
H5PY_AVAIL = False # Optional
NETWORKX_AVAIL = False # optional
PYDOT_AVAIL = False # optional
NUMEXPR_AVAIL = False # Optional
PYGRAPHVIZ_AVAIL = False # Optional
UMAP_LEARN_AVAIL = False # Optional

# Some simple options for printing genelists
NUM_ITEMS_TO_PRINT = 3 # number of items to print by default.
PRINT_LAST_ITEM = True

# size of buckets for collide() and overlap()
# If this is changed then glload will not work correctly.
bucket_size = 10000 # in bp - tested, seems a reasonable choice.

DEFAULT_DPI = 150 # not working?
draw_mode = "png"
draw_size = "medium"
draw_aspect = "normal"
valid_draw_modes = frozenset(["png", "ps", "eps", "svg", 'pdf', 'tiff'])

# -------------- Start of the new-style options:

class flat:
    block_size = 2000
    cache = 100000 # maximum number of blocks to keep in memory.
    # these are equivalent to about 800 Mb's of memory on a Fedora box
    # And 650 Mbs on Mac OSX. Never tested on Windows

def change_draw_size(size):
    # These would be better as attached to the variable? __call__() ?
    # This is usually correctly respected.
    assert size in ["small", "medium", "large", "huge"], "drawing size '%s' not found" % size
    draw_size = size

def change_draw_aspect(aspect):
    # This is often ignored though by the drawing methods themselves...
    assert size in ["normal", "square", "long"], "drawing aspect '%s' not found" % size
    draw_size = aspect

def get_interpolation_mode(filename):
    # Get a suitable interpolator for non-aliased heatmaps and hist2d's, etc.
    # config.draw_mode is no longer reliable, so I have to get it from the actual filename now:
    tmode = filename.split('.')[-1]
    if tmode in ("svg", 'pdf', 'eps', 'ps'):
        return('nearest') # Yes, it really is nearest. Otherwise it will go to something like bilinear
    return('nearest') # seems fixed in newer matplotlibs?

# -------------- set up the logger here.
# You can access it using config.log()
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)-8s: %(message)s',
                    datefmt='%m-%d %H:%M')

# use config.log. ... () to get to the logger
log = logging.getLogger('glbase3')

mpl_logger = logging.getLogger('matplotlib') # Bodge to silence the matplotlib logging
mpl_logger.setLevel(logging.WARNING)

# helpers [Should be deprecated. Call config.log.<level>() to call info]
info = log.info
warning = log.warning
debug = log.debug
error = log.error

def silence_log(): # I think a NullHandler would be better for this?
    do_logging = False
    # by pointing to empty lambdas:
    log.info = lambda x:x
    log.warning = lambda x:x
    log.error = lambda x:x
    log.debug = lambda x:x
    # Keep critical

if SILENT:
    log.setLevel(logging.CRITICAL) # not acutally silenced...

def set_log_level(level):
    """
    Change the logging level for the on-screen logger.
    the console logger is not affected by this call.
    """
    level_map = {None: logging.CRITICAL, # not correct?
        "info": logging.INFO,
        "debug": logging.DEBUG}

    assert level in level_map, "no valid level used to set the logger, valid modes are 'info' or 'debug'"

    log.setLevel(level_map[level])
    return(True)
