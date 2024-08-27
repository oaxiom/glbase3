"""
Utilities

Various utilities to support the genome scanning scripts.

Many of these predate glbase3, but are a little tricky to remove as I am not sure where
they are used (if at all).

So excuse the terrible code in places. I will deprecate occasional functions from this.

"""

import sys, os, numpy, string, csv, random, math, pickle, gzip
import scipy.stats as stats
from scipy.ndimage.filters import uniform_filter1d

from . import config

def library(args):
    """
    lovely generator...
    I'm still not sure exactly how this works...
    """
    if not args:
        yield ""
        return
    for i in args[0]:
        for tmp in library(args[1:]):
            yield i + tmp
    return

def expandDegenerateMotifs(motif):
    """
    expand a degeneratemotif into its constituant parts;
    returns a list;
    motif should be of the form:
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]
    e.g. "raangt"
    """
    mlen = len(motif)
    nm = motif.lower() # just in case;

    # scan the motif and count the no of degenerate motifs
    newlist = []
    for n in range(mlen):
        if nm[n] in ["r", "y", "k", "m", "s", "w"]: # AG
            newlist.append((2, n, nm[n])) # append a triple, with the number of flips and its location
        elif nm[n] == "n":
            newlist.append((4, n, nm[n]))
        else:
            newlist.append((1, n, nm[n]))

    if not newlist: return([motif]) # no degenerate elements in the motif

    l = []
    #print newlist
    iti(newlist, 0, None, l)

    return l
'''
def moving_average(_list, window=20, normalise=False, bAbsiscaCorrect=True):
    assert window < len(_list), "the window size for the moving average is too large"
    assert window >= 1, "the window size is too small (%s < 1)" % window

    if window == 1: # just return the original array
        return numpy.arange(0, len(_list)), _list

    return uniform_filter1d(_list, size=window, mode='reflect')[window//2:-window//2] # emulates previous version;
'''

def moving_average(a, window):
    ret = numpy.cumsum(a, dtype=float)
    ret[window:] = ret[window:] - ret[:-window]
    return ret[window:] / window

def movingAverage(_list, window=20, normalise=False, bAbsiscaCorrect=True):
    assert window < len(_list), "the window size for the moving average is too large"
    assert window >= 1, "the window size is too small (%s < 1)" % window

    if window == 1: # just return the original array
        return numpy.arange(0, len(_list)), _list

    if bAbsiscaCorrect:
        half_window_left = int(math.ceil(window / 2.0)) # correct for floating error division.
        half_window_right = int(math.floor(window / 2.0))
        x = numpy.arange(half_window_left, len(_list)-half_window_right)
    else:
        x = numpy.arange(0, len(_list)-window)

    y = []

    for n in range(half_window_left, len(_list)-half_window_right):
        score = sum(
            _list[i]
            for i in range(n - half_window_left, n + half_window_right, 1)
        )

        if normalise:
            y.append(float(score) / window)
        else:
            y.append(score)

    return (x, y)

def rc(seq):
    """
    get the reverse complemnt of seq
    """
    compdict = {'A': 'T',
                'C': 'G',
                'G': 'C',
                'T': 'A',
                'N': 'N',
                "a": "t",
                "c": "g",
                "g": "c",
                "t": "a",
                "n": "n"
                }
    tseq = [compdict[i] for i in seq] # new list
    tseq.reverse()

    return "".join(tseq)

def rc_expanded(seq):
    """
    get the reverse complemnt of seq
    this one works for the expanded alphabet;
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    [ACT] = h, [ACG] = v, [AGT] = d, [CGT] = b
    degenerate character N=[ATCG]
    """
    compdict = {'a': 't',
                'c': 'g',
                'g': 'c',
                't': 'a',
                'r': 'y',
                'y': 'r',
                'k': 'm',
                'm': 'k',
                's': 's', # palindromic
                'w': 'w', # palindromic
                'n': 'n', # all;
                "h": "d", # threes
                "v": "b",
                "d": "h",
                "b": "v"
                }
    tseq = [compdict[i] for i in seq] # new list
    tseq.reverse()

    return "".join(tseq)

def convertFASTAtoDict(filename, gzip_input=False):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    returns a list of the form [{name, seq}, ... {name, seq}]
    """
    assert os.path.isfile(filename), f"filename {filename} not found"

    openfile = gzip.open(filename, "rt") if gzip_input else open(filename, "rt")
    result = []
    for line in openfile:
        line = line.strip()
        if not line:
            continue

        if line[0] != ">": # not a FASTA block, so add this line to the sequence
            entry["seq"] += line.strip().replace('\r', '').replace("\n", "")

        if line[0] == ">": # fasta start block
            entry = {"seq": "", "name": "empty"} # make a new node
            # start recording
            entry["name"] = line.replace(">", "").strip()
            # add the old Node to the list
            if entry["name"] != "empty":
                # convert the list to a tuple
                result.append(entry) # You have to think about this one, but it works as it appends a view!
                # And so will not miss the last item!
    return result

def qcollide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return False
    if Aleft > Bright:
        return False

    if Aleft <= Bright and Aright >= Bright:
        return True  # Bright point is within A, collision

    if Aright >= Bleft and Aleft <= Bleft:
        return True # Bleft point is within A, collision.

    if Bleft <= Aright and Bright >= Aright:
        return True # Aright point is within B, collision

    if Bright >= Aleft and Bleft <= Aleft:
        return True # Aleft point is within B, collision.

    #print "unhandled!"
    return False

def FASTAToLIST(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    try:
        openfile = open(filename, "rb")
    except IOError:
        print("Error opening File:", filename)
        sys.exit()

    record = ""
    elementList = []
    entry = Node("empty")
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.replace("\n", "").replace("\r", "")

        if line[:1] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                n = entry.seq.lower()
                elementList.append(n) # and stick it on the big list
            entry = Node(line) # make a new node

    # all done;
    return elementList

# This code comes from http://www.johndcook.com/standard_deviation.html
# the point of all this complexity is to allow incremental computation of mean and std in
# a numerically stable way.
class accumulate_mean:
    def __init__(self):
        self.m_n = self.m_oldM = self.m_newM = 0

    def fields(self):
        return ["mean", "stdev"]

    def add(self, x, binWidth):
        x = float(x)/binWidth
        self.m_n += 1
        if self.m_n == 1:
            self.m_oldM = self.m_newM = x
            self.m_oldS = 0.0
        else:
            self.m_newM = self.m_oldM + (x - self.m_oldM) / self.m_n
            self.m_newS = self.m_oldS + (x - self.m_oldM) * (x - self.m_newM)
            # set up for next iteration
            self.m_oldM = self.m_newM
            self.m_oldS = self.m_newS

    def finishPosition(self):
        pass # nothing to do for mean

    def finalize(self, count):
        mean = float(self.m_newM) * self.m_n / count
        stdev = math.sqrt(self.m_newS/ (self.m_n - 1)) if self.m_n > 1 else 0.0
        self.val = {"mean":mean, "stdev":stdev}

    def value(self):
        return self.val

def transpose(list):
    """
    FOR DEPRECATION, use numpy.T
    A transpose command, rotates a matrix or equivalent by 90
    """
    try:
        rows = len(list[0])
    except (TypeError, IndexError):
        rows = 1 # probably.
    cols = len(list)

    newl = [[0 for x in range(cols)] for _ in range(rows)]
    for r in range(rows):
        for c in range(cols):
            newl[r][c] = list[c][r]
    return newl

def isPalindromic(seq):
    """
    is a sequence palindromic?
    returns True or False
    """
    return rc_expanded(seq.lower()) == seq.lower()

def isPalindrome(seq):
    """
    is a sequence a palindrome?
    """
    return rc_expanded(seq.lower()) == seq.lower()

def bin_data(array_like, bin_size):
    """
    This is an old alias, please use bin_sum_data or bin_mean_data
    """
    return [sum(array_like[i:i+bin_size]) for i in range(0, len(array_like), bin_size)]

def bin_sum_data(array_like, bin_size):
    return [sum(array_like[i:i+bin_size]) for i in range(0, len(array_like), bin_size)]

def bin_mean_data(array_like, bin_size):
    return [(sum(array_like[i:i+bin_size]) / float(len(array_like[i:i+bin_size]))) for i in range(0, len(array_like), bin_size)]

def scale_data(array_like, range=(0, 100)):
    """
    rescale the data within range
    """
    assert range[0] == 0, "sorry, only ranges 0..n are supported at present"

    scaled = numpy.zeros(range[1])
    s = len(array_like)/ float(range[1])

    for fi, f in enumerate(numpy.arange(range[0], len(array_like)-1, s)):
        val = array_like[int(math.floor(f)):int(math.floor(f+s))]
        if len(val) >= 1: # when s < 1.0 sometimes the recovered array will be empty
            scaled[fi] += numpy.average(val)
    return scaled

def kde(val_list, range=(0,1), covariance=0.02, bins=20):
    """
    kernal denstity estimation of an array

    covariance not working?
    """
    a = numpy.linspace(range[0], range[1], bins)

    def covariance_factor(self):
        return covariance

    kde = stats.gaussian_kde(val_list)
    #setattr(kde, 'covariance_factor', covariance_factor.__get__(kde, type(kde))) # Hack gaussian_kde()
    #kde._compute_covariance()

    kk = kde.evaluate(a) # resacle to get in integer range.

    return numpy.array(kk, dtype=numpy.float64)

def fold_change(c1, c2, log=2, pad=1e-6):
    """
    Calculate the fold-change of two values

    by default returns the log fold-change
    """
    try:
        c1v = c1 + pad
        c2v = c2 + pad
        if c2v > c1v:
            if log:
                return(math.log((c2v/c1v), log))
            else:
                return((c2v/c1v))
        else:
            if log:
                return(-math.log(c1v/c2v, log))
            else:
                return(-(c1v/c2v))

    except (OverflowError, ZeroDivisionError, ValueError):
        if c2v > c1v:
            config.log.error("(%.2f/%.2f) failed" % (c2v, c1v))
        else:
            config.log.error("(%.2f/%.2f) failed" % (c1v, c2v))
        raise Exception("fold_change() encountered an error, possibly the pad value is too small, or you are trying to apply fold-change to log transformed data")

def rgba_to_hex(rgba_color):
    return ('#{r:02x}{g:02x}{b:02x}'.format(r=int(rgba_color[0]*255),g=int(rgba_color[1]*255),b=int(rgba_color[2]*255)))

def hex_to_rgb(hex_str):
    return tuple((int(hex_str.lstrip('#')[i:i+2], 16)/255) for i in (0, 2, 4))

def qdeepcopy(anobject):
    return pickle.loads(pickle.dumps(anobject, -1))

def fastq(filename, gziped=False):
    """
    generator object to parse a fastQ file

    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((

    """
    oh = gzip.open(filename, "rt") if gziped else open(filename, "rU")
    name = "."
    while name != "":
        name = oh.readline().strip()
        seq = oh.readline().strip()
        strand = oh.readline().strip()
        qual = oh.readline().strip()

        yield {"name": name, "strand": strand, "seq": seq, "qual": qual}
    return

def fastqPE(filename1, filename2, gziped=True):
    """
    generator object to parse fastQ PE files
    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((
    """
    if gziped:
        oh1 = gzip.open(filename1, "rt")
        oh2 = gzip.open(filename2, "rt")
    else:
        oh1 = open(filename1, "rt")
        oh2 = open(filename2, "rt")

    name1 = "."
    while name1 != "":
        name1 = oh1.readline().strip()
        seq1 = oh1.readline().strip()
        strand1 = oh1.readline().strip()
        qual1 = oh1.readline().strip()

        name2 = oh2.readline().strip()
        seq2 = oh2.readline().strip()
        strand2 = oh2.readline().strip()
        qual2 = oh2.readline().strip()

        yield ({"name": name1, "strand": strand1, "seq": seq1, "qual": qual1},
            {"name": name2, "strand": strand2, "seq": seq2, "qual": qual2})
    return
