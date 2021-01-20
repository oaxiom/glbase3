"""
Utilities

Various utilities to support the genome scanning scripts.

Many of these predate glbase3, but are a little tricky to remove as I am not sure where
they are used (if at all).

So excuse the terrible code in places. I will deprecate occasional functions from this.

R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
degenerate character N=[ATCG]
3-fold degenerate motifs re not used like the Lander paper.

"""

import sys, os, numpy, string, csv, random, math, pickle, gzip
import scipy.stats as stats
from scipy.spatial.distance import pdist

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
        if nm[n] == "r" or nm[n] == "y" or nm[n] == "k" or nm[n] == "m" or nm[n] == "s" or nm[n] == "w": # AG
            newlist.append((2, n, nm[n])) # append a triple, with the number of flips and its location
        elif nm[n] == "n":
            newlist.append((4, n, nm[n]))
        else:
            newlist.append((1, n, nm[n]))

    if len(newlist) == 0 : return([motif]) # no degenerate elements in the motif

    l = []
    #print newlist
    iti(newlist, 0, None, l)

    return l

def iti(_fm, _fmcpos, _cl, _list): # my iterator
    """
    This is possibly the best piece of code I've ever made!
    Mainly because it's almost entirely unitelligable....
    It's some sort of iterator to to generate
    N-mers.
    I think this has been replaced by regexs now.
    fm = a triple of the form (number to flip, location, seq_at_pos)
    """
    # do some set up
    if not _cl: # also okay if _cl == False; be careful with these, as may be False, but not None
        _cl = []
        for n in range(len(_fm)):
            _cl.append("")

    # the main iterator;
    if _fmcpos >= len(_fm):
        return(False) # reached end of list, signal to stick it back on the end;
    else:
        n = _fm[_fmcpos]
        #print n
        for x in range(n[0]):
            _cl[n[1]] = osc(_cl[n[1]], n[2])
            if not iti(_fm, _fmcpos+1, _cl, _list): # each time we iterate at the end of the motif add it to the list;
                # convert the list back to a string
                _copy = string.join(_cl, '')
                _list.append(_copy)
        return True
    return True

def movingAverage(_list, window=20, normalise=False, bAbsiscaCorrect=True):
    """
    actually a sliding window
    """
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
        score = 0
        for i in range(n-half_window_left, n+half_window_right, 1):
            score += _list[i]

        if normalise:
            y.append(float(score) / window)
        else:
            y.append(score)

    return (x, y)

def cumulative_line(listIn, percent=True):
    """
    convert the array to a cumulative array

    The array is expected to be some sort of list of numbers.
    This will return the array added in a cumulative manner from 0 .. 100 %
    The returned list will always be floats.

    (If percent is left as True - the default)
    """
    s = sum(listIn)
    m = max(listIn)

    c = 0
    rc = 0
    n = []
    for i in listIn:
        c += i
        if c >= 50:
            rc -= i
        else:
            rc += i
        n.append(rc)

    if percent:
        a = (numpy.array(n) / float(s)) * 100
    else:
        a = n
    return a

def osc(last, type):
    """
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]
    """
    if type == "r":
        if last == "a": return("g")
        if last == "g": return("a")
        return "a"
    if type == "y":
        if last == "c": return("t")
        if last == "t": return("c")
        return"c"
    if type == "k":
        if last == "g": return("t")
        if last == "t": return("g")
        return "g"
    if type == "m":
        if last == "a": return("c")
        if last == "c": return("a")
        return "a"
    if type == "s":
        if last == "g": return("c")
        if last == "c": return("g")
        return "g"
    if type == "w":
        if last == "a": return("t")
        if last == "t": return("a")
        return "a"
    if type == "n":
        if last == "a": return("c")
        if last == "c": return("g")
        if last == "g": return("t")
        if last == "t": return("a")
        return "a"
    return type

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

def expandElement(_element, _preserve=True): #
    """
    _element = (string) to expand

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for left in range(4): # base iterator;
        for right in range(4):
            lib.append(dir[left]+_element+dir[right])

    return lib

def expandElementRightOnly(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in range(4):
        lib.append(_element+dir[right])

    return lib

def expandElementRightOnly_degenerate(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'
    will use a degenerate code;
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t",
           4 : "r",
           5 : "y",
           6 : "k",
           7: "m",
           8: "s",
           9: "w",
           10:"n"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in range(11):
        lib.append(_element+dir[right])

    return lib

def expandElementRightOnly_degenerate_n(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'
    will use a degenerate code;
    this one n only
    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t",
           4:"n"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in range(5):
        lib.append(_element+dir[right])

    return lib

def convertFASTAtoCSV(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    assert os.path.exists(filename), "filename %s not found" % filename

    try:
        openfile = open(filename, "rb")
        savefile = open(filename+'_out.csv', "wb")
    except IOError:
        print("Error opening File")
        sys.exit()

    writer = csv.writer(savefile)

    record = ""
    entry = Node("empty")
    for line in openfile:
        if line[0] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.strip().replace('\r\n', '') # strip out the new line WINDOWS specific!

        if line[0] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                writer.writerow([entry.name, "", "", "", entry.seq])
                del entry
            entry = Node(line) # make a new node

def convertFASTAtoDict(filename, gzip_input=False):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    returns a list of the form [{name, seq}, ... {name, seq}]
    """
    assert os.path.exists(filename), "filename %s not found" % filename

    if gzip_input:
        openfile = gzip.open(filename, "rt")
    else:
        openfile = open(filename, "rt")

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

def scanNumberOfBasePairs(fastafilehandle):
    """
    pass me a fasta file handle (already open()'d) and this will return the raw count
    or some other iterable object;
    """
    dict = {"a" : 0,
            "A" : 0,
            "c" : 1,
            "C" : 1,
            "g" : 2,
            "G" : 2,
            "t" : 3,
            "T" : 3}

    a = []
    a.append(0)
    a.append(0)
    a.append(0)
    a.append(0)

    for line in fastafilehandle:
        if line[0] != ">": # there is a more elegant way to do this...
            lcline = line.lower()
            a[dict["a"]] += lcline.count("a")
            a[dict["c"]] += lcline.count("c")
            a[dict["g"]] += lcline.count("g")
            a[dict["t"]] += lcline.count("t")
    return a

def collide(Aleft, Aright, Bleft, Bright):
    # quickest rejections first;
    if Aright < Bleft:
        return False
    if Aleft > Bright:
        return False

    if Aleft == Bleft: return 1 # I have to cheat here otherwise it returns 0 which will evaluate as False;
    if Aleft == Bright: return 1
    if Aright == Bleft: return 1
    if Aright == Bright: return 1

    if Aleft <= Bright and Aright >= Bright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return closest # Bright point is within A, thus collision

    if Aright >= Bleft and Aleft <= Bleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return closest # Bleft point is within A, thus collision.

    if Bleft <= Aright and Bright >= Aright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return closest # Aright point is within B, thus collision

    if Bright >= Aleft and Bleft <= Aleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return closest # Aleft point is within B, thus collision.

    #print "unhandled!"
    return False

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

def removeDuplicatesFromListOfDicts(list_of_dicts, key):
    """
    remove duplicates from a list of dicts based on key,
    returns the list in the format it arrived.
    only checks key. Does not check anything else.
    """
    ulist = []
    newlist = []
    dupecount = 0

    for line in list_of_dicts:
        if line[key] in ulist:
            # don't write this enty,
            dupecount +=1
        else:
            # add to ulist and write to file;
            if line[key]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[key])
            newlist.append(line)
    #print ">> Duplicates Found:", dupecount
    return newlist

def removeDuplicatesFrom2DList(_list, column_no = 3):
    """
    delete duplicates based on the column no
    returns a list
    deletes dupes in a 2D-csv like list;
    """
    ulist = []
    newlist = []
    dupecount = 0

    for line in _list:
        if line[column_no] in ulist:
            # don't write this enty,
            dupecount +=1
            print("Duplicate:%s" % (line[column_no]))
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
            newlist.append(line)
    #print ">> Duplicates Found:", dupecount
    return newlist

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
    return (elementList)

def loadTSVAsLIST(file):
    oh = open(file, "rU")
    reader = csv.reader(oh, dialect=csv.excel_tab)
    newl = []
    for line in reader:
        newl.append(line)
    return newl

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
        if self.m_n > 1:
            stdev = math.sqrt(self.m_newS/ (self.m_n - 1))
        else:
            stdev = 0.0
        self.val = {"mean":mean, "stdev":stdev}

    def value(self):
        return self.val

def mean(intList):
    try:
        return sum(intList) / float(len(intList))
    except TypeError:
        return intList # intList is probably a single int -

def std(intList):
    return math.sqrt(mean([(abs(x - mean(intList) ) ** 2) for x in intList]))

def transpose(list):
    """
    a transpose command, rotates a mtrix or equivalent by 90

    more like a transpose from R than anything else.
    """
    newl = []
    try:
        rows = len(list[0])
    except:
        rows = 1 # probably.
    cols = len(list)

    for row in range(rows):
        newl.append([0 for x in range(cols)])
    for r in range(rows):
        for c in range(cols):
            newl[r][c] = list[c][r]
    return newl

def isPalindromic(seq):
    """
    is a sequence palindromic?
    returns True or False
    """
    if rc_expanded(seq.lower()) == seq.lower():
        return True
    return False

def isPalindrome(seq):
    """
    is a sequence a palindrome?
    """
    if rc_expanded(seq.lower()) == seq.lower():
        return True
    return False

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
    # You should wrap me in a try: except:
    return pickle.loads(pickle.dumps(anobject, -1))

def fastq(filename, gziped=False):
    """
    generator object to parse a fastQ file

    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((

    """
    if gzip:
        oh = gzip.open(filename, "rt")
    else:
        oh = open(filename, "rU")

    name = "dummy"
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

    name1 = "dummy"
    while name1 != "":
        name1 = oh1.readline().strip()
        seq1 = oh1.readline().strip()
        strand1 = oh1.readline().strip()
        qual1 = oh1.readline().strip()

        name2 = oh2.readline().strip()
        seq2 = oh2.readline().strip()
        strand2 = oh2.readline().strip()
        qual2 = oh2.readline().strip()

        res = ({"name": name1, "strand": strand1, "seq": seq1, "qual": qual1},
            {"name": name2, "strand": strand2, "seq": seq2, "qual": qual2})
        yield res
    return
