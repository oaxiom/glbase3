"""
Utilities

Various utilities to support the genome scanning scripts.

MAny of these predate glbase3, but are a little tricky to remove as I am not sure where
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

#from errors import AssertionError

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

    return(l)

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
        return(True)
    return(True)

def movingAverage(listIn, window=20, normalise=False, bAbsiscaCorrect=True):
    """
    actually a sliding window
    """
    assert window < len(listIn), "the window size for the moving average is too large"
    assert window >= 1, "the window size is too small (%s < 1)" % window

    if window == 1: # just return the original array
        return(numpy.arange(0, len(listIn)), listIn)

    if bAbsiscaCorrect:
        half_window_left = int(math.ceil(window / 2.0)) # correct for floating error division.
        half_window_right = int(math.floor(window / 2.0))
        x = numpy.arange(half_window_left, len(listIn)-half_window_right)
    else:
        x = numpy.arange(0, len(listIn)-window)

    y = []

    for n in range(half_window_left, len(listIn)-half_window_right):
        score = 0
        for i in range(n-half_window_left, n+half_window_right, 1):
            score += listIn[i]

        if normalise:
            y.append(float(score) / window)
        else:
            y.append(score)

    return(x, y)

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
    return(a)

def osc(last, type):
    """
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]
    """
    if type == "r":
        if last == "a": return("g")
        if last == "g": return("a")
        return("a")
    if type == "y":
        if last == "c": return("t")
        if last == "t": return("c")
        return("c")
    if type == "k":
        if last == "g": return("t")
        if last == "t": return("g")
        return("g")
    if type == "m":
        if last == "a": return("c")
        if last == "c": return("a")
        return("a")
    if type == "s":
        if last == "g": return("c")
        if last == "c": return("g")
        return("g")
    if type == "w":
        if last == "a": return("t")
        if last == "t": return("a")
        return("a")
    if type == "n":
        if last == "a": return("c")
        if last == "c": return("g")
        if last == "g": return("t")
        if last == "t": return("a")
        return("a")
    return(type)

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

    return("".join(tseq))

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

    return("".join(tseq))

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

    return  lib

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

    return  lib

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

    return  lib

def convertCSVtoFASTA(csvfilename, outfile, sequenceCol, fastaNameCol = None):
    """
    csvfilename = the csv filename to load without the path, assumed to be in environment.mouseGenomePath
    outfile = outfile name, also in environment.mouseGenomePath
    sequenceCol (integer) = the number of the column (starts from 0)
    fastaNameCol (integer) (optional) = the column to use as a FASTA name, otherwise will default to 0..n
    ignores a row if the first column begins with "#"
    """
    ofh = open(os.path.join(sys.path[0], csvfilename), "rb")
    sfh = open(os.path.join(sys.path[0], outfile), "wb")

    csvreader = csv.reader(ofh)

    fastanameseries = 0

    for n in csvreader:
        t = n[0]
        if t[0] != "#":
            seq = n[sequenceCol]
            if fastaNameCol:
                fastname = n[fastaNameCol]
                sfh.write('>'+fastname+'\r\n')
            else: # no fasta name so use the series
                sfh.write('>'+str(fastanameseries)+'\r\n')
                fastanameseries += 1
            sfh.write(seq+'\r\n')

    ofh.close()
    sfh.close()

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
    return(result)

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
    return(a)

def removeDuplicatesFromCSV(path, csvfile, outfile, column_no = 3, bKeepEmptyCols=False):
    """
    delete duplicates based on the column no
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []

    for line in reader:
        if line[column_no] in ulist:
            # don't write this enty,
            print("Duplicate:", line[column_no])
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
                writer.writerow(line) # if I tab this in - don't keep empty rows
            else: # col is empty;
                #print "Duplicate: <empty>"
                if bKeepEmptyCols: writer.writerow(line)
    inf.close()
    outf.close()

def collide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft == Bleft: return(1) # I have to cheat here otherwise it returns 0 which will evaluate as False;
    if Aleft == Bright: return(1)
    if Aright == Bleft: return(1)
    if Aright == Bright: return(1)

    if Aleft <= Bright and Aright >= Bright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Bright point is within A, thus collision

    if Aright >= Bleft and Aleft <= Bleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Bleft point is within A, thus collision.

    if Bleft <= Aright and Bright >= Aright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Aright point is within B, thus collision

    if Bright >= Aleft and Bleft <= Aleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Aleft point is within B, thus collision.

    #print "unhandled!"
    return(False)

def qcollide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft <= Bright and Aright >= Bright:
        return(True) # Bright point is within A, collision

    if Aright >= Bleft and Aleft <= Bleft:
        return(True) # Bleft point is within A, collision.

    if Bleft <= Aright and Bright >= Aright:
        return(True) # Aright point is within B, collision

    if Bright >= Aleft and Bleft <= Aleft:
        return(True) # Aleft point is within B, collision.

    #print "unhandled!"
    return(False)

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
    return(newlist)

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
    return(newlist)

def removeDuplicatesFromCSV_2Cols(path, csvfile, outfile, column_no1 = 0, column_no2 = 1):
    """
    delete duplicates based on the column no1 and 2
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []

    for line in reader:
        if line[column_no1]+line[column_no2] in ulist:
            # don't write this enty,
            print("duplicate:", line[column_no1]+line[column_no2])
        else:
            # add to ulist and write to file;
            ulist.append(line[column_no1]+line[column_no2])
            writer.writerow(line)
    inf.close()
    outf.close()

def keepRowOnlyIfColXHasValue(path, _in, _out, _colX):
    """
    what it says
    """
    print("keepRowOnlyIfColXHasValue(",path, _in, _out, _colX,")")
    inf = open(os.path.join(path, _in), "rb")
    outf = open(os.path.join(path, _out), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)

    for line in reader:
        if len(line) > _colX:
            if line[_colX]:
                writer.writerow(line)
    inf.close()
    outf.close()

def FASTAToCSV(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    #try:
    openfile = open(filename, "rb")
    savefile = open(filename+'_out.csv', "wb")
    #except IOError:
    #    print "Error opening File: %s" % filename
    #    sys.exit()

    writer = csv.writer(savefile)

    record = ""
    entry = Node("empty")
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.replace("\r", "").replace("\n", "")

        if line[:1] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                writer.writerow([entry.name, "", "", entry.seq])
                del entry
            entry = Node(line) # make a new node

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

def repeat_mask(seq):
    return(seq.replace("a", "n").replace("c", "n").replace("g", "n").replace("t", "n"))

def loadTSVAsLIST(file):
    oh = open(file, "rU")
    reader = csv.reader(oh, dialect=csv.excel_tab)
    newl = []
    for line in reader:
        newl.append(line)
    return(newl)

def renameDuplicatesFromCSV(path, csvfile, outfile, column_no = 3, bKeepEmptyCols=False):
    """
    append _1 .. _n to duplicates based on the column no
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []
    nFound = {}

    for line in reader:
        if line[column_no] in ulist:
            # don't write this enty,
            print("Duplicate:", line[column_no])
            if bKeepEmptyCols:
                if column_no == 0:
                    writer.writerow(["%s_%s" % (line[column_no], nFound[line[column_no]])] + line[column_no+1:])
                elif column_no == len(line):
                    writer.writerow(line[:column_no] + ["%s_1" % line[column_no]])
                else:
                    writer.writerow(line[:column_no] + ["%s_1" % line[column_no]] + line[column_no+1:])
                nFound[line[column_no]] += 1
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
                # add a key to the nTimeFound dict;
                nFound[line[column_no]] = 1
                writer.writerow(line) # if I tab this in - don't keep empty rows

            else: # col is empty;
                #print "Duplicate: <empty>"
                if bKeepEmptyCols: writer.writerow(line)
    inf.close()
    outf.close()

# This code comes from http://www.johndcook.com/standard_deviation.html
# the point of all this complexity is to allow incremental computation of mean and std in
# a numerically stable way.
# Taken from ACT
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

"""
    It's hard to believe, but these custom routines below are 10x as
    fast as their numpy equivalents...
    Numpy is a dog.
"""
def mean(intList):
    try:
        return sum(intList) / float(len(intList))
    except TypeError:
        return(intList) # intList is probably a single int -

def std(intList):
    return(math.sqrt(mean([(abs(x - mean(intList) ) ** 2) for x in intList])))

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
    return(newl)

def isPalindromic(seq):
    """
    is a sequence palindromic?
    returns True or False
    """
    if rc_expanded(seq.lower()) == seq.lower():
        return(True)
    return(False)

def isPalindrome(seq):
    """
    is a sequence a palindrome?
    """
    if rc_expanded(seq.lower()) == seq.lower():
        return(True)
    return(False)

def bin_data(array_like, bin_size):
    """
    This is an old alias, please use bin_sum_data or bin_mean_data
    """
    return([sum(array_like[i:i+bin_size]) for i in range(0, len(array_like), bin_size)])

def bin_sum_data(array_like, bin_size):
    return([sum(array_like[i:i+bin_size]) for i in range(0, len(array_like), bin_size)])

def bin_mean_data(array_like, bin_size):
    return([(sum(array_like[i:i+bin_size]) / float(len(array_like[i:i+bin_size]))) for i in range(0, len(array_like), bin_size)])

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
    return(scaled)

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

    return(numpy.array(kk, dtype=numpy.float64))

# I think this is non-functional below:
class pp:
    reachability_distance = None
    done = False

def _core_distance(pt, epsilon, MinPts):
    return()

def _optics_update(N, p, Seeds, eps, Minpts):

    coredist = core_distance(p, eps, MinPts)

    for o in N:
        if (o is not processed):
            new_reach_dist = max(coredist, dist(p,o))

            if not o.reachability_distance: # o is not in Seeds
                o.reachability_distance = new_reach_dist
                Seeds.insert(o, new_reach_dist)
            else: # o in Seeds, check for improvement
                if (new_reach_dist < o.reachability_distance):
                    o.reachability_distance = new_reach_dist
                    Seeds.move_up(o, new_reach_dist)

def OPTICS(DB, eps, MinPts):

    #for p in DB:
    #   p.reachability_distance = None
    ordered = []

    for p in DB:
        if not p.done:
            N = getNeighbors(p, eps)
            p.done = True
            ordered.append(p)
            Seeds = heapq([])

            if not core_distance(p, eps, Minpts):
                _optics_update(N, p, Seeds, eps, Minpts)
                for q in Seeds:
                    Np = getNeighbors(q, eps)
                    q.done = True
                    ordered.append(q)
                    if not core_distance(q, eps, Minpts):
                        _optics_update(Np, q, Seeds, eps, Minpts)

# The below code was taken and modified from sklearn.
# It's mostly here for reference, I just use the sklearn kit.

"""
Sklearn was released under the BSD license, which is repeated here:

New BSD License

Copyright (c) 20072013 The scikit-learn developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the Scikit-learn Developers  nor the names of
     its contributors may be used to endorse or promote products
     derived from this software without specific prior written
     permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

"""
def euclidean_distances(X, Y=None):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    For efficiency reasons, the euclidean distance between a pair of row
    vector x and y is computed as::

        dist(x, y) = sqrt(dot(x, x) - 2 * dot(x, y) + dot(y, y))

    This formulation has two main advantages. First, it is computationally
    efficient when dealing with sparse data. Second, if x varies but y
    remains unchanged, then the right-most dot-product 'dot(y, y)' can be
    pre-computed.

    APH: Removed support for Sparse matrices to simplify things a bit and squared returns

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = [n_samples_1, n_features]

    Y : {array-like, sparse matrix}, shape = [n_samples_2, n_features]

    Y_norm_squared : array-like, shape = [n_samples_2], optional
        Pre-computed dot-products of vectors in Y (e.g.,
        ''(Y**2).sum(axis=1)'')

    Returns
    -------
    distances : {array, sparse matrix}, shape = [n_samples_1, n_samples_2]

    Examples
    --------
    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> X = [[0, 1], [1, 1]]
    >>> # distance between rows of X
    >>> euclidean_distances(X, X)
    array([[ 0.,  1.],
           [ 1.,  0.]])
    >>> # get distance to origin
    >>> euclidean_distances(X, [[0, 0]])
    array([[ 1.        ],
           [ 1.41421356]])
    """
    # should not need X_norm_squared because if you could precompute that as
    # well as Y, then you should just pre-compute the output and not even
    # call this function.
    #X, Y = check_pairwise_arrays(X, Y)
    if Y is None:
        X = Y = numpy.asarray(numpy.atleast_2d(X), dtype=numpy.float, order=None)
    else: # make sure float;
        X = numpy.asarray(numpy.atleast_2d(X), dtype=numpy.float, order=None)
        Y = numpy.asarray(numpy.atleast_2d(Y), dtype=numpy.float, order=None)

    XX = numpy.sum(X * X, axis=1)[:, numpy.newaxis]

    if X is Y:  # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    else:
        YY = numpy.sum(Y ** 2.0, axis=1)[numpy.newaxis, :]

    distances = numpy.dot(X, Y.T) #safe_sparse_dot(X, Y.T, dense_output=True)
    distances *= -2
    distances += XX
    distances += YY
    numpy.maximum(distances, 0, distances)

    if X is Y:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        distances.flat[::distances.shape[0] + 1] = 0.0

    return distances

def dbscan(X, eps=0.5, min_samples=5, metric='euclidean', random_state=None):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples] or [n_samples, n_features]
        Array of distances between samples, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.
    eps: float, optional
        The maximum distance between two samples for them to be considered
        as in the same neighborhood.
    min_samples: int, optional
        The number of samples in a neighborhood for a point to be considered
        as a core point.
    metric:
        Only euclidean in this variant.
    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

    Returns
    -------
    core_samples: array [n_core_samples]
        Indices of core samples.

    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.

    """
    X = numpy.asarray(X)
    n = X.shape[0]

    # If index order not given, create random order.
    random_state = numpy.random.mtrand._rand
    index_order = numpy.arange(n)
    random_state.shuffle(index_order)
    D = pdist(X, metric="euclidean")
    print(D)
    D = euclidean_distances(X)
    print(D)
    print(D.shape)
    print(numpy.sum(D, axis=0))
    print(numpy.sum(D, axis=1))
    print(numpy.max(D))
    print(numpy.min(D))
    print(numpy.histogram(D))

    # Calculate neighborhood for all samples. This leaves the original point
    # in, which needs to be considered later (i.e. point i is the
    # neighborhood of point i. While True, its useless information)
    neighborhoods = [numpy.where(x <= eps)[0] for x in D]
    #print neighborhoods

    # Initially, all samples are noise.
    labels = -numpy.ones(n)

    # A list of all core samples found.
    core_samples = []

    # label_num is the label given to the new cluster
    label_num = 0

    # Look at all samples and determine if they are core.
    # If they are then build a new cluster from them.
    for index in index_order:
        if labels[index] != -1 or len(neighborhoods[index]) < min_samples:
            # This point is already classified, or not enough for a core point.
            continue
        core_samples.append(index)
        labels[index] = label_num

        # candidates for new core samples in the cluster.
        candidates = [index]
        while len(candidates) > 0:
            new_candidates = []

            # A candidate is a core point in the current cluster that has
            # not yet been used to expand the current cluster.
            for c in candidates:
                noise = numpy.where(labels[neighborhoods[c]] == -1)[0]
                noise = neighborhoods[c][noise]
                labels[noise] = label_num
                for neighbor in noise:

                    # check if its a core point as well
                    if len(neighborhoods[neighbor]) >= min_samples:
                        # is new core point
                        new_candidates.append(neighbor)
                        core_samples.append(neighbor)

            # Update candidates for next round of cluster expansion.
            candidates = new_candidates

        # Current cluster finished.
        # Next core point found will start a new cluster.
        label_num += 1
    return core_samples, labels

# ---- End of sklearn code.

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
        raise Exception("__fold_change() encountered an error, possibly the pad value is too small, or you are trying to apply fold-change to log transformed data")

def rgba_to_hex(rgba_color):
    return ('#{r:02x}{g:02x}{b:02x}'.format(r=int(rgba_color[0]*255),g=int(rgba_color[1]*255),b=int(rgba_color[2]*255)))

def hex_to_rgb(hex_str):
    return tuple((int(hex_str.lstrip('#')[i:i+2], 16)/255) for i in (0, 2, 4))

def qdeepcopy(anobject):
    # You should wrap me in a try: except:
    return(pickle.loads(pickle.dumps(anobject, -1)))

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

def fastqPE(filename1, filename2):
    """
    generator object to parse fastQ PE files

    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((

    """
    oh1 = open(filename1, "rU")
    oh2 = open(filename2, "rU")

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
