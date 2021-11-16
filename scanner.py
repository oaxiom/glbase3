"""

**scanner.py**

motif scanner, part of glbase

Interacts strongly with fastalist
"""


"""
scanPeriodicity

scan's a list of elements and returns their positional bias.

todo:
-----
o use degenerate elements

"""

import sys, os, re, string
from numpy import zeros

# get glbase
from . import config
from .draw import draw
from .history import historyContainer
from .data import *
from .genelist import genelist

class scanner(genelist.genelist):
    """
    DNA motif scanner back scanners and other tools.
    """
    def __init__(self):
        genelist.__init__(self)

    def scanPeriodicityFASTA(degenmotif, lib, path, fastafile):
        """
        degenmotif = the original degenerate motif used to extract
        lib = list of elements of any length (within reason)
        fastafile = er... fasta file

        fasta file MUST be of the form:
        >ID
        accacactttcactasequence\r\n
        >nextID

        returns a list of (object) rObj ?
        o would really like to do a reg expression for degenerate motifs here;
        """
        print("Info: Started")
        oh = file(os.path.join(path, fastafile), "rU") # Don't cache the FASTA so can perform arbitrarily massive searches.

        linelength = None
        scanpos = 0

        regex = re.compile(string.join([regex_dict[bp] for bp in degenmotif], "")) # get regex of search term...

        n = 0
        m = 0
        f = 0

        for line in oh:
            if line[0] != ">":
                if not linelength: # I don't know the line length yet
                    linelength = len(line) # now I know it so make an array ->
                    a = zeros(linelength) # assume all subsequent line lengths are the same...

                line = line.strip("\n").strip("\r").lower()

                for rc, sequence in enumerate([line, rerverse_complement(line)]):
                    found = regex.finditer(sequence)
                    if found:
                        for match in found:
                            #temp_finds.append({"seq": sequence[match.start():match.end()], "pos": "%s:%s:%s:%s" % (n, rc, match.start(), match.end() ) })
                            a[match.start()] += 1
                            f += 1
                found_list = []

                n += 1
                if n > 1000:
                    m += 1
                    print("Done: %s000, Found: %s" % (m, f))
                    n = 0

        for w in [10, 50, 100, 500, 1000, 2000]:
            plotGraph(a, degenmotif, linelength, path, fastafile, w)

        oh.close()

    def plotGraph(a, degenmotif, linelength, path, fastafile, window):
        print("Info: Plotting Graph")
        x, ma = movingAverage(a, window)
        plot(x, ma)
        title(degenmotif+' periodicity')
        savefig(os.path.join(path, "%s_m%s_w%s_per.png" % (fastafile.split(".")[0], degenmotif, window)))
        cla()


class matrixmotif:
    def __init__(self, matrix, threshold):
        self._matrix = matrix
        self._threshold = threshold
        self._length = len(self._matrix)

    def find_matches(self, seq):
        """
            find the matches in a particular sequence.
            returns a list v
            containing [strand, position, score, and sequence]
        """
        if seq.upper().count("N") > 0:
            return(False) # reject if contains an N

        con_setpos = {"A" : 0, "C": 1, "G":2, "T": 3}
        result = [] # new list

        up_seq = seq.upper() # make sure the sequence is all CAPS
        rc_seq = self.rc(up_seq) # because I need to search the rc sequence too.

        for strand, seq in enumerate([up_seq, rc_seq]):
            for pos in range(len(up_seq)):

                if (pos + self._length) > len(seq): # fail is seq < length
                    break

                totalscore = 0.0
                for i, pwm_atpos in enumerate(self._matrix._pwm): # iterate through the pwm.
                    seq_atpos = seq[i+pos] # get sequence at position
                    totalscore += pwm_atpos[con_setpos[seq_atpos]]

                score = (totalscore - self._matrix._minscore) / (self._matrix._maxscore - self._matrix._minscore) # get norm score; 0..1

                if score >= self._threshold: # keep it.
                    if strand == 0:
                        result.append(['+', pos, score, up_seq[pos:pos+self._length]]) # 1 position short?
                    elif strand == 1:
                        result.append(['-', len(rc_seq)-pos, score, rc_seq[pos:pos+self._length]])

        if result:
            return(result)
        else:
            return(None)

    def findAllScores(self, seq):
        """
        similar to the old find_matches, but this returns a list of
        threshold scores across the length of the sequence.
        seq must be a single string.
        """

        if seq.upper().count("N") > 0:
            return(False) # reject if contains an N

        pos_seq = seq.upper()
        neg_seq = self.rc(pos_seq)
        result = []
        hits = []

        con_setpos = {"A" : 0, "C": 1, "G": 2, "T": 3}

        for p_pos in range(self._length, len(pos_seq)-self._length, 1):
            n_pos = len(pos_seq) - p_pos # work out the corresponding -ve strand

            seq_both_strands = [pos_seq[p_pos:p_pos+self._length], neg_seq[n_pos-self._length:n_pos]]

            scores = []
            for s in seq_both_strands:
                # get the threshold score for this motif.
                totalscore = 0.0
                for i, pwm_atpos in enumerate(self._matrix._pwm): # iterate through the pwm.
                    seq_atpos = s[i] # get sequence at position
                    if seq_atpos in con_setpos:
                        totalscore += pwm_atpos[con_setpos[seq_atpos]]

                scores.append((totalscore - self._matrix._minscore) / (self._matrix._maxscore - self._matrix._minscore)) # get norm score; 0..1

            if max(scores) >= self._threshold: # keep it.
                # found motif:
                print("here:", seq_both_strands, p_pos)
                hits.append(p_pos)
            result.append(max(scores))
        return({"scores": result, "hits": hits})

    def EstimateErrorRate(self, seq, steps):
        """
        Estimates the number of matches for a particular data set for the thresholds ranging from 0 to 1 in 0.05 steps
        Slow, essentially loops through the entire search 30 times.
        """
        total = []
        oldthreshold = self._threshold
        threshold = 0.75
        self._threshold = 0.0
        for n in range(steps):
            self._threshold = threshold + (n * 0.007) #(0.1/float(steps-1)))
            r = self.find_matches(seq)
            if r:
                total.append(len(self.find_matches(seq)))

        self._threshold = oldthreshold # swap back old threshold;
        return(total)

    def rc(self, seq):
        """
        get the reverse complemnt of seq
        """
        compdict = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A',
                    'N': 'N'
                    }
        tseq = [compdict[i] for i in seq] # new list
        tseq.reverse()

        return("".join(tseq))

class matrix:
    def __init__(self, matrixname, pwm):
        """
            Class constructor
        """
        self._matrixname = matrixname
        self._pwm = pwm
        self._length = len(pwm)
        self._minscore = self.min_score()
        self._maxscore = self.max_score()
        self.ConvertPFMtoPWM() # because almost always matrices are in PFM format, not PWM format.

    def value(self, position):
        """
        Return a row of the list
        """
        return(self._pwm[position])

    def ConvertPFMtoPWM(self):
        """
        Converts a pfm to a pwm
        matrices must be in pwm format for accurate searching.
        via this algorithm:

        w = log2 ( ( f + sqrt(N) * p ) / ( N + sqrt(N) ) / p )

        where:
        w - is a weight for the current nucleotide we are calculating
        f - is a number of occurences of the current nucleotide in the current column (e.g., "1" for A in column 1, "8" for C etc)
        N - total number of observations, the sum of all nucleotides occurences in a column (13 in this example)
        p - [prior] [background] frequency of the current nucleotide; this one usually defaults to 0.25 (i.e. one nucleotide out of four)

        Shamelessly borrowed from:
        http://bogdan.org.ua/2006/09/11/position-frequency-matrix-to-position-weight-matrix-pfm2pwm.html
        """
        n = 0
        w = 0.0
        bign = 0 # number of occurences in a particular column.
        newpwm = []
        for n in self._pwm:
            newbase = []
            bign = int(n[0]) + int(n[1]) + int(n[2]) + int(n[3])
            for e in n:
                w = math.log( ( e + math.sqrt(bign) * 0.25) / (bign + math.sqrt(bign)) / 0.25, 2) # log 2, not log n
                newbase.append(w)

            newpwm.append(newbase)
        self._pwm = newpwm
        self.redo_minmax() # redo min/max scores

    def __len__(self):
        """
        Method to give a valid len(self) assignment.
        """
        return(len(self._pwm)) # just pass back the length

    def min_score(self):
        """
        Calculate the minimum score for a particular pmw.
        """
        n = 0
        self._minscore = 0
        for n in self._pwm:
            self._minscore += min(n)

    def max_score(self):
        """
        Calculate the maximum score for a particular pmw.
        """
        n = 0
        self._maxscore = 0
        for n in self._pwm:
            self._maxscore += max(n)

    def redo_minmax(self):
        """
        Redo the min max score (for example, after doing pfm -> pwm)
        """
        self.min_score()
        self.max_score()

    def getName(self):
        return(self.__str__())

    def __str__(self):
        return(self._matrixname)


def BatchConvert(fasta_file, matrix, threshold = 0.883, bKeepAll=False, errorf=False, path=None):
    """
    Only switch on Error if you really mean it! Very very slow!
    """
    # join together the base dirs;

    if not path:
        path = sys.path[0]

    writerfile = open(os.path.join(path, "Matches_%s_%s.csv" % (fasta_file.split("_")[0], matrix._matrixname)), "w")
    fastafile = open(os.path.join(path, "Matches_%s_%s.fa" % (fasta_file.split("_")[0], matrix._matrixname)), "w")
    if errorf: errorfile = open(os.path.join(path, fasta_file+'_'+matrix._matrixname+'_ErrorEstimate.csv'), "w")

    reader = utils.convertFASTAtoDict(os.path.join(path, fasta_file)) # reader
    writer = csv.writer(writerfile)
    with open(os.path.join(path, "Matches_%s_%s.bed" % (fasta_file.split("_")[0], matrix._matrixname)), "w") as bed_out:
        if errorf: errorf = csv.writer(errorfile)

        mm = matrixmotif(matrix, threshold)

        li = ["name", "location", "motif_locations", "count", "motif_seq", "results", "sequence (f)"]
        writer.writerow(li)

        total = 0
        tscore = []
        posinhit = []
        listlen = len(reader)
        i = 0

        headers = frozenset("Chr")
        regex = re.compile("\([+\-]\)")

        for data in reader:
            name = data["name"]
            if name.find(">hg18") == -1: # probably my format from extractSeqLists.?
                location = name.split("_")[1]
                strand = "+"

            else: # ralfs format, this is COSMIC or liftOver format?
                #>hg18.chr8(+):67737014-67737227|hg18_3
                #chr8(+):67737014-67737227
                t = name.split(".")[1].split("|")[0]
                location = t.replace("(+)", "").replace("(-)", "")
                strand = regex.search(t).group().replace("(", "").replace(")", "")
                    # always feed the + strand from the persepctive of the motif_scanner.
            seq = data["f"] if strand == "+" else data["r"]
                    #seq = data["f"]

            if errorf: # Do the error stuff;
                print("n:",i)
                i += 1
                errorf.writerow(mm.EstimateErrorRate(seq, 30))

            else:
                res = mm.find_matches(seq) # returns: [strand, pos, score, seq]
                if res:
                    total += len(res)
                    motiflocs = []
                    motifseqs = []
                    testers = []
                    for item in res:
                        pos = int(item[1])
                        t = utils.getLocation(location)
                        hit_strand = item[0]

                        # write out a BED file:

                         # I jsut need to correct the location depending upon the feed strand.
                        if hit_strand == "+":
                            newloc = "chr%s:%s-%s" % (t["chr"], int(t["left"])+pos+1, int(t["left"])+pos+len(item[3]))
                            bed_out.write("chr%s\t%s\t%s\t0\t%s\n" % (t["chr"], int(t["left"])+pos+1, int(t["left"])+pos+len(item[3]), strand))
                        elif hit_strand == "-":
                            newloc = "chr%s:%s-%s" % (t["chr"], int(t["left"])+pos+1-len(item[3]), int(t["left"])+pos)
                            bed_out.write("chr%s\t%s\t%s\t0\t%s\n" % (t["chr"], int(t["left"])+pos+1-len(item[3]), int(t["left"])+pos, strand))

                        motiflocs.append(newloc)
                        motifseqs.append(item[3])

                        testcode = "testLocation(hg18, \"%s\", \"%s\", \"%s\")" % (newloc, item[0], item[3])
                        testers.append(testcode)
                    writer.writerow([name, location, " ".join(motiflocs), len(res), " ".join(motifseqs), res, " ".join(testers), seq])

                    # write the fasta list, discard empty finds.. .er obviously...
                    for item in res: # iterate over the list;
                        fastafile.write(">seq\n")
                        fastafile.write(item[3]+'\n\n')
                        tscore.append(item[2])
                        posinhit.append(item[1])
                else:
                    if bKeepAll:
                        writer.writerow([name, location, 0, "None", seq])

        writerfile.close()
        fastafile.close()
    if errorf: errorfile.close()

    return({"total": total, "t_score_dist": tscore, "pos_in_hit": posinhit, "listlen": listlen})


if __name__ == "__main__":
    # windows
    path = "/home/hutchinsa/genomePeriodicity/"
    """
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]"""

    #
    # multiple sox's elements: c(a/t)ttgt(c/g)-c(a/t)ttgt(c/g) = cwttgtscwttgts
    scanPeriodicityFASTA("mwttswnnnnnnnnnynkccty", None, path, "mm8_refseq.csv_l8000bp_r2500bp.fa") # soxf_obp_soxf

