"""

logo.py

tools for drawing logos from sequences, fasta files, genelists, etc.

"""



import sys, os, math
from operator import itemgetter

import numpy

from . import config, utils
from .genelist import genelist
from .errors import NotImplementedError

import matplotlib.pyplot as plot
import matplotlib.image as mpimg

class logo:
    """
    **Purpose**
        Draw sequence logos:
        http://en.wikipedia.org/wiki/Sequence_logo
        
        NOTE:
        logo.py will only deal with DNA motifs at present.
        
    **Arguments**
        sequence_data (Optional)
            the sequence data to load, can be a list, a fasta filename or a genelist-like object
            
        numpy_matrix (Optional)
            A numpy matrix:
                [[A, C, G, T] [A, C, G, T] ...]
    """
    
    def __init__(self, sequence_data=None, numpy_matrix=None):
        self.__calc_done = False
        if sequence_data:
            self.load_seq(sequence_data)
        elif numpy_matrix.any():
            self.__load_numpy(numpy_matrix)
        
    def __repr__(self):
        return("<glbase.logos.logo>")
        
    def __str__(self):
        return(str(self.data))

    def get_pfm(self):
        self.__gen_freq_table()
        return(self.freq)

    def print_pfm(self, style="normal", transpose=False):
        """
        **Purpose**
            print out the pfm in a variety of styles
            
        **Arguments**
            style (Optional, default="normal")
                valid styles are: 
                    "normal" - a tab separated pfm
                    "lol" - list of lists.
            
            transpose (Optional, default=False)
                transpose the matrix. By default matrieces are printed:
                    a []
                    c []
                    g []
                    t []
                If transpose is set to True:
                    [a, c, g, t],
                    [a, c, g, t],
                    ...
                    [a, c, g, t]
        
        **Returns**
            Nothing
        """        
        if style not in ["normal", "lol"]:
            raise NotImplementedError("%s style not found" % style)
        
        print(self.freq)
        
        if transpose:
            if style == "normal":
                raise NotImplementedError("%s, transposed, not done yet, sorry" % style)
            elif style == "lol":
                newa = []
                for k in ["a", "c", "g", "t"]:
                    newa.append([i[k] for i in self.freq])
                n = numpy.array(newa).T
                print(n)
        else:
            if style == "normal":
                for k in ["a", "c", "g", "t"]:
                    print("%s:%s" % (k, "\t".join([str(i[k]) for i in self.freq])))
            elif style == "lol":
                for k in ["a", "c", "g", "t"]:
                    print("[%s]" % (", ".join([str(i[k]) for i in self.freq])))
            
        
    def load_seq(self, sequence_data):
        """
        **Purpose**
            Load in a bunch of sequence data
            
        **Arguments**
            sequence_data (Required)
                sequenece data, can be a list, a fasta filename, or a genelist-like object with a 
                "seq" key containing suitable sequence data
        
        **Returns**
            None
        """
        assert sequence_data, "sequence data appears empty"
        
        if isinstance(sequence_data, list):
            self.__load_list(sequence_data)
        elif isinstance(sequence_data, str):
            self.__load_file(sequence_data)
        else:
            try:
                self.__load_genelist(sequence_data)
            except Exception:
                pass
        
    def __load_file(self, filename):
        # See if the file is a fasta file:
        oh = open(filename, "rU")
        for i in range(100): # take a big sample
            l = oh.readline()
            if ">" in l:
                # It's probably a FASTA
                self.__load_list(utils.FASTAToLIST)
                return(True)
    
        seql = []
        # Okay, assume it's a list of sequences
        oh = open(filename, "rU")
        for line in oh:
            seql.append(line.strip())
        self.__load_list(seql)
        return(True)
    
    def __load_genelist(self, gl_object):
        raise NotImplementedError("genelist loading for logos not implemented")
    
    def __load_list(self, list_obj):
        self.data = numpy.array([[c for c in l] for l in list_obj], dtype='|S1').T # Why do I use |S1?
        self.seq_len = len(list_obj[0])
        self.num_items = len(self.data[0]) 
        self.__gen_freq_table()

    def __load_numpy(self, numpy_array):
        self.data = numpy_array
        self.seq_len = len(numpy_array)
        self.num_items = sum(self.data[0])
        # assumes ACGT
        self.freq = []
        for bp in self.data:
            self.freq.append(dict(a=bp[0], c=bp[1], g=bp[2], t=bp[3]))
        
    def __gen_freq_table(self):
        self.freq = []
        for row in self.data:
            t = {"a": 0, "c": 0, "g": 0, "t":0}
            for c in row:
                t[c.lower()] += 1
            self.freq.append(t)
    
    def __calc(self):
        """
        Calculate the data
        Do this just before drawing.
        """
        #self.freq must be valid
       
        self.heights = []
        en = (4-1) / (2 * math.log(2) * self.num_items)
        en=0
        
        for freq in self.freq:
            height = {"a": 0, "c": 0, "g": 0, "t":0}
            rel_freq = {}
            for k in freq:
                if freq[k] == 0:
                    rel_freq[k] = 0
                else:
                    rel_freq[k] = (freq[k]/self.num_items)
            sum_freq = sum([i for i in list(rel_freq.values())])
            
            ind = -sum([rel_freq[k] * math.log(rel_freq[k], 2) for k in rel_freq if rel_freq[k] != 0])
                      
            total_height = math.log(4,2) - (ind - en)
            
            for k in freq:
                if freq[k] == 0:
                    height[k] = 0.0
                else: 
                    height[k] = rel_freq[k] * total_height
                    """
                    rel_freq = (freq[k]/self.num_items)
                    Hi = -rel_freq * math.log(rel_freq, 2)
                    Ri = 2 - (Hi + en)
                    # work out each letter height
                    height[k] = rel_freq * Ri 
                    print Hi, en, Ri, height[k]
                    """
            self.heights.append(height)            
              
        self.__calc_done = True
        
    def test_reqd_chars_avail(self, char_col_list):
        """
        The logo drawer relies upon char_X.png being available, to get around severe limitations in
        the text drawing. However, I need matplotlib to generate these pngs first and make certain 
        they are available.
        """
        for k in char_col_list:
            #if not os.path.exists("char_%s.png" % k):
                fig = plot.figure(figsize=(2,2))
                ax = fig.add_subplot(111)
                ax.set_position([0,0,1,1])
                ax.set_xlim([0,1])
                ax.set_ylim([0,1])
                ax.text(0,0,k, size=190, weight=1000, color=char_col_list[k])
                [item.set_markeredgewidth(0.0) for item in ax.xaxis.get_ticklines()]
                [item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]
                [t.set_fontsize(8) for t in ax.get_yticklabels()] # generally has to go last.
                [t.set_fontsize(8) for t in ax.get_xticklabels()]
                ax.set_frame_on(False)
                fig.savefig("char_%s.png" % k)
                plot.close(fig)
        return(True)

    def draw(self, filename):
        if not self.__calc_done:
            self.__calc()
            
        self.test_reqd_chars_avail(dict(A="orange", C="blue", G="blue", T="orange")) # The function is generic, but the rest of log is not (yet)
            
        # get chars:
        imgs = {}
        for c in ["A", "C", "G", "T"]:
            imgs[c] = mpimg.imread("char_%s.png" % c)
            
        fig = plot.figure(figsize=(7,3))
        ax = fig.add_subplot(111)
        ax.set_position([0.07,0.05,0.9,0.9])
        
        a = [i["a"] for i in self.heights]
        c = [i["c"] for i in self.heights]
        g = [i["g"] for i in self.heights]
        t = [i["t"] for i in self.heights]
        sums = [i["a"] + i["c"] + i["g"] + i["t"] for i in self.heights]
        
        bars = numpy.arange(self.seq_len)
        
        for i, b in enumerate(bars):
            sor = [dict(bp="A", val=a[i]), 
                dict(bp="C", val=c[i]), 
                dict(bp="G", val=g[i]), 
                dict(bp="T", val=t[i])]
            sor = sorted(sor, key=itemgetter("val"))
            #sor.reverse() # smallest to biggest
            ctot = 0
            for item in sor:
                if item["val"] != 0:
                    ax.imshow(imgs[item["bp"]], extent=[b, b+1, ctot, ctot+item["val"]],
                        aspect="equal", interpolation="bilinear")
                    #ax.text(b, ctot, item["bp"], size=24, va="bottom", ha="center")
                    ctot += item["val"]
        
        #ax.bar(bars, a, ec="green", fc="none")
        #ax.bar(bars, c, bottom=a, ec="blue", fc="none")
        #ax.bar(bars, g, bottom=c, ec="orange", fc="none")
        #ax.bar(bars, t, bottom=g, ec="red", fc="none")
        #ax.bar(bars, sums)
        
        ax.set_xlim([0, self.seq_len])
        ax.set_ylim([0, 2])
        ax.set_ylabel("bits")
        ax.set_xticklabels(bars+1)
        
        [item.set_markeredgewidth(0.0) for item in ax.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]
        [t.set_fontsize(8) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(8) for t in ax.get_xticklabels()]
        #ax.set_frame_on(False)
        
        fig.savefig(filename)
        plot.close(fig)
        
if __name__ == "__main__":
    """
    test suite items 
    """
    a = [
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat",
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat",
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat",
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat",
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat",
        "aaaaaa",
        "acaaac",
        "acaaag",
        "agaaat"]
    
    l = logo(a)
    l.draw("logo.png")
    print(l.data)
    print(l.freq)
    