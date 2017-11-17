"""

These are a series of undocumented routines for QC of RNA-seq, primarily FASTQ files.


"""


import sys, os, math, time, random # Don't worry, random is just for db name
from numpy import array
import numpy
import matplotlib
import pylab as plot
import sqlite3

def rnaseqqc(filename, tmp_dir="/tmp/"):
    """
    **Purpose**
        Perform some common QC on high throughput sequencing data.
    
    **Arguments**
        filename
            A fastq file of sequence data
            
        tmp_dir (Optional, default="/tmp/")
            Override the location to store tmp data in
            
    """
    print("Started '%s'" % filename)
    oh = open(filename, "rU")
    
    res = None
    num_reads = 0
    # set-up sql db
    conn = sqlite3.connect(os.path.join(tmp_dir, "temp_db-%s-%s.sql" % (time.time(), hex(100000000+random.randint(0, 10000000)))))
    curs = conn.cursor()
    curs.execute("CREATE TABLE main (tabs TEXT PRIMARY KEY)")
    n = 0
    m = 0
    h = 0 
    max_q_score = 42 # maximum accepted quality score
    
    while 1:
        try: # Test still in file and on top of correct entry
            if oh.readline()[0] == "@":
                pass
        except IndexError:
            break
    
        seq = oh.readline().strip() # sequnece
        oh.readline() # strand (always +)
        qual = oh.readline().strip("\n")
    
        num_reads += 1  
    
        try: 
            len(res) # sample a typical entry to get the length of the read
        except TypeError: # Avoid numpy truth testing sillyness
            # This works as len(None) = TypeError, whereas len(numpy) == len(numpy) i.e. not None
            # If you try len(array([], numpy.int64)) then it raises an exception that truth testing doesn't make sense
            # on Numpy arrays, I only want __nonzero__() dumbasses!!! Stupid eh?
    
            # Only gets run once. 
            res = array([0 for x in range(len(qual))], numpy.int64)
            bp_dist = [{"A": 0, "C": 0, "G": 0, "T": 0} for x in range(len(qual))]
            res_std = numpy.zeros((len(qual), max_q_score))
    
        for i, c in enumerate(qual): # Collect the aggregate scores
            #print res
            score = ord(c) - 64 # post 1.3 scores
            res[i] += score
            res_std[i, score] += 1
        
        for i, bp in enumerate(seq): # Collect bp bias.
            if bp.upper() in ["A", "C", "G", "T"]:
                bp_dist[i][bp.upper()] += 1
    
        # see if table is already available
        result = curs.execute("SELECT * FROM main WHERE tabs=?", (seq[0:4],))
        result = result.fetchone()
        if not result:
            curs.execute("CREATE TABLE tab_%s (seq TEXT, count INT)" % (seq[0:4],))
            curs.execute("INSERT INTO main VALUES (?)", (seq[0:4], ))
        
        # check whether to increment or insert:
        result = curs.execute("SELECT * FROM tab_%s WHERE seq=?" % seq[0:4], (seq,))
        
        result = result.fetchone()
        if result: # increment existing value
            #curs.execute("UPDATE tab_%s SET count=count+1 WHERE seq=?" % seq[0:4], (seq,))
            curs.execute("UPDATE tab_%s SET count=? WHERE seq=?" % seq[0:4], (result[1]+1, seq))
        else:
            curs.execute("INSERT INTO tab_%s VALUES (?, 1)" % seq[0:4], (seq, ))
        
        h += 1
        if h > 1000:
            conn.commit() # Stop too much caching
            h = 0
        
        n += 1
        if n > 100000:
            m += 1
            print("%se5 reads done" % m)
            n = 0
            #break
        
    oh.close()
    conn.commit()
    
    frac_res = [x / num_reads for x in res]
    print("File: %s" % filename)
    print("Number of FASTQ entries: %s" % num_reads)

    # Draw some graphs:
    filename_base = "_".join(os.path.split(filename)[1].split(".")[:-1])
    fig = plot.figure()
    ax = fig.add_subplot(111)
    for bp in ["A", "C", "G", "T"]:
        ax.plot([(bp_dist[i][bp] / num_reads) * 100.0 for i in range(len(bp_dist))], label=bp)
    ax.legend()
    ax.set_title("%s - bp bias" % filename_base)
    ax.set_xlabel("base pair")
    ax.set_ylabel("Percent")
    ax.set_ylim([0,40])
    fig.savefig("%s_bpdist.png" % filename_base)
    
    # Work out standard error.
    err = [(numpy.std(i)) for i in res_std]
    
    fig = plot.figure()
    ax = fig.add_subplot(111)
    ax.pcolor(res_std.T, antialiased=False)
    ax.set_title("%s - Phred Quality Scores" % filename_base)
    ax.set_xlabel("base pair")
    ax.set_ylabel("Phred Quality Score")
    ax.set_ylim([0,max_q_score])
    fig.savefig("%s_qual_score.png" % filename_base)
    
    # and finally the uniquness plot.
            
    freq = numpy.zeros(500, dtype=numpy.int64)
    x_axis = numpy.arange(500)
    
    # sql get
    for table in conn.execute("SELECT * FROM main"):
        for read in conn.execute("SELECT * FROM tab_%s" % table):
            if int(read[1]) >= 499:
                freq[499] += 1
            else:
                freq[int(read[1])] += 1
        # count unqs and dupes
    
    fig = plot.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x_axis, freq, s=20, c="red", alpha=0.8, edgecolor="none")
    ax.set_yscale("log", basey=10)
    ax.set_ylim([0.1, max(freq)*10])
    ax.set_xlim([-10, 510])
    ax.set_xlabel("Number of occurences of read")
    ax.set_ylabel("Number of reads (Log10)")
    ax.set_title("%s - Phred Quality Scores" % filename_base)
    ax.text(300, 10000, "Uniques:")
    ax.text(360, 10000, "%.1f%%" % ((freq[1]/ sum(freq))*100))
    ax.text(420, 10000, "%s" % int(freq[1]))
    
    ax.text(300, 100, "Dupes:")
    ax.text(360, 100, "%.1f%%" % ((sum(freq[2:])/ sum(freq))*100))
    ax.text(420, 100, "%s" % int(sum(freq[2:])))
    
    fig.savefig("%s_repeat_reads.png" % filename_base)
    
    

    