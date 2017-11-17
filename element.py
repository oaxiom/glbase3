"""

element.py

lists and tools to deal with motifs (words).

pwms, pwm and element need to be merged into a single 'motif.py'...

Part of glbase.

"""

import config

import re, random, sys, os, string, numpy, copy

import utils
from genelist import genelist as new_gl
from errors import AssertionError
from draw import draw
from numpy import zeros

# R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT],
# [ACT] = h, [ACG] = v, [AGT] = d, [CGT] = b
regex_dict = {
    "a" : "a",
    "c" : "c",
    "g" : "g",
    "t" : "t",
    "r" : "[ag]", # twos
    "y" : "[ct]",
    "k" : "[gt]",
    "m" : "[ac]",
    "s" : "[gc]",
    "w" : "[at]",
    "h" : "[act]", # threes
    "v" : "[acg]",
    "d" : "[agt]",
    "b" : "[cgt]",
    "n" : "[acgt]" # four
}

class motif:
    """
    **Purpose**
        a motif for word-searching in DNA elements
        
    **Arguments**
        name 
            name of the motif
            
        sequence
            the sequence in a degenerate DNA alphabet
            
        You can make an empty element with::
        
            e = motif("empty_name")
        
        Here are two examples::
        
            e = motif("Sox2", "ywttswn")
            e = motif("Oct4", "atgywnww")
        
        See pwm and pwms part of glbase for position weight matrices
    """
    def __init__(self, name=None, sequence="", **kargs):
        assert sequence, "no specified sequence for the motif"
        assert name, "you must give your motif a name"

        self.name = name
        self.draw = draw(self)
        self.seq = sequence
        self.palindromic = utils.isPalindromic(sequence)

        if "scramble" in kargs and kargs["scramble"]:   
            self._scrambleMotif()

    def isPalindromic(self):
        return(self.palindromic)

    def isPalindrome(self):
        return(self.palindromic)

    def getRegEx(self):
        """
        return the regex.
        """
        self.re = re.compile("".join([regex_dict[bp] for bp in self.seq.lower()]))
        return(self.re)

    def _scrambleMotif(self, number=10):
        """
        (Internal)
        generate scrambled motifs.

        Be careful with this, it does not work the same as in fexcom.
        It does not scramble the motif in-place, instead it returns a set of
        'regexd' motifs.
        """
        res = []
        for n in xrange(number):
            newSeq = oldSeq = self.seq
            while newSeq == oldSeq: # stops the new motif == old motif.
                q = []
                for item in self.seq:
                    q.append(item)

                new = []
                while q:
                    i = random.randint(0, len(q)-1)
                    new.append(q[i])
                    q.remove(q[i])

                newSeq = "".join(new)
            res.append(re.compile(newSeq))
        return(res)

    def __str__(self):
        return("motif name:%s sequence:%s" % (self.name, self.seq))

    def __len__(self):
        return(len(self.seq))

    def scan_sequence(self, seq=None, both_strands=True, mismatches=0):
        """
        **Purpose**
            scan the word across a single sequence of DNA
            
        **Arguments**
            seq (Required)
                the DNA sequence (as a string) to be scanned
                
            both_strands (Optional, default=True)
                search both strands of the sequence?
                
            mismatches (Optional, default=0)
                number of base pair mismatches to tolerate. 
        
        **Returns**
            A dictionary, containing the keys:
                num_motifs_found: The number of motifs discovered
                locations: The locations (relative to the 5' end of the seq)
                strands: List of strands the motif comes from
                sequences: A list of sequences discovered
            
        """
        assert seq, "You must provide a sequence"

        if both_strands:
            seqs = {"+": seq.lower(), "-": utils.rc(seq.lower())}
        else:
            seqs = {"+": seq.lower()}
        
        if mismatches == 0:
            return(self._scan_no_mismatch(seqs))
        else:
            return(self._scan_seqs_mismatches(seqs, mismatches))

    def _scan_no_mismatch(self, seqs=None):
        '''
        Search path for no mismatches
        
        This could be done a lot faster I think...
        '''
        seq_len = len(seqs['+'])
        local_pos = []
        found_seqs = []
        strands = []
        num_motifs = 0
        
        motif = self.getRegEx()
        
        for strand in seqs: 
            i = motif.finditer(seqs[strand])
            for m in i:
                if strand == "+":
                    local_pos.append(m.start()) 
                    found_seqs.append(seqs['+'][m.start():m.start()+len(self)])
                    strands.append(strand)
                    
                elif strand == "-":
                    st = abs(seq_len - m.end())
                    local_pos.append(st) # convert to get the 5' most
                    found_seqs.append(utils.rc(seqs['+'][st:st+len(self)]))
                    strands.append(strand)
                    
                num_motifs += 1
        return({"num_motifs_found": num_motifs, "locations": local_pos, "sequences": found_seqs, "strands": strands})

    def _getMismatchRegEx(self, mismatches):
        """
        return a list of regexs that take into account the
        number of degenerate base pairs.
        motif should be a degenerate acgtykrmn type sequence string.
        """
        rel = [regex_dict[bp] for bp in self.seq.lower()] # make the base regex to modify below.
        seq = self.seq.lower()

        rels = []
        
        if mismatches:
            pos = 0
            while pos < len(seq):
                if seq[pos] == "n":
                    pos += 1 # already a degenerate bp.
                else:
                    old = rel[pos]
                    rel[pos] = regex_dict["n"]
                    rels.append(re.compile("".join(rel)))
                    #print "".join(rel)
                    rel[pos] = old # replace the old bp.
                    pos += 1 # move on a bp.

        self.re = rels
        return(rels)

    def _scan_seqs_mismatches(self, seqs, mis_matches=1):
        '''
        Search path for mismatch searching
        
        At the moment this only reports the number of matched motifs
        
        '''
        
        assert mis_matches < 2, "more than 1 bp degenerate in the back search is not supported. Forcing 1 bp mismatch"
        assert mis_matches > 0, 'At least 1 bp mismatch required'

        rels = self._getMismatchRegEx(mis_matches)

        temp_finds = []

        for regex in rels:
            for strand in seqs:
                found = regex.finditer(seqs[strand])

                if found:
                    for match in found:
                        # Positions are currently wrong as they are strand relative?
                        temp_finds.append({"seq": seqs[strand][match.start():match.end()], "pos": "%s:%s:%s" % (strand, match.start(), match.end()) })
        
        found_list = []
        # prune duplicate locations;
        unique = []
        for item in temp_finds: # keep only the first element.
            if not item["pos"] in unique:
                found_list.append(item["seq"])
                unique.append(item["pos"])

        return({"num_motifs_found": len(unique), "locations": None, "sequences": found_list, "strands": None})
                
    def scan_sequences(self, seq_data=None, both_strands=True, keep_found_only=False):
        """
        **Purpose**
            Scan over a set of FASTA sequences and add a new key <name of motif> with
            the result either "Found" or None
            
        **Arguments**
            seq_data (Required)
                the DNA sequence, loaded from a FASTA or a genelist containing a "seq"
                key
                
            both_strands (Optional, default=True)
                search both strands of the sequence?
        
            keep_found_only (Optional, default=False)
                If set to True, then items in the genelist that do not contain the motif are
                not added to the returned list (returned list contains only DNA sequences
                that contain a motif)
        
        **Returns**
            A new genelist copy of seq_data with an additional new key named after the motifs name
            and indicating if the motif was found or not
        """
        assert seq_data, "You must provide sequence data"
        
        newl = []
        for item in seq_data:               
            res = self.scan_sequence(item["seq"], both_strands=both_strands)
                            
            if res["num_motifs_found"]:
                newi = copy.deepcopy(item)
                item[self.name] = "Found"
                item["%s_seqs" % self.name] = ", ".join(res["sequences"])
                newl.append(item)
            elif not keep_found_only:
                newi = copy.deepcopy(item)
                item[self.name] = "None"
                item["%s_seqs" % self.name] = "None"
                newl.append(item)
                    
        newgl = seq_data.shallowcopy()
        newgl.load_list(newl)
        return(newgl)                

    def count_sequences_with_motif(self, seq_data=None, both_strands=True, **kargs):
        """
        **Purpose**
            Return some statistics about the numbers and percent of motifs found in each list.
            
            Returns some information in a dict:
                {"num_motifs": <The total number of motifs found>,
                "with_motif": <The number of items in genelist with 1 or more motif>.
                "percent": <The percent of the genelist/fasta with 1 or more motif>}
                
        **Arguments**
            seq_data (Required)
                the DNA sequence, loaded from a FASTA or a genelist containing a "seq"
                key
                
            both_strands (Optional, default=True)
                search both strands of the sequence?
        """
        assert seq_data, "You must provide sequence data"
        
        results = {"num_motifs": 0, "with_motif": 0}
        newl = []
        for item in seq_data:               
            res = self.scan_sequence(item["seq"], both_strands=both_strands)

            if res["num_motifs_found"]:
                results["num_motifs"] += res["num_motifs_found"]
                results["with_motif"] += 1
        
        results["percent"] = (results["with_motif"] / float(len(seq_data))) * 100.0
        return(results)

    def getMotifMatches(self, genelist=None, both_strands=True, return_num_motifs_only=False):
        """
        **Purpose**
        
        scan across the genelist (must include sequence data) with this element and then 
        return a pileup spanning the length of the fasta, with a score for each time the element
        occurs.
        
        The fasta sequence lengths must be the same across all entries in the genelist.
        
        **Arguments**
            
        genelist (Required)
            the genelist to scan. Must contain sequence data in "seq", "sequence" or "f".
            If you load a FASTA file::
            
                f = genelist(..., format=format.fasta)
            
            then the sequence is already stored in "seq"
            
        both_strands (Optional, default=True)
            search both strands by default.
            
        return_num_motifs_only (Optional, default=False)
        
            By default getMotifMatches() returns a numpy array. 
            
            If this is set to True it returns a single integer of the number of motifs found
            
        **Returns**
        
        A numpy array going from 0 .. n bp in length. Each location there is a motif the score at that
        position will be +1. No motifs are 0.
        
        THe motif score extends along the length of the motif, so for example if your motif is 8 bp long then
        the data will end up looking like this:
        
        0,0,0,1,1,1,1,1,1,1,1,0,0,0.
        
        
        """
        assert genelist, "you must provide a genelist"
        
        # get the seq key:
        if "seq" in genelist.linearData[0]:
            seq_key = "seq"
        elif "f" in genelist.linearData[0]:
            seq_key = "f"
        elif "sequence" in genelist.linearData[0]:
            seq_key = "sequence"
        else:
            raise AssertionError, "the genelist does not appear to have a sequence attached" # fake an assert        
        
        seq_len = len(genelist[0][seq_key])
        pileup = zeros(seq_len)
        motif = self.getRegEx()
        local_pos = []
        num_motifs = 0

        for item in genelist.linearData:
            found = 0
            
            if both_strands:
                seq = {"+": item[seq_key].lower(), "-": utils.rc(item[seq_key].lower())}
            else:
                seq = {"+": item[seq_key].lower()}
                
            for strand in seq: 
                i = motif.finditer(seq[strand])
                if i:
                    tm = 0 # need to be inited, as an empty i will not enter the for loop
                    for m in i:
                        if strand == "+":
                            local_pos.append(m.start()) 
                        elif strand == "-":
                            local_pos.append(abs(seq_len - m.end())) # convert to get the 5' most
                        num_motifs += 1

        for pos in local_pos:
            for p in range(pos, pos+len(self)):
                pileup[p] += 1
        
        if return_num_motifs_only:
            return(num_motifs)
        return(pileup)

    def scanMotifFrequency(self, genelist=None, filename=None, sum_randoms=True, centre=True, 
        random_fastas=None, normalise=True, **kargs):
        """
        **Purpose**
    
            scan a set of sequences for the presence of this particular DNA word.
            score each fasta sequence for a match for a particular motif
            then output the result as a movingAverage graph. See also getMotifMatches()
    
            NOTE:
                All Fasta sequences must be of the same length, if not it will produce an error

        **Arguments**

            genelist (Required)
                genelist with sequence data attached, must be in a "seq", "sequence"
                "f". If you load a fasta, the sequence is loaded into "seq".
                You can send a list of genelists as well, and these will all be plotted
                on the graph.
    
            random_fastas (Optional, default=None)
                a single or list of fastas genelist to act as a random background.
    
            filename
                filename of the resulting image file.
    
            window (Optional, default = 20% of the first list in 'genelist')
                the moving window average size in % of the total size of the fasta
    
            normalise (Optional, default=True)
                normalise the motif frequency to the size of the list
    
            sum_randoms (default=False)
                sum the random lists together, default is to show a grey line for each random
                list supplied.
            
            centre (Optional, default=True)
                centre the x axis (ie. for a span of bp from 0->400 bp, move the axis to
                -200 -> 200 bp instead, so that 0 is the centre of the sequence.
    
        **Result**
            returns a dictionary containing the data and the labels and a figure saved to filename
        """
        # draw a moving average plot
        assert filename, "you must specify a filename to save the graph to"
        assert genelist, "you must provide a genelist"

        # turn genelist into a vanilla list if it is not already.
        if not isinstance(genelist, list):
            genelist = [genelist]

        if random_fastas and not isinstance(random_fastas, list):
            random_fastas = [random_fastas]

        motif = self.getRegEx()
        __bWarningDoneAlready = False
        fasta_lists = genelist

        # get the seq key:
        if "seq" in genelist[0][0]:
            seq_key = "seq"
        elif "f" in genelist[0][0]:
            seq_key = "f"
        elif "sequence" in genelist[0][0]:
            seq_key = "sequence"
        else:
            raise AssertionError, "the genelist does not appear to have a sequence attached" # fake an assert

        motif_result = []
        labels = []
        local_pos = []
        pos_results = []
        seq_len = len(genelist[0][0][seq_key])

        windowf = 0.2
        window = int(seq_len * windowf) # 20% of the list
        if "window" in kargs and kargs["window"]:
            window = int(kargs["window"])
            windowf = window/100.0

        pileups = []
        # Do scan
        for fasta in fasta_lists:           
            pileup = self.getMotifMatches(genelist=fasta)
            if normalise:
                pileup /= len(fasta)
            pileups.append(pileup)
            labels.append(fasta.name)
        
        rand_pileups = []
        if random_fastas:
            labels.append("Background")
            for fasta in random_fastas:
                pileup = self.getMotifMatches(genelist=fasta)
                if normalise:
                    pileup /= len(fasta)
                rand_pileups.append(pileup)
                
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        for i, p in enumerate(pileups):
            x, ma = utils.movingAverage(p, len(p)*windowf, False)
            if centre:
                x = numpy.arange(len(ma)) - (len(ma)//2)
            ax.plot(x, ma, label=labels[i])
        
        if random_fastas:
            for i, p in enumerate(rand_pileups):
                x, ma = utils.movingAverage(p, len(p)*windowf, False)
                if centre:
                    x = numpy.arange(len(ma)) - (len(ma)//2)
                if i == 0:
                    ax.plot(x, ma, color="grey", alpha=0.5, label="Background")
                else:
                    ax.plot(x, ma, color="grey", alpha=0.5)
        
        ax.legend()
        ax.set_xlabel("Position around FASTA (bp)")
        if normalise:
            ax.set_ylabel("Normalised motif frequency (arbitrary units)")
        else:
            ax.set_ylabel("Motif frequency (counts)")
        
        #if centre:
        #    xlabs = [
        
        self.draw.do_common_args(ax, **kargs)
        actual_filename = self.draw.savefigure(fig, filename)

        #config.log.info("Found: %s (%.1f%%) of '%s' motifs in fasta list '%s'" % (sum(motif_result[0]), (sum(motif_result[0])/float(len(fasta_lists[0])))*100, self.seq, fasta_lists[0].name))
        config.log.info("Saved figure '%s'" % actual_filename)
        return({"data": pileups, "labels": labels, "background": rand_pileups})
