"""

pwms.py

An amalgamation of a set of pwms.

"""

import sys, os, csv, random, string
import config, utils

from numpy import zeros, array, mean, std
from draw import draw
from genelist import genelist
from errors import AssertionError, NotImplementedError
from progress import progressbar
from location import location
#from history import historyContainer

import pwm

class pwms(genelist):
    """
    **Purpose**
        A collection of pwm items with a few specialised
        methods for the analysis of multiple motifs.

    **Arguments**
        filename (Required)
            some sort of file or path that I can make sense of as a pwm.

        name (Optional)
            The name of the collection of pwms.

        format (Required)
            The format of the pwm file.
            Valid formats currently supported are:
            "uniprobe"
            "transfac"
            "jaspar_matrix_only"
            "ACGT_rows"
            "Jolma_csv"
            "HOMER"

            to add, but currently unsupported formats:
            "jaspar"
    """
    def __init__(self, filename=None, name="None", format=None, **kargs):
        genelist.__init__(self, filename=None, name=name, format=format, **kargs) # no filename, get an empty genelist
        assert filename, "filename argument not specified"
        assert format, "no format specifier provided"
        assert os.path.exists(filename), "file '%s' not found" % filename

        self.linearData = []
        self.name = name
        self.draw = draw(self)
        #self._history = historyContainer(self) # a list of historyItem's

        self.format = format # this is used if saving

        if format == "uniprobe":
            self.__load_uniprobe(filename)
        elif format == "jaspar":
            raise NotImplementedError
        elif format == "jaspar_matrix_only":
            self.__load_jaspar_matrix_only(filename)
        elif format == "transfac":
            self.__load_transfac(filename)
        elif format == "ACGT_rows":
            self.__load_ACGT_rows(filename)
        elif format == "Jolma_csv":
            self.__load_jolma_csv(filename)
        elif format == 'HOMER':
            self.__load_HOMER(filename)
        else: 
            config.log.critical("'%s' format is unknown" % format)
        
        config.log.info("loaded '%s' with %s pwms" % (filename, len(self)))

        self._optimiseData()    

    def __load_HOMER(self, filename):
        '''
        load HOMER format, basically load_ACGT_rows with enhanced header splitting
        
        >AAAGCATTGA	14-AAAGCATTGA,BestGuess:PB0197.1_Zfp105_2/Jaspar(0.625)	10.296193	-829.609933	0	T:187.0(2.20%),B:4.8(0.01%),P:1e-360
        0.700	0.100	0.100	0.100
        0.700	0.100	0.100	0.100
        0.700	0.100	0.100	0.100
        0.100	0.100	0.700	0.100
        0.100	0.700	0.100	0.100
        0.700	0.100	0.100	0.100
        0.100	0.100	0.100	0.700
        0.100	0.100	0.100	0.700
        0.100	0.100	0.700	0.100
        0.700	0.100	0.100	0.100
        '''
        oh = open(filename, 'rU')
        
        pwm_store = None
        
        for line in oh:
            if line and line[0] != "#":
                if line[0] == ">":
                    if pwm_store: # save the finished pwm
                        self.linearData.append({"name": name, 'motif_seq': motif_seq, 
                            'p': p, 'score': score,
                            "pwm": pwm.pwm(name=name, pwm_matrix=pwm_store, isPFM=False)})
                    # setup a new pwm:
                    pwm_store = []
                    tt = line.strip().replace(">", "").split('\t')
                    #print tt
                    name = tt[1].split(':')[1].split('/')[0]
                    motif_seq = tt[0]
                    score = tt[2]
                    p = tt[5].split(':')[-1]
                else:
                    pwm_store.append(line.strip().split())
                    
        # store the last motif in the file
        if pwm_store: # save the finished pwm
            self.linearData.append({"name": name, 'motif_seq': motif_seq, 
                'p': p, 'score': score,
                "pwm": pwm.pwm(name=name, pwm_matrix=pwm_store, isPFM=False)})
        oh.close()   
        return()     

    def __load_jolma_csv(self, filename):
        """
        Load in a csv of Jolma's supp table S3
        """
        oh = open(filename, "rU")
        
        name = oh.readline().strip().split(",")[0]
        
        names = [name] # Jolma names are not unique
        
        while name:
            A = [int(i) for i in oh.readline().strip().split(",")[1:] if i]
            C = [int(i) for i in oh.readline().strip().split(",")[1:] if i]
            G = [int(i) for i in oh.readline().strip().split(",")[1:] if i]
            T = [int(i) for i in oh.readline().strip().split(",")[1:] if i]
        
            data = [A, C, G, T]
                
            pfm = array(data)
            pfm = pfm.T
        
            self.linearData.append({"name": name, 
                "pwm": pwm.pwm(name=name, pwm_matrix=pfm, isPFM=False)})
            
            oriname = oh.readline().strip().split(",")[0] # get the next name
            name = oriname
            n = 1
            while name in names:
                n += 1
                name = "%s_%s" % (oriname, n)
            names.append(name)
        
        oh.close()

    def __load_uniprobe(self, path):
        """
        Load the uniprobe pwm data set.

        To get the uniprobe data set, go to uniprobe>Downloads and
        download 'all data'. You should have a file called "all_data.zip"

        unzip that to give a folder:
        "ALL_DATA/"

        this function will then step through the data and load all of it in.
        """
        post_process = False
        
        for subdir, dirs, files in os.walk(path):
            for file in files:
                # get rid of the extension and use this to append the motif.
                if (file[0] != ".") and ("readme" not in file): # Mac OSX specific malarky.
                    name = file.split(".")[0]

                    oh = open(os.path.join(subdir, file), "rU")
                    
                    data = []
                    
                    for line in oh:
                        if "#" not in line:
                            l = line.replace("\n", "").split("\t")[1:]
                            if len(l) > 3:
                                l = [float(i) for i in l]
                                data.append(l)
                            if "Probability matrix" in line:
                            	post_process = True
                            	# THere's a change in format 
                            	# later, and they give 4 motifs
                            	# for the same TF (different representations)
                            	# we only want the last one 9the probabilty matrix
						
                    # convert to numpy.
                    pfm = array(data)
                    
                    if post_process:
						# keep only the last four columns
						pfm = pfm[14:18]
						post_process = False
                    
                    pfm = pfm.T

                    self.linearData.append({"name": name, 
                        "pwm": pwm.pwm(name=name, pwm_matrix=pfm, isPFM=False)})

    def __load_jaspar_matrix_only(self, filename):
        """
        Load the JASPAR data, this accepts the "matrix_only" option from the JASPAR website
        
        the file should be something like "matrix_only.txt" and should look like this:
        >MA0002.1 RUNX1
        A  [10 12  4  1  2  2  0  0  0  8 13 ]
        C  [ 2  2  7  1  0  8  0  0  1  2  2 ]
        G  [ 3  1  1  0 23  0 26 26  0  0  4 ]
        T  [11 11 14 24  1 16  0  0 25 16  7 ]
        """
        oh = open(filename, "rU")
        
        done = False
        while not done:
            line = oh.readline()
            if ">" in line:
                name = line.replace(">", "").replace("\n", "").split(" ")[1]
                id = line.replace(">", "").split(" ")[0]
                A = [float(x) for x in oh.readline().strip("ACGT[]\n ").replace("   ", " ").replace("  ", " ").replace("    ", " ").split(" ")]
                C = [float(x) for x in oh.readline().strip("ACGT[]\n ").replace("   ", " ").replace("  ", " ").replace("    ", " ").split(" ")]
                G = [float(x) for x in oh.readline().strip("ACGT[]\n ").replace("   ", " ").replace("  ", " ").replace("    ", " ").split(" ")]
                T = [float(x) for x in oh.readline().strip("ACGT[]\n ").replace("   ", " ").replace("  ", " ").replace("    ", " ").split(" ")]
                data = [A, C, G, T]
                #print data
                
                pfm = array(data)
                pfm = pfm.T
                
                self.linearData.append({"name": name, "pwm": pwm.pwm(name=name, pwm_matrix=pfm, isPFM=False), "id": id})
            else:
                done = True # no empty lines in the JASPAR format
                
    def __load_ACGT_rows(self, filename):
        """
        Load a ACGT in rows arrangement (e.g. output from fexcom):
        
            >ap2af_-1bp_gcm1f
            22	678	612	21
            27	1238	31	37
            24	1261	11	37
            378	498	4	453
            294	430	592	17
            419	373	522	19
            518	12	776	27
            22	46	1234	31
            310	571	442	10
            377	28	514	414
            43	1078	85	127
            32	1252	18	31
            576	705	7	45
            324	45	613	351
            27	1218	76	12
            1193	67	24	49
            25	32	839	437
            347	52	672	262
        
        """
        
        oh = open(filename, "rU")
        
        pwm_store = None
        
        for line in oh:
            if line and line[0] != "#":
                if line[0] == ">":
                    if pwm_store: # save the finished pwm
                        self.linearData.append({"name": name, "pwm": pwm.pwm(name=name, pwm_matrix=pwm_store, isPFM=False)})
                    # setup a new pwm:
                    pwm_store = []
                    name = line.strip().replace(">", "")
                else:
                    pwm_store.append(line.strip().split())
        oh.close()   
                
    def __load_transfac(self, filename):
        """
        Load in transfac db. In this case transfac comes as a single file,
        so you need to send me that file.
        """
        oh = open(filename, "rU")

        newpwm = {"desc":""} # temp store.
        data = []
        degen = []

        for line in oh:
            two_chars = line[0:2]

            # process easily read labels:
            if two_chars == "ID":
                newpwm["ID"] = line[2:].replace("\n", "")
            elif two_chars == "NA":
                newpwm["name"] = line[2:].replace("\n", "")
            elif two_chars == "DE":
                newpwm["desc"] = line[2:].replace("\n", "")
            elif two_chars == "//":
                # finish the pwm, clean it up and generate a new one.
                pfm = array(data)
                self.linearData.append({"id": newpwm["ID"], "pwm": pwm.pwm(name=newpwm["name"], pwm_matrix=pfm, isPFM=False),
                    "degen": "".join(degen), "name": newpwm["name"], "desc": newpwm["desc"]}) # don't convert to pwm
                newpwm = {"desc": ""}
                data = []
                degen = []

            # try to grab the motif:
            try:
                p = int(line[:2])
                while p > len(data):
                    data.append([])
                while p > len(degen):
                    degen.append("")

                t = line.replace("\n", "").split(" ")[1:]
                newt = [float(i) for i in t[:-1] if i != ""]

                data[p-1] = newt
                degen[p-1] = t[-1]
            except ValueError:
                pass

        return(None)

    def __repr__(self):
        return("glbase.pwms")

    def scan_sequence(self, loc, sequence, features, filename=None,
        merge_strands=True, summary_file=None,
        z_score=4.0, **kargs):
        """
        **Purpose**

        **Arguments**
            sequence
                The DNA sequence to scan.

            features
                send a dictionary, derived from a genelist loaded from
                say refseq/ucsc and this will know what to do with it.

            filename (Required)
                The filename to save the image files to.

            summary_filename (Required)
                filename to save a summary of the interesting hits to,
                saves a tsv file.

            merge_strands (Optional, default=True)
                merge the score strands. At the moment only this method
                is supported and separate scores per strand are not supported.

            z_score (Optional, default = 4.0)
                The z-score to use to call an interesting peak.
                By default it is set to 4 standard deviations
                away from the mean.

        **Returns**
            saves an image file to filename.
        """
        result = []

        # cooerce the location
        loc = location(loc=loc)

        p = progressbar(len(self))
        for pi, pwm in enumerate(self.linearData):
            data = pwm["pwm"].scan_sequence(sequence, merge_strands)

            # get the significant hits
            m = mean(data)
            s = std(data)
            # collect interesting results
            z_reqd = m+(s*z_score)
            for i, v in enumerate(data):
                if v > z_reqd:
                    result.append({"loc": location(chr=loc["chr"], left=loc["left"]+i, right=loc["left"]+i+len(pwm["pwm"])),
                        "score": v, "strand": "?", "seq": sequence[i:i+len(pwm["pwm"])], "Z-score": abs(m-v)/s,
                        "name": pwm["name"], "local_location": i})

            real_filename = "%s_%s.png" % (self.name, pwm["name"])
            self.draw._plot_and_histogram(filename=real_filename, data=data,
                figsize=(20,5), title=pwm["name"], ylims=(0,1), xlims=(0, len(data)),
                loc=loc, genomic_features=features)

            del data # clean up the array, some really big searches kill
                     # the 'puter.

            p.update(pi)
        config.log.info("Saved %s image files" % len(self))

        oh = open(summary_file, "w")
        oh.write("motif\tloc\tscore\tstrand\tseq\tZ-score\n")
        for item in result:
            oh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (item["name"], str(item["loc"]), item["score"],
                item["strand"], item["seq"], item["Z-score"]))
        oh.close()
        config.log.info("Saved a summary file to '%s'" % summary_file)

        # save a summary label image.
        sum = []
        for item in result:
            sum.append({"pos": (item["__local_location"], item["Z-score"]), "label": item["name"]})
            
        real_filename = self.draw._labeled_figure(data=sum, axissize=(len(sequence), 0), ylim=(4,5.5),
            figsize=(20,5), genomic_features=features, filename=summary_file,
            loc=loc)
            
        config.log.info("Saved a summary image file to '%s'" % real_filename)

        return(real_filename)

    def saveCSV(self, filename=None, **kargs):
        """
        **Purpose**
            save the pwm list as a csv.
            This doens't work completely - the pwm matrix is not correctly saved.

        **Arguments**
            see genelist.saveCSV()
        """
        if "key_order" in kargs and kargs["key_order"]:
            key_order = kargs["key_order"]
            if "pwm" in key_order:
                key_order.remove("pwm")
        else:
            key_order = [k for k in self.linearData[0] if k != "pwm"]
        key_order.append("pwm")

        genelist.saveCSV(self, filename=filename, key_order=key_order, **kargs) # inherit

    def saveFASTA(self, filename=None, **kargs):
        """
        **Purpose**
            save the pwm list as a FASTA-style entry
            all other keys are stuffed into the > line

        **Arguments**
            filename (REquired)
                the filename to save to
                
            key_order (Optional, default=<'name' first, all other keys random>)
                key_order for stuffing the FASTA header line
        """
        
        if "key_order" in kargs and kargs["key_order"]:
            key_order = kargs["key_order"]
        else:
            ks = self.keys()
            ks.remove('name')
            ks.remove('pwm')
            if ks:
                key_order = ks

        oh = open(filename, 'w')
        for pwm in self.linearData:
            if key_order:
                oh.write('>%s\t%s\n' % (pwm['name'], '\t'.join([pwm[k] for k in key_order])))
            else:
                oh.write('>%s\n' % (pwm['name'], ))
                
            #print pwm['pwm'].get_matrix()
            for bp in pwm['pwm'].get_matrix():
                oh.write('%s\n' % ('\t'.join([str(i) for i in bp])), )
        oh.close()

    def save(self, filename=None, mode="binary"):
        """
        **Purpose**
            Save the matrix motifs in a set format.
            
        **Arguments**
            filename (Required)
                the filename to save the data to.
                
            mode (Optional, default="binary")
                the mode to save the file as. Defualts to the glbase binary.
                modes:
                    binary - glb binary file
                    cisfinder - a mode compatible with cisfinder
        **Returns**
            none
        """
        assert filename, "filename not specified"
        
        if mode == "binary":    
            base_genelist.save(self, filename)
        elif mode == "cisfinder":
            self.__save_cisfinder(filename)
        config.log.info("Saved '%s'" % filename)
        return(None)
        
    def __save_cisfinder(self, filename):
        oh = open(filename, "w")
        
        for p in self:
            m = p["pwm"].get_matrix()
            
            oh.write(">%s\tNNNNNNN\tNNNNNNNN\t0\t0\t0\t0\n" % p["name"])
            
            for i, row in enumerate(m):
                if self.format == "uniprobe":
                    oh.write("%s\t%s\n" % (i, "\t".join([str(int(x*1000)) for x in row])))
                else:
                    oh.write("%s\t%s\n" % (i, "\t".join([str(int(x)) for x in row])))
                
            oh.write("\n")
        
        oh.close()