"""

glbase class to perform functions on fastq files.

Generally glbase is not efficient at this and other tools may be preferable. 

"""

import sys, os

from . import config
from .genelist import genelist
from .errors import AssertionError, NotSupportedError, DelayedListError
from .delayedlist import delayedlist
# Todo: Move rnaseqqc into fastq

class fastq(delayedlist):
    """
    reimplementation of delayedlist
    """
    def __init__(self, filename=None, qual_format=False, **kargs):
        """
        **Process the fastq file**        
        
        **Arguments**
            filename (Required)
                The filename to load
            
            qual_format (Required)
                Expects one of phred33 or phred64
                
        """
        genelist.__init__(self) # no kargs though. I want the empty version.

        assert filename, "No Filename"
        assert os.path.exists(filename), "%s not found" % (filename)
        assert qual_format, "You must specify the qual format"
        assert qual_format in ("phred33", "phred64"), "qual format '%s' not recognised"
        
        self.__phred_warning = False
        
        self.filename = filename
        self.name = os.path.split(self.filename)[1]
        self.filehandle = None
        self.qual_format = qual_format
        
        self._optimiseData()
        
    def __repr__(self):
        return("glbase.fastq")

    def getSequences(self):
        raise NotSupportedError("getSequences() not supported for fastq file")

    def _optimiseData(self):
        """
        The optimesData() __iter__ cycle is a bit funky inside delayedlists. Basically 
        _optimiseData opens the file and gets it ready for iter and resets the file between
        operations.
        
        From the user perspective they should run:
        
        fq = fastq(filename)
        for item in fq:
            print item
        
        fq.splitPE() 
        """
        if self.filehandle:
            self.filehandle.close()
        self.filehandle = open(self.filename, "rU")       
            
    def __iter__(self):
        """
        (Override)
        
        iter is now locked to a single fastq format.
        """
        id = self.filehandle.readline().strip()
        if not id:
            self._optimiseData()
            return
        seq = self.filehandle.readline().strip()
        strand = self.filehandle.readline().strip()
        qual = self.__conv_quals(self.filehandle.readline().strip())
        
        d = {"id": id, "seq": seq, "strand": strand, "qual": qual}
        
        yield d 

    def __conv_quals(self, raw_quals):
        """
        convert the quals into phred scores for internal representation
        """
        
        if self.qual_format == "phred33":
            qq = [ord(i)-33 for i in raw_quals]
        elif self.qual_format == "phred64":
            qq = [ord(i)-64 for i in raw_quals]
        
        if not self.__phred_warning:
            if True in [i>40 for i in qq]:
                config.log.warning("Phred Quality > 40, are you sure you selected the correct quality format?")
                self.__phred_warning = True
            elif True in [i<0 for i in qq]:
                config.log.warning("Phred Quality < 0, are you sure you selected the correct quality format?")
                self.__phred_warning = True
            
        return(qq)        

    def __phred33_str(self, qual):
        """ convert the quals to Phred33 """
        return "".join(chr(i+33) for i in qual)
        
    def splitPE(self, file1, file2):
        """
        **Purpose**
            split a paired-end file into two separate files.
            
            It's for when you get an interleaved file like this:
            
            @HWI-ST507:76:A81MKNABXX:7:1:5257:1711:1
            CTCCTAGAGGGAAATATGGAGCAATTACATATTGTTCTCTAGGA
            +
            *625229@@@@@@@@@@@@7;@@;@@@@@;@@@@=@7@;@7@##
            @HWI-ST507:76:A81MKNABXX:7:1:5257:1711:2
            CAGGAGGGTCTGTGGTAGAAGGCTGTTACATACATAATAAA
            +
            HHHHHHHHCHBEEEE9EEBEEEDCECFBFBCB>?ACC>C##
            
            Where the paired-end are sequentially ordered. 
            
            Note that there is only light checking of the paired reads and it assumes that
            each 
            
            Also note that by default this tool ALWAYS writes the quals in Phred33 format.
            
        **Arguments**
            file1, file2 (Required)
                names of the files for the two 
                
        """
        config.log.info("Started '%s'" % self.name)
        __unpaired_warning = False

        h1 = open(file1, "w")
        with open(file2, "w") as h2:
            paired = 0
            unpaired = 0

            read1 = True
            while True:
                try:
                    read1 = next(self.__iter__())
                    read2 = next(self.__iter__())
                except StopIteration:
                    break

                if read1["id"].split(":")[0:6] == read2["id"].split(":")[0:6]:
                    paired += 1
                    h1.write("%s\n" % read1["id"])
                    h1.write("%s\n" % read1["seq"])
                    h1.write("%s\n" % read1["strand"])
                    h1.write("%s\n" % self.__phred33_str(read1["qual"]))

                    h2.write("%s\n" % read2["id"])
                    h2.write("%s\n" % read2["seq"])
                    h2.write("%s\n" % read2["strand"])
                    h2.write("%s\n" % self.__phred33_str(read2["qual"]))
                else:
                    unpaired += 1
                    if __unpaired_warning:
                        config.log.warning("Unpaired reads found, there is a high chance that splitPE will incorrectly format the output!")
                        __unpaired_warning = True

            h1.close()
        config.log.info("splitPE(): Correctly paired: %s, unpaired: %s" % (paired, unpaired))