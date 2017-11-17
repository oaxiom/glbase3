"""

realtime.py

part of glbase

Used to tidy up qRT-PCR data.

This is essentially deprecated now in favour of the new realtime functionality 
(which is basically a tool wrapper around expression())

Supported platforms and file formats:
-------------------------------------
. copy and paste from SDS 2.2.2 and 2.4 (tested 2010)
. biomark (not really tested though since 2009)

"""

import csv, math, numpy

import config

class assay:
    def __init__(self, _sampleName, _primerName, _output = None):
        self.sample = _sampleName
        self.primer = _primerName
        self.sampleRef = None
        self.primerRef = None

        self.Ct = []
        self.Ct_call = []
        # derived numbers;
        self.log10ddCt = []
        self.log2ddCt = []
        self.percentddCt = []

    def __str__(self):
        return("sample: %s, primer: %s, ct: '%s'" % (self.sample, self.primer, self.Ct))

    def addCt(self, _Ct, _Call = None):
        self.Ct.append(_Ct)
        if _Call:
            if _Call == "Pass" or _Call == 1 or _Call == "True" or (_Call == False):
                self.Ct_call.append(True)
            elif _Call == "Fail" or _Call == 0 or _Call == "False" or (not _Call):
                self.Ct_call.append(False)
        else:
            self.Ct_call.append(True) # otherwise assume all Ct's are good;

    def setSampleReference(self, a_refSample):
        self.sampleRef = a_refSample

    def setPrimerReference(self, a_refPrimer):
        self.primerRef = a_refPrimer

    def calc_pass1(self):
        # calculate the first delta and weed out the bad calls from Ct_call
        self.tp1 = []
        self.fp = []
        for i in xrange(len(self.Ct)):
            for t in xrange(len(self.primerRef.Ct)):
                if (self.Ct_call[i]) and (self.primerRef.Ct_call[t]):
                    self.tp1.append(self.Ct[i] - self.primerRef.Ct[t])
                    self.fp.append(self.Ct[i])

    def calc_pass2(self):
        # multiply and calculate all ddCts against one another;
        # do primer ref first;

        pow_ddCt = []

        # now normalise against the sample
        for i in self.tp1:
            for t in self.sampleRef.tp1:
                pow_ddCt.append(math.pow(2, -i) / math.pow(2, -t))

        pow_foldChange = []
        for i in self.fp:
            for t in self.sampleRef.fp:
                pow_foldChange.append(math.pow(2, -i) / math.pow(2, -t))

        self.log10ddCt = []
        self.log2ddCt = []

        # now the power;
        if pow_ddCt:
            for i in pow_ddCt:
                self.log10ddCt.append(math.log10(i))
                self.log2ddCt.append(math.log(i, 2))
        else:
            pass

        if len(self.log10ddCt) == 0: # no values calculated
            self.ddcts = []
            self.mean_log10ddCt = 0
            self.err_log10ddCt = 0
            self.mean_log2ddCt = 0
            self.err_log2ddCt = 0
            self.mean_percentddCt = 0
            self.err_percentddCt = 0
            self.mean_foldOverControl = 0
        else:
            self.ddcts = self.tp1
            self.mean_dCt = numpy.mean(self.tp1)
            self.mean_log10ddCt = numpy.mean(self.log10ddCt)
            self.err_log10ddCt = numpy.std(self.log10ddCt) / math.sqrt(len(self.Ct)) # divide by the sqrt of the original number of independent samples (ie 3)
            self.mean_log2ddCt = numpy.mean(self.log2ddCt)
            self.err_log2ddCt = numpy.std(self.log2ddCt) / math.sqrt(len(self.Ct)) # divide by the sqrt of the original number of independent samples (ie 3)
            self.mean_percentddCt = numpy.mean(pow_ddCt)*100
            self.err_percentddCt = numpy.std([x*100 for x in pow_ddCt]) / math.sqrt(len(self.Ct))
            self.mean_foldOverControl = numpy.mean(pow_foldChange)
            self.err_foldOverControl = numpy.std(pow_foldChange) / math.sqrt(len(self.Ct))


    # debug printers;
    def _pprint(self):
        self.sprint()
        if self.primerRef:
            print ">>>>>primerRef"
            self.primerRef.rprint()
        else:
            print ">>>>>primerRef\r\nNone"
        if self.sampleRef:
            print ">>>>>sampleRef"
            self.sampleRef.rprint()
        else:
            print ">>>>>sampleRef\r\nNone"
        print "\r\n"

    def _sprint(self):
        print self.sample, "\r\n",self.primer, "\r\n",  self.Ct, "\r\n", self.pow_ddCt, "\r\n", self.ddCt, "\r\n", self.mean, "\r\n",
    def _rprint(self):
        print self.sample, "\r\n",self.primer, "\r\n",  self.Ct, "\r\n", self.mean, "\r\n",

class realtime2:
    def __init__(self, filename=None, format=None, **kargs):
        """
        **Purpose**
            initialise and load in qRT-PCR data
        """
        # set up the lists;
        self.assayList = []
        self.primerList = []
        self.sampleList = []
        self.fileFormat = None
        self.cRefSampleName = None
        self.cRefPrimerName = None

        if format == "ABI_SDS":
            self.importSDSCopyPaste(filename)
        elif format == "ABI_SDS2.4":
            self.importSDSCopyPaste24(filename)
        elif format == "gene_row_sam_col":
            self.import_gene_row_sam_col(filename)

    def importSDSCopyPaste(self, filename):
        """
        **Purpose**
            import a SDS/ABI file that was copy and pasted
            from the view in the bottom left corner of the SDS software

        **Arguments**
            filename (Required)
                the filename to load
        """
        oh = open(filename, "rU")
        reader = csv.reader(oh, dialect=csv.excel_tab)

        for line in reader:
            try:
                sample_name = line[1] # catch other column missing wierdness
                probe_name = line[2]
                ct = float(line[4]) # catch undetermined's
            except:
                ct = None

            if sample_name and probe_name and ct: # only keep ones with value
                n = self.contains(probe_name, sample_name)
                if n:
                    n.addCt(ct, True)
                else: # make a new assay
                    a = assay(sample_name, probe_name) # not present so make a new element
                    a.addCt(ct, True)
                    self.assayList.append(a)

                # is this a new sample type?
                if sample_name not in self.sampleList:
                    self.sampleList.append(sample_name)
                # is this a new primer type?
                if probe_name not in self.primerList:
                    self.primerList.append(probe_name)
        oh.close()
        return(True)

    def importSDSCopyPaste24(self, filename):
        """
        **Purpose**
            import a SDS/ABI file that was copy and pasted from version 2.4
            from the view in the bottom left corner of the SDS software

        **Arguments**
            filename (Required)
                the filename to load
        """
        oh = open(filename, "rU")
        reader = csv.reader(oh, dialect=csv.excel_tab)

        for line in oh:
            line = line.split("\t") # csv borks on the unicode.
            try:
                sample_name = line[4] # catch other column missing wierdness
                probe_name = line[5]
                ct = float(line[7]) # catch undetermined's
            except Exception:
                ct = None

            if sample_name and probe_name and ct: # only keep ones with value
                n = self.contains(probe_name, sample_name)
                if n:
                    n.addCt(ct, True)
                else: # make a new assay
                    a = assay(sample_name, probe_name) # not present so make a new element
                    a.addCt(ct, True)
                    self.assayList.append(a)

                # is this a new sample type?
                if sample_name not in self.sampleList:
                    self.sampleList.append(sample_name)
                # is this a new primer type?
                if probe_name not in self.primerList:
                    self.primerList.append(probe_name)
        oh.close()
        return(True)

    def importBioMarkCSVFilter(self, _csvfile):
        fh = open(_csvfile, "rb")
        csvr = csv.reader(fh)
        # read the header row;
        # start at row 10 and read in the three lines
        for n in xrange(10):
            line = csvr.next()
        headers = []
        # then read in three lines and store them in a list of lists
        for n in xrange(2):
            line = csvr.next()
            headers.append(line)
        # squash the three lines into 1
        nheaders = []
        for n in xrange(len(headers[0])):
            a = headers[0]
            b = headers[1]
            nheaders.append(a[n]+b[n])

        #print nheaders

        for n in xrange(len(nheaders)): # find the column numbers for all of the elements;
            item = nheaders[n].lower()

            if item == "samplename": cSamName = n
            elif (item == "fam-mgbname") or (item == "fam-tamraname"): cPrimerName = n
            elif item == "ctvalue": cCt = n
            elif item == "ctcall": cCtCall = n

        # and now iterate through the rest of the file and load in the data;
        for row in csvr: # this continues from where the last read left off
            n = self.contains(row[cPrimerName], row[cSamName])
            if n:
                n.addCt(float(row[cCt]) , row[cCtCall]) # element already present, so just add the Ct
            else:
                a = assay(row[cSamName], row[cPrimerName]) # not present so make a new element
                a.addCt(float(row[cCt]), row[cCtCall])
                self.assayList.append(a)

            # is this a new sample type?
            if not self.sampleList.count(row[cSamName]):
                self.sampleList.append(row[cSamName])
            # is this a new primer type?
            if not self.primerList.count(row[cPrimerName]):
                self.primerList.append(row[cPrimerName])

        fh.close()

    def importGenericFilter(self, _csvfile):
        print "> Generic File Import"
        fh = open(_csvfile, "rb")
        csvr = csv.reader(fh)
        # read the header row;
        line = csvr.next()
        cCtCall = None
        for n in xrange(len(line)): # find the column numbers for all of the elements;
            item = line[n].lower()

            if item == "sample_name": cSamName = n
            elif item == "primer_name": cPrimerName = n
            elif item == "ct": cCt = n
            elif item == "ct_call": cCtCall = n

        fh.close()

        self.importByColumnNumber(_csvfile, 1, cSamName, cPrimerName, cCt, cCtCall)

    def importByColumnNumber(self, csvfile, _numRowsToSkip, _SamColNum, _PrimerColNum, _CtColNum, tsv=False, **kargs):
        """
        (Internal)
        This should probably be internal
        Supply the known column numbers for the particular categories.
        only CtCallColNum can be None
        """
        fh = open(csvfile, "rU")
        if tsv:
            csvr = csv.reader(fh, dialect=csv.excel_tab)
        else:
            csvr = csv.reader(fh)

        lineno = 0
        # and now iterate through the rest of the file and load in the data
        for row in csvr:
            lineno += 1
            if lineno > _numRowsToSkip:
                try:
                    entry = self.contains(row[_PrimerColNum], row[_SamColNum]) # check if pre-exisiting entry
                    if entry:
                        entry.addCt(float(row[_CtColNum]))
                    else: # not present so make a new element
                        Ct = float(row[_CtColNum]) # I do this first as if it fails the Ct is empty and so it will fail to add a superfluous entry.
                        a = assay(row[_SamColNum], row[_PrimerColNum])
                        a.addCt(float(row[_CtColNum]))
                        self.assayList.append(a)
                    
                    sample_name = row[_SamColNum]
                    probe_name = row[_PrimerColNum]
                    
                    # is this a new sample type?
                    if sample_name not in self.sampleList:
                        self.sampleList.append(sample_name)
                    # is this a new primer type?
                    if probe_name not in self.primerList:
                        self.primerList.append(probe_name)
                        
                except ValueError, IndexError:
                    # Ct column is empty and this is an empty row, skip it.
                    pass

        fh.close()

    def import_gene_row_sam_col(self, filename):
        """
        Import files where the gene is in the row and the columns are samples, e.g.:
        
        NA	MEF	D9_Oct4	d15_Oct4	d30_Oct4
        Actin	13.84773636	15.46720409	15.36977069	16.74893824
        Nestin	20.92989349	21.5210886	21.4988575	24.12586117
        pax6	31.56213379	30.98410511	27.82558918	27.57571411
        Sox1	29.75990963	24.99845028	24.59917641	27.45584869
        mT	32.67443275	31.92605591	25.98482418	30.10064888
        Eomes	29.95026398	29.28219509	25.44676971	26.60316944
        Lhx1	32.42358112	30.61198425	24.54306221	30.56797695
        
        """
        oh = open(filename, "rU")
        
        header = oh.readline().strip().split("\t")[1:]
        
        for line in oh:
            tt = line.strip().split("\t")
            probe_name = tt[0]
                       
            cts = [float(i) for i in tt[1:]]
            
            print probe_name, cts
            
            for i, ct in enumerate(cts):
                sample_name = header[i]
                print probe_name, ct, sample_name
                n = self.contains(probe_name, sample_name)
                if n:
                    n.addCt(ct, True)
                else: # make a new assay
                    a = assay(sample_name, probe_name) # not present so make a new element
                    a.addCt(ct, True)
                    self.assayList.append(a)
                                        
                    # is this a new sample type?
                    if sample_name not in self.sampleList:
                        self.sampleList.append(sample_name)
                    # is this a new primer type?
                    if probe_name not in self.primerList:
                        self.primerList.append(probe_name)
        oh.close()       

    def calculate(self):
        # these to be userset
        found_PrimerRef = False
        found_sampleRef = False

        self.refList = []
        # Am I a reference? If so add me to a new list
        for item in self.assayList:
            if item.primer == self.cRefPrimerName:
                self.refList.append(item)
                found_PrimerRef = True
                #print "Add Primer:", item.primer
            elif item.sample == self.cRefSampleName:
                self.refList.append(item)
                #print "Add Sample:", item.sample
                found_sampleRef = True

        if not found_sampleRef: # with the GUI, these *should* be impossible;
            print "> Error: No reference Sample found:", self.cRefSampleName
            return(False)
        if not found_PrimerRef:
            print "> Error: No reference Primer found: ", self.cRefPrimerName
            return(False)

        # now pair up the assays with their controls;
        for i in self.assayList:
            for t in self.refList:
                if not t == i:
                    if (t.primer == self.cRefPrimerName) and (t.sample == i.sample):
                        # the corresponding reference primer;
                        i.setPrimerReference(t)
                    if (t.primer == i.primer) and (t.sample == self.cRefSampleName):
                        i.setSampleReference(t)


        # I have to make a special case for the reference samples:
        # this assumes only one reference
        for t in self.refList:
            for i in self.refList:
                if not t == i:
                    if t.primer == self.cRefPrimerName:
                        # the reference is itself in this instance
                        t.setPrimerReference(t)
                    t.setSampleReference(t)

        for i in self.assayList:
            i.calc_pass1()
        # This is a bit weird, but I need to complete two loops.
        for i in self.assayList:
            i.calc_pass2() # requires two passes;

        config.log.info("Calculated Ct data")
        return(True)

    def contains(self, primerName, sampleName):
        for item in self.assayList:
            if primerName == item.primer:
                if sampleName == item.sample:
                    return(item)
        return(False)

    def write_output_data(self, csvfile):
        """
        **Purpose**
            write the output data to a csvfile.

        **Arguments**
            csvfile
                the filename to save the data to.

        **Returns**
            None and a csv file containing
        """
        try:
            self.assayList[0].mean_log10ddCt
        except:
            self.calculate()

        f = open(csvfile, "w")
        csvw = csv.writer(f, dialect=csv.excel_tab)

        line = ["Sample", "Primer", "ddCts", "mean_log10ddCt", "err_log10ddCt", "mean_log2ddCt", "err_log2ddCt", "mean_percentddCt", "err_percentddCt"]
        csvw.writerow(line)

        for i in self.assayList:
            line = [i.sample, i.primer, i.ddcts, i.mean_log10ddCt, i.err_log10ddCt, i.mean_log2ddCt, i.err_log2ddCt, i.mean_percentddCt, i.err_percentddCt]
            csvw.writerow(line)

        f.close()
        config.log.info("Saved file containing generated data '%s'" % csvfile)

    def detectFileFormat(self, _filename):
        """
        **Purpose**
            This will try to guess the file format based on a small
            set of criteria.
            To be honest this method does not work particularly well and
            its use is not recommended

        **Arguments**
            filename
                the filename to attempt to load

        **Returns**
            returns the type of file it thinks it finds.
        """
        ret = FF_UNKNOWN
        fh = open(_filename, "rb")

        # read the first few bytes from the file;
        s = fh.readline()
        if s.count("well") or s.count("Well"):
            ret = FF_GENERIC
        elif s.count("Chip Run Filename"):
            ret = FF_BIOMARK
        self.fileFormat = ret
        fh.close()

        if self.fileFormat == FF_GENERIC:
            self.importGenericFilter(_filename)
        elif self.fileFormat == FF_BIOMARK:
            self.importBioMarkCSVFilter(_filename)
        else:
            return(False)

        return(ret)

    def set_reference_probe(self, probe_name):
        """
        **Purpose**
            set the reference probe by probe_name

        **Arguments**
            probe_name
                the name of the probe, must be present in the realtime data

        **Returns**
            True
        """
        assert probe_name in self.primerList, "probe not found in this data set"
        self.cRefPrimerName = probe_name
        config.log.info("Set Probe Reference to '%s'" % probe_name)
        return(True)

    def set_reference_sample(self, sample_name):
        """
        **Purpose**
            set the reference sample by sample_name

        **Arguments**
            samples_name
                the name of the sample, must be present in the realtime data

        **Returns**
            True
        """
        assert sample_name in self.sampleList, "probe not found in this data set"
        self.cRefSampleName = sample_name
        config.log.info("Set Sample Reference to '%s'" % sample_name)
        return(True)

    def __str__(self):
        """
        (Override)
        print a nice summary of the data
        """
        ret = []
        ret.append("Realtime data:")
        ret.append("Primers: %s" % (", ".join(self.primerList)))
        ret.append("Samples: %s" % (", ".join(self.sampleList)))
        return("\n".join(ret))

    def print_all_data(self):
        """
        **Purpose**
            print all of the data contained in the realtime set

        **Arguments**
            None
        """
        print "Data:"
        print self.__str__()
        print
        for item in self.assayList:
            print item
