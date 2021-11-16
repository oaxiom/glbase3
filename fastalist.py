"""

a fasta file of sequences.

a gene list where the sequence is stored under "fasta" key


fasta list may be obselete?

"""

import sys, os

from . import utils

from .genelist import genelist

class fastalist(genelist):
    def __init__(self, **kargs):
        genelist.__init__(self)

        if "filename" in kargs:
            self.importFASTA(kargs["filename"])

    def importFASTA(self, filename=None, **kargs):
        """
        import a fasta file.
        """
        if (not filename) or (not os.path.exists(filename)):
            print("Error: fasta filename is not valid")
            return(False)

        try:
            self.linearData = utils.convertFASTAtoDict()
        except:
            print("Error: Not a valid FASTA file")
            return(False)
        # sequence is stored under the f and r tags...
        # see if I can make some tags from the FASTA header:
        self._optimiseData()
        return(True)

    def exportFASTA(self, filename=None, **kargs):
        assert filename
        raise NotImplementedError

    def scanMotif(self, motif_word=None, **kargs):
        """
        get the freq of motifs in motif_word.

        must be a degenerate word, using the IUPAC naming system.

        returns a list with the number of hits per fasta
        """
        if not motif_word:
            print("Error: scanMotif() requires a motif word")
            return(False)

        try:
            rel = [regex_dict[bp] for bp in self.seq]
            regex = re.compile(string.join(rel, ""))
        except:
            print("Error: your 'motif_word' is not valid")
            return(False)

        res = []
        for item in self.linearData:
            for seq in [item["f"], item["r"]]:
                matches = regex.findall(seq)
                res.append(len(matches))

        # draw figures for this run



