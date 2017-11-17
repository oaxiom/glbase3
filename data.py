"""

**data.py**

data stuff, part of glbase

Contains various data that doens't belong in helpers.py, flags.py or opt.py.
And probably shouldn't be exported outside of glbase

These are constants, although they can be changed, in case, for example,
the value of pi changes.
"""

# This file must be clean of imports

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

ignorekeys = frozenset( # these are functional tags - so I should ignore them.
    ["dialect", 
    "duplicates_key",
    "skiplines", 
    "debug", 
    "special", 
    "skiptill", 
    "force_tsv",
    "gtf_decorators", 
    "endwith", 
    "__description__",
    "commentlines", 
    "keepifxin", 
    '__column_must_be_used',
    '__ignore_empty_columns'
    ]) 

typical_headers = frozenset(["chipseq_loc", "loc", "chr", "#",
"Gene Name", "", "GenBank", "RefSeq",
"Systematic", "mm8.refGene.chrom", "mm8", "loc", 'chromosome',
"mm9.refGene.chrom", "mm9",
"======================================================================", # stupid sissrs format garbage.
"=====================================================================", # stupid sissrs format garbage.
"======================================================================="] # stupid sissrs format garbage.
) # typical header labels;

positive_strand_labels = frozenset(["+", "1", "f", "F", 1])
negative_strand_labels = frozenset(["-", "0", "r", "R", -1, 0, "-1"])
