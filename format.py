"""

format specifiers. Part of glbase

This is now the approved way to get at the format specifiers:

format=format.bed

TODO:
Implement all the fileformats from::

http://genome.ucsc.edu/FAQ/FAQformat#format5.1

See below for the catalogue of file formats

"""

import csv, re, copy

import gzip as gzipfile
from .format_container import fc
from .helpers import lst_find

# ------------------- 'sniffer' formats - tell glbase to guess the file format.

default = {"sniffer": 0} # the default, it loads your file based on the heading labels in the csv.
sniffer = default # alternative name
sniffer_tsv = {"sniffer": 0, "force_tsv": True} # alternative name

# ------------------- standard formats:
bed = fc(name="bed",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",
        "name": 3, "score": 4, "strand": 5,
        "force_tsv": True, "skiplines": -1},
    description= "The 'normal' 6 column definition of a BED file, containing location, name, score and strand")

minimal_bed = fc(name="minimal_bed",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])", "force_tsv": True,
        "skiplines": -1},
    description= "A minimal BED file contains just the location data, in columns 0, 1 and 2.")

full_bed = fc(name="full_bed", # The full formal definition of a BED file.
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",
        "name": 3, "score": 4, "strand": 5, "thickStart": 6, "thickEnd": 7,
        "itemRgb": 8, "blockCount": 9, "blockSizes": 10, "blockStarts": 11,
        "force_tsv": True, "skiplines": -1, "commentlines": "#"},
    description= "This is the full 12-column definition of a BED file")

# Would be useful to add a "optional" and "required" for bed files?

bed_no_strand = minimal_bed # Old name, kept for compatability

psl = fc(name="psl",
    format={"skiplines": -1, "force_tsv": True,
        "matches": 0, "misMatches": 1, "repMatches": 2, "nCount": 3,
        "qNumInsert": 4, "qBaseInsert": 5, "tNumInsert": 6, "tBaseInsert": 7,
        "strand": 8, "qName": 9, "qSize": 10, "qStart": 11, "qEnd": 12,
        "tName": 13, "tSize": 14, "tStart": 15, "tEnd": 16,"blockCount": 17,
        "blockSizes": 18, "qStarts": 19, "tStarts": 20},
    description="PSL files, as produced by BLAT")

fasta = fc(name="fasta",
    description="FASTA sequence file",
    format={"special": "fasta"})

gtf = fc(name="gtf",
    description="GTF, gene transfer file format.",
    format={"feature_type": 1, "feature": 2, "gtf_decorators": 8, "commentlines": "#",
        "loc": "location(chr=column[0], left=column[3], right=column[4])",
        "strand": 6, "skiplines": -1, "force_tsv": True})

snp = fc(name="snp",
    format=dict(bin=0, name=4, score=5, strand=6, refNCBI=7, refUCSC=8, observed=9,
        molType=10, cclass=11, valid=12, avHet=13, avHetSE=14, func=15, locType=16, weight=17,
        loc={"code": "location(chr=column[1], left=column[2], right=column[3])"}, force_tsv=True),
    description="snp.txt file format reader")

pgsnp = fc(name="pgsnp",
    format={"loc": "location(chrom=column[0], left=column[1], right=column[2]",
        "name": 3, "alleleCount": 4, "alleleFreq": 5,
        "alleleScores": 6, "skiplines": -1, "forse_tsv": True},
    description="Personal Genome SNP file format")

encode_rna_expn = fc("encode_rna_expn",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",
        "name": 3, "score": 4, "strand": 5,
        "level": 6, "signif": 7, "score2": 8,
        "force_tsv": True, "skiplines": -1},
    description="ENCODE RNA elements: BED6 + 3 scores format")

encode_broadpeak = fc("encode_broadpeak",
    format={"loc": "location(chr=column[0], left=column[1], right=column[2])",
        "name": 3, "score": 4, "strand": 5,
        "signalValue": 6, "pValue": 7, "qValue": 8,
        "force_tsv": True, "skiplines": -1},
    description="ENCODE broadPeak: Broad Peaks (or Regions) format")

# --------------------- motif discovery formats
fimo_out = fc(name="fimo_out",
    description="Load in the fimo.txt file output by FIMO, part of the MEME suite.",
    format={"force_tsv": True, "pattern-name": 0, "name": 1, "score": 4,
        "p-value": 5, "q-value": 6, "sequence": 7,
        "loc": "location(chr=column[1], left=min(column[2], column[3]), right=max(column[2], column[3]))"})

homer_known_motifs = fc(name="homer_known_motifs",
    description='HOMER knownMotifs.txt format',
    format={'name': 0,
        'motf_seq': 1,
        'q': 4,
        'p': 2,
        '-log10P': 3,
        "force_tsv": True,
        })

# --------------------- peak-discovery tool file-types
macs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 5,
    "skiptill": "chr", "force_tsv": True, "fold_change": 7}

macs_summit = {"loc": {"code": "location(chr=column[0], left=int(column[1])+int(column[4]), right=int(column[1])+int(column[4]))"},
    "fold_change": 7, "tag_count": 5, "fdr": 8,
    "skiptill": "chr",
    "force_tsv": True}

ccat_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"}, "force_tsv": True,
    "tag_height": 4, "fold": 6, "skiplines": -1}

ccat_output_csv = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"},
    "tag_height": 4, "fold": 6, "fold_change": {"code": "barSplitter(3)[2]"}, "skiplines": -1}

# load in SISSRS file.
sissrs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 3,
    "fold": 4, "p-value": 5, "force_tsv": True, "skiplines": 57}

homer_peaks = fc(name="homer_peaks",
    format={"loc": {"code": "location(chr=column[1], left=column[2], right=column[3])"},
    "force_tsv": True, "commentlines": "#", "tagcount": 5},
    description="file format output by HOMER findPeaks")

dfilter_bed = fc(name="dfilter_bed",
    format={"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
        "p_value": 5, "tag_count": 6, "force_tsv": True, "skiplines": 1},
    description="file format output by DFilter")

# --------------------- next-gen sequencing formats:
# This is not a full implementation of the sam file specification.
# but it will accept a sam file as outputted by tophat.
sam_tophat = fc(name="sam_tophat",
    format={"name": 0, "loc": {"code": "location(chr=column[2], left=column[3], right=column[3])"},
        "mapq": 3, "seq": 9, "force_tsv": True, "skiplines": -1},
    description="SAM file. Note this only implements SAM files as output by tophat")

# lst_find() is in helpers.py
sam_tophat_xs = fc(name="sam_tophat_xs",
    format={"name": 0, "loc": {"code": "location(chr=column[2], left=column[3], right=column[3])"},
        "mapq": 4, "seq": 9, "force_tsv": True, "skiplines": -1, "keepifxin": "XS:",
        "strand": {"code": "re.sub(r'[A-Z]*[:]*', '', column[lst_find(column, lambda x: 'XS:' in x)])"}},
    description="Tophat SAM file. This uses the XS tag for the strand, for strand-specific RNA-seq")

sam = fc(name="sam",
    format={"qname": 0, "flag": 1, "cigar": 5,
        "loc": {"code": "location(chr=column[2], left=column[3], right=column[3])"},
        "mapq": 3, "seq": 9, "force_tsv": True, "skiplines": -1, "commentlines": "@"},
    description="SAM file. This is a partial implementation")

exporttxt_loc_only = fc(name="exporttxt_loc_only",
    format={"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "force_tsv": True}, # export.txt file (output from the illumina pipeline), but only loads the location and strand.
    description="Load in the export.txt file as output by the Illumina HTS pipeline.\n\t\tThis variant loads only the read location")

exporttxt_all = fc(name="exporttxt_all",
    format={"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13, "seq": 6, "quality_score": 7,
        "force_tsv": True}, # export.txt file (output from the illumina pipeline), but only loads the location and strand.
    description="Load in the export.txt file as output by the Illumina HTS pipeline.\n\t\tThis variant loads the location, sequence, strand and quality score")

bowtie = fc(name="bowtie",
    format={"loc": "location(chr=column[2], left=column[3], right=column[3])", "strand": 1,
        "seq": 4, "quality_scores": 5, "force_tsv": True},
    description="Loads most of the relevant data from the default bowtie output")

bowtie_loc_only = fc(name="bowtie_loc_only",
    format={"loc": "location(chr=column[2], left=column[3], right=column[3])", "strand": 1,
        "force_tsv": True},
    description="Loads the genomic location and strand from the default bowtie output")

# --------------------- miscellaneous file formats
ann_list = {"loc": 4, "strand": 6, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1}# output format from my annotation script
peak_loc = {"loc": 4, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1} # A list of peak locations, annotated?
genespring_out = {"refseq": 4, "entrez": 8, "fold_change": 1, "array_systematic_name": 0, "force_tsv": True}# A list out of GeneSpring
array_simplified = {"refseq": 10, "array_systematic_name": 0, "entrez": 8}
illumina_anotations = {"array_systematic_name":0, "refseq": 3, "entrez": 7}

# --------------------- snpXXX.txt file
snp_txt = dict(bin=0, name=4, score=5, strand=6, refNCBI=7, refUCSC=8, observed=9,
    molType=10, cclass=11, valid=12, avHet=13, avHetSE=14, func=15, locType=16, weight=17,
    loc={"code": "location(chr=column[1], left=column[2], right=column[3])"}, force_tsv=True)

ncbi_gwas = fc(name="ncbi_gwas",
    format={'force_tsv': True, 'loc': 'location(chr=column[11], left=int(column[12]), right=int(column[12])+1)',
    'snp_name': 20, 'associated_gene': 13, '__column_must_be_used': 11},
    description="Export from the NCBI GWAS Catalog. Will only load data that has a valid genomic location")

# --------------------- microarray-like file formats
array_tsv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False,
        "force_tsv": True}

array_csv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False}

mm8_refgene = fc(name="mm8_refgene",
    format={"loc": {"code": "location(chr=column[0], left=column[2], right=column[3])"},
        "strand": 1, "name": 4, "description": 5, "dialect" : csv.excel_tab,
        "refseq": 6, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"},
        },
    description="The mm8 refGene table downloaded from UCSC Genome Browser")

# --------------------- UCSC/Ensembl annotations
mm9_refgene = fc(name="mm9_refgene",
    format={"loc": "location(chr=column[2], left=column[4], right=column[5])",
        "strand": 3, "name": 12, "force_tsv": True,
        "refseq": 1, "tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8
        }, # description is lost from the mm9 table?
    description="The mm9 refGene table downloaded from UCSC Genome Browser")

mm10_refgene = fc(name="mm10_refgene",
    format={"loc": "location(chr=column[2], left=column[4], right=column[5])",
        "strand": 3, "name": 12, "force_tsv": True,
        "refseq": 1, "tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8
        }, # description is lost from the mm9 table?
    description="The mm10 refGene table downloaded from UCSC Genome Browser")

ensembl = {"loc": {"code": "location(chr=column[2], left=column[4], right=column[5])"},
    "tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
    "ensmbl": 1, "name": 12, "exon_count": 8,
    "force_tsv": True, "skiplines": 0}

hg19_refgene = fc(name="hg19_refgene",
    format = {
        "loc": {"code": "location(chr=column[2], left=int(column[4]), right=int(column[5]))"},
        "strand": 3,
        "name": 12,
        "refseq": 1,
        #"tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"}, # Broken in this version of glbase
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8,
        "exonStarts": {"code": "[int(x) for x in column[9].strip(\",\").split(\",\")]"},
        "exonEnds": {"code": "[int(x) for x in column[10].strip(\",\").split(\",\")]"},
        "force_tsv" : True},
    description="The hg19 refGene table downloaded from UCSC Genome Browser")

# hg18 default refseq export.
hg18_refseq = fc(name="hg18_refseq",
    format={"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
        "strand": 5, "force_tsv": True,
        "refseq": 3, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}},
    description="The Hg18 RefSeq gene table as downloaded from the UCSC Genome Browser")

# --------------------- miscellaneous
homer_annotated = {"loc": "location(chr=column[0], left=column[1], right=column[2])",
    "peak_id":3,"motif_score":4,"strand": 5,"motif_seq":6,"summit_dist":7,
    "summit_height":8,"nearest_gene":9,"TSS_dist":10,"annotation":11,
    "force_tsv": True,"skiplines": -1}

MACS_combined = {"loc": "location(chr=column[0], left=column[1], right=column[2])","peak_id":3,
    "summit_height":4,"p-value_score": 5,"number_tags":6,"fold_enrichment":7,"FDR":8, "force_tsv": True,
    "skiplines": -1}

# --------------------- RNA-seq data

# 1.2.0 version:
rsem_gene = fc(name="rsem_gene",
    format={"ensg": 0, "tpm": 8, "tpm_lo": 10, "tpm_hi": 11,
    "transcripts": 1, "force_tsv": True},
    description="The *.gene.* file format output by RSEM >=1.2.0")

# --------------------- BLAST and related formats

blast_tabular = fc(name="blast_tabular",
    format={'queryId': 0, 'subjectId': 1, 'percIdentity': 2, 'alnLength': 3,
        'mismatchCount': 4, 'gapOpenCount': 5, 'queryStart': 6, 'queryEnd': 7,
        'subjectStart': 8, 'subjectEnd': 9, 'eVal': 10, 'bitScore': 11,
        "force_tsv": True, "skiplines": -1},
    description="The default BLAST tabular format output")

# --------------------- Interproscan and HMMER

hmmer_tbl = fc(name="hmmer_tbl",
    description="hmmsearch/hmmscan --tbl_out loader",
    format={"special": "hmmer_tbl"})


def _load_hmmer_tbl(filename, gzip=False):
    """
    # Unbelievably stupid format for hmmer:
    Load the hmmer tbl_out table
    """
    oh = gzipfile.open(filename, "rt") if gzip else open(filename, "rt")
    res = []
    for line in oh:
        if "#" not in line:
            ll = line.split()

            name = ll[18:]
            gene = "NA"
            # split the name up into k:v pairs
            for item in name:
                if ":" in item and "gene:" in item:
                    gene = item.split(":")[1]
            try:
                res.append({"peptide": ll[0], "dom_acc": ll[3], "dom_name": ll[2], "score": float(ll[4]),
                    "gene": gene})
            except ValueError:
                res.append({"peptide": ll[0], "dom_acc": ll[3], "dom_name": ll[2], "score": ll[4],
                    "gene": gene})
    return res

hmmer_domtbl = fc(name="hmmer_domtbl",
    description="hmmsearch/hmmscan --domtbl_out loader",
    format={"special": "hmmer_domtbl"})

def _load_hmmer_domtbl(filename, gzip=False):
    """
    # Irritating format for hmmer:
    Load the hmmer domtblout format table
    """
    oh = gzipfile.open(filename, "rt") if gzip else open(filename, "rt")
    res = []
    for line in oh:
        if "#" in line:
            continue

        ll = line.split()

        row = {"peptide": ll[0],
            "dom_acc": ll[4], "dom_name": ll[3],
            'dom_loc': (int(ll[17]), int(ll[18])),
            'e': float(ll[6]),
            'tlen': int(ll[2]),
            'qlen': int(ll[5]),
            "gene": ll[22],
            'name': ll[22]}

        #try:
        #    row[
        #except ValueError:

        res.append(row)
    return res


# --------------------- GO outputs

go_GREAT_shown = fc(name="go_GREAT_shown",
    description="GO shown-* tables from GREAT",
    format={"commentlines": "#", "name": 0,  "qvalue": 3, "force_tsv": True}
    )

go_DAVID = fc(name="go_DAVID",
    description="GO table saved from DAVID 'functional Chart' view",
    format={"force_tsv": True, "ontology": 0, 'count': 2, "name": 1, "qvalue": 11, "pvalue": 4}
    )

go_GOseq = fc(name="go_GOSeq",
    description="GO table by GOSeq",
    format={"force_tsv": True, "ontology": 7, 'count': 4, "name": {'code': 'cat_columns(column[1], column[6])'}, "qvalue": 2, "pvalue": 2} # no q-value?
    )

# --------------------- class container

class fccatalogue():
    def __init__(self, formats):
        # put in dict by name

        self.formats = {}
        for item in formats:
            self.formats[item.name] = item

    def __str__(self):
        print("Found %s formats" % len(self.formats))
        keys = list(self.formats.keys()) # report in alphabetical order
        keys.sort()
        a = [
            "format.{name:24} - {desc:5}".format(
                name=k, desc=self.formats[k].description
            )
            for k in keys
        ]

        return("\n".join(a))

    def __iter__(self):
        for k in self.formats:
            yield self.formats[k]

    def __len__(self):
        return(len(self.formats))

    def find(self, value):
        """
        **Purpose**
            find a particular motif format on value

        **Arguments**
            value
                searches through the list of formats and returns possible matches, searches the name
                and description to find relevant formats.
        """

        a = [
            "format.{name:22} - {desc:6}".format(
                name=key, desc=self.formats[key].description
            )
            for key in self.formats
            if value.lower() in self.formats[key].name.lower()
            or value in self.formats[key].description.lower()
        ]

        if a:
            print("\n".join(a))
        else:
            print("None found")

catalogue = fccatalogue([
    fimo_out, homer_known_motifs,
    fasta, gtf, bed, full_bed, minimal_bed,
    exporttxt_loc_only, exporttxt_all,
    mm8_refgene, mm9_refgene, mm10_refgene, hg18_refseq, hg19_refgene,
    snp, pgsnp, psl, encode_rna_expn, encode_broadpeak,
    sam_tophat, sam_tophat_xs, sam,
    blast_tabular,
    rsem_gene, bowtie_loc_only,
    bowtie, homer_peaks,
    hmmer_tbl, hmmer_domtbl, #HMMER
    dfilter_bed,
    # GO lists:
    go_GREAT_shown, go_DAVID,
    ncbi_gwas,
    ])

if __name__ == "__main__":
    print(catalogue)
    print()
    print("Find:")
    catalogue.find("bed")

