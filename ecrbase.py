
from . import config

import sys, os

from .genelist import genelist as genelist_object
from .location import location

class tfbs_iter:
    def __init__(self, filename):
        self.oh = open(filename)

    def __iter__(self):
        """
        This will output a tfbs_iter object
        that behaves as a list of dicts
        Each yield produces a dict with these keys:
            loc - the location in genomic coordinates
            motifs - all of the motifs
            motif_pos - a dict of the form {"tf_name": "chr1:100-200"}
        """
        tfbs = {}
        loc = None
        for line in self.oh:
            if line[0] == ">": # moved onto a new loc.
                motifs = [i.lstrip("  ").split(" ")[0] for i in tfbs]
                pos = [i.lstrip("  ").split(" ")[-1:] for i in tfbs]

                # need to process position offsets here:
                newpos = {}

                for i, m in enumerate(motifs):
                    pp = abs(int(pos[i][0]))
                    if int(pos[i][0]) > 0:
                        newpos[m] = location(chr=loc["chr"], left=loc["left"] + pp, right=loc["left"] + pp+8) # don't know length of pwm
                    else:
                        newpos[m] = location(chr=loc["chr"], left=loc["left"] + pp-8, right=loc["left"] + pp) # don't know length of pwm

                if loc:
                    yield {"loc": loc, "motifs": motifs, "motif_pos": newpos}

                # and now set the new loc.
                loc = location(loc=line.replace("> ", ""))
                tfbs = []
            else:
                tfbs.append(line.replace("\n", "").replace("\r", ""))

        self.oh.close()

        yield {"loc": loc, "motifs": [i.lstrip("  ").split(" ")[0] for i in tfbs], "motif_pos": newpos} # release a final set
        return # raise StopIteration


class ecrbase:
    def __init__(self, path, cache=False, **kargs):
        """
        **Purpose**
            a few quick utils to put into glbase.
            Make it nice and glbase-like at a later date.

        **Arguments**
            path
                path to the relevant ecrbase files

            cache (Optional, default=False)
                cache the ecrbase data
                This can be very slow depending upon which
                file you are trying to cache into memory.
                On computers <4 Gb of memory this may well fail
                and fail pretty badly too.

                Files are only cached when they are asked for.
        """
        assert os.path.exists(path)
        self.cache = cache
        self.path = path
        self.__cached_items = {}
        config.log.info("Setup ecrbase object")

    def __repr__(self):
        return("<glbase.ecrbase>")

    def __str__(self):
        return(self.__repr__())

    def __cache_alignment_file(self, filename):
        """
        (Internal)

        This will cache an alignment file into memory
        """
        if filename in self.__cached_items:
            return(True)

        config.log.info("Caching '%s'..." % filename)
        oh = open(os.path.join(self.path, filename), "rU")
        store = []
        for line in oh:
            # process the ecrbase file:
            t = line.replace("\n", "").split("-vs- ")
            tt = t[0].split(" ")
            store.append(location(loc=tt[0])) # the mouse conserved region.
        oh.close()
        config.log.info("Finished Caching...")

        self.__cached_items[filename] = store

    def __cache_tfbs_file(self, erc_tfbs_filename):
        """
        (Internal)

        This will cache a tfbs file
        """
        oh = open(erc_tfbs_filename, "rU")

        config.log.info("Caching TFBS File '%s'" % erc_tfbs_filename)
        tfbs = []
        done = 0
        loc = None
        for line in oh:
            if line[0] == ">": # moved onto a new loc.
                motifs = [i.lstrip("  ").split(" ")[0] for i in tfbs]
                #pos = [i.lstrip("  ").split(" ")[-1:] for i in tfbs]

                if loc:
                    tfbs_file.append({"loc": loc, "motifs": motifs})

                # and now set the new loc.
                loc = location(loc=line.replace("> ", ""))
                tfbs = []
            else:
                tfbs.append(line.replace("\n", "").replace("\r", ""))
        oh.close()
        config.log.info("Finished Caching...")
        self.__cached_items[filename] = tfbs_file
        return(tfbs_file)

    def map_ecr_locs(self, genelist, ecrbase_organism_pair, key_name="in_conserved_block?"):
        """
        Adds a new key, with "conserved" yes or no
        """
        if self.cache:
            self.__cache_item("ecrs.%s.txt" % ecrbase_organism_pair)
        else:
            raise NotImplementedError

        for rr, region in enumerate(self.__cached_items["ecrs.%s.txt" % ecrbase_organism_pair]):
            for ii, gene in enumerate(genelist):
                if not key_name in gene:
                    gene[key_name] = "no"

                if gene["window"].qcollide(region):
                    gene[key_name] = "yes"
        return(genelist)

    def map_ecr_tfbs_locs(self, genelist, ecrbase_organism_pair, motif_key="motif",
        version="102", key_name="in_conserved_block?", bed_filename=None, **kargs):
        """
        **Purpose**
            Retrieves all the pwms in the genelist
            the genelist must have a valid 'loc' key

        **Arguments**
            genelist (Required)
                A genelist-like object or valid iterable

            ecrbase_organism_pair (Required)
                the species names to compare. For example:
                "mm9monDom4"
                "mm9hg18"
                etc.

            version (Optional, default=102)
                the version of the tfbs file to use
                The current version is 102 (5/7/2010)

            motif_key

            key_name

            bed_filename (Optional)
                save the TFBSs as a bed file that can be uploaded to the UCSC
                genome browser.

        **Returns**
            A new genelist with a new key "motif" which contains a list
            for the pwms and a second new key "locs" containing the genomic
            coordinates
        """
        if self.cache:
            store = self.__cache_tfbs_file("ecrs.%s.txt" % ecrbase_organism_pair)
        else:
            store = tfbs_iter(os.path.join(self.path, "tfbs_ecrs.%s.v%s.txt" % (ecrbase_organism_pair, version)))

        res = genelist_object()

        config.log.info("Searching for TFBS in '%s'" % ("tfbs_ecrs.%s.v%s.txt" % (ecrbase_organism_pair, version)))
        for rr, item in enumerate(store):
            if item["motifs"]: # only search if there are motifs:
                for ii, gene in enumerate(genelist):
                    if gene["loc"].qcollide(item["loc"]):
                        res.append(item) # THe slow one here is probably okay for the vast majority of applications.

        if bed_filename:
            name = "TFBS %s" % ecrbase_organism_pair
            desc = "TFBS from ecrbase file %s" % "tfbs_ecrs.%s.v%s.txt" % (ecrbase_organism_pair, version)
            oh = open(bed_filename, "w")
            oh.write("track name='%s' description='%s' visibility=2 color=0,128,255, itemRgb='On'\n" % (name, desc))
            for r in res:
                for tf in r["motif_pos"]:
                    oh.write("chr%s\t%s\t%s\t%s\t0\n" % (r["motif_pos"][tf]["chr"], r["motif_pos"][tf]["left"], r["motif_pos"][tf]["right"], str(tf)))
            oh.close()
            config.log.info("Saved bed to '%s'" % bed_filename)

        config.log.info("Finished...")
        return(res)


