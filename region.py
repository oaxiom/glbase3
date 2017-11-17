
import sys, os, csv, re
from . import config, utils

import matplotlib.pyplot as plot

from .genelist import genelist
from .flags import *
from .errors import AssertionError
from .location import location
from .base_genelist import _base_genelist
from .draw import draw

class region(_base_genelist):
    def __init__(self, genome=None, sequence=None, loc=None, name=None, **kargs):
        """
        **Purpose**
            a class that contains data and descriptions of discreet genomic
            regions.

        **Arguments**
            sequence
                the sequence data

            loc (Required)
                the genomic coordinates of the region

            assembly (Optional)
                the genome assembly identifier (e.g. mm8, Hg18, monDom4, NCBI37 etc)

        """
        assert loc, "a genomic location must be provided to region, got '%s'" % loc

        self.name = name
        if not self.name: # no valid name, guess one.
            self.name = "region-%s_%s_%s" % (loc["chr"], loc["left"], loc["right"])

        self.sequence = sequence
        self.loc = location(loc=loc)
        self._draw = draw(self)

        if "assembly" in kargs and kargs["assembly"]:
            self.assembly = kargs["assembly"]

        if genome: # if genome sent, try to get some annotations:
            self.sequence = genome.getSequence(loc=self.loc)
            self.sequence_lower = self.sequence.lower()
            self.features = genome.getFeatures(loc=self.loc)
            self.assembly = genome.name
            self.__genome = genome

    def __repr__(self):
        return("<base.region>")

    def mark_primer_pairs(self, primer_pair_list, filename=None, bed_filename=None, **kargs):
        """
        **Purpose**
            Draw a figure containing the locations of the primer pairs

        **Arguments**
            primer_pair_list
                a genelist containing primer pair sets (sequence), in the form:
                {"left": "acac...acac", "right": "acc...cac", "name": "PrimerA"}

            filename
                The name to save the image file to.

            bed_filename (Optional, default=None)
                If a valid filename, output a bed format file for loading into UCSC
                browser.

        **Returns and Effects**
            Returns None

            Writes to the log any primer pairs not found
            Will throw an error if more than one match for the primer pair
                is found.
            Saves a bed file if bed_filename is valid. This contains
            the locations of the primers and their name. This can be uploaded
            directly to the UCSC genome browser.
            Saves an image file if
        """
        assert self.sequence, "This region does not have sequence data"
        assert filename and self.features, "This region does not have any gene features - the image cannot be drawn"

        primers = []
        for primer in primer_pair_list:
            primer_details = {"name": primer["name"]}

            left = re.compile(primer["left"], re.IGNORECASE)
            # find matches:
            left_matches = left.finditer(self.sequence)
            if left_matches:
                for i, m in enumerate(left_matches):
                    primer_details["left"] = (m.start(), m.end())
                    if i >= 1:
                        raise AssertionError("Found more than one matching sequence! Discard primer: %s, (%s)" % (primer["name"], primer["left"]))

            right = re.compile(utils.rc(primer["right"]), re.IGNORECASE)
            #right = re.compile(primer["right"], re.IGNORECASE)
            right_matches = right.finditer(self.sequence)
            if right_matches:
                for i, m in enumerate(right_matches):
                    primer_details["right"] = (m.start(), m.end())
                    if i >= 1:
                        raise AssertionError("Found more than one matching sequence! Discard primer: %s, (%s)" % (primer["name"], primer["left"]))
            primers.append(primer_details)

        fig = plot.figure(dpi=config.DEFAULT_DPI, figsize=(14,3))
        #position_genomic = [0.02, 0.05, 0.96, 0.1]
        ax1 = fig.add_subplot(111)
        ax1.set_position([0.02, 0.02, 0.96, 0.96])
        #ax1.set_position(position_genomic)
        self._draw._genome_segment(ax1, self.loc, self.features)
        for p in primers:
            ax1.text(p["left"][0], 5.6, p["name"], size=6, rotation="vertical")
            ax1.arrow(p["left"][0], 5.5, p["right"][1] - p["left"][0], 0, alpha=1, width=0.02)

        if bed_filename:
            oh = open(bed_filename, "w")
            oh.write("track name='%s' description='%s' visibility=2 color=0,128,255, itemRgb='On'\n" % (self.name, self.name))
            for p in primers:
                oh.write("chr%s\t%s\t%s\t%s\t0\n" % (self.loc["chr"], self.loc["left"] + p["left"][0], self.loc["left"] + p["right"][1], p["name"]))
            oh.close()

        actual_filename = self._draw._saveFigure(fig, filename)
        config.log.info("Saved '%s' image file" % actual_filename)
        return(None)
