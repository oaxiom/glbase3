"""
Framework for GO analysis in glbase
"""

import sys, os, math
import numpy
import matplotlib.pyplot as plot
from matplotlib.lines import Line2D
import matplotlib.lines as mpllines

from . import config
from .genelist import genelist
from .location import location
from .draw import draw

# Available formats
# go_DAVID
# go_GREAT_shown
# go_GOseq

# Bind the GO data into glbase?
# If so make the data available here:

class glgo(genelist):
    def __init__(self, filename=None,  **kargs):
        assert "loadable_list" not in kargs, "glgo() does not support loadable_lists"
        assert "format" in kargs and kargs["format"], "glgo requries a 'format' argument"
        assert "qvalue" in kargs["format"] or 'pvalue' in kargs["format"], "the format must contain a 'qvalue' or 'pvalue' key"

        genelist.__init__(self, filename=filename, **kargs)

        # Tidy up '~' symbols in GO terms from typical DAVID/R output
        if 'name' in self.linearData[0]:
            for item in self.linearData:
                if not isinstance(item['name'], location): # weird but true. Occasionally appears as an enriched term
                    item['name'] = item['name'].replace("~", ":")
            self._optimiseData()

    def bind_GO_set(self, filename):
        """
        **Purpose**
            bind a GO set (or list) of data into glgo and make it available for testing

        **Arguments**
            filename
                the filename for the GO set.
        """
        pass

    def barh(self,
        filename:str,
        ontology=None,
        maximum=None,
        pad:float = 1e-105,
        color:str = 'lightgrey',
        key:str ="qvalue",
        label_fontsize:int = 6,
        mark_significance:bool = True,
        qvalue_limit:float = 0.1,
        bot_pad:float = 0.1,
        **kargs):
        """
        **Purpose**
            Draw a barchart for some ontology

        **Arguments**
            filename (Required)
                THe filename to save the image to

            ontology (Optional)
                The name of the ontology class to use. If not specified then use all of the data
                in the list.

            maximum (Optional)
                the maximum number of GO items to plot. If not specified then use all of the items
                in the list

            qvalue_limit (Optional, default=0.1)
                The significance limit for this data

            color (Optional, default='lightgrey')
                color for thew bars

            pad (Optional, default=1e-105)
                By default glgo will -log10 transform the pvalues/qvalues.
                Occasionally some GO discovery programs will output a p/q-value
                of 0 or NA. Pad allows these values to be used. But also sets a maximum
                size of the bar chart.

                Warning:
                If your GO program has higher precision than pad (1E-105) then items with a p/q-value
                of 0 will be set to pad and so may appear in the middle of your list and
                not at the top, as you would expect.

            label_fontsize (Optional, default=&)
                The font size of the labels

            mark_significance (Optional, default=True)
                mark lines of significance at 0.05 and 0.01 p-value/q-value.

        **Returns**
            None
        """
        self.sort(key)

        if ontology:
            data = self.get(key="ontology", value=ontology, mode="greedy") # get only for this ontology
        else:
            data = self.linearData # take all data

        bar_width = 0.7
        half_bar_width = bar_width / 2.0

        if not data: # Check that we got some entries
            config.log.warning(f"No items found for ontology '{ontology}'")
            return None # Don't throw an error just warn and ignore.

        if maximum:
            data = data[0:maximum]
        data.reverse()

        vals = [i[key] for i in data] # data ends up unfolded.
        labs = [i["name"] for i in data]

        if qvalue_limit:
            newl = []
            newv = []
            for l, v in zip(labs, vals):
                if v < qvalue_limit:
                    newl.append(l)
                    newv.append(v)
            labs = newl
            vals = [-math.log10(i+pad) for i in newv]
            print(labs, vals)
        else:
            vals = [-math.log10(i[key]+pad) for i in data] # now padded, etc

        if "figsize" not in kargs:
            kargs["figsize"] = (7, 6) # Can only accommodate a set number of results;

        mmrelative_hei = 0.1+(0.022*len(data)) # 0.022 is space per item;

        fig = self.draw.getfigure(**kargs)
        fig.subplots_adjust(left=0.82, right=0.98, top=mmrelative_hei, bottom=bot_pad)
        ax = fig.add_subplot(111)

        if "position" in kargs and kargs["position"] and len(kargs["position"]) == 4:
            ax.set_position(kargs["position"])

        x = numpy.arange(0, len(vals))
        ax.barh(x, vals, bar_width, ec="none", color=color)
        ax.set_yticklabels(labs)
        ax.set_yticks(x) # This is fixed in new matplotlib
        ax.set_frame_on(False)

        if ontology:
            ax.set_title(ontology)
        ax.set_ylim([-bar_width, len(vals)-1+bar_width]) # num_items + bar_width
        ax.set_xlim([0, max(vals + [-math.log10(0.05)])])
        if mark_significance:
            #ax.axvline(-math.log10(0.05), ls="-", color="black")
            ax.axvline(-math.log10(0.01), ls=":", color="red", lw=0.5)

        [t.set_fontsize(label_fontsize) for t in ax.get_yticklabels()]
        [t.set_fontsize(label_fontsize) for t in ax.get_xticklabels()]

        ax.tick_params(top=False, bottom=True, left=False, right=False)

        ax.tick_params(axis='both', size=6)
        ax.set_xlabel(key)


        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info(f"barh(): saved '{actual_filename}'")

        return None
