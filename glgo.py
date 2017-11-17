
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
    
    def barh(self, filename, ontology=None, maximum=None, pad=1e-105, key="qvalue", 
        label_fontsize=7, mark_significance=True, tight_layout=True, qvalue_limit=0.1, **kargs):
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
                
            tight_layout (Optional, default=True)
                use tight_layout to maximise space usage, but the location of the y axis 
                will shift around to make space for the labels and may not be in the same
                position across multiple graphs. 
                
                If tight_layout is False though then the lables may be truncated on the left hand edge.
                Use the 'position' argument to modify the position of the graph
                
        **Returns**
            None
        """
        self.sort(key)
        if ontology:
            data = self.get(key="ontology", value=ontology, mode="greedy") # get only for this ontology
        else:
            data = self.linearData # take all data
        
        bar_width = 0.5
        half_bar_width = bar_width / 2.0
        
        if not data: # Check that we got some entries
            config.log.warning("No items found for ontology '%s'" % ontology)
            # Don't throw an error just warn and ignore. 
            return(None)
        
        if maximum:
            data = data[0:maximum]
        data.reverse()
        
        #print data
        
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
            kargs["figsize"] = (8, 5)
        
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
    
        x = numpy.arange(0, len(vals))
        ax.barh(x, vals, bar_width, ec="none", color="grey")
        ax.set_yticklabels(labs)
        ax.set_yticks(x+half_bar_width)
        ax.set_frame_on(False)
        if ontology:
            ax.set_title(ontology)
        ax.set_ylim([-half_bar_width, len(vals)-1+bar_width+half_bar_width]) # num_items + bar_width
        ax.set_xlim([0, max(vals + [-math.log10(0.05)])])
        if mark_significance:
            ax.axvline(-math.log10(0.05), ls="-", color="black")
            ax.axvline(-math.log10(0.01), ls=":", color="grey")
        [t.set_fontsize(label_fontsize) for t in ax.get_yticklabels()]
        [t.set_fontsize(8) for t in ax.get_xticklabels()]
        ax.tick_params(top=False, bottom=True, left=False, right=False)
        #ax.axes.get_yaxis().set_visible(False)
        for line in ax.get_xticklines():
            line.set_marker(mpllines.TICKDOWN)
        #[item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]
        #[item.set_markeredgewidth(0.0) for item in ax.xaxis.get_ticklines()]
        
        # put back x axis:
        xs = ax.get_xaxis().get_view_interval()
        ymin = ax.get_yaxis().get_view_interval()[0]
        ax.add_artist(Line2D(xs, (ymin, ymin), color='black', linewidth=2))
        
        if tight_layout:
            fig.tight_layout()
        elif "position" in kargs and kargs["position"] and len(kargs["position"]) == 4:
            ax.set_position(kargs["position"])
        else:
            ax.set_position([0.7, 0.05, 0.25, 0.90]) # My favourite position. Uhm...
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("barh(): saved '%s'" % actual_filename)
        return(None)
        