'''

hic.py

Analysis for HiC data.

'''

import pickle, numpy, math

import matplotlib.cm as cm
from scipy import ndimage
from scipy.signal import argrelextrema

from operator import itemgetter
from . import config
from .draw import draw
from .progress import progressbar
from .location import location

if config.H5PY_AVAIL:
    import h5py
else:
    raise AssertionError('Asked for a hic, but h5py is not avaialble')

# I'm not using the base_genelist class, so, you need to add in defs as needed.
class hic:
    def __init__(self, filename=None, name='', new=False, inter_chrom_only=True):
        """
        **Purpose**
            store for HiC data.
            
        **Arguments**
            filename (Required)
                hic is an on-disk array. So, you need to provide 
                a filename to save the data into.
                
                I suggest a '.hiccy' extension, so as not to confuse with JuiceBox's .hic
        
            new (Optional, default=False)
                set this to True if you want to make a new hiccy
        
            name (Optional)
                a name to use as a label, etc
                
            inter_chrom_only (Required, default=True)
                At the moment only inter-chromsomal contacts are recorded.
                
        """
        self.readonly = True
        # These need to also be packaged into the hiccy
        if new:
            self.readonly = False
            #try:
            self.hdf5_handle = h5py.File(filename, 'w')  # This should be x to stop overwriting an exisiting file
            
            # h5py hierarchy (Made at load_time)
            # Things stored:
            # str name
            # matrix # A numpy style matrix
            # chrom_edges slice locations for all the chromsome names
            self.hdf5_handle.attrs['name'] = name
            self.hdf5_handle.attrs['inter_chrom_only'] = inter_chrom_only
            self.all_chrom_names = [] # Made into a set later
            self.draw = draw()              
        else: # old
            self.hdf5_handle = h5py.File(filename, 'r')
            # TODO: fetch it back out as a set:
                    
            dat = self.hdf5_handle['all_chrom_names'].value 
            dat = [i[0] for i in dat]# flatten
            self.all_chrom_names = [n.decode("ascii", "ignore") for n in dat]

            # Rescue bin_lookup from disk:
            # save the bin data as an emulated dict:
            self.bin_lookup_by_binID = {}
            self.bin_lookup_by_chrom = {}
            for chrom in self.all_chrom_names:
                flat_bin = self.hdf5_handle['bin_lookup/chrom_%s/bins' % chrom].value 
                self.bin_lookup_by_chrom[chrom] = [(row[0], row[1], row[2], row[3]) for row in flat_bin]
                
                for row in flat_bin:
                    self.bin_lookup_by_binID[row[0]] = (chrom, row[1], row[2], row[3]) # row[3] == newBinID

            #print(self.bin_data)

            self.mats = {}
            for chrom in self.all_chrom_names:
                self.mats[chrom] = self.hdf5_handle['matrix_%s/mat' % chrom]
                
            self.draw = draw()
            config.log.info('Bound "%s" Hic file' % filename)
        return(None)
    
    def visit(self):
        self.hdf5_handle.visit(print)

    def __len__(self):
        return(self.num_bins)

    def __getitem__(self, index):
        """
        Confers:
        
        a = expn["condition_name"]
        
        and inherits normal genelist slicing behaviour
        """
        return(self.hdf5_handle.attrs[index]) 

    def keys(self): 
        return(self.hdf5_handle.keys()) 

    def _optimiseData(self):
        pass
    
    def close(self):
        self.hdf5_handle.flush()
        self.hdf5_handle.close()

    def __find_binID_spans(self, loc):
        """
        (Internal)
        
        Fing the binIDs for the span loc. 
        
        Returns a tuple containing the local binIDs for these coordiantes, 
        and the new updated genomic coordinates taken from the binspans
        
        """
        assert self.readonly, 'To __find_binID_spans, this must be read-only, set new=False' 
        
        # These are in global bins.
        binLeft = 1e50
        binRight = -1
        locLeft = 0
        locRight = 0        
        mostLeft = 1e50
        mostRight = -1
        
        for bin in self.bin_lookup_by_chrom[loc['chr']]:
            # bin = (binID, left, right)
            if bin[2] > loc['left'] and bin[0] < binLeft:
                binLeft = bin[0]
                locLeft = bin[1]
            if loc['right'] > bin[1] and bin[0] > binRight:
                binRight = bin[0] 
                locRight = bin[2]
            
            # To work out the local bin positions: 
            if bin[0] < mostLeft:
                mostLeft = bin[0]
            if bin[0] > mostRight:
                mostRight = bin[0]
                mostRightLoc = bin[2]
        
        # in case locations are asked for off the edge of the chromsome:
        if binLeft < 0:
            binLeft = 0
            locLeft = 0
        if binRight > mostRight:
            binRight = mostRight
            locLeft = mostRightLoc
        
        newloc = location(chr=loc['chr'], left=locLeft+1, right=locRight-1) # the +1/-1 stops the bins from overlapping by 1bp, and solves a lot of problems
        
        binLeft = binLeft - mostLeft
        binRight = binRight - mostLeft
        mostRight = self.bin_lookup_by_binID[mostRight][3] # Convert to newBinID:
        mostLeft = self.bin_lookup_by_binID[mostLeft][3]
        assert binRight-binLeft > 2, 'the genome view (loc) is too small, and contains < 2 bins'
        
        return(binLeft, binRight, newloc, mostLeft, mostRight)
        
    def load_hicpro_matrix(self, matrix_filename, bed_file):
        """
        **Purpose**
            Load a .matrix file output from HiC-Pro
            
            bin1, bin2, norm_score

        **Arguments**
            filename (Required)
                the Filename to load.
        """
        assert not self.readonly, 'To load, this must be read-only, set new=True'
        
        # First go through the bed file and get all the bins, chroms and sizes.
        bins = []
        self.all_chrom_names = []
        oh = open(bed_file, 'r')
        for lin in oh:
            # format = chr, left, right, bin#
            lin = lin.strip().split('\t')
            
            chr = lin[0].replace('chr', '')
            try:
                chr = int(chr)
            except ValueError:
                pass # It's chrX etc.
                
            bins.append((chr, int(lin[1]), int(lin[2]), int(lin[3])))
            if chr not in self.all_chrom_names:
                self.all_chrom_names.append(chr)
                
        oh.close()        
        self.all_chrom_names = set(self.all_chrom_names)
        self.hdf5_handle.attrs['num_bins'] = len(bins)
        
        config.log.info('Found %s bins' % self['num_bins'])
        # the bins are not sorted, 
        bin_order = sorted(bins, key=lambda v: (isinstance(v[0], str), v[0], v[1])) # Sort chroms by str, then int
        bin_lookup = {} # python 3: no need OrderedDict
        
        for newbinID, bin in enumerate(bin_order):
            # Basically, bin_lookup[oldbinID] = (chrom, left, right, newbinID)
            bin_lookup[bin[3]] = (bin[0], bin[1], bin[2], newbinID) # left, right, bin#

        # Get the edges of the chorms, and so the matrix sizes
        self.chrom_edges = {}
        for chrom in self.all_chrom_names:
            self.chrom_edges[chrom] = [1e20,-1]
            
            for bin in bin_lookup:
                if bin_lookup[bin][0] != chrom:
                    continue
                    
                # fill in:
                if bin_lookup[bin][3] < self.chrom_edges[chrom][0]:
                    self.chrom_edges[chrom][0] = bin_lookup[bin][3]
                if bin_lookup[bin][3] > self.chrom_edges[chrom][1]:
                    self.chrom_edges[chrom][1] = bin_lookup[bin][3]
        # Have to then convert them all to the new bins:
        # Although not strictly necessary in the matrix-based, I preserve this fiddly step so that
        # when I implement the inter-chrom interactions it is easier.
        matrices = {}
        bins = {} # The bin -> matrix convertor
        for chrom in self.all_chrom_names:            
            # Now I know how big the arrays need to be, and all of the chromosomes:
            size =  self.chrom_edges[chrom][1]-self.chrom_edges[chrom][0]
            matrices[chrom] = numpy.zeros((size+1, size+1), dtype='float32')
            #bins[k] = 
        # TODO: Support inter-chromosomal contacts with a sparse array:
        
        oh = open(matrix_filename, 'r')
        p = progressbar(self['num_bins'] * self['num_bins']) # expected maximum
        for idx, lin in enumerate(oh):
            lin = lin.strip().split('\t')
            # First check the two bins are on the same chromosome:
            bin1 = bin_lookup[int(lin[0])]
            bin2 = bin_lookup[int(lin[1])]
            
            if bin1[0] == bin2[0]:
                chrom = bin1[0]
                chrom_edge = self.chrom_edges[chrom][0]
                x = bin1[3]-chrom_edge
                y = bin2[3]-chrom_edge
                
                matrices[chrom][x,y] = float(lin[2])
                matrices[chrom][y,x] = float(lin[2])
            else:
                # TODO: Support for interchrom with a sparse array:
                pass
                
            p.update(idx)
        p.update(self['num_bins'] * self['num_bins']) # the above is unlikely to make it to 100%, so fix the progressbar.
        oh.close()

        # save the bin data as an emulated dict:
        grp = self.hdf5_handle.create_group('bin_lookup')
        for chrom in self.all_chrom_names: 
            grp = self.hdf5_handle.create_group('bin_lookup/chrom_%s' % chrom)
            flat_bin = []
            
            for oldbinID in bin_lookup:
                if bin_lookup[oldbinID][0] == chrom:
                    flat_bin.append([oldbinID, bin_lookup[oldbinID][1], bin_lookup[oldbinID][2], bin_lookup[oldbinID][3]]) # i.e. oldId, left, right, newID
            flat_bin = numpy.array(flat_bin)
            grp.create_dataset('bins', (flat_bin.shape), dtype=int, data=flat_bin)
        
            #self.hdf5_handle.create_dataset('bin_lookup', (len(to_store), 1), dtype='S10', data=to_store)
        
        # A sparse array would be better. The matrix is only around ~40% full in 100k, and
        # likely this gets much lower as the resolution drops...
        config.log.info('Loaded matrix %s*%s' % (self['num_bins'], self['num_bins']))
        for chrom in self.all_chrom_names:
            config.log.info('Added chrom=%s to table' % chrom)
            grp = self.hdf5_handle.create_group('matrix_%s' % chrom)
            grp.create_dataset('mat', matrices[chrom].shape, dtype=numpy.float32, data=matrices[chrom])
        config.log.info('Saved all matrices to hdf5')
        
        dat = [str(n).encode("ascii", "ignore") for n in self.all_chrom_names]
        self.hdf5_handle.create_dataset('all_chrom_names', (len(self.all_chrom_names), 1), 'S10', dat)
                
        return(True)
        
    def heatmap(self, filename, chr=None, loc=None,
        bracket=None, colour_map=cm.inferno_r, **kargs):
        """
        **Purpose** 
            Draw an interaction heatmap
            
        **Arguments**
            filename (Required)
                filename to save the image to
                
            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr, 
                or you can set a specific location with loc
            
            aspect (Optional, default='square')
                image aspect for the heatmap
        
            colbar_label (Optional, default='log2(Density)')
                Label for the colour bar.
            
            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.
                
            bracket (Optional, default=None)
                clip the data within the ranges [low, high]
                
        **Returns**
            None, and a file in filename.               
            
        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')
        
        if chr:            
            data = self.mats[str(chr).replace('chr', '')]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)
            
            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, _, _ = self.__find_binID_spans(loc)
            
            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]
        else:
            data = self.matrix # use the whole lot
        
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs: 
            kargs["colbar_label"] = "log2(Density)"
                        
        fig = self.draw.getfigure(**kargs)
        
        # positions of the items in the plot:
        heatmap_location  = [0.05,   0.01,   0.90,   0.90]
        scalebar_location = [0.05,   0.97,   0.90,   0.02]

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            data[data == -numpy.inf] = 0
            vmin = data.min()
            vmax = data.max()

        # ---------------- (heatmap) -----------------------
        ax3 = fig.add_subplot(121)
        
        ax3.set_position(heatmap_location) # must be done early for imshow
        hm = ax3.imshow(data, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]], 
            interpolation=config.get_interpolation_mode()) 

        #ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)
        ax3.set_xlim([0,data.shape[1]])
        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")
        
        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]        
        
        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("heatmap: Saved '%s'" % actual_filename)       
    
    def tri_plot(self, filename, chr=None, loc=None, bracket=None, colour_map=cm.inferno_r, 
        **kargs):
        """
        **Purpose**
            Draw the side-on half triangle plot. 
            
        **Arguments**
            filename (Required)
                filename to save the image to.

            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr, 
                or you can set a specific location with loc
                
            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]
                
        **Returns**
            None and a file in filename.
        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')
        
        if chr:            
            data = self.mats[str(chr).replace('chr', '')]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)
            
            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, _ = self.__find_binID_spans(loc)
                        
            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]
        else:
            data = self.matrix # use the whole lot

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            vmin = data.min()
            vmax = data.max()

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs: 
            kargs["colbar_label"] = "log2(Density)"

        heatmap_location =  [0.05,   0.01,   0.90,   0.80]
        scalebar_location = [0.05,  0.97,   0.90,   0.02]
            
        if 'figsize' not in kargs:
            kargs['figsize'] = (8,3)
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        dst = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False, cval=0)
        
        ax.set_position(heatmap_location)
        hm = ax.imshow(dst, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]], 
            interpolation=config.get_interpolation_mode()) 
        
        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([data.shape[1]/2,(data.shape[1]/2)+data.shape[1]/10]) # get the middle start position, then about 1/10 of the way up
        #print([dst.shape[1]/2,data.shape[1]])
        #ax.set_yticklabels("")
        #ax.set_xticklabels("")
        
        ax.yaxis.tick_right()
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]     
            
        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("tri_plot: Saved '%s'" % actual_filename)         
        
        return()

    def tri_plot_and_freq_plot(self, filename, expn=None, expn_cond_name=None, 
        chr=None, loc=None,
        bracket=None, colour_map=cm.inferno_r, 
        **kargs):
        """
        **Purpose**
            Draw the side-on half triangle plot. and a frequency plot below
            
        **Arguments**
            filename (Required)
                filename to save the image to.
                
            expn (Required)
                an expression object, with conditions to use in the frequency plot.
                
            expn_cond_name

            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr, 
                or you can set a specific location with loc
                
            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]
                
        **Returns**
            None and a file in filename.
        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')
        
        if chr:            
            data = self.mats[str(chr).replace('chr', '')]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)
            
            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, mostLeft, mostRight = self.__find_binID_spans(loc)            
            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]

            # I just assume the bins match
            this_chrom = [0] * (localRight-localLeft+1)       
            cindex = expn.getConditionNames().index(expn_cond_name)
            for i in expn:
                if i['loc']['chr'] == loc['chr']:
                    # Take the old BinID and convert it to the new binID:
                    # If inside this part of the chromosome:
                    if loc.qcollide(i['loc']):
                        local_bin_num = (self.bin_lookup_by_binID[i['bin#']][3] - mostLeft) - localLeft
                        this_chrom[local_bin_num] = i['conditions'][cindex]
            plot_y = this_chrom
            plot_x = numpy.arange(0, len(plot_y))
            
        else: 
            raise(AssertionError, 'Not implemented here! :(')
        '''
        if chr:
            # get a specific chromosome:
            bin_top = self.chrom_edges[chr.replace('chr', '')][0]
            bin_bot = self.chrom_edges[chr.replace('chr', '')][1]
            
            data = self.matrix[bin_top:bin_bot, bin_top:bin_bot]
            scores, borders = self.score_insulation(data) # Supposed to be done on non-log transformed
            
            this_chrom = [0] * (bin_bot - bin_top + 2) 
            
            # I just assume the bins match
            cindex = expn.getConditionNames().index(expn_cond_name)
            for i in expn:
                if i['loc']['chr'] == chr.replace('chr', ''):
                    local_bin_num = i['bin#'] - bin_top                     
                    this_chrom[local_bin_num] = i['conditions'][cindex]
            plot_y = this_chrom
            plot_x = numpy.arange(0, len(plot_y))
        else:
            data = self.matrix # use the whole lot
            # can't be done at the moment as I resorted the bins
        '''
    
        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            with numpy.errstate(divide='ignore'): 
                data = numpy.log2(data)
            vmin = data.min()
            vmax = data.max()

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs: 
            kargs["colbar_label"] = "log2(Density)"

        scalebar_location = [0.05,  0.97,   0.90,   0.02]
        heatmap_location  = [0.05,  0.50,   0.90,   0.40]
        tadtri_location   = [0.05,  0.40,   0.90,   0.10]
        tadborder_location= [0.05,  0.35,   0.90,   0.05] # Bashes into;
        freq_location     = [0.05,  0.01,   0.90,   0.28]

        if 'figsize' not in kargs:
            kargs['figsize'] = (8,5)
            
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(512)
        
        dst = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False, cval=0)
        
        ax.set_position(heatmap_location)
        hm = ax.imshow(dst, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]], 
            interpolation=config.get_interpolation_mode()) 
        
        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([data.shape[1]/2,(data.shape[1]+data.shape[1]/10)]) # the /10 determines how much it gets stretched.
        ax.set_yticklabels("")
        ax.set_xticklabels("")
        
        ax.yaxis.tick_right()
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(511)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        t = cb.set_label(kargs["colbar_label"], fontsize=6)
        #t.set_fontsize(6)
        cb.ax.tick_params(labelsize=5)
        #[label.set_fontsize(5) for label in ax0.get_xticklabels()]     
        '''
        # TAD boundaries
        ax = fig.add_subplot(513)
        ax.set_position(tadborder_location)
        bar = [0] * (len(plot_y)+1)
        for b in borders:
            bar[b] = 1
            
        ax.imshow(numpy.array((bar,bar)), cmap=cm.Greys, vmin=0, vmax=1, aspect="auto",
            origin='lower', extent=[0, len(bar), 0, 2], 
            interpolation=config.get_interpolation_mode()) 
            
        #ax.plot(plot_x, plot_y)     
        ax.set_xlim([0, len(plot_y)])
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")
        
        # TAD triangle plot
        ax1 = fig.add_subplot(514)
        ax1.set_position(tadtri_location)    
        # zip up border_edges:
        bo = zip(borders[:-1], borders[1:-1])        
        xx = []
        yy = []
        for b in bo:
            midH = (b[1] - b[0])/2
            mid = midH + b[0]
            
            xx.append(b[0])
            yy.append(0)
            
            xx.append(mid)
            yy.append(midH)
        # need to fill in rightmost point:
        xx.append(borders[-1])
        yy.append(0)
        
        ax1.plot(xx, yy)
        ax1.tick_params(top=False, bottom=False, left=False, right=False)
        ax1.set_xlim([0, len(plot_y)])
        ax1.set_yticklabels("")
        ax1.set_xticklabels("")
        '''
        # freq plot
        ax1 = fig.add_subplot(515)
        ax1.set_position(freq_location)
        #ax1.set_frame_on(False)
        # get a chrom if required:
        ax1.plot(plot_x, plot_y)     
        ax1.set_xlim([0, len(plot_y)])
        ax1.tick_params(top=False, bottom=False, left=True, right=False)
        ax1.set_xticklabels("")
        
        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("tri_plot: Saved '%s'" % actual_filename)         
        
        return()
        
    def density_plot(self, filename, vmin=0, vmax=50000, **kargs):
        """
        **Purpose**
            Draw a density histogram of the matrix scores.
            
        **Arguments**
            filename (Required)
                filenmame to save the image to
                
        **Returns**
            None and image in filename
        """
        assert self.readonly, 'must be readonly to draw a density_plot. set new=False'
        assert filename, "density_plot: you must specify a filename"       
        
        fig = self.draw.getfigure(**kargs)
        ax1 = fig.add_subplot(121) # All
        ax2 = fig.add_subplot(122) # inter-chromosome
        
        max = 0
        for i in self.matrix: # this will return each row
            if i.max() > max: 
                max = i.max()
                
        print(max)
        data = [0] * (int(max)+1)
        for r in self.matrix:
            for i in r:
                if int(i) > 0: # There is a lot of 0
                    data[int(i)] += 1
                
        ax1.plot(numpy.arange(0, len(data)), data)
        
        l = numpy.log2(numpy.array(data))
        ax2.plot(numpy.arange(0, len(data)), l)
        
        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("density_plot: Saved '%s'" % actual_filename)         
        
        return()
        