"""

Mass Spectrometry processing genelist-like object

"""


import sys, os, gzip, csv
from statistics import mean
import numpy

import matplotlib.colors as colors
from collections.abc import Iterable

from . import config, utils
from .base_expression import base_expression
from .draw import draw
from .progress import progressbar
from .errors import AssertionError, ArgumentError
from .utils import fold_change

class massspec(base_expression):
    supported_formats = {'maxquant_txt'}
    valid_call_modes = {'coip'}
    valid_datasets = {'intensities', 'call', 'fold_change'}

    def __init__(self, 
        filenames=None, 
        format=None, 
        gzip:bool = False, 
        mq_use_lfq:bool = False,
        name = None
        ):
        '''
        **Purpose**
            Handler and processer for Mass Spec data
            
            Note that replicate sample names will be merged by taking the average 
        
        **Arguments**
            filenames (Required)
                list of filenames to process
                
            format (Required)
                One of {}
                
            gzip (Optional)
                Input file is gzipped
                
            name (Optional)
                name for the masspec object
            
            # Options for MaxQuant:
            mq_use_lfq (Optional, default=False)
                MaxQuant gives two intensity scores, LFQ intensity estimates the amount of protein based
                on a pool of proteins that don't change much. Whilst intensity is the raw intensity.
                LFQ is preferred for whole cell MS, whilst without is probably better for 
                enrichment-based MS, e.g. Co-IP. Use this switch to use lfq or not.
            
        '''.format(self.supported_formats)
        assert isinstance(filenames, list), 'filenames must be a list'
        assert format in self.supported_formats
        
        self.name = name
        
        if format == 'maxquant_txt':
            peps, sample_names = self.__maxquant_txt(filenames, gzip, mq_use_lfq=mq_use_lfq)
        else:
            raise NotImplementedError
        
        # Average samples with the same names;
        new_linear_data = []
        for peptide, vals in peps.items():
            new_peptide = {'name': peptide, 'PIDs': vals['pids']}
            intensities = []
            peptide_counts = []
            for s in sample_names:
                intensities.append(mean(vals['intensities'][s]) if vals['intensities'][s] else 0)    
                peptide_counts.append(mean(vals['peptide_counts'][s]) if vals['peptide_counts'][s] else 0)
                
            new_peptide['intensities'] = intensities
            new_peptide['peptide_counts'] = peptide_counts
            new_peptide['call'] = [None] * len(peptide_counts)
            new_peptide['fold_change'] = ['-'] * len(peptide_counts)
            new_linear_data.append(new_peptide)
        
        # Pack to a genelist-type object
        
        self._conditions = sample_names
        self.linearData = new_linear_data
        self._optimiseData()
        
        self.draw = draw()
        
        self.called = False # call is valid;
        self.fold_change = False
        config.log.info('Loaded {filenames}')

    def __str__(self):
        out = base_expression.__str__(self)
        out += f'\nCalled: {self.called}\nFold change done: {self.fold_change}'
        return out

    def __repr__(self):
        return "<massspec class>"

    '''
    def sort(self): # To implement
        raise NotImplementedError
    '''

    def find(self): # To implement
        raise NotImplementedError

    def from_pandas(self): # Not to implement
        raise NotImplementedError

    def _optimiseData(self):
        '''
        Needs a custom one
        '''
        # unpack the intensity data based on the order in peps;
        self._intensities = [i['intensities'] for i in self.linearData]
        self._peptide_counts = [i['peptide_counts'] for i in self.linearData]
        self._call = [i['call'] for i in self.linearData]
        self._call = [i['fold_change'] for i in self.linearData]
        
        self.qkeyfind = {}
        for index, item in enumerate(self.linearData):
            for key in item:
                if key not in self.qkeyfind:
                    self.qkeyfind[key] = {}

                try:
                    if item[key] not in self.qkeyfind[key]:
                        self.qkeyfind[key][item[key]] = []
                    self.qkeyfind[key][item[key]].append(index)
                except TypeError:
                    # The item in unhashable and cannot be added to the qkeyfind
                    # This should be pretty rare if not impossible.
                    pass

    def saveCSV(self, filename=None, interleave_errors=True, no_header=False, **kargs):
        """
        A CSV version of saveTSV(), see saveTSV() for syntax
        """
        self.saveTSV(filename=filename, tsv=False, interleave_errors=True, no_header=False, **kargs)
        config.log.info(f"saveCSV: Saved '{filename}'")
        
        return None

    def saveTSV(self, filename=None, tsv=True, **kargs):
        """
        (Override)
        **Purpose**
            Save the massspec data as a tsv file
            This is a little different from the normal genelist.saveTSV()
            as I want to make certain that the various numeric data 
            is written in a sensible manner at
            the end of the TSV.
            
            In general TSVs are not designed for loading back into glbase3, please use the 
            massspec.save() and glload() for persistent objects.
            
        **Arguments**
            filename
                The filename (with a valid path) to save the file to.

        **Returns**
            returns None
        """
        self._save_TSV_CSV(filename=filename, tsv=tsv, no_header=False, **kargs)
        config.log.info(f"saveTSV: Saved '{filename}'")
        
        return None

    def _save_TSV_CSV(self, filename=None, tsv=True, no_header=False, **kargs):
        """
        Internal unified saveCSV/TSV for expression objects
        """        
        valig_args = ["filename", "tsv", "key_order", "no_header"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.saveCSV, k)

        assert filename, "you must specify a filename"

        with open(os.path.realpath(filename), "wt") as oh:
            writer = csv.writer(oh, dialect=csv.excel_tab) if tsv else csv.writer(oh)
            array_data_keys_to_skip = ("intensities", "peptide_counts", "call", 'fold_change')
            all_keys = list(self.keys())
            write_keys = []
            if "key_order" in kargs:
                write_keys = kargs["key_order"]
                # now add in any missing keys to the right side of the list:
                write_keys += [k for k in all_keys if k not in write_keys and k not in array_data_keys_to_skip]
            else:
                # just select them all:
                write_keys = [k for k in all_keys if k not in array_data_keys]
                        
            if not no_header:
                title_row = [k for k in write_keys if k in all_keys]
                title_row += [f'Intensity {i}/M' for i in self.getConditionNames()]
                title_row += [f'Peptide_counts {i}' for i in self.getConditionNames()]
                if self.called: title_row += [f'Call {i}' for i in self.getConditionNames()]
                if self.fold_change: title_row += [f'FC {i}' for i in self.getConditionNames()]
                writer.writerow(title_row)

            for data in self.linearData:
                line = [data[k] for k in write_keys if k in data]
                line += data['intensities'] + data['peptide_counts']
                if self.called: line += [str(b) for b in data['call']]
                if self.fold_change: line += [str(b) for b in data['fold_change']]
                writer.writerow(line)
                
        return None
        
    def __maxquant_txt(self, filenames, gzip, mq_use_lfq):
        peps = {}
        sample_names = []
    
        if mq_use_lfq:
            intensity_search_key = 'Intensity '
        else:
            intensity_search_key = 'LFQ intensity '

        # scan all the filenames to get all the sample_names;
        for filename in filenames:
            if gzip:
                oh1 = gzip.open(filename, 'rt')
            else:
                oh1 = open(filename, 'rt')
            
            for hit in oh1:
                if 'Protein IDs' in hit:
                    hit = hit.strip().split('\t')
                    for lab in hit:
                        if 'Razor + unique peptides ' in lab:
                            sample_names.append(lab.replace('Razor + unique peptides ', ''))
                    break
            oh1.close()
    
        for filename in filenames:
            if gzip:
                oh1 = gzip.open(filename, 'rt')
            else:
                oh1 = open(filename, 'rt')
            
            # need to get the sample_names:
            
            for hit in oh1:
                if 'Protein IDs' in hit:
                    # gets all the keys;
                    intensity_keys = {}
                    hit = hit.strip().split('\t')
                    for idx, lab in enumerate(hit):
                        if lab.startswith(intensity_search_key): 
                            intensity_keys[lab.replace(intensity_search_key, '')] = idx
                            
                    unique_count_keys = {}
                    for idx, lab in enumerate(hit):
                        if 'Razor + unique peptides ' in lab:
                            unique_count_keys[lab.replace('Razor + unique peptides ', '')] = idx
                    continue
                    
                hit = hit.strip().split('\t')
                pep_name = hit[6]
                pids = hit[1]

                if not hit[6]:
                    continue

                if pep_name not in peps:
                    peps[pep_name] = {
                        'pids': pids,
                        'intensities': {k: [] for k in sample_names},
                        'peptide_counts': {k: [] for k in sample_names},
                        'call': {k: None for k in sample_names},
                        'fold_change': {k: None for k in sample_names},
                        }
                
                for sample_name in sample_names:
                    if sample_name not in intensity_keys:
                        # From a previous filename
                        continue 
                    intensity = hit[intensity_keys[sample_name]]
                    peps[pep_name]['intensities'][sample_name].append(float(intensity)/1e6)
                    unq_peps = hit[unique_count_keys[sample_name]]
                    peps[pep_name]['peptide_counts'][sample_name].append(int(unq_peps))
            oh1.close()
                    
        return peps, sample_names
        
    def filter(self, 
        mode, 
        minimum_unique_peptides,
        minimum_intensity
        ):
        '''
        **Purpose**
            Call a peptide, and filter out likely false positives based on some scheme and
            thresholds
            
            NOTE: This is an INPLACE method that will remove not-called
            matches
            
        **Purpose**
            mode (Required)
                one of {}
                
                schemes are:
                    coip: Strategy for coip MS, filters out typical contaminants found in CoIPMS
                    
            minimum_unique_peptides (Required)
                minimum unique peptides in any one condition to keep peptide.
                
            minimum_intensity (Required)
                minimum intensity (in /1e6) required in any one condition
                
        **Returns**
            number of peptides deleted
                    
        '''.format(self.valid_call_modes)
        assert mode in self.valid_call_modes, f'mode "{mode}" not found in {self.valid_call_modes}'
        
        starting_len = len(self)
        
        if mode == 'coip':
            newl = self.__filter_coipms(minimum_unique_peptides=minimum_unique_peptides, minimum_intensity=minimum_intensity)
        
        if not newl:
            raise AssertionError('All peptides have been filtered out, list is empty')
        
        self.linearData = newl
        self._optimiseData()
        
        config.log.info(f'filter: Trimmed {starting_len - len(self)} possible false positives and peptides not meeting thresholds')
        config.log.info(f'filter: List now {len(self)} peptides long')
        return starting_len - len(self)
        
    def __filter_coipms(self, minimum_unique_peptides, minimum_intensity):
        newl = []
    
        for pep in self.linearData:
            pep_name = pep['name']
            # Simple filter for common contaminants
            if 'Rpl' in pep_name: continue
            elif 'Rps' in pep_name: continue
            elif 'Tub' in pep_name: continue
            elif 'Gapdh' in pep_name: continue
            elif 'Act' in pep_name: continue
            elif 'Myh' in pep_name: continue
            elif 'Ighg' in pep_name: continue
            elif 'Iglv' in pep_name: continue
            elif 'Col' in pep_name: continue
            elif 'Golga' in pep_name: continue
            elif 'Kif' in pep_name: continue
            elif 'Myl' in pep_name: continue
            elif 'Krt' in pep_name: continue
            elif 'Eif' in pep_name: continue
            elif 'Vim' in pep_name: continue
            elif 'Atp' in pep_name: continue
            elif 'Igkv' in pep_name: continue 
            elif 'Ighv' in pep_name: continue 
        
            if True not in [i>minimum_unique_peptides for i in pep['peptide_counts']]:
                continue
 
            if True not in [i>minimum_intensity for i in pep['intensities']]:
                continue
  
            newl.append(pep)
        
        return newl

    def call(self, expt_scheme, intensity_fold_change_threshold=1.0, cull_empty_calls=True):
        '''
        **Purpose**
            Call a peptide in a specific condition.
            based on an enrichment threshold.
            This method is primarily aimed at Co-IP MS. For whole proteome MS it's 
            better advised to use something like DESeq2 or equivalent.
            
            NOTE: This is an INPLACE method that modifies the underlying list
            
        **Arguments**
            expt_scheme (Required)
                This calls a peptide based on enrichment metric 
                derived from intensity, basically the fold-change over the control. 
                Consequently it needs to know the matching control to perform the 
                threshold comparison. 
                
                The experimental scheme should be in this form:
                
                {
                'ctrl_condition_name1': ['condition1', 'condition2' ... 'conditionN'],
                'ctrl_condition_name2': ['condition1', 'condition2' ... 'conditionN'],
                ...
                'ctrl_condition_nameN': ['condition1', 'condition2' ... 'conditionN']
                }
            
            intensity_fold_change_threshold (Optional, default=1.0)
                fold-change increase intensity from ctrl to target required to call
                a peptide.
            
            cull_empty_calls (Optional, default=True)
                If not called in any condition, then delete the peptide.
            
        **Returns**
            None
        '''
        # check all the conditions are actually present
        for ctrl in expt_scheme:
            assert ctrl in self._conditions, f'Control data {ctrl} not found in this massspec data'
            for expt in expt_scheme[ctrl]:
                assert expt in self._conditions, f'{expt} not found in this massspec data'
        
        cond_index_lookup = {c: i for i, c in enumerate(self._conditions)}

        for pep in self.linearData:
            pep_name = pep['name']

            for ctrl in expt_scheme:
                ctrl_index = cond_index_lookup[ctrl]
                
                for ip in expt_scheme[ctrl]:
                    ip_index = cond_index_lookup[ip]
                    fc = fold_change(pep['intensities'][ip_index], pep['intensities'][ctrl_index], 0.1)
                    pep['fold_change'][ip_index] = fc
                    pep['fold_change'][ctrl_index] = 0.0

                    if fc >= intensity_fold_change_threshold:
                        pep['call'][ip_index] = True
                    '''
                    if pep['peptide_counts'][ctrl] >= 0:
                        # High stringency path
                        if fc >= intensity_fold_change_threshold:
                            pep['call'] = True
                            
                    elif peps[pep_name]['peptide_counts'][ctrl] == 0:
                        # Low stringency path
                        pep['call'] = True
                    '''

        if cull_empty_calls:
            #print(self.linearData)
            #print([pep['call'] for pep in self.linearData])
            newl = [pep for pep in self.linearData if True in pep['call']]
            self.linearData = newl
        
        self._optimiseData()
        
        self.called = True
        self.fold_change = True
        config.log.info('call:')

    def sliceConditions(self,
        conditions:Iterable=None,
        **kargs):
        """
        **Purpose**

            return a copy of the expression-data, but only containing
            the condition names specified in conditions

            Note that you can also use this method to change the order of the conditions.
            Just slice all of the keys in the order you want them to occur in.

            Additionally, you can use this to replicate a key, e.g.

            gl = gl.sliceConditions(["cond1", "cond2", "cond1"])

            will now give you an expression set with two 'cond1' conditions

        **Arguments**
            conditions (Required)
                A list, or other iterable of condition names to extract
                from the expression-data. Every condition name must be present
                on the expression-data.

        **Result**

            A new expression-data object with the same settings as the original,
            but containing only the expression-data conditions specified in
            the 'conditions' argument.
        """
        # TODO: Implement
        pass

    def heatmap(self,
        dataset:str,
        filename:str,
        row_label_key:str ="name",
        row_color_threshold=None,
        optimal_ordering=True,
        dpi:int =300,
        _draw_supplied_cell_labels=None,
        cmap=None,
        bracket=None,
        **kargs):
        """
        **Purpose**

            draw a simple heatmap of the current massspec-data using the specified dataset.

        **Arguments**
            filename (Required)
                the filename of the image to save. depending upon the current
                drawing settings it will save either a png (default) svg or eps.
                
            dataset (Required)
                massspec objects have several datasets that could be plotted.
                You must specify one of:
                
                {}

            bracket (Optional, default = no bracketing performed)
                bracket the data within a certain range of values.
                For example to bracket to 0 .. 1 you would use the syntax::

                    result = array.heatmap(filename="ma.png", bracket=[0,1])

                Or for something like log2 normalised array data::

                    result = array.heatmap(filename="ma.png", bracket=[-2,2])

                "Bracket' chops off the edges of the data, using this logic::

                    if value > high_bracket then value := high_bracket
                    if value < low_bracket then value := low_bracket

                See normal for a method that modifies the data.

            row_label_key (Optional, default="name")
                A key in your genelist to use to label the rows. Examples would be gene names accesion
                numbers or something else.

            normal (Optional, default = no normalising)
                Unimplemented

            row_cluster (Optional, default = True)
                cluster the rows? True or False

            row_color_threshold (Optional, default=None)
                color_threshold to color the rows clustering dendrogram.

                See also scipy.hierarchy.dendrogram

            col_cluster (Optional, default = True)
                cluster the column conditions, True or False

            log (Optional, defualt=False, True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

            row_tree (Optional, default=False)
                provide your own tree for drawing. Should be a valid dendrogram tree.
                Probably the output from tree(), although you could role your own with Scipy.
                row_labels and the data
                will be rearranged based on the tree, so don't rearrange the data yourself.

            col_tree (Optional, default=False)
                provide your own tree for ordering the data by. See row_tree for details.
                This one is applied to the columns.

            highlights (Optional, default=None)
                sometimes the row_labels will be suppressed as there is too many labels on the plot.
                But you still want to highlight a few specific genes/rows on the plot.
                Send a list to highlights that matches entries in the row_names.

            digitize (Optional, default=False)
                change the colourmap (either supplied in cmap or the default) into a 'digitized' version
                that has large blocks of colours, defined by the number you send to discretize.
                In other words, place the expression values into the number of 'digitized' bins

            cmap (Optional, default=matplotlib.cm.RdBu)
                colour map for the heatmaps. Use something like this:

                import matplotlib.cm as cm

                gl.heatmap(..., cmap=cm.afmhot)

            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0

            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1

            row_font_size (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements.

            col_font_size (Optional, default=8)
                the size of the column labels (in points)

            heat_wid (Optional, default=0.25)
                The width of the heatmap panel. The image goes from 0..1 and the left most
                side of the heatmap begins at 0.3 (making the heatmap span from 0.3 -> 0.55).
                You can expand or shrink this value depending wether you want it a bit larger
                or smaller.

            heat_hei (Optional, default=0.85)
                The height of the heatmap. Heatmap runs from 0.1 to heat_hei, with a maximum of 0.9 (i.e. a total of 1.0)
                value is a fraction of the entire figure size.

            colbar_label (Optional, default="expression")
                the label to place beneath the colour scale bar

            grid (Optional, default=False)
                draw a grid around each cell in the heatmap.

            draw_numbers (Optional, default=False)
                draw the values of the heatmaps in each cell see also draw_numbers_threshold

            draw_numbers_threshold (Optional, default=-9e14)
                draw the values in the cell if > draw_numbers_threshold

            draw_numbers_fmt (Optional, default= '{{:.1f}}')
                string formatting for the displayed values

                You can also send arbitrary text here, (for example, if you wanted to
                mark significane with a '*' then you could set draw_numbers_fmt='*').

            draw_numbers_font_size (Optional, default=6)
                the font size for the numbers in each cell

            _draw_supplied_cell_labels (Optional, default=False)
                semi-undocumented function to draw text in each cell.

                Please provide a 2D list, with the same dimensions as the heatmap, and this text
                will be drawn in each cell. Useful for tings like drawing a heatmap of expression
                and then overlaying p-values on top of all significant cells.

            imshow (Optional, default=False)
                Embed the heatmap as an image inside a vector file. (Uses matplotlib imshow
                to draw the heatmap part of the figure. Allows very large matrices to
                be saved as an svg, with the heatmap part as a raster image and all other elements
                as vectors).

            sample_label_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the condition names

            col_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the condition names

            row_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the row names.

                Note that unclustered data goes from the bottom to the top!

            optimal_ordering (Optional, default=True)
                An improved ordering for the tree, at some computational and memory cost.
                Can be trouble on very large heatmaps

                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html

        **Result**
            saves an image to the 'filename' location and

            returns a dictionary containing the real_filename and re-ordered labels after clustering (if any).

            returns the 'actual filename' that really gets saved (glbase
            will modify e.g. a '.png' ending to '.svg' or '.eps' etc. depending
            upon the current setings).
        """.format(self.valid_datasets)
        
        # checks for option here please.
        assert filename, "heatmap: you must specify a filename"
        assert row_label_key in list(self.keys()), 'row_label_key "%s" not found in this genelist' % row_label_key
        assert dataset in self.valid_datasets, f'dataset {dataset} is not found in {self.valid_datasets}'
        if dataset == 'call': assert self.called, 'Asking for the call dataset, but this massspec has not been called, run massspec.call()'
        if dataset == 'fold_change': assert self.fold_change, 'Asking for the fold_change dataset, but this massspec does not have a fold_change, run massspec.call()'
        

        if dataset == 'intensities':
            data = numpy.array([p['intensities'] for p in self.linearData]).T
            if not cmap: cmap = 'copper'
            if not bracket: bracket = [0, 10]
            
        elif dataset == 'call':
            def booler(b):
                if b is None: return 0.1
                elif b is True: return 1
            data = numpy.array([[booler(i) for i in p['call']] for p in self.linearData]).T
            cmap1 = colors.LinearSegmentedColormap.from_list("grey-or", ["lightgrey", "orange"])
            if not cmap: cmap = cmap1
            if not bracket: bracket = [0, 1]
            
        elif dataset == 'fold_change':
            data = numpy.array([p['fold_change'] for p in self.linearData]).T
            if not cmap: cmap = 'RdBu_r'
            if not bracket: bracket = [-2, 2]

        if "log" in kargs:
            data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])

        # convert it into the serialisedArrayDataDict that _heatmap() expects.
        newdata = {}
        for index, name in enumerate(self._conditions):
            newdata[name] = data[index] # get the particular column

        res = self.draw.heatmap(data=newdata,
            cmap=cmap,
            bracket=bracket,
            row_names=self[row_label_key],
            col_names=self.getConditionNames(),
            filename=filename,
            row_color_threshold=row_color_threshold,
            optimal_ordering=optimal_ordering, dpi=dpi,
            _draw_supplied_cell_labels=_draw_supplied_cell_labels,
            **kargs)

        config.log.info("heatmap: Saved %s" % res["real_filename"])
        return res

    def correlation_heatmap(self, 
        dataset, 
        filename:str, 
        label_key:str =None, 
        mode:str ="r", 
        aspect="square", 
        bracket=(0,1),
        optimal_ordering:bool =True,
        **kargs):
        """
        **Purpose**
            Plot a heatmap of the (R, R^2, Pearson or Spearman) correlation for all pairs of
            samples in the expression object.

            Note that R and R^2 the calculation is very fast but pearson and spearman are quite a bit slower.

        **Arguments**
            filename (Required)
                the filename to save the image to.

            dataset (Required)
                massspec objects have several datasets that could be plotted.
                You must specify one of:
                
                {}

            label_key (Optional, but required if axis=genes or rows)
                You must specify a label key if the axis is rows or genes

            mode (Optional, default="r")

                by default the R (Coefficient of determination) is squared. Set to 'r' for
                Coefficient of determination value.
                use 'pearson' for a Pearson correlation score and 'spearman' for a
                Spearman-ranked score.

            bracket (optional, default=(0,1))
                bracket the heatmap by these values.

            Other heatmap options that will work:
                heat_hei
                heat_wid
                ...

        **Returns**
            A dict containing::
                "data": 2D numpy array containing the correlation scores (depending upon the mode)
                "labels": the labels along the top and bottom of the array, sorted according to the
                clustering (if any)

            Note that the arrays in 'data' are not clustered and will not be in the same order as the figure.
        """
        assert mode in ("r", "r2", "pearson", "spearman"), f"'{mode}' mode not found"
        assert filename, "correlation_heatmap: you must specify a filename"
        assert dataset in self.valid_datasets, f'correlation_heatmap: dataset {dataset} is not found in {self.valid_datasets}'
        
        labels = self._conditions
        if dataset == 'intensities':
            data_table = numpy.array([p['intensities'] for p in self.linearData]).T
            
        elif dataset == 'call':
            assert self.called, 'correlation_heatmap: Asking for the call dataset, but this massspec has not been called, run massspec.call()'
            def booler(b):
                if b is None: return 0.1
                elif b is True: return 1
            data_table = numpy.array([[booler(i) for i in p['call']] for p in self.linearData]).T
            
        elif dataset == 'fold_change':
            assert self.fold_change, 'correlation_heatmap: Asking for the fold_change dataset, but this massspec does not have a fold_change, run massspec.call()'
            data_table = numpy.array([p['fold_change'] for p in self.linearData]).T

        if mode in ("r", "r2"):
            arr = numpy.corrcoef(data_table) # PMCC or little r,
            if mode == "r2":
                arr *= arr

        else:
            arr = numpy.zeros((len(labels), len(labels)))
            ps = numpy.zeros((len(labels), len(labels)))
            p = progressbar(len(labels))
            for ic1, c1 in enumerate(labels):
                for ic2, c2 in enumerate(labels):
                    if ic1 == ic2: # stop identical labels mucking it up.
                        arr[ic1, ic2] = 1.0
                    else:
                        x = data_table[ic1,:]
                        y = data_table[ic2,:]

                        if mode == "pearson":
                            arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.pearsonr(x, y)
                        elif mode == "spearman":
                            arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.spearmanr(x, y)
                p.update(ic1)

        square = True
        if "heat_hei" in kargs or "heat_wid" in kargs:
            square=False

        results = self.draw.heatmap(filename=filename, data=arr, square=square,
            bracket=bracket, aspect=aspect, row_names=labels, col_names=labels,
            colbar_label=f"Correlation ({mode})", optimal_ordering=optimal_ordering, **kargs)
        config.log.info("correlation_heatmap: Saved '%s'" % results["real_filename"])
        
        return {"data": results["reordered_data"], "labels": results["reordered_cols"]}

    def volcano(self):
        pass
    
    def network(self, filename, expt_to_bait=None):
        """
        **Purpose**
            Draw a connected network PPI-style plot.
            
            Designed for Co-IP MS experiments
            
            Requires networkx
            
        **Arguments**
            filename (Required)
                filename to save the image of the network
                
            expt_to_bait (Required)
               A dict with the bait for each experiment. Use '-' for controls.
               
               
                
        **Returns**
            The networkx network
        """
        assert filename, 'filename is required'
        assert self.called, 'network can only be used if the data has been called'
        assert expt_to_bait, 'Requires a expt_to_bait argument'
        assert isinstance(expt_to_bait, dict), 'Must be a dict'
        
        import networkx as nx
        
        g = nx.Graph()
        
        # Nodes
        #for pep in self.linearData:
        #    g.add_node(pep['name'])

        co_interactors = []

        #Edges         
        for pep in self.linearData:
            for ip, call in zip(self._conditions, pep['call']):
                try:
                    if expt_to_bait[ip] == '-':
                       continue
                except KeyError:
                    config.log.warning(f'Treating {ip} as a control sample, skipping;') 
                
                if call:
                    g.add_edge(expt_to_bait[ip], pep['name'])
                    co_interactors.append(pep['name'])

        # Remove isolated nodes;
        #pos = nx.nx.kamada_kawai_layout(g)
        #pos = nx.circular_layout(g, scale=0.2)
        pos = nx.spring_layout(g, k=1, iterations=1000, seed=1234)

        fig = plot.figure(figsize=[8,8])
        plot.axis("off")
        
        ax = fig.add_subplot(111)
        ax.set_position([0,0, 1,1])
        nx.draw_networkx_edges(g, pos=pos, alpha=0.1)
        nx.draw_networkx_nodes(g, pos=pos, nodelist=co_interactors, alpha=0.6, node_size=100, node_color="tab:red")
        nx.draw_networkx_nodes(g, pos=pos, nodelist=baits, alpha=0.9, node_size=100, node_color="tab:orange")
        nx.draw_networkx_labels(g, pos=pos, font_size=7, alpha=0.9)
        fig.savefig(filename)
        plot.close(fig)

        return g