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
from .expression import expression
from .draw import draw
from .progress import progressbar
from .errors import AssertionError, ArgumentError
from .utils import fold_change

def booler(b):
    r = []
    for i in b:
        if i is None: r.append(0.1)
        elif i is True: r.append(1)
        else: raise AssertionErorr(f"booler failed to collapse data on {i}")
    return r

class massspec(base_expression):
    supported_formats = {'maxquant_txt'}
    valid_call_modes = {'coip'}
    valid_datasets = {'intensities', 'call', 'fold_change'}
    supported_species = {
        'Hs': {'RPL', 'RPS', 'TUB', 'GAPDH', 'ACT',
            'MRPS', 'sm-', 'MYH', 'FLN', 'MYO', 'APOA1',
            'MYL',
            'IGHG', 'IGLV', 'IGHV', 'IGKV',
            'COL',
            'KRT', 'EIF', 'ATP'},
        'Mm': {'Rpl', 'Rps', 'Tub', 'Gapdh', 'Act', 'Myh',
            'Ighg', 'Iglv', 'Col', 'Golga', 'Kif', 'Myl',
            'Krt', 'Eif', 'Vim', 'Atp', 'Igkv', 'Ighv'},
        }

    def __init__(self,
        filenames=None,
        format=None,
        gzip:bool = False,
        mq_use_lfq:bool = False,
        mq_split_ambiguous_names:bool = True,
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

            mq_split_ambiguous_names (Optional, default=True)
                MaxQuant reports ambiguous peptide matches with PIDs and names like this:

                Gene1;Gene2;Gene3, PID1; PID2; PID3

                If mq_use_lfq is True then each ambiguous match is treated as a full match
                (generally more useful, as map() will correctly work), if False then preserve each
                matching protein entry as is.

        '''.format(self.supported_formats)
        assert isinstance(filenames, list), 'filenames must be a list'
        assert format in self.supported_formats

        self.name = name

        if format == 'maxquant_txt':
            peps, sample_names = self.__maxquant_txt(filenames, gzip, mq_use_lfq=mq_use_lfq, mq_split_ambiguous_names=mq_split_ambiguous_names)
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

    def __maxquant_txt(self, filenames, gzip, mq_use_lfq, mq_split_ambiguous_names):
        peps = {}
        sample_names = []

        if mq_use_lfq:
            intensity_search_key = 'LFQ intensity '
        else:
            intensity_search_key = 'Intensity '

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

        # Now load the peptide matches:
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

                if not hit[6]: # Unidentified, just skip it;
                    continue

                if mq_split_ambiguous_names and ';' in pep_name:
                    pep_name = pep_name.split(';')
                    pids = pids.split(';')
                else:
                    pep_name = [pep_name,]
                    pids = [pids, ]

                for name, pid in zip(pep_name, pids):
                    if name not in peps:
                        peps[name] = {
                            'pids': pid,
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
                        peps[name]['intensities'][sample_name].append(float(intensity)/1e6)
                        unq_peps = hit[unique_count_keys[sample_name]]
                        peps[name]['peptide_counts'][sample_name].append(int(unq_peps))

            oh1.close()

        return peps, sample_names

    def filter(self,
        mode,
        minimum_unique_peptides,
        minimum_intensity,
        species,
        remove_fusions=True,
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

            species (Required)
                Hs : Human
                Mm : Mouse

            remove_fusions (Optional, default=True)
                remove fusion genes (i.e. with a '-' character)

        **Returns**
            number of peptides deleted

        '''.format(self.valid_call_modes)
        assert mode in self.valid_call_modes, f'mode "{mode}" not found in {self.valid_call_modes}'
        assert species in self.supported_species, f"species {species} not found in {self.supported_species}"

        starting_len = len(self)

        if mode == 'coip':
            newl = self.__filter_coipms(minimum_unique_peptides=minimum_unique_peptides,
                minimum_intensity=minimum_intensity,
                species=species,
                remove_fusions=remove_fusions)

        if not newl:
            raise AssertionError('All peptides have been filtered out, list is empty')

        self.linearData = newl
        self._optimiseData()

        config.log.info(f'filter: Trimmed {starting_len - len(self)} possible false positives and peptides not meeting thresholds')
        config.log.info(f'filter: List now {len(self)} peptides long')
        return starting_len - len(self)

    def __filter_coipms(self, minimum_unique_peptides, minimum_intensity, species, remove_fusions):
        newl = []

        for pep in self.linearData:
            pep_name = pep['name']
            # Simple filter for common contaminants
            if remove_fusions and '-' in pep_name: continue

            if True in [pep_name.startswith(i) for i in self.supported_species[species]]: continue

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
        start_len = len(self)
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
        config.log.info(f'call: {len(self)} peptides met threshold, {start_len-len(self)} peptides removed')

    def sliceConditions(self,
        conditions:Iterable=None,
        **kargs):
        """
        **Purpose**

            return a copy of the expression-data, but only containing
            the condition names specified in conditions

            Note that you can also use this method to change the order of the conditions.
            Just slice all of the keys in the order you want them to occur in.

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
        assert conditions, "sliceConditions: You must specify a list of conditions to keep"
        assert not isinstance(conditions, str), "sliceConditions: You must specify an iterable of conditions to keep"
        assert isinstance(conditions, (tuple, list, set, Iterable)), "sliceConditions: You must specify an iterable of conditions to keep"
        conditions = list(conditions) # Some weird bugs if not a list;
        for item in conditions:
            assert item in self._conditions, f"sliceConditions: '{item}' condition not found in this expression data"

        newgl = self.deepcopy()
        newl = []

        # Massspec objects don't use numpy, so do the hard way:
        idxs_to_keep = [self._conditions.index(c) for c in conditions]
        for pep in self.linearData:
            newi = dict(pep)
            newi['intensities'] = [pep['intensities'][idx] for idx in idxs_to_keep]
            newi['peptide_counts'] = [pep['peptide_counts'][idx] for idx in idxs_to_keep]
            newi['call'] = [pep['call'][idx] for idx in idxs_to_keep]
            newi['fold_change'] = [pep['fold_change'][idx] for idx in idxs_to_keep]

            newl.append(newi)

        newgl.linearData = newl
        newgl._conditions = conditions
        newgl._optimiseData()

        config.log.info(f'sliceConditions: sliced for {len(newgl[0]["call"])} conditions')
        return newgl

    def merge_replicates(self, merge_data):
        """
        **Purpose**
            merge replicate MS samples.

        **Arguments**
            merge_data (Required)
                merge dict of lists:
                    {
                    'condition1_name': [cond1, cond2', ... 'condN'],
                    'condition2_name': [cond1, cond2', ... 'condN'],
                    etc
                    }

        **Returns**
            None
        """
        assert merge_data, 'You must specify merge_data'
        assert isinstance(merge_data, dict), 'merge_data must be a dict'

        newe = []
        all_conds = self._conditions
        done = set([])

        all_reps = set(sum([merge_data[x] for x in merge_data], []))
        samples_by_keepname = {}
        for k in merge_data:
            for cond_name in merge_data[k]:
                samples_by_keepname[cond_name] = k

        # Possible problem if names are not unique;
        assert len(samples_by_keepname) == len(all_reps), 'merge_data condition names are not unique!'
        assert len(sum([merge_data[x] for x in merge_data], [])) == len(all_reps), 'merge_data condition names are not unique!'

        if False in [c in self._conditions for c in all_reps]:
            missing_conds = [c for c in all_reps if c not in self._conditions]
            raise AssertionError("merge_replicates: '{}' condition names not found".format(", ".join(sorted(missing_conds))))

        new_serialisedArrayDataDict = {}
        errors = {}
        new_condition_name_list = sorted(list(merge_data.keys()))

        newl = self.deepcopy()
        # Blank the metric keys;
        for pep in newl:
            pep['intensities'] = []
            pep['peptide_counts'] = []
            pep['call'] = []
            pep['fold_change'] = []

        # max the merges;
        for new_merged_condition_name in new_condition_name_list:
            idxs_to_keep = [self._conditions.index(c) for c in merge_data[new_merged_condition_name]]

            for newi, original in zip(newl, self.linearData):
                newi['intensities'].append(max([original['intensities'][idx] for idx in idxs_to_keep]))
                newi['peptide_counts'].append(max([original['peptide_counts'][idx] for idx in idxs_to_keep]))
                newi['call'].append(True if True in [original['call'][idx] for idx in idxs_to_keep] else None)
                newi['fold_change'].append(max([original['fold_change'][idx] for idx in idxs_to_keep]))

        newl._conditions = new_condition_name_list
        newl._optimiseData()

        config.log.info(f"merge_replicates: Started with {len(self._conditions)} conditions, ended with {len(newl[0]['intensities'])}")
        return newl

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

            colbar_label (Optional, default=same as dataset argument)
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

        import matplotlib.colors as colors

        if dataset == 'intensities':
            data = numpy.array([p['intensities'] for p in self.linearData]).T
            if not cmap: cmap = colors.LinearSegmentedColormap.from_list("grey-or", ["lightgrey", "orange", 'red', 'purple'])
            if not bracket: bracket = [0, 20]

        elif dataset == 'call':
            data = numpy.array([booler(p['call']) for p in self.linearData]).T
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

        if 'colbar_label' in kargs and kargs['colbar_label']:
            pass
        else:
            kargs['colbar_label'] = dataset.capitalize()

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
            data_table += 0.01
            numpy.log2(data_table)

        elif dataset == 'call':
            assert self.called, 'correlation_heatmap: Asking for the call dataset, but this massspec has not been called, run massspec.call()'
            data_table = numpy.array([booler(p['call']) for p in self.linearData]).T

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

    def network(self,
        filename:str,
        expt_to_bait=None,
        figsize=[8,8],
        fontsize=6,
        **kargs):
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

            fontsize (Optional, default=6)
                fontsize for the nodes;

        **Returns**
            The networkx network
        """
        assert filename, 'filename is required'
        assert self.called, 'network can only be used if the data has been called'
        assert expt_to_bait, 'Requires a expt_to_bait argument'
        assert isinstance(expt_to_bait, dict), 'Must be a dict'

        import networkx as nx

        g = nx.Graph()
        __already_warned_using_as_bait = set([])

        # Nodes
        #for pep in self.linearData:
        #    g.add_node(pep['name'])

        co_interactors = []

        baits = [b for b in expt_to_bait.values() if b != '-']

        #Edges
        for pep in self.linearData:
            for ip, call in zip(self._conditions, pep['call']):
                try:
                    if expt_to_bait[ip] == '-':
                       continue

                except KeyError:
                    if ip not in __already_warned_using_as_bait:
                        __already_warned_using_as_bait.add(ip)
                        config.log.warning(f'Treating "{ip}" as a control sample. Add {{"{ip}": "-"}} to expt_to_bait to remove this warning')

                if call:
                    if expt_to_bait[ip] == pep['name']:
                        continue
                    g.add_edge(expt_to_bait[ip], pep['name'])
                    co_interactors.append(pep['name'])

        # Check all the baits made it into the network
        baits = [b for b in baits if b in g]

        config.log.info('Built network')
        # Remove isolated nodes;
        pos = nx.spring_layout(g, k=1, iterations=200, seed=1234)
        config.log.info('Layout done')

        fig = self.draw.getfigure(figsize=figsize)

        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.set_position([0,0, 1,1])
        nx.draw_networkx_edges(g, pos=pos, alpha=0.1)
        nx.draw_networkx_nodes(g, pos=pos, nodelist=set(co_interactors), alpha=0.4, node_size=200, linewidths=0, node_color="tab:red")
        nx.draw_networkx_nodes(g, pos=pos, nodelist=baits, alpha=0.9, node_size=400, linewidths=0, node_color="tab:orange")
        nx.draw_networkx_labels(g, pos=pos, font_size=fontsize, alpha=0.9)

        self.draw.do_common_args(ax, **kargs)
        real_filename = self.draw.savefigure(fig, filename)

        return g

    def venn(self, filename, *conditions, **kargs):
        """
        **Purpose**
            Perform a Venn diagram between 2-5 condition names

        **Arguments**
            filename (Required)
                filename to save the Venn to

            *conditions (Required)
                A list between 2-5 conditions names (must be in this massspec object)
                to use for the Venn overlap.

            **kargs (Optional)
                mainly to modify the figure. See glbase3/draw.py/do_common_args()

        **Returns**
            Sublists containing the various overlaps, contained in a dict
            Simple enough in the 2-venn, gets very complex in the 4-venn.
        """
        assert self.called, 'massspec must be called for this method to work'
        assert filename, 'you must provide a filename'
        assert (len(conditions) >= 2 and len(conditions) <= 5), 'conditions not >= 2 or <= 5, which are the only supported Venn sizes'
        assert False not in [c in self._conditions for c in conditions], 'A condition was not found in this massspec object'
        assert len([i['name'] for i in self.linearData]) == len(set([i['name'] for i in self.linearData])), '"name" key is not unqiue!'

        if len(conditions) == 2:
            Aidx = self._conditions.index(conditions[0])
            Bidx = self._conditions.index(conditions[1])
            Acalls = set([i['name'] for i in self.linearData if i['call'][Aidx]])
            Bcalls = set([i['name'] for i in self.linearData if i['call'][Bidx]])
            ABcalls = Acalls & Bcalls
            A = len(Acalls)
            B = len(Bcalls) # overlap subtracted in venn2
            AB = len(ABcalls)

            self.draw.venn2(A, B, AB, conditions[0], conditions[1], filename, **kargs)

            # Subtract overlap here so the returns are correct;
            Acalls -= ABcalls
            Bcalls -= ABcalls

            # get return versions:
            Aret = self.deepcopy()
            Aret.linearData = [pep for pep in Aret if pep['name'] in Acalls]
            Aret._optimiseData()

            Bret = self.deepcopy()
            Bret.linearData = [pep for pep in Bret if pep['name'] in Bcalls]
            Bret._optimiseData()

            ABret = self.deepcopy()
            ABret.linearData = [pep for pep in ABret if pep['name'] in ABcalls]
            ABret._optimiseData()

            return {'left': Aret, 'overlap': ABret, 'right': Bret}

        if len(conditions) == 3:
            raise NotImplementedError('Venn3 not implemented')

        if len(conditions) == 4:
            raise NotImplementedError('Venn4 not implemented')

        if len(conditions) == 5:
            raise NotImplementedError('Venn5 not implemented')

        return None # Should not reach

    def tree(self, mode="conditions",
        dataset=None,
        row_label_key=None,
        filename=None, row_name_key=None,
        cluster_mode="euclidean", color_threshold=None, label_size=6, cut=False,
        radial=False, optimal_ordering=True,
        **kargs):
        """
        **Purpose**
            Draw a hierarchical clustered tree of either the 'conditions' or 'rows'

        **Arguments**
            filename (Required)
                filename to save an image to.

            mode (Optional, default="conditions")
                cluster either the "conditions" or the "rows" or "genes". (rows are usually genes)

            dataset (Required)
                massspec objects have several datasets that could be plotted.
                You must specify one of:

                'call', 'intensities', 'fold-change'

            row_name_key (Optional, default=None)
                A key to use for the row_names, if mode == "row"

            cluster_mode (Optional, default="euclidean")
                A metric to cluster the data by.

            color_threshold (Optional, default=None)
                By default tree() uses the Scipy/MATLAB default colours to signify
                links with similar distances. Set to -1 to change the tree to all blue.

                see scipy.cluster.hierarch.dendrogram

            label_size (Optional, default=7)
                The size of the text attached to the leaf labels.

            cut (Optional, default=False)
                the cut threshold on the tree, to select groups and return the cut.
                Will draw a line on the dendrgram tree indicating its position.
                The value should range from 0..1.0

                Note that if cut is used then row_name_key must also be valid

                Also note that the order of the clusters returned is random.

            bracket (Optional, default=False)
                The min and max to bracket the data. This is mainly here for compatibility with heatmap()
                So that the row/col trees in heatmap and here exactly match

            optimal_ordering (Optional, default=True)
                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html

        **Returns**
            if cut is False:
                The Tree in a dictionary {"dendrogram": ..., "linkage": ..., "distance": ...}
            if cut is True:
                A list of genelists, containing all of the items from the cut.

        """
        assert filename, "heatmap: you must specify a filename"
        assert row_label_key in list(self.keys()), f'row_label_key "{row_label_key}" not found in this genelist'
        assert dataset in self.valid_datasets, f'dataset {dataset} is not found in {self.valid_datasets}'
        if dataset == 'call': assert self.called, 'Asking for the call dataset, but this massspec has not been called, run massspec.call()'
        if dataset == 'fold_change': assert self.fold_change, 'Asking for the fold_change dataset, but this massspec does not have a fold_change, run massspec.call()'

        if dataset == 'call':
            _data = [booler(i['calls']) for i in self.linearData]
            _data = numpy.array(_data).T # Transpose to get the expected shape for expression.tree()
        elif dataset == 'intensities':
            _data = [i['intensities'] for i in self.linearData]
            _data = numpy.array(_data).T # Transpose to get the expected shape for expression.tree()
        else:
            raise AssertionError(f'No method implemented for {dataset}')

        ret = expression.tree(self, mode=mode,
            _data=_data, # send out own custom matrix
            filename=filename, row_name_key=row_name_key,
            cluster_mode="euclidean", color_threshold=color_threshold, label_size=label_size, cut=cut,
            radial=radial, optimal_ordering=optimal_ordering,
            **kargs)

        return ret
