
import gzip
from statistics import mean

from . import config, utils
from .base_expression import base_expression
from .draw import draw
from .progress import progressbar
from .errors import AssertionError, ArgumentError

class massspec(base_expression):
    supported_formats = {'maxquant_txt'}

    def __init__(self, 
        filenames=None, 
        format=None, 
        gzip:bool = False, 
        mq_use_lfq:bool = False
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
            
            # Options for MaxQuant:
            mq_use_lfq (Optional, default=False)
                MaxQuant gives two intensity scores, LFQ intensity estimates the amount of protein based
                on a pool of proteins that don't change much. Whilst intensity is the raw intensity.
                LFQ is preferred for whole cell MS, whilst without is probably better for 
                enrichment-based MS, e.g. Co-IP. Use this switch to use lfq or not.
            
            gzip (Optional)
                Input file is gzipped
            
        '''.format(self.supported_formats)
        assert isinstance(filenames, list), 'filenames must be a list'
        assert format in self.supported_formats
        
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
            new_linear_data.append(new_peptide)
        
        # Pack to a genelist-type object
        
        self._conditions = sample_names
        self.linearData = new_linear_data
        self._optimiseData()
        
        self.called = False # call is valid;

    def __repr__(self):
        return "<massspec class>"

    def sort(self): # To implement
        raise NotImplementedError

    def find(self): # To implement
        raise NotImplementedError

    def from_pandas(self): # Not to implement
        raise NotImplementedError

    def _optimiseData(self):
        '''
        Needs a custom one
        '''
        print(self.linearData)
        # unpack the intensity data based on the order in peps;
        self._intensities = [i['intensities'] for i in self.linearData]
        self._peptide_counts = [i['peptide_counts'] for i in self.linearData]
        self._call = [i['call'] for i in self.linearData]

    def saveCSV(self, filename=None, interleave_errors=True, no_header=False, no_col1_header=False, **kargs):
        """
        A CSV version of saveTSV(), see saveTSV() for syntax
        """
        self.saveTSV(filename=filename, tsv=False, interleave_errors=True, no_header=False, no_col1_header=False, **kargs)
        config.log.info(f"saveCSV(): Saved '{filename}'")
        
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
        self._save_TSV_CSV(filename=filename, tsv=tsv, no_header=False, no_col1_header=False, **kargs)
        config.log.info(f"saveTSV(): Saved '{filename}'")
        
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
            array_data_keys = ("conditions", "err", "cv_err")

            write_keys = []
            if "key_order" in kargs:
                write_keys = kargs["key_order"]
                # now add in any missing keys to the right side of the list:
                for item in list(self.keys()):
                    if item not in write_keys and item not in array_data_keys: # But omit the array_data_keys
                        write_keys.append(item)
            else:
                        # just select them all:
                write_keys = [k for k in list(self.keys()) if k not in array_data_keys]

            if "err" in list(self.keys()):
                if interleave_errors:
                    conds = ["mean_%s" % c for c in self.getConditionNames()]
                    errs = ["err_%s" % c for c in self.getConditionNames()]
                    paired = [val for pair in zip(conds, errs) for val in pair]

                    if not no_header:
                        title_row = [k for k in write_keys if k in list(self.keys())]
                        writer.writerow(title_row + paired)

                    for data in self.linearData:
                        line = [data[k] for k in write_keys if k in data]

                        interleaved_data = [val for pair in zip(data["conditions"], data["err"]) for val in pair] # I never understand how these work, but what the hell.

                        writer.writerow(line + interleaved_data)# conditions go last.
                else:
                    if not no_header:
                        title_row = [k for k in write_keys in k in list(self.keys())]
                        writer.writerow(write_keys + self.getConditionNames() + ["err_%s" % i for i in self.getConditionNames()])

                    for data in self.linearData:
                        line = [data[k] for k in write_keys if k in data]
                        writer.writerow(line + data["conditions"] + data["err"])# conditions go last.
            else: # no error, very easy:
                if not no_header:
                    title_row = [k for k in write_keys if k in list(self.keys())]
                    if no_col1_header:
                        title_row[0] = ""
                    writer.writerow(title_row + self.getConditionNames())

                for data in self.linearData:
                    line = [data[k] for k in write_keys if k in data]
                    writer.writerow(line + data["conditions"])# conditions go last.
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
        
    def call_matches(self, mode=None):
        '''
        '''
        
        self.called = True