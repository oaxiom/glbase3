
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
            
        '''.format(supported_formats)
        assert isinstance(filenames, list), 'filenames must be a list'
        assert format in supported_formats
        
        if format == 'maxquant_txt':
            peps = self.__maxquant_txt(filename, gzip, mq_use_lfq=mq_use_lfq)
        
        # Pack to a genelist-type object
    
        self.linearData = [{'name': k, 'PIDs': peps[k]['pids']} for k in peps]
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
        # unpack the intensity data based on the order in peps;
        self._intensities = [i['intensities'] for i n self.linearData]
        self._peptide_counts = [i['peptide_counts'] for i n self.linearData]
        self._call = [i['call'] for i n self.linearData]
        
    def __maxquant_txt(self, filename, gzip, mq_use_lfq):
        peps = {}
        sample_names = None
    
        if mq_use_lfq:
            intensity_search_key = 'Intensity '
        else:
            intensity_search_key = 'LFQ intensity '
    
        for filename in filenames:
            if gzip:
                oh1 = open(filename, 'rt')
            else:
                oh1 = open(filename, 'rt')
            
            # need to get the sample_names:
            
            for protein_match in oh1:
                if 'Protein IDs' in hit:
                    intensity_keys = {}
                    hit = protein_match.strip().split('\t')
                    for idx, lab in enumerate(hit):
                        if lab.startswith(intensity_search_key): 
                            intensity_keys[lab.replace(intensity_search_key, '')] = idx
                            
                    #print(intensity_keys)

                    unique_count_keys = {}
                    for idx, lab in enumerate(hit):
                        if 'Razor + unique peptides ' in lab:
                            unique_count_keys[lab.replace('Razor + unique peptides ', '')] = idx

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
                    intensity = hit[intensity_keys[sample_name]]
                    peps[pep_name]['intensities'][sample_name].append(float(intensity)/1e6)
                    unq_peps = hit[unique_count_keys[sample_name]]
                    peps[pep_name]['peptide_counts'][sample_name].append(int(unq_peps))
                    
        return peps, sample_names
        
    def call_matches(self, mode=None):
        '''
        '''
        
        self.called = True