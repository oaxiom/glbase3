"""Classes and functions for working with datasets."""


import re
import copy
from itertools import groupby

import numpy as N

from .util import *
from . import discretizer
from . import config

#
# Module parameters
#
_pfilename = config.StringParameter(
    'data.filename',
    'File to read data from.',
    config.fileexists(),
)

_ptext = config.StringParameter(
    'data.text',
    'The text of a dataset included in config file.',
    default=''
)

_pdiscretize = config.IntParameter(
    'data.discretize',
    'Number of bins used to discretize data. Specify 0 to indicate that '+\
    'data should not be discretized.',
    default=0
)

#
# Exceptions
#
class ParsingError(Exception): 
    """Error encountered while parsing an ill-formed datafile."""
    pass

class IncorrectArityError(Exception):
    """Error encountered when the datafile speifies an incorrect variable arity.

    If variable arity is specified, it should be greater than the number of
    unique observation values for the variable.

    """
    
    def __init__(self, errors):
        self.errors = errors

    def __repr__(self):
        msg = "Incorrect arity specified for some variables.\n"
        for v,uniquevals in errors:
            msg += "Variable %s has arity of %d but %d unique values.\n" % \
                   (v.name, v.arity, uniquevals)

class ClassVariableError(Exception):
    """Error with a class variable."""
    msg = """Data for class variables must include only the labels specified in
    the variable annotation."""


#
# Variables and Samples
#
class Annotation(object):
    """Additional information about a sample or variable."""

    def __init__(self, name, *args):
        # *args is for subclasses
        self.name = str(name)

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,  self.name)

class Sample(Annotation):
    """Additional information about a sample."""
    pass 

class Variable(Annotation): 
    """Additional information about a variable."""
    arity = -1

class ContinuousVariable(Variable): 
    """A variable from a continuous domain."""
    def __init__(self, name, param):
        self.name = str(name)

class DiscreteVariable(Variable):
    """A variable from a discrete domain."""
    def __init__(self, name, param):
        self.name = str(name)
        self.arity = int(param)

class ClassVariable(DiscreteVariable):
    """A labeled, discrete variable."""
    def __init__(self, name, param):
        self.name = str(name)
        self.labels = [l.strip() for l in param.split(',')]
        self.label2int = {l: i for i,l in enumerate(self.labels)}
        self.arity = len(self.labels)

#
# Main class for dataset
#
class Dataset(object):
    def __init__(self, observations, missing=None, interventions=None, 
                 variables=None, samples=None, skip_stats=False):
        """
        Create a pebl Dataset instance.

        A Dataset consists of the following data structures which are all
        numpy.ndarray instances:

        * observations: a 2D matrix of observed values. 
            - dimension 1 is over samples, dimension 2 is over variables.
            - observations[i,j] is the observed value for jth variable in the ith
              sample.

        * missing: a 2D binary mask for missing values
            - missing[i,j] = 1 IFF observation[i,j] is missing
        
        * interventions: a 2D binary mask for interventions
            - interventions[i,j] = 1 IFF the jth variable was intervened upon in
              the ith sample.
        
        * variables,samples: 1D array of variable or sample annotations
        
        This class provides a few public methods to manipulate datasets; one can
        also use numpy functions/methods directly.

        Required/Default values:

             * The only required argument is observations (a 2D numpy array).
             * If missing or interventions are not specified, they are assumed to
               be all zeros (no missing values and no interventions).
             * If variables or samples are not specified, appropriate Variable or
               Sample annotations are created with only the name attribute.

        Note:
            If you alter Dataset.interventions or Dataset.missing, you must
            call Dataset._calc_stats(). This is a terrible hack but it speeds
            up pebl when used with datasets without interventions or missing
            values (a common case).

        """

        self.observations = observations
        self.missing = missing
        self.interventions = interventions
        self.variables = variables
        self.samples = samples
        
        obs = observations
        if missing is None:
            self.missing = N.zeros(obs.shape, dtype=bool)
        if interventions is None:
            self.interventions = N.zeros(obs.shape, dtype=bool)
        if variables is None:
            self.variables = N.array([Variable(str(i)) for i in range(obs.shape[1])])
            self._guess_arities()
        if samples is None:
            self.samples = N.array([Sample(str(i)) for i in range(obs.shape[0])])

        if not skip_stats:
            self._calc_stats()

    # 
    # public methods
    # 
    def subset(self, variables=None, samples=None):
        """Returns a subset of the dataset (and metadata).
        
        Specify the variables and samples for creating a subset of the data.
        variables and samples should be a list of ids. If not specified, it is
        assumed to be all variables or samples. 

        Some examples:
        
            - d.subset([3], [4])
            - d.subset([3,1,2])
            - d.subset(samples=[5,2,7,1])
        
        Note: order matters! d.subset([3,1,2]) != d.subset([1,2,3])

        """

        variables = variables if variables is not None else list(range(self.variables.size))
        samples = samples if samples is not None else list(range(self.samples.size))
        skip_stats = not (self.has_interventions or self.has_missing)
        d = Dataset(
            self.observations[N.ix_(samples,variables)],
            self.missing[N.ix_(samples,variables)],
            self.interventions[N.ix_(samples,variables)],
            self.variables[variables],
            self.samples[samples],
            skip_stats = skip_stats
        )
        
        # if self does not have interventions or missing, the subset can't.
        if skip_stats:
            d._has_interventions = False
            d._has_missing = False

        return d

    def _subset_ni_fast(self, variables):
        ds = _FastDataset.__new__(_FastDataset)

        if not self.has_interventions:
            ds.observations = self.observations[:,variables]
            ds.samples = self.samples
        else:
            samples = N.where(self.interventions[:,variables[0]] == False)[0] 
            ds.observations = self.observations[samples][:,variables]
            ds.samples = self.samples[samples]

        ds.variables = self.variables[variables]
        return ds

    # TODO: test
    def subset_byname(self, variables=None, samples=None):
        """Returns a subset of the dataset (and metadata).

        Same as Dataset.subset() except that variables and samples can be
        specified by their names.  
        
        Some examples:

            - d.subset(variables=['shh', 'genex'])
            - s.subset(samples=["control%d" % i for i in xrange(10)])

        """

        vardict = {v.name: i for i,v in enumerate(self.variables)}
        sampledict = {s.name: i for i,s in enumerate(self.samples)}

        # if name not found, we let the KeyError be raised
        variables = [vardict[v] for v in variables] if variables else variables
        samples = [sampledict[s] for s in samples] if samples else samples

        return self.subset(variables, samples)

    def discretize(self, includevars=None, excludevars=[], numbins=3):
        """Discretize (or bin) the data in-place.

        This method is just an alias for pebl.discretizer.maximum_entropy_discretizer()
        See the module documentation for pebl.discretizer for more information.

        """
        self.original_observations = self.observations.copy()
        self = discretizer.maximum_entropy_discretize(
           self, 
           includevars, excludevars, 
           numbins
        )

    def tofile(self, filename, *args, **kwargs):
        """Write the data and metadata to file in a tab-delimited format."""
        
        with file(filename, 'w') as f:
            f.write(self.tostring(*args, **kwargs))

    def tostring(self, linesep='\n', variable_header=True, sample_header=True):
        """Return the data and metadata as a string in a tab-delimited format.
        
        If variable_header is True, include variable names and type.
        If sample_header is True, include sample names.
        Both are True by default.

        """

        def dataitem(row, col):
            val = "X" if self.missing[row,col] else str(self.observations[row,col])
            val += "!" if self.interventions[row,col] else ''
            return val
        
        def variable(v):
            name = v.name

            if isinstance(v, ClassVariable):
                return "%s,class(%s)" % (name, ','.join(v.labels))    
            elif isinstance(v, DiscreteVariable):
                return "%s,discrete(%d)" % (name, v.arity)
            elif isinstance(v, ContinuousVariable):
                return "%s,continuous" % name
            else:
                return v.name

        # ---------------------------------------------------------------------

        # python strings are immutable, so string concatenation is expensive!
        # preferred way is to make list of lines, then use one join.
        lines = []

        # add variable annotations
        if sample_header:
            lines.append("\t".join([variable(v) for v in self.variables]))
        
        # format data
        nrows,ncols = self.shape
        d = [[dataitem(r,c) for c in range(ncols)] for r in range(nrows)]
        
        # add sample names if we have them
        if sample_header and hasattr(self.samples[0], 'name'):
            d = [[s.name] + row for row,s in zip(d,self.samples)]

        # add data to lines
        lines.extend(["\t".join(row) for row in d])
        
        return linesep.join(lines)


    #
    # public properties
    #
    @property
    def shape(self):
        """The shape of the dataset as (number of samples, number of variables)."""
        return self.observations.shape

    @property
    def has_interventions(self):
        """Whether the dataset has any interventions."""
        if not hasattr(self, '_has_interventions'):
            self._has_interventions = self.interventions.any()

        return self._has_interventions

    @property
    def has_missing(self):
        """Whether the dataset has any missing values."""
        if not hasattr(self, '_has_missing'):
            self._has_missing = self.missing.any()

        return self._has_missing

    #
    # private methods/properties
    #
    def _calc_stats(self):
        self._has_interventions = self.interventions.any()
        self._has_missing = self.missing.any()
    
    def _guess_arities(self):
        """Guesses variable arity by counting the number of unique observations."""

        for col,var in enumerate(self.variables):
            var.arity = N.unique(self.observations[:,col]).size
            var.__class__ = DiscreteVariable


    def check_arities(self):
        """Checks whether the specified airty >= number of unique observations.

        The check is only performed for discrete variables.

        If this check fails, the CPT and other data structures would fail.
        So, we should raise error while loading the data. Fail Early and Explicitly!

        """
        
        errors = [] 
        for col,v in enumerate(self.variables):
            if isinstance(v, DiscreteVariable):
                uniquevals = N.unique(self.observations[:,col]).size
                if v.arity < uniquevals:
                    errors.append((v, uniquevals))

        if errors:
            raise IncorrectArityError(errors)


class _FastDataset(Dataset):
    """A version of the Dataset class created by the _subset_ni_fast method.

    The Dataset._subset_ni_fast method creates a quick and dirty subset that
    skips many steps. It's a private method used by the evaluator module. Do
    not use this unless you know what you're doing.  
    
    """
    pass

def merge(datasets, axis=None):
    """Merges multiple datasets.

    datasets should be a list of Dataset objects.
    axis should be either 'variables' or 'samples' and determines how the
    datasets are merged.  
    
    """

    if axis == 'variables':
        variables = N.hstack(tuple(d.variables for d in datasets))
        samples = datasets[0].samples
        stacker = N.hstack
    else:
        samples = N.hstack(tuple(d.samples for d in datasets))
        variables = datasets[0].variables
        stacker = N.vstack

    missing = stacker(tuple(d.missing for d in datasets))
    interventions = stacker(tuple(d.interventions for d in datasets))
    observations = stacker(tuple(d.observations for d in datasets))

    return Dataset(observations, missing, interventions, variables, samples)


